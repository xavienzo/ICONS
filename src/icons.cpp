// Core compiled routines for ICONS.
//
// Two things live here because they dominate the run time and are awkward to
// express efficiently in R:
//
//   1. greedy peeling with l0 shrinkage, and the adaptive loop around it.
//      The thresholded matrix is stored once in compressed sparse row form,
//      so a peel costs O(nnz + m log m) rather than O(m^2), and the adaptive
//      loop reuses that one structure instead of copying a shrinking
//      submatrix on every round.
//
//   2. the block statistics that drive SCFA.  Both the community sums and the
//      per-variable sums of squares come out of a single pass over the data,
//      which is what lets the estimators avoid ever forming a p-by-p matrix.
//
// Plain Rcpp only: no BLAS/LAPACK and no Fortran, so the package builds
// anywhere Rcpp does.

#include <Rcpp.h>
#include <vector>
#include <queue>
#include <cmath>
#include <limits>

using namespace Rcpp;

namespace {

// Compressed sparse row storage of the thresholded matrix.  Entries strictly
// below `threshold` are dropped, matching W[W < threshold] <- 0, and the
// diagonal is always dropped.
struct Csr {
  std::vector<int> ptr;
  std::vector<int> idx;
  std::vector<double> val;
};

Csr build_csr(const NumericMatrix& W, double threshold) {
  const int n = W.nrow();
  Csr g;
  g.ptr.assign(n + 1, 0);

  for (int i = 0; i < n; ++i) {
    int cnt = 0;
    for (int j = 0; j < n; ++j) {
      if (j == i) continue;
      double w = W(i, j);
      if (!ISNAN(w) && w >= threshold && w != 0.0) ++cnt;
    }
    g.ptr[i + 1] = g.ptr[i] + cnt;
  }

  g.idx.resize(g.ptr[n]);
  g.val.resize(g.ptr[n]);

  for (int i = 0; i < n; ++i) {
    int at = g.ptr[i];
    for (int j = 0; j < n; ++j) {
      if (j == i) continue;
      double w = W(i, j);
      if (!ISNAN(w) && w >= threshold && w != 0.0) {
        g.idx[at] = j;
        g.val[at] = w;
        ++at;
      }
    }
  }
  return g;
}

// Heap entry for the lazy-deletion priority queue.  Ordered by degree, ties
// broken on the smaller node index so the result matches which.min() in R and
// is reproducible across platforms.
struct Item {
  double deg;
  int v;
};

struct ItemGreater {
  bool operator()(const Item& a, const Item& b) const {
    if (a.deg != b.deg) return a.deg > b.deg;
    return a.v > b.v;
  }
};

// One greedy peeling pass restricted to `active`.
//
// Repeatedly strips the minimum-degree node and scores the surviving set with
// the generalised density  w(S) / |S|^(2 * lambda)  of Tsourakakis et al.
// (2013) as used by Chen et al. (2024).  `keep` is filled with the nodes of
// the best-scoring set in reverse removal order (densest core first) and
// `drop` with the nodes peeled away before it, both as original node indices.
//
// Returns the best score, or NA_REAL when `active` has fewer than two nodes.
double peel_active(const Csr& g,
                   const std::vector<int>& active,
                   double lambda,
                   std::vector<int>& keep,
                   std::vector<int>& drop) {
  const int m = static_cast<int>(active.size());
  keep.clear();
  drop.clear();

  if (m == 0) return NA_REAL;
  if (m == 1) {
    keep.push_back(active[0]);
    return NA_REAL;
  }

  const int n = static_cast<int>(g.ptr.size()) - 1;
  std::vector<char> alive(n, 0);
  std::vector<double> deg(n, 0.0);
  for (int a = 0; a < m; ++a) alive[active[a]] = 1;

  double total = 0.0;  // total edge weight inside the surviving set
  for (int a = 0; a < m; ++a) {
    const int v = active[a];
    double d = 0.0;
    for (int e = g.ptr[v]; e < g.ptr[v + 1]; ++e) {
      if (alive[g.idx[e]]) d += g.val[e];
    }
    deg[v] = d;
    total += d;
  }
  total *= 0.5;

  std::priority_queue<Item, std::vector<Item>, ItemGreater> pq;
  for (int a = 0; a < m; ++a) pq.push(Item{deg[active[a]], active[a]});

  std::vector<int> order;
  order.reserve(m);

  double best = -std::numeric_limits<double>::infinity();
  int best_cut = -1;  // number of nodes removed at the best score

  for (int step = 0; step < m; ++step) {
    int v = -1;
    while (!pq.empty()) {
      const Item it = pq.top();
      pq.pop();
      // Skip entries made stale by a later degree update.
      if (alive[it.v] && it.deg == deg[it.v]) { v = it.v; break; }
    }
    if (v < 0) break;  // unreachable, but keeps the loop total

    alive[v] = 0;
    order.push_back(v);
    total -= deg[v];

    for (int e = g.ptr[v]; e < g.ptr[v + 1]; ++e) {
      const int u = g.idx[e];
      if (!alive[u]) continue;
      deg[u] -= g.val[e];
      pq.push(Item{deg[u], u});
    }

    const int remaining = m - step - 1;
    if (remaining >= 1) {
      const double score = total / std::pow(static_cast<double>(remaining), 2.0 * lambda);
      if (score > best) {          // strict: first maximum wins, as which.max does
        best = score;
        best_cut = step + 1;
      }
    }
  }

  if (best_cut < 0) {  // every score was degenerate; keep the whole set
    for (int a = static_cast<int>(order.size()) - 1; a >= 0; --a) keep.push_back(order[a]);
    return NA_REAL;
  }

  drop.assign(order.begin(), order.begin() + best_cut);
  for (int a = static_cast<int>(order.size()) - 1; a >= best_cut; --a) keep.push_back(order[a]);
  return best;
}

}  // namespace

// Single greedy peeling pass (compiled)
//
// @param W Numeric adjacency matrix.
// @param lambda Shrinkage exponent.
// @param threshold Edge threshold applied before peeling.
// @return List with `keep`, `drop` (1-based node indices) and `score`.
// @keywords internal
// [[Rcpp::export(.peel_cpp)]]
List peel_cpp(NumericMatrix W, double lambda, double threshold) {
  const int n = W.nrow();
  const Csr g = build_csr(W, threshold);

  std::vector<int> active(n);
  for (int i = 0; i < n; ++i) active[i] = i;

  std::vector<int> keep, drop;
  const double score = peel_active(g, active, lambda, keep, drop);

  IntegerVector rkeep(keep.size()), rdrop(drop.size());
  for (size_t i = 0; i < keep.size(); ++i) rkeep[i] = keep[i] + 1;
  for (size_t i = 0; i < drop.size(); ++i) rdrop[i] = drop[i] + 1;

  return List::create(_["keep"] = rkeep, _["drop"] = rdrop, _["score"] = score);
}

// Adaptive dense subgraph extraction (compiled)
//
// @param W Numeric adjacency matrix.
// @param threshold Edge threshold.
// @param lambda Shrinkage exponent.
// @param min_size Smallest community that is retained.
// @param max_k Stop after this many communities; `<= 0` means no limit.
// @return List with `order`, `sizes` and `scores`.
// @keywords internal
// [[Rcpp::export(.detect_cpp)]]
List detect_cpp(NumericMatrix W, double threshold, double lambda,
                int min_size, int max_k) {
  const int n = W.nrow();
  const Csr g = build_csr(W, threshold);

  std::vector<int> active(n);
  for (int i = 0; i < n; ++i) active[i] = i;

  std::vector<int> order;
  order.reserve(n);
  std::vector<int> sizes;
  std::vector<double> scores;
  std::vector<int> keep, drop;

  while (static_cast<int>(active.size()) >= min_size) {
    if (max_k > 0 && static_cast<int>(sizes.size()) >= max_k) break;

    Rcpp::checkUserInterrupt();
    const double score = peel_active(g, active, lambda, keep, drop);

    // Stop once the densest remaining subgraph is too small to be a factor;
    // whatever is left becomes the singleton set.
    if (static_cast<int>(keep.size()) < min_size) break;
    // Nothing was peeled off: no progress is possible, so stop.
    if (drop.empty()) {
      order.insert(order.end(), keep.begin(), keep.end());
      sizes.push_back(static_cast<int>(keep.size()));
      scores.push_back(score);
      active.clear();
      break;
    }

    order.insert(order.end(), keep.begin(), keep.end());
    sizes.push_back(static_cast<int>(keep.size()));
    scores.push_back(score);
    active = drop;
  }

  // Everything not assigned to a community, in ascending node order.
  std::vector<char> taken(n, 0);
  for (size_t i = 0; i < order.size(); ++i) taken[order[i]] = 1;
  const int n_single = n - static_cast<int>(order.size());
  for (int i = 0; i < n; ++i) if (!taken[i]) order.push_back(i);

  IntegerVector rorder(n);
  for (int i = 0; i < n; ++i) rorder[i] = order[i] + 1;

  return List::create(
    _["order"] = rorder,
    _["sizes"] = wrap(sizes),
    _["scores"] = wrap(scores),
    _["n_singletons"] = n_single);
}

// Block statistics for SCFA (compiled)
//
// One pass over the data returning the centred within-community sums
// (`T`, n by K) and the centred per-variable sums of squares (`ss`, length p).
// Everything the closed-form estimators need follows from these, which is why
// no p-by-p matrix is ever built.
//
// @param X Numeric data matrix, observations in rows.
// @param center Numeric vector of column means.
// @param g Integer vector of 1-based community labels, `0` for unassigned.
// @param K Number of communities.
// @return List with `T` and `ss`.
// @keywords internal
// [[Rcpp::export(.block_stats_cpp)]]
List block_stats_cpp(NumericMatrix X, NumericVector center, IntegerVector g, int K) {
  const int n = X.nrow();
  const int p = X.ncol();

  NumericMatrix Tk(n, K);
  NumericVector ss(p);

  for (int j = 0; j < p; ++j) {
    const double mu = center[j];
    const int k = g[j] - 1;  // -1 when unassigned
    double s2 = 0.0;
    if (k >= 0 && k < K) {
      for (int i = 0; i < n; ++i) {
        const double z = X(i, j) - mu;
        s2 += z * z;
        Tk(i, k) += z;
      }
    } else {
      for (int i = 0; i < n; ++i) {
        const double z = X(i, j) - mu;
        s2 += z * z;
      }
    }
    ss[j] = s2;
  }

  return List::create(_["T"] = Tk, _["ss"] = ss);
}
