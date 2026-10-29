"""Depth-neutral stage-to-stage change: thin the deeper stage to the shallower stage's depth before fitting.

Usage (numpy / scipy / anndata only):
    import anndata as ad
    from depth_neutral import recover_counts, to_log, depth_neutral_delta
    prev, donor = ad.read_h5ad("E8.5_RNA.h5ad"), ad.read_h5ad("E9.5_RNA.h5ad")
    delta, ratio = depth_neutral_delta(prev.X, donor.X, prev.obs["celltype"].values, donor.obs["celltype"].values)
    # delta: {celltype: per-gene change in log1p(CP10K), measured at the donor stage's depth}
"""
import numpy as np, scipy.sparse as sp


def recover_counts(Xlog):
    """log1p(CP10K) -> integer UMI counts (per cell, 1 UMI = the smallest non-zero expm1 value)."""
    M = sp.csr_matrix(Xlog, dtype=np.float64, copy=True)
    M.data = np.expm1(M.data); M.data[M.data < 1e-9] = 0; M.eliminate_zeros()
    nnz = np.diff(M.indptr); assert (nnz > 0).all(), "empty cell"
    unit = np.minimum.reduceat(M.data, M.indptr[:-1])
    c = M.data / unit[np.repeat(np.arange(M.shape[0]), nnz)]
    err = float(np.abs(c - np.round(c)).max()); M.data = np.round(c)
    return M, err


def to_log(C, target=1e4):
    C = sp.csr_matrix(C); s = np.asarray(C.sum(1)).ravel()
    L = sp.csr_matrix(sp.diags(target / np.maximum(s, 1e-12)) @ C, dtype=np.float64); L.data = np.log1p(L.data)
    return L


def thin(C, ratio, seed=0):
    """Binomial thinning of integer counts to a fraction `ratio` of their depth."""
    rng = np.random.default_rng(seed); T = sp.csr_matrix(C, copy=True)
    T.data = rng.binomial(T.data.astype(np.int64), ratio).astype(np.float64); T.eliminate_zeros()
    return T


def depth_neutral_delta(X_prev, X_donor, lab_prev, lab_donor, seed=0, min_cells=20):
    """Per-cell-type change donor - previous in log1p(CP10K), with the previous stage first thinned to the
    donor stage's median depth when it is deeper. Returns ({celltype: delta vector}, depth_ratio)."""
    C_prev, _ = recover_counts(X_prev); C_donor, _ = recover_counts(X_donor)
    ratio = float(np.median(np.asarray(C_donor.sum(1)).ravel()) / np.median(np.asarray(C_prev.sum(1)).ravel()))
    if ratio < 1.0: C_prev = thin(C_prev, ratio, seed)
    Lp, Ld = to_log(C_prev), to_log(C_donor)
    lab_prev, lab_donor = np.asarray(lab_prev).astype(str), np.asarray(lab_donor).astype(str)
    out = {}
    for ct in np.intersect1d(np.unique(lab_prev), np.unique(lab_donor)):
        ip, idn = np.where(lab_prev == ct)[0], np.where(lab_donor == ct)[0]
        if len(ip) >= min_cells and len(idn) >= min_cells:
            out[ct] = np.asarray(Ld[idn].mean(0)).ravel() - np.asarray(Lp[ip].mean(0)).ravel()
    return out, ratio
