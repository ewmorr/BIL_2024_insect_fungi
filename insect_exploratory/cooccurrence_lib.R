##############################################################################
# Shared engine for the presence/absence co-occurrence screens
# (target_species_ophiostomatales_cooccurrence.R,
#  scolytinae_ophiostomatales_cooccurrence.R).
#
# jaccard_grid_test() computes, for every beetle x fungal-taxon pair in a
# grid, the observed Jaccard similarity J = a / (a + b + c) plus a restricted
# permutation null and a Fisher's-exact cross-check. It does NOT compare
# beetle-vs-beetle or fungus-vs-fungus pairs.
#
#   beetle_pa   n x B  0/1 matrix (samples x beetle taxa)
#   fungal_pa   n x K  0/1 matrix (samples x fungal taxa), SAME row order
#   perm_mat    n_perm x n index matrix from permute::shuffleSet() -- each row
#               a within-trap permutation of the sample order; only the beetle
#               vector is permuted, the fungal vector is held observed
#   grid_name   label copied into the $grid column
#
# Returns one row per (beetle, fungal_unit): the 2x2 counts (a/b/c/d),
# jaccard, the null mean/sd and z_score = (J - mean_null)/sd_null, one-sided
# p_pos (co-occurrence) and p_neg (avoidance), two-sided p_perm =
# 2*min(p_pos,p_neg), Fisher p, and BH q-values (q_perm, q_fisher) computed
# WITHIN the grid.
##############################################################################

jaccard_grid_test <- function(beetle_pa, fungal_pa, perm_mat, grid_name) {
  stopifnot(nrow(beetle_pa) == nrow(fungal_pa),
            ncol(perm_mat) == nrow(beetle_pa))
  n      <- nrow(fungal_pa)
  n_perm <- nrow(perm_mat)
  fn     <- colSums(fungal_pa)                     # fungus presence counts (K)
  K      <- ncol(fungal_pa)
  rows   <- list()

  for (b in colnames(beetle_pa)) {
    bvec <- beetle_pa[, b]
    bn   <- sum(bvec)

    a_obs <- as.vector(crossprod(bvec, fungal_pa))            # K
    j_obs <- a_obs / (bn + fn - a_obs)

    # permuted beetle vector -> shared-count for every pair at once (n_perm x K)
    bperm <- bvec[perm_mat]; dim(bperm) <- dim(perm_mat)
    a_p   <- bperm %*% fungal_pa
    j_p   <- a_p / sweep(-a_p, 2, bn + fn, "+")

    p_pos <- (1 + colSums(sweep(j_p, 2, j_obs, ">="))) / (n_perm + 1)
    p_neg <- (1 + colSums(sweep(j_p, 2, j_obs, "<="))) / (n_perm + 1)
    exp_j <- colMeans(j_p)
    sd_j  <- apply(j_p, 2, sd)

    fish_p <- vapply(seq_len(K), function(k) {
      a <- a_obs[k]; bb <- bn - a; cc <- fn[k] - a; dd <- n - a - bb - cc
      fisher.test(matrix(c(a, bb, cc, dd), 2))$p.value
    }, numeric(1))

    rows[[b]] <- tibble::tibble(
      grid        = grid_name,
      beetle      = b,
      fungal_unit = colnames(fungal_pa),
      n_samples   = n,
      beetle_n    = bn,
      fungus_n    = fn,
      shared_n    = a_obs,
      beetle_only = bn - a_obs,
      fungus_only = fn - a_obs,
      both_absent = n - bn - fn + a_obs,
      jaccard     = j_obs,
      exp_jaccard_null = exp_j,
      sd_jaccard_null  = sd_j,
      z_score     = ifelse(sd_j > 0, (j_obs - exp_j) / sd_j, NA_real_),
      p_pos = p_pos, p_neg = p_neg,
      p_perm   = pmin(1, 2 * pmin(p_pos, p_neg)),
      p_fisher = fish_p,
      direction = ifelse(j_obs >= exp_j, "co-occurring", "avoiding")
    )
  }

  dplyr::bind_rows(rows) |>
    dplyr::mutate(q_perm   = p.adjust(p_perm,   method = "BH"),
                  q_fisher = p.adjust(p_fisher, method = "BH")) |>
    dplyr::relocate(q_perm,   .after = p_perm) |>
    dplyr::relocate(q_fisher, .after = p_fisher) |>
    dplyr::arrange(beetle, dplyr::desc(jaccard))
}
