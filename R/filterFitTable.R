filterFitTable <- function(table = NULL,
                           metadata = NULL,
                           filter_underpowered = NULL,
                           filter_homozygous = NULL,
                           af_cutoff = NULL,
                           ranks = 10) {
  fitTable <- dplyr::left_join(table, metadata, by = "SAMPLE_ID") %>%
    dplyr::select(-file) %>%
    dplyr::relocate(PATIENT_ID, .after = SAMPLE_ID) %>%
    dplyr::relocate(TP53freq, smooth, .after = expected_TP53_AF) %>%
    dplyr::relocate(precPloidy, precPurity, .after = purity)

  ## Apply hard filters
  ##  filter under powered fits when config variable is TRUE
  if (filter_underpowered) {
    fitTable <- fitTable %>%
      dplyr::filter(powered == 1)
  }
  ## filter high prop homozygous loss when config variable is TRUE
  if (filter_homozygous) {
    fitTable <- fitTable %>%
      dplyr::filter(homozygousLoss <= homozygous_prop)
  }

  # standard filtering
  filtered_results <- fitTable %>%
    dplyr::group_by(SAMPLE_ID, ploidy) %>%
    dplyr::mutate(rank_clonality = dplyr::min_rank(clonality)) %>% #rank clonality within a unique ploidy state
    dplyr::filter(rank_clonality == 1) %>% #select ploidy with the lowest clonality within a unique ploidy state
    dplyr::group_by(SAMPLE_ID) %>%
    dplyr::top_n(-ranks, wt = clonality) %>% # select top 10 ploidy states with the lowest clonality values
    dplyr::mutate(rank_clonality = dplyr::min_rank(clonality)) %>% # rank by clonality within a sample across ploidy in top 10
    # retain samples without TP53 mutations and where expected and observed TP53freq <=0.15
    dplyr::filter(is.na(TP53freq) |
                    dplyr::near(expected_TP53_AF, TP53freq, tol = af_cutoff)) %>%
    dplyr::arrange(PATIENT_ID, SAMPLE_ID)

  pruned_results <- filtered_results %>%
    dplyr::arrange(SAMPLE_ID, ploidy) %>%
    dplyr::group_by(SAMPLE_ID) %>%
    dplyr::mutate(pl_diff = abs(ploidy - dplyr::lag(ploidy))) %>% #, pu_diff = abs(purity - dplyr::lag(purity)) not used
    dplyr::mutate(new_state_n = dplyr::row_number() == 1 |
                    pl_diff > 0.3) %>%
    dplyr::mutate(new_state = cumsum(new_state_n)) %>%
    dplyr::group_by(SAMPLE_ID, new_state) %>%
    dplyr::filter(rank_clonality == min(rank_clonality)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(use = rep(NA, times = nrow(.)),
                  notes = rep(NA, times = nrow(.))) %>%
    dplyr::arrange(ploidy)

  return(list(filtered = filtered_results, pruned = pruned_results))
}
