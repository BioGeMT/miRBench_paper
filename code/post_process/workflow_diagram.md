```mermaid
graph TD
  %% Define files
  M[Manakov]
  H[Hejret]
  K[Klimentova]
  
  %% Process 0: pp_0_filter_and_deduplicate
  pp_0[pp_0_filter_and_deduplicate]
  M_filtered[Manakov.filtered.deduplicated.tsv]
  H_filtered[Hejret.filtered.deduplicated.tsv]
  K_filtered[Klimentova.filtered.deduplicated.tsv]
  
  %% Process 1: pp_1_exclude_mirna_families
  pp_1[pp_1_exclude_mirna_families]
  M_excluded[Manakov.filtered.deduplicated.excluded.tsv]
  M_remaining[Manakov.filtered.deduplicated.remaining.tsv]
  
  %% Process 2: pp_2_make_negatives
  pp_2[pp_2_make_negatives]
  M_neg[Manakov.filtered.deduplicated.remaining.gene_clusters_added.mirfam_sorted.negatives.tsv]
  H_neg[Hejret.filtered.deduplicated.gene_clusters_added.mirfam_sorted.negatives.tsv]
  K_neg[Klimentova.filtered.deduplicated.gene_clusters_added.mirfam_sorted.negatives.tsv]
  M_excluded_neg[Manakov.filtered.deduplicated.excluded.gene_clusters_added.mirfam_sorted.negatives.tsv]
  
  %% Process 3: pp_3_train_test_splits
  pp_3[pp_3_train_test_splits]
  M_neg_train[Manakov.filtered.deduplicated.remaining.gene_clusters_added.mirfam_sorted.negatives.train.tsv]
  M_neg_test[Manakov.filtered.deduplicated.remaining.gene_clusters_added.mirfam_sorted.negatives.test.tsv]
  H_neg_train[Hejret.filtered.deduplicated.gene_clusters_added.mirfam_sorted.negatives.train.tsv]
  H_neg_test[Hejret.filtered.deduplicated.gene_clusters_added.mirfam_sorted.negatives.test.tsv]
  
  %% Process 4: pp_4_add_conservation
  pp_4[pp_4_add_conservation]
  M_train_conserv[Manakov.filtered.deduplicated.remaining.gene_clusters_added.mirfam_sorted.negatives.train.conservation.tsv]
  M_test_conserv[Manakov.filtered.deduplicated.remaining.gene_clusters_added.mirfam_sorted.negatives.test.conservation.tsv]
  H_train_conserv[Hejret.filtered.deduplicated.gene_clusters_added.mirfam_sorted.negatives.train.conservation.tsv]
  H_test_conserv[Hejret.filtered.deduplicated.gene_clusters_added.mirfam_sorted.negatives.test.conservation.tsv]
  K_neg_conserv[Klimentova.filtered.deduplicated.gene_clusters_added.mirfam_sorted.negatives.conservation.tsv]
  M_excluded_conserv[Manakov.filtered.deduplicated.excluded.gene_clusters_added.mirfam_sorted.negatives.conservation.tsv]
  
  %% Connections
  M --> pp_0
  H --> pp_0
  K --> pp_0
  pp_0 --> M_filtered
  pp_0 --> H_filtered
  pp_0 --> K_filtered
  M_filtered --> pp_1
  pp_1 --> M_excluded
  pp_1 --> M_remaining
  M_remaining --> pp_2
  M_excluded --> pp_2
  H_filtered --> pp_2
  K_filtered --> pp_2
  pp_2 --> M_neg
  pp_2 --> H_neg
  pp_2 --> K_neg
  pp_2 --> M_excluded_neg
  M_neg --> pp_3
  H_neg --> pp_3
  pp_3 --> M_neg_train
  pp_3 --> M_neg_test
  pp_3 --> H_neg_train
  pp_3 --> H_neg_test
  M_neg_train --> pp_4
  M_neg_test --> pp_4
  H_neg_train --> pp_4
  H_neg_test --> pp_4
  K_neg --> pp_4
  M_excluded_neg --> pp_4
  pp_4 --> M_train_conserv
  pp_4 --> M_test_conserv
  pp_4 --> H_train_conserv
  pp_4 --> H_test_conserv
  pp_4 --> K_neg_conserv
  pp_4 --> M_excluded_conserv
