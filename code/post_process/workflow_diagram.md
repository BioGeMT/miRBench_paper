```mermaid
graph LR
  %% Define files
  M[AGO2_eCLIP_Manakov22_full_dataset.tsv]
  
  %% Process 0: pp_0_filter_and_deduplicate
  pp_0[pp_0_filter_and_deduplicate]
  M_filtered[AGO2_eCLIP_Manakov22.filtered.deduplicated.tsv]
  
  %% Process 1: pp_1_exclude_mirna_families
  pp_1[pp_1_exclude_mirna_families]
  M_excluded[AGO2_eCLIP_Manakov22.filtered.deduplicated.excluded.tsv]
  M_remaining[AGO2_eCLIP_Manakov22.filtered.deduplicated.remaining.tsv]

  %% Process 1a: pp_1a_add_seeds_and_filter_interactions
  pp_1a[pp_1a_add_seeds_and_filter_interactions]
  M_excluded_canonical[AGO2_eCLIP_Manakov22.filtered.deduplicated.excluded.canonical.tsv]
  M_excluded_noncanonical[AGO2_eCLIP_Manakov22.filtered.deduplicated.excluded.noncanonical.tsv]
  M_excluded_nonseed[AGO2_eCLIP_Manakov22.filtered.deduplicated.excluded.nonseed.tsv]
  M_remaining_canonical[AGO2_eCLIP_Manakov22.filtered.deduplicated.remaining.canonical.tsv]
  M_remaining_noncanonical[AGO2_eCLIP_Manakov22.filtered.deduplicated.remaining.noncanonical.tsv]
  M_remaining_nonseed[AGO2_eCLIP_Manakov22.filtered.deduplicated.remaining.nonseed.tsv]
  
  %% Process 2: pp_2_make_negatives
  pp_2[pp_2_make_negatives]
  M_excluded_canonical_neg[AGO2_eCLIP_Manakov22.filtered.deduplicated.excluded.canonical.gene_clusters_added.mirfam_sorted.negatives.tsv]
  M_excluded_noncanonical_neg[AGO2_eCLIP_Manakov22.filtered.deduplicated.excluded.noncanonical.gene_clusters_added.mirfam_sorted.negatives.tsv]
  M_excluded_nonseed_neg[AGO2_eCLIP_Manakov22.filtered.deduplicated.excluded.nonseed.gene_clusters_added.mirfam_sorted.negatives.tsv]
  M_remaining_canonical_neg[AGO2_eCLIP_Manakov22.filtered.deduplicated.remaining.canonical.gene_clusters_added.mirfam_sorted.negatives.tsv]
  M_remaining_noncanonical_neg[AGO2_eCLIP_Manakov22.filtered.deduplicated.remaining.noncanonical.gene_clusters_added.mirfam_sorted.negatives.tsv]
  M_remaining_nonseed_neg[AGO2_eCLIP_Manakov22.filtered.deduplicated.remaining.nonseed.gene_clusters_added.mirfam_sorted.negatives.tsv]

    %% Process 3: pp_3_train_test_splits
  pp_3[pp_3_train_test_splits]
  M_remaining_canonical_neg_train[AGO2_eCLIP_Manakov22.filtered.deduplicated.remaining.canonical.gene_clusters_added.mirfam_sorted.negatives.train.tsv]
  M_remaining_canonical_neg_test[AGO2_eCLIP_Manakov22.filtered.deduplicated.remaining.canonical.gene_clusters_added.mirfam_sorted.negatives.test.tsv]
  M_remaining_noncanonical_neg_train[AGO2_eCLIP_Manakov22.filtered.deduplicated.remaining.noncanonical.gene_clusters_added.mirfam_sorted.negatives.train.tsv]
  M_remaining_noncanonical_neg_test[AGO2_eCLIP_Manakov22.filtered.deduplicated.remaining.noncanonical.gene_clusters_added.mirfam_sorted.negatives.test.tsv]
  M_remaining_nonseed_neg_train[AGO2_eCLIP_Manakov22.filtered.deduplicated.remaining.nonseed.gene_clusters_added.mirfam_sorted.negatives.train.tsv]
  M_remaining_nonseed_neg_test[AGO2_eCLIP_Manakov22.filtered.deduplicated.remaining.nonseed.gene_clusters_added.mirfam_sorted.negatives.test.tsv]

  %% Process 4: pp_4_drop_test_col
  pp_4[pp_4_drop_test_col]
  M_remaining_canonical_neg_train_drop[AGO2_eCLIP_Manakov22.filtered.deduplicated.remaining.canonical.gene_clusters_added.mirfam_sorted.negatives.train.drop_test_col.tsv]
  M_remaining_canonical_neg_test_drop[AGO2_eCLIP_Manakov22.filtered.deduplicated.remaining.canonical.gene_clusters_added.mirfam_sorted.negatives.test.drop_test_col.tsv]
  M_remaining_noncanonical_neg_train_drop[AGO2_eCLIP_Manakov22.filtered.deduplicated.remaining.noncanonical.gene_clusters_added.mirfam_sorted.negatives.train.drop_test_col.tsv]
  M_remaining_noncanonical_neg_test_drop[AGO2_eCLIP_Manakov22.filtered.deduplicated.remaining.noncanonical.gene_clusters_added.mirfam_sorted.negatives.test.drop_test_col.tsv]
  M_remaining_nonseed_neg_train_drop[AGO2_eCLIP_Manakov22.filtered.deduplicated.remaining.nonseed.gene_clusters_added.mirfam_sorted.negatives.train.drop_test_col.tsv]
  M_remaining_nonseed_neg_test_drop[AGO2_eCLIP_Manakov22.filtered.deduplicated.remaining.nonseed.gene_clusters_added.mirfam_sorted.negatives.test.drop_test_col.tsv]
  %% Process 4: pp_4_drop_test_col
  pp_4[pp_4_drop_test_col]
  M_excluded_canonical_neg_drop[AGO2_eCLIP_Manakov22.filtered.deduplicated.excluded.canonical.gene_clusters_added.mirfam_sorted.negatives.drop_test_col.tsv]
  M_excluded_noncanonical_neg_drop[AGO2_eCLIP_Manakov22.filtered.deduplicated.excluded.noncanonical.gene_clusters_added.mirfam_sorted.negatives.drop_test_col.tsv]
  M_excluded_nonseed_neg_drop[AGO2_eCLIP_Manakov22.filtered.deduplicated.excluded.nonseed.gene_clusters_added.mirfam_sorted.negatives.drop_test_col.tsv]

  %% Process 5: pp_5_add_conservation
  pp_5[pp_5_add_conservation]
  M_remaining_canonical_neg_train_conserv[AGO2_eCLIP_Manakov22.filtered.deduplicated.remaining.canonical.gene_clusters_added.mirfam_sorted.negatives.train.drop_test_col.conservation.tsv]
  M_remaining_canonical_neg_test_conserv[AGO2_eCLIP_Manakov22.filtered.deduplicated.remaining.canonical.gene_clusters_added.mirfam_sorted.negatives.test.drop_test_col.conservation.tsv]
  M_remaining_noncanonical_neg_train_conserv[AGO2_eCLIP_Manakov22.filtered.deduplicated.remaining.noncanonical.gene_clusters_added.mirfam_sorted.negatives.train.drop_test_col.conservation.tsv]
  M_remaining_noncanonical_neg_test_conserv[AGO2_eCLIP_Manakov22.filtered.deduplicated.remaining.noncanonical.gene_clusters_added.mirfam_sorted.negatives.test.drop_test_col.conservation.tsv]
  M_remaining_nonseed_neg_train_conserv[AGO2_eCLIP_Manakov22.filtered.deduplicated.remaining.nonseed.gene_clusters_added.mirfam_sorted.negatives.train.drop_test_col.conservation.tsv]
  M_remaining_nonseed_neg_test_conserv[AGO2_eCLIP_Manakov22.filtered.deduplicated.remaining.nonseed.gene_clusters_added.mirfam_sorted.negatives.test.drop_test_col.conservation.tsv]
  M_excluded_canonical_neg_conserv[AGO2_eCLIP_Manakov22.filtered.deduplicated.excluded.canonical.gene_clusters_added.mirfam_sorted.negatives.drop_test_col.conservation.tsv]
  M_excluded_noncanonical_neg_conserv[AGO2_eCLIP_Manakov22.filtered.deduplicated.excluded.noncanonical.gene_clusters_added.mirfam_sorted.negatives.drop_test_col.conservation.tsv]
  M_excluded_nonseed_neg_conserv[AGO2_eCLIP_Manakov22.filtered.deduplicated.excluded.nonseed.gene_clusters_added.mirfam_sorted.negatives.drop_test_col.conservation.tsv]

  %% Connections
  M --> pp_0
  pp_0 --> M_filtered
  M_filtered --> pp_1
  pp_1 --> M_excluded
  pp_1 --> M_remaining
  M_excluded --> pp_1a
  M_remaining --> pp_1a
  pp_1a --> M_excluded_canonical
  pp_1a --> M_excluded_noncanonical
  pp_1a --> M_excluded_nonseed
  pp_1a --> M_remaining_canonical
  pp_1a --> M_remaining_noncanonical
  pp_1a --> M_remaining_nonseed
  M_excluded_canonical --> pp_2
  M_excluded_noncanonical --> pp_2
  M_excluded_nonseed --> pp_2
  M_remaining_canonical --> pp_2
  M_remaining_noncanonical --> pp_2
  M_remaining_nonseed --> pp_2
  pp_2 --> M_excluded_canonical_neg
  pp_2 --> M_excluded_noncanonical_neg
  pp_2 --> M_excluded_nonseed_neg
  pp_2 --> M_remaining_canonical_neg
  pp_2 --> M_remaining_noncanonical_neg
  pp_2 --> M_remaining_nonseed_neg
  M_remaining_canonical_neg --> pp_3
  M_remaining_noncanonical_neg --> pp_3
  M_remaining_nonseed_neg --> pp_3
  pp_3 --> M_remaining_canonical_neg_train
  pp_3 --> M_remaining_canonical_neg_test
  pp_3 --> M_remaining_noncanonical_neg_train
  pp_3 --> M_remaining_noncanonical_neg_test
  pp_3 --> M_remaining_nonseed_neg_train
  pp_3 --> M_remaining_nonseed_neg_test
  M_remaining_canonical_neg_train --> pp_4
  M_remaining_canonical_neg_test --> pp_4
  M_remaining_noncanonical_neg_train --> pp_4
  M_remaining_noncanonical_neg_test --> pp_4
  M_remaining_nonseed_neg_train --> pp_4
  M_remaining_nonseed_neg_test --> pp_4
  M_excluded_canonical_neg --> pp_4
  M_excluded_noncanonical_neg --> pp_4
  M_excluded_nonseed_neg --> pp_4
  pp_4 --> M_remaining_canonical_neg_train_drop
  pp_4 --> M_remaining_canonical_neg_test_drop
  pp_4 --> M_remaining_noncanonical_neg_train_drop
  pp_4 --> M_remaining_noncanonical_neg_test_drop
  pp_4 --> M_remaining_nonseed_neg_train_drop
  pp_4 --> M_remaining_nonseed_neg_test_drop
  pp_4 --> M_excluded_canonical_neg_drop
  pp_4 --> M_excluded_noncanonical_neg_drop
  pp_4 --> M_excluded_nonseed_neg_drop
  M_remaining_canonical_neg_train_drop --> pp_5
  M_remaining_canonical_neg_test_drop --> pp_5
  M_remaining_noncanonical_neg_train_drop --> pp_5
  M_remaining_noncanonical_neg_test_drop --> pp_5
  M_remaining_nonseed_neg_train_drop --> pp_5
  M_remaining_nonseed_neg_test_drop --> pp_5
  M_excluded_canonical_neg_drop --> pp_5
  M_excluded_noncanonical_neg_drop --> pp_5
  M_excluded_nonseed_neg_drop --> pp_5
  pp_5 --> M_remaining_canonical_neg_train_conserv
  pp_5 --> M_remaining_canonical_neg_test_conserv
  pp_5 --> M_remaining_noncanonical_neg_train_conserv
  pp_5 --> M_remaining_noncanonical_neg_test_conserv
  pp_5 --> M_remaining_nonseed_neg_train_conserv
  pp_5 --> M_remaining_nonseed_neg_test_conserv
  pp_5 --> M_excluded_canonical_neg_conserv
  pp_5 --> M_excluded_noncanonical_neg_conserv
  pp_5 --> M_excluded_nonseed_neg_conserv

  %% Define file and process styles
  classDef file fill:#E0F7FA,stroke:#00796B,stroke-width:2px,color:black;
  classDef process fill:#FFE0B2,stroke:#E65100,stroke-width:2px,color:black;

  %% Apply styles
  class M,M_filtered,M_excluded,M_remaining,M_excluded_canonical,M_excluded_noncanonical,M_excluded_nonseed,M_remaining_canonical,M_remaining_noncanonical,M_remaining_nonseed,M_excluded_canonical_neg,M_excluded_noncanonical_neg,M_excluded_nonseed_neg,M_remaining_canonical_neg,M_remaining_noncanonical_neg,M_remaining_nonseed_neg,M_remaining_canonical_neg_train,M_remaining_canonical_neg_test,M_remaining_noncanonical_neg_train,M_remaining_noncanonical_neg_test,M_remaining_nonseed_neg_train,M_remaining_nonseed_neg_test,M_excluded_canonical_neg_drop,M_excluded_noncanonical_neg_drop,M_excluded_nonseed_neg_drop,M_remaining_canonical_neg_train_drop,M_remaining_canonical_neg_test_drop,M_remaining_noncanonical_neg_train_drop,M_remaining_noncanonical_neg_test_drop,M_remaining_nonseed_neg_train_drop,M_remaining_nonseed_neg_test_drop,M_excluded_canonical_neg_conserv,M_excluded_noncanonical_neg_conserv,M_excluded_nonseed_neg_conserv,M_remaining_canonical_neg_train_conserv,M_remaining_canonical_neg_test_conserv,M_remaining_noncanonical_neg_train_conserv,M_remaining_noncanonical_neg_test_conserv,M_remaining_nonseed_neg_train_conserv,M_remaining_nonseed_neg_test_conserv file;
  class pp_0,pp_1,pp_1a,pp_2,pp_3,pp_4,pp_5 process;

