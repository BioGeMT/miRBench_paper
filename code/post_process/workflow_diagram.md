```mermaid
graph LR
  %% Define files
  M[AGO2_eCLIP_Manakov22_full_dataset.tsv]
  H[AGO2_CLASH_Hejret23_full_dataset.tsv]
  K[AGO2_eCLIP_Klimentova22_full_dataset.tsv]
  
  %% Process 0: pp_0_filter_and_deduplicate
  pp_0[pp_0_filter_and_deduplicate]
  M_filtered[AGO2_eCLIP_Manakov22.filtered.deduplicated.tsv]
  H_filtered[AGO2_CLASH_Hejret23_full_dataset.filtered.deduplicated.tsv]
  K_filtered[AGO2_eCLIP_Klimentova22_full_dataset.filtered.deduplicated.tsv]
  
  %% Process 1: pp_1_exclude_mirna_families
  pp_1[pp_1_exclude_mirna_families]
  M_excluded[AGO2_eCLIP_Manakov22.filtered.deduplicated.excluded.tsv]
  M_remaining[AGO2_eCLIP_Manakov22.filtered.deduplicated.remaining.tsv]
  
  %% Process 2: pp_2_make_negatives
  pp_2[pp_2_make_negatives]
  M_neg[AGO2_eCLIP_Manakov22.filtered.deduplicated.remaining.gene_clusters_added.mirfam_sorted.negatives.tsv]
  H_neg[AGO2_CLASH_Hejret23_full_dataset.filtered.deduplicated.gene_clusters_added.mirfam_sorted.negatives.tsv]
  K_neg[AGO2_eCLIP_Klimentova22_full_dataset.filtered.deduplicated.gene_clusters_added.mirfam_sorted.negatives.tsv]
  M_excluded_neg[AGO2_eCLIP_Manakov22.filtered.deduplicated.excluded.gene_clusters_added.mirfam_sorted.negatives.tsv]
  
  %% Process 3: pp_3_train_test_splits
  pp_3[pp_3_train_test_splits]
  M_neg_train[AGO2_eCLIP_Manakov22.filtered.deduplicated.remaining.gene_clusters_added.mirfam_sorted.negatives.train.tsv]
  M_neg_test[AGO2_eCLIP_Manakov22.filtered.deduplicated.remaining.gene_clusters_added.mirfam_sorted.negatives.test.tsv]
  H_neg_train[AGO2_CLASH_Hejret23_full_dataset.filtered.deduplicated.gene_clusters_added.mirfam_sorted.negatives.train.tsv]
  H_neg_test[AGO2_CLASH_Hejret23_full_dataset.filtered.deduplicated.gene_clusters_added.mirfam_sorted.negatives.test.tsv]

    %% Process 4: pp_4_drop_test_col
  pp_4[pp_4_drop_test_col]
  M_train_drop[AGO2_eCLIP_Manakov22.filtered.deduplicated.remaining.gene_clusters_added.mirfam_sorted.negatives.train.drop_test_col.tsv]
  M_test_drop[AGO2_eCLIP_Manakov22.filtered.deduplicated.remaining.gene_clusters_added.mirfam_sorted.negatives.test.drop_test_col.tsv]
  H_train_drop[AGO2_CLASH_Hejret23_full_dataset.filtered.deduplicated.gene_clusters_added.mirfam_sorted.negatives.train.drop_test_col.tsv]
  H_test_drop[AGO2_CLASH_Hejret23_full_dataset.filtered.deduplicated.gene_clusters_added.mirfam_sorted.negatives.test.drop_test_col.tsv]
  K_neg_drop[AGO2_eCLIP_Klimentova22_full_dataset.filtered.deduplicated.gene_clusters_added.mirfam_sorted.negatives.drop_test_col.tsv]
  M_excluded_drop[AGO2_eCLIP_Manakov22.filtered.deduplicated.excluded.gene_clusters_added.mirfam_sorted.negatives.drop_test_col.tsv]

  %% Process 5: pp_5_add_conservation
  pp_5[pp_5_add_conservation]
  M_train_conserv[AGO2_eCLIP_Manakov22.filtered.deduplicated.remaining.gene_clusters_added.mirfam_sorted.negatives.train.drop_test_col.conservation.tsv]
  M_test_conserv[AGO2_eCLIP_Manakov22.filtered.deduplicated.remaining.gene_clusters_added.mirfam_sorted.negatives.test.drop_test_col.conservation.tsv]
  H_train_conserv[AGO2_CLASH_Hejret23_full_dataset.filtered.deduplicated.gene_clusters_added.mirfam_sorted.negatives.train.drop_test_col.conservation.tsv]
  H_test_conserv[AGO2_CLASH_Hejret23_full_dataset.filtered.deduplicated.gene_clusters_added.mirfam_sorted.negatives.test.drop_test_col.conservation.tsv]
  K_neg_conserv[AGO2_eCLIP_Klimentova22_full_dataset.filtered.deduplicated.gene_clusters_added.mirfam_sorted.negatives.drop_test_col.conservation.tsv]
  M_excluded_conserv[AGO2_eCLIP_Manakov22.filtered.deduplicated.excluded.gene_clusters_added.mirfam_sorted.negatives.drop_test_col.conservation.tsv]

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
  pp_4 --> M_train_drop
  pp_4 --> M_test_drop
  pp_4 --> H_train_drop
  pp_4 --> H_test_drop
  pp_4 --> K_neg_drop
  pp_4 --> M_excluded_drop
  M_train_drop --> pp_5
  M_test_drop --> pp_5
  H_train_drop --> pp_5
  H_test_drop --> pp_5
  K_neg_drop --> pp_5
  M_excluded_drop --> pp_5
  pp_5 --> M_train_conserv
  pp_5 --> M_test_conserv
  pp_5 --> H_train_conserv
  pp_5 --> H_test_conserv
  pp_5 --> K_neg_conserv
  pp_5 --> M_excluded_conserv

  %% Define file and process styles
  classDef file fill:#E0F7FA,stroke:#00796B,stroke-width:2px,color:black;
  classDef process fill:#FFE0B2,stroke:#E65100,stroke-width:2px,color:black;

  %% Apply styles
  class M,H,K,M_filtered,H_filtered,K_filtered,M_excluded,M_remaining,M_neg,H_neg,K_neg,M_excluded_neg,M_neg_train,M_neg_test,H_neg_train,H_neg_test,M_train_drop,M_test_drop,H_train_drop,H_test_drop,K_neg_drop,M_excluded_drop,M_train_conserv,M_test_conserv,H_train_conserv,H_test_conserv,K_neg_conserv,M_excluded_conserv file;
  class pp_0,pp_1,pp_2,pp_3,pp_4,pp_5 process;