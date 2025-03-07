library(msigdbr)
library(dplyr)

# Define pathways
H <- c("HALLMARK_APOPTOSIS", "HALLMARK_COAGULATION", "HALLMARK_COMPLEMENT", 
              "HALLMARK_DNA_REPAIR", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
              "HALLMARK_FATTY_ACID_METABOLISM", "HALLMARK_HYPOXIA", "HALLMARK_INFLAMMATORY_RESPONSE", 
              "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE",
              "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_P53_PATHWAY", 
              "HALLMARK_WNT_BETA_CATENIN_SIGNALING", "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY")

K <- c("KEGG_APOPTOSIS", "KEGG_CELL_ADHESION_MOLECULES_CAMS", "KEGG_CHEMOKINE_SIGNALING_PATHWAY",
          "KEGG_ECM_RECEPTOR_INTERACTION", "KEGG_INSULIN_SIGNALING_PATHWAY", 
          "KEGG_JAK_STAT_SIGNALING_PATHWAY", "KEGG_MAPK_SIGNALING_PATHWAY", "KEGG_MELANOGENESIS", 
          "KEGG_MELANOMA", "KEGG_MTOR_SIGNALING_PATHWAY", "KEGG_NOTCH_SIGNALING_PATHWAY",
          "KEGG_OXIDATIVE_PHOSPHORYLATION", "KEGG_SPHINGOLIPID_METABOLISM")

R <- c("REACTOME_ADAPTIVE_IMMUNE_SYSTEM", "REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM", 
              "REACTOME_INNATE_IMMUNE_SYSTEM", "REACTOME_LECTIN_PATHWAY_OF_COMPLEMENT_ACTIVATION",
              "REACTOME_TERMINAL_PATHWAY_OF_COMPLEMENT", "REACTOME_COMPLEMENT_CASCADE", 
              "REACTOME_ALTERNATIVE_COMPLEMENT_ACTIVATION")

# Fetch gene sets
pathway_geneset <- bind_rows(
  msigdbr("Homo sapiens", "H") %>% filter(gs_name %in% H),
  msigdbr("Homo sapiens", "C2", "CP:KEGG") %>% filter(gs_name %in% K),
  msigdbr("Homo sapiens", "C2", "CP:REACTOME") %>% filter(gs_name %in% R)
) %>%
  select(gs_name, gene_symbol) %>%
  distinct() %>%
  mutate(
    gs_name = recode(gs_name,
      "HALLMARK_APOPTOSIS" = "H_Apoptosis", "HALLMARK_COAGULATION" = "H_Coagulation",
      "HALLMARK_COMPLEMENT" = "H_Complement_cascade", "HALLMARK_DNA_REPAIR" = "H_DNA_repair",
      "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION" = "H_EMT", "HALLMARK_FATTY_ACID_METABOLISM" = "H_FA_metabolism",
      "HALLMARK_HYPOXIA" = "Hypoxia", "HALLMARK_INFLAMMATORY_RESPONSE" = "H_Inflammatory_response",
      "HALLMARK_INTERFERON_ALPHA_RESPONSE" = "H_INFA_response", "HALLMARK_INTERFERON_GAMMA_RESPONSE" = "H_INFG_response",
      "HALLMARK_OXIDATIVE_PHOSPHORYLATION" = "H_Oxidative_phosphorylation", "HALLMARK_P53_PATHWAY" = "H_p53",
      "HALLMARK_WNT_BETA_CATENIN_SIGNALING" = "H_WNTB_catenin_signaling", "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY" = "H_ROS",
      "KEGG_APOPTOSIS" = "K_Apoptosis", "KEGG_CELL_ADHESION_MOLECULES_CAMS" = "K_CAMs",
      "KEGG_CHEMOKINE_SIGNALING_PATHWAY" = "K_Chemokine_signaling", "KEGG_ECM_RECEPTOR_INTERACTION" = "K_ECM_receptor",
      "KEGG_INSULIN_SIGNALING_PATHWAY" = "K_Insulin_signaling", "KEGG_JAK_STAT_SIGNALING_PATHWAY" = "K_JAK_STAT_signaling",
      "KEGG_MAPK_SIGNALING_PATHWAY" = "K_MAPK_signaling", "KEGG_MELANOGENESIS" = "K_Melanogenesis",
      "KEGG_MELANOMA" = "K_Melanoma", "KEGG_MTOR_SIGNALING_PATHWAY" = "K_mTOR_signaling",
      "KEGG_NOTCH_SIGNALING_PATHWAY" = "K_Notch_signaling", "KEGG_OXIDATIVE_PHOSPHORYLATION" = "K_Oxidative_phosphorylation",
      "KEGG_SPHINGOLIPID_METABOLISM" = "K_Sphingolipid_metabolism",
      "REACTOME_ADAPTIVE_IMMUNE_SYSTEM" = "R_Adaptive_immune_response",
      "REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM" = "R_Cytokine_signaling",
      "REACTOME_INNATE_IMMUNE_SYSTEM" = "R_Innate_immune_response",
      "REACTOME_LECTIN_PATHWAY_OF_COMPLEMENT_ACTIVATION" = "R_Lectin_pathway",
      "REACTOME_TERMINAL_PATHWAY_OF_COMPLEMENT" = "R_Terminal_pathway",
      "REACTOME_COMPLEMENT_CASCADE" = "R_Complement_cascade",
      "REACTOME_ALTERNATIVE_COMPLEMENT_ACTIVATION" = "R_Alternative_pathway"
    ),
    category = case_when(
      gs_name %in% c("H_FA_metabolism", "K_FA_metabolism", "H_Oxidative_phosphorylation", "K_Oxidative_phosphorylation", "K_Sphingolipid_metabolism") ~ "Metabolic_Pathways",
      gs_name %in% c("H_Complement_cascade", "H_Inflammatory_response", "H_INFA_response", "H_INFG_response", "R_Adaptive_immune_response", "R_Cytokine_signaling", "R_Innate_immune_response", "R_Lectin_pathway", "R_Alternative_pathway") ~ "Immune_Regulation_Pathways",
      gs_name %in% c("H_WNTB_catenin_signaling", "K_Chemokine_signaling", "K_Insulin_signaling", "K_JAK_STAT_signaling", "K_MAPK_signaling", "K_mTOR_signaling", "K_Notch_signaling") ~ "Signaling_Cascade",
      gs_name %in% c("H_Apoptosis", "H_DNA_repair", "H_p53") ~ "Cell_Cycle_and_Apoptosis",
      gs_name %in% c("H_Coagulation", "H_EMT", "Hypoxia", "K_CAMs", "K_ECM_receptor", "K_Melanogenesis", "K_Melanoma") ~ "Cell_differentiation"
    ))

write.csv(pathway_geneset, "data/input/pathway_geneset.csv", row.names = FALSE)
