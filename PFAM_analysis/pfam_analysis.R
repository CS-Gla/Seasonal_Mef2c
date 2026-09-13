library(dplyr)
library(readr)
library(stringr)
library(purrr)
library(tidyr)
library(ggplot2)
library(rstudioapi)

setwd(dirname(getActiveDocumentContext()$path))

  read_domtblout_all <- function(file, species) {
    lines <- readr::read_lines(file)
    lines <- lines[!grepl("^#", lines)]
    if (length(lines) == 0) return(NULL)
    
    df <- readr::read_table2(
      paste(lines, collapse="\n"),
      col_names = c(
        "target_name","accession","tlen",
        "query_name","qacc","qlen",
        "full_evalue","full_score","full_bias",
        "dom_evalue","dom_score","dom_bias",
        "hmm_from","hmm_to","ali_from","ali_to",
        "env_from","env_to","acc","description"
      ),
      col_types = readr::cols(.default="c")
    )
    
    df %>%
      mutate(
        pfam = target_name,
        gene_id = tolower(query_name),
        species = species,
        full_evalue = as.numeric(full_evalue),
        full_score  = as.numeric(full_score)
      ) %>%
      select(species, gene_id, pfam, full_evalue, full_score)
  }
  
  extract_best_hits <- function(df) {
    df %>%
      group_by(species, gene_id) %>%
      arrange(full_evalue, desc(full_score)) %>%   # best first
      slice(1) %>%                                  # keep only top hit
      ungroup()
  }
  
  pfam_files <- list(
    psun = "pfam/psun.domtblout",
    ovis = "pfam/ovis.domtblout",
    athal = "pfam/athal.domtblout",
    cjap = "pfam/cjap.domtblout",
    taest = "pfam/taest.domtblout"
  )
  
  pfam_all <- purrr::imap_dfr(pfam_files, read_domtblout_all)

  
  pfam_best <- pfam_all %>%
    select(species, gene_id, pfam, full_evalue, full_score) %>%
    group_by(species, gene_id) %>%
    arrange(full_evalue, desc(full_score), .by_group = TRUE) %>%
    slice_head(n = 1) %>%
    ungroup()
  
  
  #No longer works. Why?
  #pfam_best <- extract_best_hits(pfam_all)

  pfam_species_presence <- pfam_best %>%
    group_by(pfam) %>%
    summarise(
      species_present = list(unique(species)),
      n_species_present = n_distinct(species),
      .groups = "drop"
    )
  
  all_species <- unique(pfam_best$species)
  n_species_total <- length(all_species)
  
  
  pfam_universal <- pfam_species_presence %>%
    filter(n_species_present == n_species_total)
  
  pfams_in_all_species <- pfam_universal$pfam
  
  # define read function
  read_expression <- function(path, species) {
    expr <- read.table(
      path,
      header = TRUE,
      sep = ",",
      check.names = FALSE,
      comment.char = "",
      quote = "\""
    )
    
    names(expr)[1] <- "gene_id"
    
    expr %>%
      mutate(species = species)
  }
  
  
  
  #list species expression files
  expr_files <- list(
    psun = "expression_pfam/psun.csv",
    ovis = "expression_pfam/ovis.csv",
    athal = "expression_pfam/athal.csv",
    cjap = "expression_pfam/cjap.csv",
    taest = "expression_pfam/taest.csv"
  )
  
  #execute read function
  expr_all <- imap(expr_files, read_expression) %>%
    bind_rows()
  
 sig_files <- list(
   psun = "sig/psun.csv",
   ovis = "sig/ovis.csv",
   athal = "sig/athal.csv",
   cjap = "sig/cjap.csv",
   taest = "sig/taest.csv"
  )
  
  sig_all <- imap(sig_files, read_expression) %>%
   bind_rows()
  
# expr_all <- expr_all %>%
 #   filter(gene_id %in% sig_all$gene_id)
  
  
  #list cols
  cols <- setdiff(colnames(expr_all), c("gene_id", "species"))
  
  conditions <- tibble(
    sample = cols,
    condition = case_when(
      grepl("^LD", cols, ignore.case=TRUE) ~ "LD",
      grepl("^SD", cols, ignore.case=TRUE)   ~ "SD",
      TRUE ~ "other"
    ),
    replicate = as.numeric(gsub(".*_(\\d+)$", "\\1", cols))
  )
  
  expr_all$gene_id <- tolower(expr_all$gene_id)
  
  
  #load phylostatiography
  #age <- read.table("gene_ages.tsv", sep="\t", header=TRUE, check.names=FALSE)
  
 # age <- na.omit(age)
  
  #age$gene_id <- tolower(age$gene_id)
  
  expr_pfam <- expr_all %>%
    left_join(pfam_best, by=c("gene_id", "species")) %>%
    filter(!is.na(pfam))
  
  
  ld_cols <- grep("^LD_", colnames(expr_pfam), value = TRUE)
  sd_cols <- grep("^SD_", colnames(expr_pfam), value = TRUE)
  all_cols <- c(ld_cols, sd_cols)
  

  
  expr_pfam_long <- expr_pfam %>%
    pivot_longer(
      cols = all_of(all_cols),
      names_to = "sample",
      values_to = "expression"
    ) %>%
    mutate(
      condition = ifelse(grepl("^LD_", sample), "LD",
                         ifelse(grepl("^SD_", sample), "SD", NA))
    )
  
  
  expr_pfam_z_long <- expr_pfam_long %>%
    group_by(gene_id, pfam, species) %>%
    mutate(
      z = (expression - mean(expression, na.rm = TRUE)) /
        sd(expression, na.rm = TRUE)
    ) %>%
    ungroup()

  expr_pfam_z_long <- na.omit(expr_pfam_z_long)
  
  
  expr_pfam_z_cond <- expr_pfam_z_long %>%
    group_by(gene_id, pfam, species, condition) %>%
    summarise(
      z_mean = mean(z)
    )
    
    expr_pfam_z_diff <- expr_pfam_z_cond %>%
      pivot_wider(names_from = condition, values_from = z_mean) %>%
      mutate(
        z_diff_abs = abs(LD - SD)
      )
  
    pfam_max_zdiff <- expr_pfam_z_diff %>%
      group_by(pfam, species) %>%
      summarise(max_z_diff = max(z_diff_abs, na.rm = TRUE), .groups = "drop")
    
    #pfam_expression[pfam_expression == "NULL"] <- NA
    #pfam_expression <- na.omit(pfam_expression)
    #pfam_expression$age <- sapply(pfam_expression$age, function(x) if(length(x)==0) NA else x)
    #pfam_expression$age <- as.numeric(pfam_expression$age)
    #pfam_expression <- na.omit(pfam_expression)
    
    
    
   # TEI_species <- pfam_expression %>%
   #   group_by(species, condition) %>%
   #   summarise(
   #      TEI = sum(expression * age) / sum(expression),
   #     .groups="drop"
   #    )
    
   #  TEI_species


    expr_pfam_z <- expr_pfam %>%
      pivot_longer(
        cols = all_of(cols),
        names_to = "sample",
        values_to = "expression"
      ) %>%
      group_by(pfam, species) %>%
      mutate(
        z = (expression - mean(expression, na.rm=TRUE)) /
          sd(expression, na.rm=TRUE)
      ) %>%
      ungroup()
    
    expr_pfam_z <- expr_pfam_z %>% left_join(conditions, by = "sample")
    
    #expr_pfam_z <- expr_pfam_z %>%
    #  left_join(age, by=c("gene_id")) %>%
     # filter(!is.na(pfam))
    
  
    expr_pfam_z <- expr_pfam_z[, -10]
    
    expr_pfam_z <- na.omit(expr_pfam_z)
    
    pfam_expression_z <- expr_pfam_z %>%
      group_by(pfam,species, condition) %>%
      summarise(
        z = sum(z),
        .groups = "drop"
      )
    
    pfam_expression_z <- pfam_expression_z %>%
      filter(pfam %in% pfams_in_all_species)
    
    
    
    plot_pfam <- function(pfam_id, data = pfam_expression_z) {
      df <- data %>% filter(pfam == pfam_id)
      
      ggplot(df, aes(x = condition, y = z, fill = species)) +
        geom_col(position = position_dodge(width = 0.7)) +
        theme_bw(base_size = 14) +
        labs(
          title = paste("Expression of", pfam_id),
          y = "Summed expression",
          x = "Condition"
        )
    }
    
    plot_pfam("SRF-TF")
    
    
    
    
    
    
    pfam_delta <- pfam_expression_z %>%
      pivot_wider(names_from = condition, values_from = z) %>%
      mutate(
        logFC = abs(LD - SD)
      )
    
    
    #filter for results
    pfam_delta <- pfam_delta %>%
      filter(abs(LD) >= 1,
             abs(SD) >= 1,
             logFC > 2)
    
    collapsed_delta <- pfam_delta %>%
      group_by(pfam) %>%
      summarise(
        species = list(sort(unique(species))),
        .groups = "drop"
      )
    
    target_species <- sort(c("ovis", "athal", "psun"))
    
    filtered_delta <- collapsed_delta %>%
      filter(purrr::map_lgl(species, ~ identical(.x, target_species)))
    
    
    ggplot(pfam_delta, aes(x = age, y = logFC, color = species)) +
      geom_point(alpha = 0.6) +
      theme_bw(base_size = 14) +
      labs(
        title = "pfam response vs evolutionary age",
        x = "Node age (phylo)",
        y = "z-score diff LD V SD"
      )

    
    library(ggVennDiagram)
    
    pfam_sets <- pfam_delta %>%
      distinct(pfam, species) %>%
      group_by(species) %>%
      summarise(pfams = list(pfam), .groups = "drop") %>%
      tibble::deframe()
    
    ggVennDiagram(pfam_sets, label_alpha = 0)
    
    
    
    top_pfam <- pfam_delta %>%
      group_by(species) %>%
      slice_max(abs(logFC), n = 10)
    
    ggplot(top_pfam, aes(x = reorder(pfam, logFC), y = logFC, fill = species)) +
      geom_col() +
      facet_wrap(~ species, scales = "free") +
      coord_flip() +
      theme_bw(base_size = 14) +
      labs(
        title = "Top 10 most-responsive orthogroups",
        x = "Orthogroup",
        y = "log2(LD/SD)"
      )

    
    
    
    
expression_srf <- pfam_expression_z[pfam_expression_z$pfam == 'SRF-TF',]
        
    
        
