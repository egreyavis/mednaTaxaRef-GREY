library(rentrez)
library(Biostrings)
library(data.table)

########### ACCESSION FINDING AND DOWNLOADING FUNCTIONS #######################
# Function to search Genbank for accessions by species and locus
search_species_accessions <- function(species_name, locus_searchterm, delay = 0.35) {
  
  search_name <- paste0(species_name, "[ORGN]")
  mito_term <- paste(search_name, "AND mitochondrion[TITL] AND complete genome[TITL]")
  target_term <- paste(search_name, locus_searchterm)
  
  result <- list(
    n_mitogenome = 0,
    ids_mitogenome = NA_character_,
    n_target = 0,
    ids_target = NA_character_
  )
  
  mitogenomes <- tryCatch({
    entrez_search(db = "nucleotide", term = mito_term, retmax = 9999)
  }, error = function(e) NULL)
  
  if (!is.null(mitogenomes)) {
    result$n_mitogenome <- mitogenomes$count
    if (mitogenomes$count > 0) {
      result$ids_mitogenome <- paste(mitogenomes$ids, collapse = "|")
    }
  }
  
  Sys.sleep(delay)
  
  targets <- tryCatch({
    entrez_search(db = "nucleotide", term = target_term, retmax = 9999)
  }, error = function(e) NULL)
  
  if (!is.null(targets)) {
    result$n_target <- targets$count
    if (targets$count > 0) {
      result$ids_target <- paste(targets$ids, collapse = "|")
    }
  }
  
  return(result)
}


########### PROCESS MISSING CLADES ####################
# Function to get species list for a clade using NCBI Taxonomy
get_clade_species <- function(clade_id, max_species = 1000) {
  
  tryCatch({
    # Search for all species under this taxon
    search_result <- entrez_search(
      db = "taxonomy",
      term = paste0("txid", clade_id, "[Subtree] AND species[Rank]"),
      retmax = max_species
    )
    
    if (length(search_result$ids) == 0) {
      return(NULL)
    }
    
    return(search_result$ids)
    
  }, error = function(e) {
    warning(paste("Failed to get species for clade", clade_id, ":", e$message))
    return(NULL)
  })
}

# Function to search for sequences (mitogenome first, then target)
search_species_sequences <- function(species_id, target_locus_searchterm) {
  
  tax_query <- paste0("txid", species_id, "[Organism]")
  
  result <- list(
    species_id = species_id,
    mito_id = NA,
    target_id = NA,
    found = FALSE
  )
  
  # Try mitogenome first
  mitogenomes <- tryCatch({
    entrez_search(
      db = "nucleotide",
      term = paste(tax_query, "AND mitochondrion[TITL] AND complete genome[TITL]"),
      retmax = 9999
    )
  }, error = function(e) NULL)
  
  Sys.sleep(0.35)
  
  if (!is.null(mitogenomes) && length(mitogenomes$ids) > 0) {
    result$mito_id <- sample(mitogenomes$ids, 1)
    result$found <- TRUE
    return(result)
  }
  
  # If no mitogenome, try target locus
  targets <- tryCatch({
    entrez_search(
      db = "nucleotide",
      term = paste(tax_query, target_locus_searchterm),
      retmax = 9999
    )
  }, error = function(e) NULL)
  
  Sys.sleep(0.35)
  
  if (!is.null(targets) && length(targets$ids) > 0) {
    result$target_id <- sample(targets$ids, 1)
    result$found <- TRUE
  }
  
  return(result)
}

# Function to get taxonomy for species IDs
get_taxonomy_info <- function(species_ids) {
  
  if (length(species_ids) == 0) {
    return(NULL)
  }
  
  # Fetch taxonomy summaries
  tax_summaries <- tryCatch({
    entrez_summary(db = "taxonomy", id = species_ids)
  }, error = function(e) {
    warning("Failed to fetch taxonomy summaries")
    return(NULL)
  })
  
  if (is.null(tax_summaries)) {
    return(NULL)
  }
  
  # Parse taxonomy for each species
  tax_list <- lapply(species_ids, function(id) {
    
    summary <- if (length(species_ids) == 1) {
      tax_summaries
    } else {
      tax_summaries[[as.character(id)]]
    }
    
    # Extract lineage information
    lineage <- summary$lineage
    sci_name <- summary$scientificname
    tax_id <- summary$taxid
    
    # Parse lineage to get family and genus
    lineage_parts <- strsplit(lineage, "; ")[[1]]
    
    list(
      species_id = id,
      scientific_name = sci_name,
      lineage = lineage,
      lineage_parts = lineage_parts
    )
  })
  
  return(tax_list)
}

# Function to extract taxonomic ranks from lineage
extract_taxonomy <- function(species_id) {
  
  # Fetch full classification
  classification <- tryCatch({
    entrez_fetch(db = "taxonomy", id = species_id, rettype = "xml")
  }, error = function(e) {
    return(NULL)
  })
  
  if (is.null(classification)) {
    return(list(family = NA, genus = NA, species = NA))
  }
  
  # Parse XML to extract ranks (simple text parsing)
  lines <- strsplit(classification, "\n")[[1]]
  
  family <- NA
  genus <- NA
  species <- NA
  
  # Look for rank tags
  for (line in lines) {
    if (grepl("<Rank>family</Rank>", line)) {
      # Get scientific name from nearby lines
      sci_line_idx <- which(lines == line)
      for (j in (sci_line_idx - 5):(sci_line_idx + 5)) {
        if (grepl("<ScientificName>", lines[j])) {
          family <- gsub(".*<ScientificName>([^<]+)</ScientificName>.*", "\\1", lines[j])
          family <- trimws(family)
          break
        }
      }
    }
    if (grepl("<Rank>genus</Rank>", line)) {
      sci_line_idx <- which(lines == line)
      for (j in (sci_line_idx - 5):(sci_line_idx + 5)) {
        if (grepl("<ScientificName>", lines[j])) {
          genus <- gsub(".*<ScientificName>([^<]+)</ScientificName>.*", "\\1", lines[j])
          genus <- trimws(genus)
          break
        }
      }
    }
    if (grepl("<Rank>species</Rank>", line)) {
      sci_line_idx <- which(lines == line)
      for (j in (sci_line_idx - 5):(sci_line_idx + 5)) {
        if (grepl("<ScientificName>", lines[j])) {
          species <- gsub(".*<ScientificName>([^<]+)</ScientificName>.*", "\\1", lines[j])
          species <- trimws(species)
          break
        }
      }
    }
  }
  
  return(list(family = family, genus = genus, species = species))
}

# Fetch accession sequences for species
scrape_species_accessions <- function(species_row, locus, output_folder, max_seqs = 100) {
  
  species_name <- species_row$search_name
  n_target <- species_row$n_target
  ids_target <- species_row$ids_target
  
  # Skip if no targets
  if (n_target == 0 || is.na(ids_target)) {
    return(NULL)
  }
  
  # Parse IDs
  ids <- unlist(strsplit(ids_target, split = "\\|"))
  
  # Subsample if > max_seqs
  if (length(ids) > max_seqs) {
    ids <- sample(ids, max_seqs)
  }
  
  # Fetch sequences from GenBank
  seqs_target <- tryCatch({
    entrez_fetch(db = "nuccore", id = ids, rettype = "fasta")
  }, error = function(e) {
    warning(paste("Failed to fetch sequences for", species_name, ":", e$message))
    return(NULL)
  })
  
  if (is.null(seqs_target)) {
    return(NULL)
  }
  
  # Write and read back sequences
  fasta_file <- file.path(output_folder, paste(species_name, paste0(locus, ".fasta")))
  write(seqs_target, fasta_file)
  
  fasta_target <- tryCatch({
    readDNAStringSet(fasta_file, format = "fasta")
  }, error = function(e) {
    warning(paste("Failed to read fasta for", species_name))
    return(NULL)
  })
  
  if (is.null(fasta_target) || length(fasta_target) == 0) {
    return(NULL)
  }
  
  # Fetch accession numbers
  seqs_target_accessions <- tryCatch({
    entrez_fetch(db = "nuccore", id = ids, rettype = "acc")
  }, error = function(e) {
    warning(paste("Failed to fetch accessions for", species_name))
    return(NULL)
  })
  
  if (is.null(seqs_target_accessions)) {
    return(NULL)
  }
  
  # Parse and format results
  seq_header <- names(fasta_target)
  sequence <- as.character(fasta_target)
  seq_accession <- unlist(strsplit(seqs_target_accessions, split = "\n"))
  
  # Create data frame
  result_df <- data.frame(
    seq_header = seq_header,
    sequence = sequence,
    seq_accession = seq_accession,
    type = "accession",
    species = species_name,
    stringsAsFactors = FALSE
  )
  
  return(result_df)
}

# Main processing function: Process all clades
process_missing_clades <- function(clades_missing, 
                                   target_locus_searchterm, 
                                   rank_column = clade_level,  # Can be "order", "family", "class", etc.
                                   max_species_per_clade = max_species_per_clade) {
  
  setDT(clades_missing)
  
  cat("Processing", nrow(clades_missing), "missing clades at rank:", rank_column, "\n")
  
  # Initialize results
  all_results <- list()
  
  for (i in 1:nrow(clades_missing)) {
    cat("\rProcessing clade", i, "of", nrow(clades_missing))
    
    # Get the current clade row (THIS IS THE KEY CHANGE)
    clade_row <- clades_missing[i, ]
    
    # Extract clade ID from the rank column (format: "Name_ID")
    clade_info <- clade_row[[rank_column]]  # Changed from get() to [[]]
    clade_id <- strsplit(clade_info, "_")[[1]][2]
    
    if (is.na(clade_id)) {
      warning(paste("No clade ID found for row", i))
      next
    }
    
    # Get species list for this clade
    species_ids <- get_clade_species(clade_id, max_species = 1000)
    
    if (is.null(species_ids) || length(species_ids) == 0) {
      warning(paste("No species found for clade", clade_id))
      next
    }
    
    Sys.sleep(0.35)
    
    # Randomize species order
    species_ids <- sample(species_ids)
    
    # Search for sequences until we find max_species_per_clade or run out
    finds <- 0
    max_to_search <- min(length(species_ids), max_species_per_clade * 10)  # Search up to 10x the target
    
    for (j in 1:max_to_search) {
      if (finds >= max_species_per_clade) break
      
      result <- search_species_sequences(species_ids[j], target_locus_searchterm)
      
      if (result$found) {
        # Get taxonomy info for this species
        tax_info <- extract_taxonomy(species_ids[j])
        
        Sys.sleep(0.35)
        
        # Create result row (USING clade_row INSTEAD OF clades_missing[i, column])
        result_row <- data.table(
          superkingdom = clade_row$superkingdom,
          kingdom = clade_row$kingdom,
          phylum = clade_row$phylum,
          class = clade_row$class,
          order = clade_row$order,
          family = ifelse(!is.na(tax_info$family), 
                          paste0(tax_info$family, "_", species_ids[j]), 
                          NA),
          genus = ifelse(!is.na(tax_info$genus), 
                         paste0(tax_info$genus, "_", species_ids[j]), 
                         NA),
          species = ifelse(!is.na(tax_info$species), 
                           paste0(tax_info$species, "_", species_ids[j]), 
                           NA),
          species_id = species_ids[j],
          ids_mitogenome = ifelse(!is.na(result$mito_id), result$mito_id, NA),
          ids_target = ifelse(!is.na(result$target_id), result$target_id, NA),
          n_mitogenome = "clade representative",
          n_target = "clade representative",
          tax_query = paste0("txid", species_ids[j])
        )
        
        all_results[[length(all_results) + 1]] <- result_row
        finds <- finds + 1
      }
    }
  }
  
  cat("\n")
  
  # Combine all results
  if (length(all_results) == 0) {
    warning("No sequences found for any clades")
    return(NULL)
  }
  
  clade_seqs <- rbindlist(all_results, fill = TRUE)
  
  return(clade_seqs)
}

############# MITOGENOME SCRAPING FUNCTIONS #######################
# Fetch GenBank record and extract sequence for a specific feature
extract_mito_feature <- function(accession, target_synonyms, is_gene = TRUE) {
  
  result <- list(
    header = paste("Unparsed mitochondrion", accession),
    sequence = NA,
    accession = accession,
    status = "scrape"
  )
  
  tryCatch({
    # Fetch GenBank record in gbwithparts format (includes features and sequence)
    gb_text <- entrez_fetch(db = "nuccore", 
                           id = accession,
                           rettype = "gbwithparts", 
                           retmode = "text")
    
    # Parse the record
    parsed <- parse_genbank_feature(gb_text, target_synonyms, is_gene)
    
    if (!is.null(parsed$sequence)) {
      result$header <- paste(parsed$feature_name, "mitochondrion", accession)
      result$sequence <- as.character(parsed$sequence)
    }
    
  }, error = function(e) {
    warning(paste("Error processing", accession, ":", e$message))
  })
  
  return(result)
}

# Parse GenBank text and extract feature sequence
parse_genbank_feature <- function(gb_text, target_synonyms, is_gene = TRUE) {
  
  lines <- strsplit(gb_text, "\n")[[1]]
  
  # Determine feature type to search for
  feature_types <- if (is_gene) {
    c("gene", "CDS")
  } else {
    c("rRNA", "tRNA", "misc_RNA")
  }
  
  # Find the target feature
  feature_info <- find_target_feature(lines, target_synonyms, feature_types)
  
  if (is.null(feature_info)) {
    return(list(sequence = NULL, feature_name = NULL))
  }
  
  # Extract the full sequence from ORIGIN section
  origin_start <- which(grepl("^ORIGIN", lines))
  if (length(origin_start) == 0) {
    return(list(sequence = NULL, feature_name = NULL))
  }
  
  # Get sequence lines (after ORIGIN, before //)
  seq_lines <- lines[(origin_start + 1):length(lines)]
  seq_lines <- seq_lines[!grepl("^//", seq_lines)]
  
  # Clean and concatenate sequence
  full_sequence <- paste(gsub("[^acgtACGT]", "", seq_lines), collapse = "")
  full_sequence <- DNAString(full_sequence)
  
  # Extract the target subsequence
  target_seq <- subseq(full_sequence, 
                       start = feature_info$start, 
                       end = feature_info$end)
  
  # Handle complement if needed
  if (feature_info$complement) {
    target_seq <- reverseComplement(target_seq)
  }
  
  return(list(
    sequence = target_seq,
    feature_name = feature_info$name
  ))
}

# Find target feature in GenBank lines
find_target_feature <- function(lines, target_synonyms, feature_types) {
  
  for (i in seq_along(lines)) {
    line <- lines[i]
    
    # Check if this is a feature of interest
    is_target_feature <- any(sapply(feature_types, function(ft) {
      grepl(paste0("^\\s{5}", ft, "\\s+"), line)
    }))
    
    if (!is_target_feature) next
    
    # Extract location from first line
    location_info <- extract_location(line)
    
    # Look ahead for gene/product name
    feature_block <- lines[i:min(i + 20, length(lines))]
    gene_name <- extract_gene_name(feature_block, target_synonyms)
    
    if (!is.na(gene_name)) {
      return(list(
        name = gene_name,
        start = location_info$start,
        end = location_info$end,
        complement = location_info$complement
      ))
    }
  }
  
  return(NULL)
}

# Extract location coordinates
extract_location <- function(line) {
  complement <- grepl("complement", line)
  
  # Extract coordinate range - handle various formats
  coords <- regmatches(line, regexpr("\\d+\\.\\.\\d+", line))
  
  if (length(coords) == 0) {
    return(list(start = NA, end = NA, complement = FALSE))
  }
  
  parts <- as.numeric(strsplit(coords, "\\.\\.")[[1]])
  
  list(
    start = parts[1],
    end = parts[2],
    complement = complement
  )
}

# Extract gene name from feature block
extract_gene_name <- function(feature_block, target_synonyms) {
  
  # Look for /gene= or /product=
  gene_line <- grep('/gene=|/product=', feature_block, value = TRUE)
  
  if (length(gene_line) == 0) return(NA)
  
  # Extract the quoted value
  gene_match <- regexpr('"[^"]+"', gene_line[1])
  if (gene_match[1] == -1) return(NA)
  
  gene_value <- regmatches(gene_line[1], gene_match)
  gene_value <- gsub('"', '', gene_value)
  
  # Check if it matches any synonym
  matches <- sapply(target_synonyms, function(syn) {
    grepl(syn, gene_value, ignore.case = TRUE)
  })
  
  if (any(matches)) {
    return(gene_value)
  }
  
  return(NA)
}

# Main processing function - base R approach
process_species_mitogenomes <- function(species_row, target_synonyms, is_gene = TRUE, max_mitos = 20) {
  
  species_name <- species_row$search_name
  n_mitogenomes <- species_row$n_mitogenome
  
  cat("\rProcessing:", species_name, "with", n_mitogenomes, "mitogenomes")
  
  # Skip if no mitogenomes
  if (n_mitogenomes == 0) {
    return(NULL)
  }
  
  # Parse mitogenome IDs
  mito_ids <- unlist(strsplit(species_row$ids_mitogenome, split = "\\|"))
  
  # Subsample if > max_mitos
  if (length(mito_ids) > max_mitos) {
    mito_ids <- sample(mito_ids, max_mitos)
  }
  
  # Fetch accessions
  mito_accessions <- tryCatch({
    acc <- entrez_fetch(mito_ids, db = "nuccore", rettype = "acc")
    unlist(strsplit(acc, split = "\n"))
  }, error = function(e) {
    warning(paste("Failed to fetch accessions for", species_name))
    return(NULL)
  })
  
  if (is.null(mito_accessions)) {
    return(NULL)
  }
  
  # Process each accession with rate limiting
  results_list <- vector("list", length(mito_accessions))
  
  for (i in seq_along(mito_accessions)) {
    Sys.sleep(0.35)  # Rate limiting
    
    extracted <- extract_mito_feature(mito_accessions[i], target_synonyms, is_gene)
    
    results_list[[i]] <- data.frame(
      header = extracted$header,
      sequence = extracted$sequence,
      accession = extracted$accession,
      source = extracted$status,
      species = species_name,
      stringsAsFactors = FALSE
    )
  }
  
  # Combine results
  do.call(rbind, results_list)
}


