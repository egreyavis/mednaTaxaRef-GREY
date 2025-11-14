library(rentrez)
library(Biostrings)
library(data.table)

################# SPECIES LIST PROCESSING
# Helper function to add taxid to taxonomy
add_taxids_to_taxonomy <- function(species_df) {
  
  setDT(species_df)
  
  # Function to get taxid for a rank
  get_rank_taxid <- function(i, rank_name) {
    path <- species_df$classification_path[i]
    path_ranks <- species_df$classification_path_ranks[i]
    path_ids <- species_df$classification_path_ids[i]
    
    names <- strsplit(path, "\\|")[[1]]
    ranks <- strsplit(path_ranks, "\\|")[[1]]
    ids <- strsplit(path_ids, "\\|")[[1]]
    
    rank_idx <- which(ranks == rank_name)
    
    if (length(rank_idx) > 0) {
      return(ids[rank_idx[1]])
    }
    return(NA_character_)
  }
  
  # Create new data.table with correct number of rows
  result <- data.table(
    superkingdom = rep(NA_character_, nrow(species_df)),
    kingdom = rep(NA_character_, nrow(species_df)),
    phylum = rep(NA_character_, nrow(species_df)),
    class = rep(NA_character_, nrow(species_df)),
    order = rep(NA_character_, nrow(species_df)),
    family = rep(NA_character_, nrow(species_df)),
    genus = rep(NA_character_, nrow(species_df)),
    species = rep(NA_character_, nrow(species_df))
  )
  
  for (i in 1:nrow(species_df)) {
    
    # For each taxonomic rank
    for (rank in c("superkingdom", "kingdom", "phylum", "class", "order", "family", "genus")) {
      name_val <- species_df[[rank]][i]
      
      if (!is.na(name_val) && name_val != "na" && name_val != "") {
        taxid <- get_rank_taxid(i, rank)
        if (!is.na(taxid)) {
          result[i, (rank) := paste0(name_val, "_", taxid)]
        }
      }
    }
    
    # Species uses taxon_id directly
    species_name <- species_df$species[i]
    if (!is.na(species_name) && species_name != "na" && species_name != "") {
      result[i, species := paste0(species_name, "_", species_df$taxon_id[i])]
    }
  }
  
  return(result)
}

# Backfill taxonomy ranks that are NA
backfill_missing_ranks <- function(taxonomy_df) {
  
  setDT(taxonomy_df)
  
  rank_hierarchy <- c("superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species")
  
  for (i in 1:nrow(taxonomy_df)) {
    for (j in 1:(length(rank_hierarchy) - 1)) {
      current_rank <- rank_hierarchy[j]
      next_rank <- rank_hierarchy[j + 1]
      
      # If current rank is NA but next rank exists
      current_val <- taxonomy_df[i, get(current_rank)]
      next_val <- taxonomy_df[i, get(next_rank)]
      
      if ((is.na(current_val) || current_val == "") && !is.na(next_val) && next_val != "") {
        taxonomy_df[i, (current_rank) := paste0(next_rank, "_", next_val)]
      }
    }
  }
  
  return(taxonomy_df)
}

########### ACCESSION FINDING AND FETCHING FUNCTIONS #######################
# Function to search Genbank for accessions by species and locus
search_species_accessions <- function(species_row, locus_searchterm, delay = 0.35) {
  
  # Extract taxid from species field (format: "Species name_taxid")
  species_info <- species_row$species
  
  # Split to get taxid
  taxid <- strsplit(species_info, "_")[[1]]
  taxid <- taxid[length(taxid)]  # Get last element in case species name has underscores
  
  if (is.na(taxid) || taxid == "" || taxid == "na") {
    warning(paste("No valid taxid found for", species_info))
    return(list(
      n_mitogenome = 0,
      ids_mitogenome = NA_character_,
      n_target = 0,
      ids_target = NA_character_
    ))
  }
  
  # Format search terms using taxid
  search_term <- paste0("txid", taxid, "[Organism]")
  mito_term <- paste(search_term, "AND mitochondrion[TITL] AND complete genome[TITL]")
  target_term <- paste(search_term, locus_searchterm)
  
  result <- list(
    n_mitogenome = 0,
    ids_mitogenome = NA_character_,
    n_target = 0,
    ids_target = NA_character_
  )
  
  # Search for mitogenomes
  mitogenomes <- tryCatch({
    entrez_search(db = "nucleotide", term = mito_term, retmax = 9999)
  }, error = function(e) {
    warning(paste("Mitogenome search failed for taxid", taxid, ":", e$message))
    return(NULL)
  })
  
  if (!is.null(mitogenomes)) {
    result$n_mitogenome <- mitogenomes$count
    if (mitogenomes$count > 0) {
      result$ids_mitogenome <- paste(mitogenomes$ids, collapse = "|")
    }
  }
  
  Sys.sleep(delay)
  
  # Search for target locus accessions
  targets <- tryCatch({
    entrez_search(db = "nucleotide", term = target_term, retmax = 9999)
  }, error = function(e) {
    warning(paste("Target search failed for taxid", taxid, ":", e$message))
    return(NULL)
  })
  
  if (!is.null(targets)) {
    result$n_target <- targets$count
    if (targets$count > 0) {
      result$ids_target <- paste(targets$ids, collapse = "|")
    }
  }
  
  return(result)
}

# Function to fetch species accessions with failure tracking
fetch_species_accessions <- function(species_row, locus, output_folder, max_seqs = 100, max_retries = 2) {
  
  species_name <- species_row$species
  n_target <- species_row$n_target
  ids_target <- species_row$ids_target
  
  # Skip if no targets
  if (is.na(n_target) || n_target == 0 || is.na(ids_target)) {
    return(list(data = NULL, status = "no_targets"))
  }
  
  # Parse IDs
  ids <- unlist(strsplit(ids_target, split = "\\|"))
  
  # Subsample if > max_seqs
  if (length(ids) > max_seqs) {
    ids <- sample(ids, max_seqs)
  }
  
  # Fetch sequences with retry logic
  seqs_target <- NULL
  for (attempt in 1:max_retries) {
    seqs_target <- tryCatch({
      entrez_fetch(db = "nuccore", id = ids, rettype = "fasta")
    }, error = function(e) {
      if (attempt < max_retries) {
        message(paste("Attempt", attempt, "failed for", species_name, "- retrying..."))
        Sys.sleep(2)  # Wait before retry
      }
      return(NULL)
    })
    
    if (!is.null(seqs_target)) break
  }
  
  if (is.null(seqs_target)) {
    return(list(data = NULL, status = "fetch_failed", species = species_name))
  }
  
  # Write and read back sequences
  fasta_file <- file.path(output_folder, paste(species_name, paste0(locus, ".fasta")))
  
  tryCatch({
    write(seqs_target, fasta_file)
  }, error = function(e) {
    return(list(data = NULL, status = "write_failed", species = species_name))
  })
  
  fasta_target <- tryCatch({
    readDNAStringSet(fasta_file, format = "fasta")
  }, error = function(e) {
    return(list(data = NULL, status = "read_failed", species = species_name))
  })
  
  if (is.null(fasta_target) || length(fasta_target) == 0) {
    return(list(data = NULL, status = "empty_fasta", species = species_name))
  }
  
  # Fetch accession numbers with retry
  seqs_target_accessions <- NULL
  for (attempt in 1:max_retries) {
    seqs_target_accessions <- tryCatch({
      entrez_fetch(db = "nuccore", id = ids, rettype = "acc")
    }, error = function(e) {
      if (attempt < max_retries) {
        message(paste("Attempt", attempt, "failed fetching accessions for", species_name, "- retrying..."))
        Sys.sleep(2)
      }
      return(NULL)
    })
    
    if (!is.null(seqs_target_accessions)) break
  }
  
  if (is.null(seqs_target_accessions)) {
    return(list(data = NULL, status = "accession_fetch_failed", species = species_name))
  }
  
  # Parse results
  seq_header <- names(fasta_target)
  sequence <- as.character(fasta_target)
  seq_accession <- unlist(strsplit(seqs_target_accessions, split = "\n"))
  
  # Check that lengths match
  n_seqs <- length(seq_header)
  n_accs <- length(seq_accession)
  
  if (n_seqs != n_accs) {
    warning(paste(
      "Mismatch for", species_name, "- Sequences:", n_seqs, "Accessions:", n_accs
    ))
    
    min_len <- min(n_seqs, n_accs)
    
    if (min_len == 0) {
      return(list(data = NULL, status = "length_mismatch", species = species_name))
    }
    
    seq_header <- seq_header[1:min_len]
    sequence <- sequence[1:min_len]
    seq_accession <- seq_accession[1:min_len]
  }
  
  # Create data frame
  result_df <- data.frame(
    seq_header = seq_header,
    sequence = sequence,
    seq_accession = seq_accession,
    type = rep("accession", length(seq_header)),
    species = rep(species_name, length(seq_header)),
    stringsAsFactors = FALSE
  )
  
  return(list(data = result_df, status = "success", species = species_name))
}

########### ADD MISSING CLADE MITOGENOMES AND ACCESSIONS (optional) ###########
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
      warning(paste("No species found for clade", clade_id))
      return(NULL)
    }
    
    cat(paste(" [Found", length(search_result$ids), "species]"))
    return(search_result$ids)
    
  }, error = function(e) {
    warning(paste("Failed to get species for clade", clade_id, ":", e$message))
    return(NULL)
  })
}

# Improved function to extract taxonomy with better XML parsing
extract_taxonomy <- function(species_id) {
  
  # Fetch full classification using efetch
  classification <- tryCatch({
    entrez_fetch(db = "taxonomy", id = species_id, rettype = "xml")
  }, error = function(e) {
    warning(paste("Failed to fetch taxonomy for", species_id))
    return(NULL)
  })
  
  if (is.null(classification)) {
    return(list(family = NA, family_id = NA, genus = NA, genus_id = NA, 
                species = NA, species_id = species_id))
  }
  
  # Initialize results
  family <- NA
  family_id <- NA
  genus <- NA
  genus_id <- NA
  species <- NA
  
  # Split into lines for parsing
  lines <- strsplit(classification, "\n")[[1]]
  
  # Find each Lineage entry with Rank and extract name and TaxId
  i <- 1
  while (i <= length(lines)) {
    
    # Look for family rank
    if (grepl("<Rank>family</Rank>", lines[i])) {
      # Search backwards and forwards for ScientificName and TaxId within same LineageEx block
      start_idx <- max(1, i - 15)
      end_idx <- min(length(lines), i + 5)
      
      for (j in start_idx:end_idx) {
        if (grepl("<ScientificName>", lines[j]) && is.na(family)) {
          family <- sub(".*<ScientificName>([^<]+)</ScientificName>.*", "\\1", lines[j])
          family <- trimws(family)
        }
        if (grepl("<TaxId>", lines[j]) && is.na(family_id) && j < i) {
          family_id <- sub(".*<TaxId>([^<]+)</TaxId>.*", "\\1", lines[j])
          family_id <- trimws(family_id)
        }
      }
    }
    
    # Look for genus rank
    if (grepl("<Rank>genus</Rank>", lines[i])) {
      start_idx <- max(1, i - 15)
      end_idx <- min(length(lines), i + 5)
      
      for (j in start_idx:end_idx) {
        if (grepl("<ScientificName>", lines[j]) && is.na(genus)) {
          genus <- sub(".*<ScientificName>([^<]+)</ScientificName>.*", "\\1", lines[j])
          genus <- trimws(genus)
        }
        if (grepl("<TaxId>", lines[j]) && is.na(genus_id) && j < i) {
          genus_id <- sub(".*<TaxId>([^<]+)</TaxId>.*", "\\1", lines[j])
          genus_id <- trimws(genus_id)
        }
      }
    }
    
    # Look for species rank - the main ScientificName at the top
    if (grepl("<Rank>species</Rank>", lines[i]) && is.na(species)) {
      start_idx <- max(1, i - 15)
      
      for (j in start_idx:i) {
        if (grepl("<ScientificName>", lines[j])) {
          species <- sub(".*<ScientificName>([^<]+)</ScientificName>.*", "\\1", lines[j])
          species <- trimws(species)
          break
        }
      }
    }
    
    i <- i + 1
  }
  
  # If species is still NA, try to get it from the main taxonomy record
  if (is.na(species)) {
    for (line in lines) {
      if (grepl("^\\s*<ScientificName>", line)) {
        species <- sub(".*<ScientificName>([^<]+)</ScientificName>.*", "\\1", line)
        species <- trimws(species)
        break
      }
    }
  }
  
  return(list(
    family = family,
    family_id = family_id,
    genus = genus,
    genus_id = genus_id,
    species = species,
    species_id = species_id
  ))
}

# Function to search for sequences (mitogenome first, then target)
search_species_sequences <- function(species_id, locus_searchterm) {
  
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
      term = paste(tax_query, locus_searchterm),
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

# Maine process missing clades function
process_missing_clades <- function(clades_missing, 
                                   locus_searchterm,
                                   rank_column = "order",
                                   max_species_per_clade = 3) {
  
  setDT(clades_missing)
  
  cat("Processing", nrow(clades_missing), "missing clades at rank:", rank_column, "\n")
  
  all_results <- list()
  
  for (i in 1:nrow(clades_missing)) {
    cat("\nProcessing clade", i, "of", nrow(clades_missing), ":")
    
    clade_row <- clades_missing[i, ]
    
    # Extract clade ID
    clade_info <- clade_row[[rank_column]]
    clade_parts <- strsplit(clade_info, "_")[[1]]
    clade_id <- clade_parts[length(clade_parts)]
    clade_name <- paste(clade_parts[-length(clade_parts)], collapse = "_")
    
    cat(" ", clade_name, "(", clade_id, ")")
    
    if (is.na(clade_id) || clade_id == "") {
      cat(" - No clade ID\n")
      next
    }
    
    # Get species list
    species_ids <- get_clade_species(clade_id, max_species = 1000)
    
    if (is.null(species_ids) || length(species_ids) == 0) {
      cat(" - No species found\n")
      next
    }
    
    Sys.sleep(0.35)
    
    # Randomize and search
    species_ids <- sample(species_ids)
    finds <- 0
    max_to_search <- min(length(species_ids), max_species_per_clade * 10)
    
    cat("\n  Searching up to", max_to_search, "species...")
    
    for (j in 1:max_to_search) {
      if (finds >= max_species_per_clade) break
      
      result <- search_species_sequences(species_ids[j], locus_searchterm)
      
      if (result$found) {
        cat("\n  Found sequences for species", species_ids[j])
        
        # Get taxonomy
        tax_info <- extract_taxonomy(species_ids[j])
        
        Sys.sleep(0.35)
        
        # Create result row
        result_row <- data.table(
          superkingdom = clade_row$superkingdom,
          kingdom = clade_row$kingdom,
          phylum = clade_row$phylum,
          class = clade_row$class,
          order = clade_row$order,
          family = ifelse(!is.na(tax_info$family) && !is.na(tax_info$family_id),
                          paste0(tax_info$family, "_", tax_info$family_id),
                          NA_character_),
          genus = ifelse(!is.na(tax_info$genus) && !is.na(tax_info$genus_id),
                         paste0(tax_info$genus, "_", tax_info$genus_id),
                         NA_character_),
          species = ifelse(!is.na(tax_info$species),
                           paste0(tax_info$species, "_", species_ids[j]),
                           NA_character_),
          species_id = as.character(species_ids[j]),
          ids_mitogenome = ifelse(!is.na(result$mito_id), as.character(result$mito_id), NA_character_),
          ids_target = ifelse(!is.na(result$target_id), as.character(result$target_id), NA_character_),
          n_mitogenome = "clade representative",
          n_target = "clade representative",
          tax_query = paste0("txid", species_ids[j])
        )
        
        all_results[[length(all_results) + 1]] <- result_row
        finds <- finds + 1
      }
    }
    
    cat("\n  Total found for", clade_name, ":", finds, "\n")
  }
  
  cat("\n\nComplete!\n")
  
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
  
  # Check if it matches any synonym (ONE AT A TIME)
  for (syn in target_synonyms) {
    if (grepl(syn, gene_value, ignore.case = TRUE)) {
      return(gene_value)
    }
  }
  
  return(NA)
}

# Function to process species mitogenomes
process_species_mitogenomes <- function(species_row, target_synonyms, is_gene = TRUE, max_mitos = 20) {
  
  species_name <- species_row$species
  n_mitogenomes <- as.numeric(species_row$n_mitogenome)
  ids_mitogenome <- species_row$ids_mitogenome
  
  cat("\rProcessing:", species_name, "with", n_mitogenomes, "mitogenomes")
  
  # Skip if no mitogenomes
  if (is.na(n_mitogenomes) || n_mitogenomes == 0 || is.na(ids_mitogenome)) {
    return(NULL)
  }
  
  # Parse mitogenome IDs
  mito_ids <- unlist(strsplit(ids_mitogenome, split = "\\|"))
  
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
  
  if (is.null(mito_accessions) || length(mito_accessions) == 0) {
    return(NULL)
  }
  
  # Process each accession with rate limiting
  results_list <- vector("list", length(mito_accessions))
  
  for (i in seq_along(mito_accessions)) {
    Sys.sleep(0.35)  # Rate limiting
    
    extracted <- extract_mito_feature(mito_accessions[i], target_synonyms, is_gene)
    
    # Only create data frame if sequence was found
    if (!is.na(extracted$sequence)) {
      results_list[[i]] <- data.frame(
        seq_header = extracted$header,
        sequence = extracted$sequence,
        seq_accession = extracted$accession,
        type = extracted$status,
        species = species_name,
        stringsAsFactors = FALSE
      )
    }
  }
  
  # Remove NULL entries
  results_list <- results_list[!sapply(results_list, is.null)]
  
  if (length(results_list) == 0) {
    return(NULL)
  }
  
  # Combine results
  do.call(rbind, results_list)
}


