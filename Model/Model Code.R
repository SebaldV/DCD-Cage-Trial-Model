#### Model description ####

# This program simulates individual gene drive crosses.
set.seed(123)

#   ____________________________________________________________________________
#   Functions                                                              ####

##  ............................................................................
##  Parameter pre-sets                                                      ####
Default_Nuc_Rates <- function(Nuc_File_Name, project_root) {
  
  default_file_location <- if (file.exists("Parameter Files/Default_Nuc.csv")) {"Parameter Files/Default_Nuc.csv"} else {file.path(project_root,"Default_Nuc.csv")}
  
  ## ----------  READ FILES WITHOUT COERCING TYPES ------------------------------
  Default_nuc_rates_table <- read.csv(
    default_file_location,
    header           = FALSE,
    col.names        = c("Parameter", "Value"),
    colClasses       = c("character", "character"),   # keep as character
    stringsAsFactors = FALSE
  )
  
  if (file.exists(file.path("Parameter Files", Nuc_File_Name))) {
  nuc_rates_table <- read.csv(
    file.path("Parameter Files", Nuc_File_Name),
    header           = FALSE,
    col.names        = c("Parameter", "Value"),
    colClasses       = c("character", "character"),
    stringsAsFactors = FALSE
  )
  } else {
    nuc_rates_table <- Default_nuc_rates_table
    message(paste0("Using: ",default_file_location))
  }
  
  nuc_rates_table <- full_join(Default_nuc_rates_table,
                               nuc_rates_table,
                               by     = "Parameter",
                               suffix = c("_default", "")) %>% 
    mutate(
      Value = ifelse(is.na(Value) | Value == "", Value_default, Value)  # patch holes
    ) %>% 
    select(Parameter, Value) %>% 
    arrange(match(Parameter, Default_nuc_rates_table$Parameter))        # optional: keep default order
  
   Nuc_File <- as.list(as.vector(nuc_rates_table[[2]]))
   names(Nuc_File) <- as.vector(nuc_rates_table[[1]])
   
   numeric_variables <- c(
     "Germline_Cutrate_Male",
     "Germline_Cutrate_Female",
     "Somatic_Cutrate_Male",
     "Somatic_Cutrate_Female",
     "Germline_HDR_Male",
     "Germline_HDR_Female",
     "Deposition_Cutrate_From_Male_Parent",
     "Deposition_Cutrate_From_Female_Parent",
     "Deposition_HDR"
   )
   
   if (!all(numeric_variables %in% names(Nuc_File))) {
     stop("Nuc file is missing parameters: ",
          paste(numeric_variables[!numeric_variables %in% names(Nuc_File)],
                collapse = ", "))
   }
   
   Nuc_File[numeric_variables] <- as.numeric(unlist(Nuc_File[numeric_variables]))
   
   logical_variables <- c(
     "gRNA_Deposition",
     "Shadow_Drive_Cas9_Deposition_Only",
     "Homozygous_Wildtype_No_HDR"
   )
   
   
   if (!all(logical_variables %in% names(Nuc_File))) {
     stop("Nuc file is missing parameters: ",
          paste(logical_variables[!logical_variables %in% names(Nuc_File)],
                collapse = ", "))
   }
   
   Nuc_File[logical_variables] <- as.logical(unlist(Nuc_File[logical_variables]))
   

   
   ## ---------------- 2.  FLAG UNDECLARED COLUMNS -------------------------------
   expected_cols   <- c(numeric_variables, logical_variables)
   unexpected_cols <- setdiff(names(Nuc_File), expected_cols)
   
   if (length(unexpected_cols) > 0) {
     message(
       "The following columns in Nuc_File are NOT in numeric_variables ",
       "or logical_variables: ",
       paste(unexpected_cols, collapse = ", ")
     )
   }
   
  # Here we check if there is deposition
  # if there is, we store two variables instead of one
  # these two variables ensure that the number of genotypes is
  # doubled in alpha (with and without deposition)
  # if there is no deposition we only need create a single copy of each genotype
  if (Nuc_File$Deposition_Cutrate_From_Female_Parent > 0) {
    Nuc_File$F_Dep <-
      c(1, 0) # By default there is no deposition. c(1,0) Means there is
  } else if (Nuc_File$Deposition_Cutrate_From_Female_Parent == 0) {
    Nuc_File$F_Dep <-
      0 # By default there is no deposition. c(1,0) Means there is
  } else {
    stop("Female deposition cut rate not 0 or >0")
  }
  
  if (Nuc_File$Deposition_Cutrate_From_Male_Parent > 0) {
    Nuc_File$M_Dep <-
      c(1, 0) # By default there is no deposition. c(1,0) Means there is
  } else if (Nuc_File$Deposition_Cutrate_From_Male_Parent == 0) {
    Nuc_File$M_Dep <-
      0 # By default there is no deposition. c(1,0) Means there is
  } else {
    stop("Male deposition cut rate not 0 or >0")
  }
   
    
  return(
    Nuc_File
  )
}

# here we set up the genotype frequency and fitness
# parameters
Default_Geno <- function(Allele_Fit_File_Name, project_root) {
  
  default_file_location <- if (file.exists("Parameter Files/Default_Alleles.csv")) {
    "Parameter Files/Default_Alleles.csv"
  } else {
    file.path(project_root,"Default_Alleles.csv")
  }

  ## --- 1.  Read defaults ----
  defaults <- read_csv(
    default_file_location,
    col_names = c("Parameter", "Default", paste0("X", 3:(ncol(read_csv(default_file_location, show_col_types = FALSE, col_names=FALSE, n_max=1))-2))), # Adjust X names if needed
    col_types = cols(.default = col_character()),
    show_col_types = FALSE,
    skip = 0, # Assuming first row is headers if col_names is used effectively
    na = c("", "NA")
  )
  
  params_file_path <- file.path("Parameter Files", Allele_Fit_File_Name)
  params_loaded_from_file <- file.exists(params_file_path)

  if (params_loaded_from_file) {
    params <- read_csv(
      params_file_path,
      col_names = FALSE, # Will be X1, X2 ...
      na = c("", "NA"),
      col_types = cols(.default = col_character())
    )
    # If params file has headers, use read_csv(..., col_names = TRUE) or let it infer
    # For this code, assuming params first column is Parameter, others are Xn style
    # If params has actual headers like "Parameter", "A", "A", ... for Locus row:
    # params <- read_csv(params_file_path, na = c("", "NA"), col_types = cols(.default = col_character()))
    # names(params)[1] <- "Parameter" # Ensure first column is Parameter
  } else {
    params <- defaults
    message(paste0("Using: ",default_file_location))
  }
  
  # Ensure the first column, whether from loaded params or defaults, is named "Parameter"
  names(params)[1] <- "Parameter"
  
  ## --- 3.  Abort early if param file contains unknown rows ----
  unknown <- setdiff(params$Parameter, defaults$Parameter)
  if (length(unknown) > 0) {
    stop("Unknown parameter(s) in parameter file: ",
         paste0(Allele_Fit_File_Name, "/", unknown, collapse = ", "))
  }
  
  # --- 4. New Patching Logic ---
  
  # --- 4a. Prepare `defaults_long` for lookup ---
  defaults_locus_row_full <- defaults[defaults$Parameter == "Locus", ]
  defaults_alleles_row_full <- defaults[defaults$Parameter == "Alleles", ]
  
  # Get actual data column names from defaults (excluding Parameter and Default)
  default_allele_data_cols <- names(defaults)[!names(defaults) %in% c("Parameter", "Default")]
  
  default_col_map_list <- lapply(default_allele_data_cols, function(col_name) {
    # Ensure that the column actually exists in the Locus/Allele rows
    def_locus_val <- if (col_name %in% names(defaults_locus_row_full)) defaults_locus_row_full[[col_name]] else NA_character_
    def_allele_val <- if (col_name %in% names(defaults_alleles_row_full)) defaults_alleles_row_full[[col_name]] else NA_character_
    
    tibble(
      def_col_name_original = col_name,
      def_locus = def_locus_val,
      def_allele = def_allele_val
    )
  })
  default_col_map <- bind_rows(default_col_map_list) %>% filter(!is.na(def_locus)) # Keep only validly mapped cols
  
  defaults_for_pivot <- defaults %>%
    filter(!Parameter %in% c("Locus", "Alleles"))
  
  if (nrow(defaults_for_pivot) > 0 && nrow(default_col_map) > 0) {
    defaults_long <- defaults_for_pivot %>%
      pivot_longer(
        cols = all_of(intersect(default_allele_data_cols, default_col_map$def_col_name_original)), # only pivot existing & mapped cols
        names_to = "def_col_name_original",
        values_to = "Value_default_specific"
      ) %>%
      left_join(default_col_map, by = "def_col_name_original") %>%
      select(Parameter, def_locus, def_allele, Value_default_specific)
  } else {
    defaults_long <- tibble(
      Parameter = character(), def_locus = character(), def_allele = character(),
      Value_default_specific = character()
    )
  }
  
  generic_defaults <- defaults %>% select(Parameter, Default)
  
  # --- 4b. Define `find_default_value` function ---
  find_default_value <- function(p_name, t_locus, t_allele, def_long_data, gen_defs_data) {
    # Try direct match for Locus and Allele
    val <- def_long_data %>%
      filter(Parameter == p_name, def_locus == t_locus, def_allele == t_allele) %>%
      pull(Value_default_specific)
    
    if (length(val) == 1 && !is.na(val) && val != "") {
      return(val[1]) # val could be vector if multiple matches, though unlikely with good data
    }
    
    # If specific allele not found, and target allele is not T, r1, or r2, try 'D' allele
    if (!t_allele %in% c("T", "r1", "r2")) {
      val_d <- def_long_data %>%
        filter(Parameter == p_name, def_locus == t_locus, def_allele == "D") %>%
        pull(Value_default_specific)
      if (length(val_d) == 1 && !is.na(val_d) && val_d != "") {
        return(val_d[1])
      }
    }
    
    # Fallback to generic default for that parameter
    gen_val <- gen_defs_data %>%
      filter(Parameter == p_name) %>%
      pull(Default)
    if (length(gen_val) == 1 && !is.na(gen_val) && gen_val != "") {
      return(gen_val[1])
    }
    
    return(NA_character_) # Should ideally not be reached if defaults are comprehensive
  }
  
  # --- 4c. Construct the final `patched` dataframe ---
  # Determine the source of allele column structure and their Locus/Allele definitions
  if (params_loaded_from_file) {
    source_col_structure <- params
    # Ensure Locus/Allele rows exist in params, if not, this will error or give NAs.
    # It's assumed params file is well-formed if it exists.
    source_locus_row <- params[params$Parameter == "Locus", ]
    source_alleles_row <- params[params$Parameter == "Alleles", ]
    final_allele_col_names <- names(params)[names(params) != "Parameter"]
  } else { # params is a copy of defaults
    source_col_structure <- defaults
    source_locus_row <- defaults[defaults$Parameter == "Locus", ]
    source_alleles_row <- defaults[defaults$Parameter == "Alleles", ]
    final_allele_col_names <- names(defaults)[!names(defaults) %in% c("Parameter", "Default")]
  }
  
  # Initialize result_df with all parameters from defaults and allele columns
  result_df <- tibble(Parameter = defaults$Parameter)
  for (col_name in final_allele_col_names) {
    result_df[[col_name]] <- NA_character_
  }
  
  # Populate Locus and Alleles rows in result_df directly
  idx_locus <- match("Locus", result_df$Parameter)
  idx_alleles <- match("Alleles", result_df$Parameter)
  
  if (length(final_allele_col_names) > 0) { # Only if there are allele columns
    for (col_name in final_allele_col_names) {
      if (!is.na(idx_locus) && col_name %in% names(source_locus_row)) {
        result_df[idx_locus, col_name] <- source_locus_row[[col_name]]
      }
      if (!is.na(idx_alleles) && col_name %in% names(source_alleles_row)) {
        result_df[idx_alleles, col_name] <- source_alleles_row[[col_name]]
      }
    }
  }
  
  
  # Fill data cells for other parameters
  data_parameters <- result_df$Parameter[!result_df$Parameter %in% c("Locus", "Alleles")]
  
  if (length(final_allele_col_names) > 0) { # Only if there are allele columns
    for (p_name in data_parameters) {
      current_param_row_in_params <- if(params_loaded_from_file && p_name %in% params$Parameter) params[params$Parameter == p_name, ] else NULL
      current_param_row_in_defaults <- if(p_name %in% defaults$Parameter) defaults[defaults$Parameter == p_name, ] else NULL
      
      for (col_name in final_allele_col_names) {
        target_locus <- if(col_name %in% names(source_locus_row)) source_locus_row[[col_name]] else NA_character_
        target_allele <- if(col_name %in% names(source_alleles_row)) source_alleles_row[[col_name]] else NA_character_
        
        if(is.na(target_locus) || is.na(target_allele)) { # Column not properly defined by Locus/Allele rows
          # Fill with generic default for the parameter if possible, or NA
          val_from_generic_default <- generic_defaults %>% filter(Parameter == p_name) %>% pull(Default)
          result_df[result_df$Parameter == p_name, col_name] <- if(length(val_from_generic_default)==1) val_from_generic_default else NA_character_
          next # Move to next column
        }
        
        param_value <- NA_character_
        # Try to get value from the original params file (if it was loaded)
        if (params_loaded_from_file && !is.null(current_param_row_in_params) && col_name %in% names(current_param_row_in_params)) {
          param_value <- current_param_row_in_params[[col_name]]
        } else if (!params_loaded_from_file) { 
          # If params was not loaded from file, it means params IS defaults.
          # So, the "original" value comes from the defaults table itself for this cell.
          if (!is.null(current_param_row_in_defaults) && col_name %in% names(current_param_row_in_defaults)) {
            param_value <- current_param_row_in_defaults[[col_name]]
          }
        }
        
        if (!is.na(param_value) && param_value != "") {
          result_df[result_df$Parameter == p_name, col_name] <- param_value
        } else {
          # Value is missing or empty, need to find default using the new function
          default_val_to_use <- find_default_value(p_name, target_locus, target_allele, defaults_long, generic_defaults)
          result_df[result_df$Parameter == p_name, col_name] <- default_val_to_use
        }
      }
    }
  }
  
  patched <- result_df %>%
    arrange(match(Parameter, defaults$Parameter)) # Ensure order is same as defaults$Parameter
  
  
  ## --- 5.  Force correct types (everything read as character so far) ----
  numeric_variables <- c(
    "Carrier_Fitness_Males", "Carrier_Fitness_Females", "Carrier_Sex_Ratio_Males",
    "Carrier_Sex_Ratio_Females", "Multiplicative_Fitness_Males", "Multiplicative_Fitness_Females",
    "nonFunctional_NHEJ_Fraction_Locus", "HDR_Penalty_Locus", "Somatic_Cutting_Sensitivity_Locus",
    "Linkage_Disequilibrium",
    paste0("Homo_Phenotype", rep(1:10, each = 3), c("_Carrier_Fitness", "_Homozygous_Fitness", "_Sex_Ratio"))
  )
  
  logical_variables <- c(
    paste0("Homo_Phenotype", rep(1:10, each = 2), c("_Male", "_Female"))
  )
  
  # Create Alleles_File (transposed data)
  if (ncol(patched) > 1) {
    Alleles_File <- as.data.frame(t(patched[, -1, drop = FALSE]))
  } else {
    # Handle case where patched might only have the "Parameter" column (no allele columns)
    Alleles_File <- data.frame(matrix(ncol = nrow(patched), nrow = 0))
  }
  names(Alleles_File) <- patched$Parameter
  
  # Apply type conversions
  # Identify which columns in Alleles_File (now representing parameters) need conversion
  cols_to_convert_numeric <- intersect(names(Alleles_File), numeric_variables)
  cols_to_convert_logical <- intersect(names(Alleles_File), logical_variables)
  
  if (length(cols_to_convert_numeric) > 0) {
    Alleles_File[cols_to_convert_numeric] <- lapply(Alleles_File[, cols_to_convert_numeric, drop = FALSE], as.numeric)
  }
  if (length(cols_to_convert_logical) > 0) {
    Alleles_File[cols_to_convert_logical] <- lapply(Alleles_File[, cols_to_convert_logical, drop = FALSE], as.logical)
  }
  
  string_variables <- c("Locus", "Alleles", "Fitness_Applied_At_Stage_Locus")
  
  ## --- 2.  FLAG UNDECLARED COLUMNS (from original script, section was misnumbered as 2) ----
  expected_cols <- c(numeric_variables, logical_variables, string_variables)
  unexpected_cols <- setdiff(names(Alleles_File), expected_cols)
  
  if (length(unexpected_cols) > 0) {
    message(
      "The following columns in Alleles_File are NOT in numeric_variables, ",
      "logical_variables, or string_variables: ",
      paste(unexpected_cols, collapse = ", ")
    )
  }
  
 
  if ("A" %in% Alleles_File$Locus) {
    A <- Alleles_File[Alleles_File$Locus == "A", ]
  } else {
    stop("No A Locus Allele")
  }
  
  if ("B" %in% Alleles_File$Locus) {
    B <- Alleles_File[Alleles_File$Locus == "B", ]
  } else {
    B <- Alleles_File[Alleles_File$Locus == "A" & Alleles_File$Allele == "T", ]
    B[numeric_variables] <- 1
    B[logical_variables] <- FALSE
    # B[c("Carrier_Sex_Ratio", "Homo_Phenotype1_Sex_Ratio", "Homo_Phenotype2_Sex_Ratio")] <-  0.5
  }
  
    if ("C" %in% Alleles_File$Locus) {
    C <- Alleles_File[Alleles_File$Locus == "C", ]
  } else {
    C <- Alleles_File[Alleles_File$Locus == "A" & Alleles_File$Allele == "T", ]
    C[numeric_variables] <- 1
    C[logical_variables] <- FALSE
    # C[c("Carrier_Sex_Ratio", "Homo_Phenotype1_Sex_Ratio", "Homo_Phenotype2_Sex_Ratio")] <-  0.5
  }

  return(list(A = A,
              B = B,
              C = C))
}

Default_Geno_Freq <- function(Alpha, Freq, A, B, C, type, value) {
  Freq <- Freq / sum(Freq)
  Alpha$Freq <- 0

  if (value > 1) {
    stop("New genotype frequency set to over 100%")
  }
  
  Sep <- FALSE
  Custom <- FALSE
  running_sum <- 0

  
  if (type == "Hom_Comb") {
    Homozygotes_Freq <- value
    Hetrozygotes_Freq <- 0
    
  } else if (type == "Het_Comb") {
    Homozygotes_Freq <- 0
    Hetrozygotes_Freq <- value
    
  } else if (type == "Het_Comb_SC") {
    Homozygotes_Freq <- 0
    Hetrozygotes_Freq <- value
    
    # Here we change the heterozygote drive frequency
    Alpha[(
      Alpha$A_M == A$Allele[1] & Alpha$A_F == c("T") &
        Alpha$B_M == B$Allele[1] & Alpha$B_F == c("T") &
        Alpha$C_M == C$Allele[1] & Alpha$C_F == c("T") &
        Alpha$grna_dep_M == levels(Alpha$grna_dep_M)[1] &
        Alpha$grna_dep_F == levels(Alpha$grna_dep_F)[1] &  
        Alpha$Dep_M == levels(Alpha$Dep_M)[1] &
        Alpha$Dep_F == 0 & 
        Alpha$SD == 1
    ) # Note that deposition is set to be only paternal
    , length(Alpha[1,])] <- Hetrozygotes_Freq
    
    Custom <- TRUE
    running_sum <-  value
    
  } else if (type == "Het_Comb_SC_2") {
    Homozygotes_Freq <- 0
    Hetrozygotes_Freq <- value
    
    # Here we change the heterozygote drive frequency
    Alpha[(
      Alpha$A_M == A$Allele[2] & Alpha$A_F == c("T") &
        Alpha$B_M == c("T") & Alpha$B_F == c("T") &
        Alpha$C_M == c("T") & Alpha$C_F == c("T") &
        Alpha$grna_dep_M == levels(Alpha$grna_dep_M)[1] &
        Alpha$grna_dep_F == levels(Alpha$grna_dep_F)[1] &  
        Alpha$Dep_M == levels(Alpha$Dep_M)[1] &
        Alpha$Dep_F == 0 & 
        Alpha$SD == 1
    ) # Note that deposition is set to be only paternal
    , length(Alpha[1,])] <- Hetrozygotes_Freq
    
    Custom <- TRUE
    running_sum <-  value
    
  } else if (type == "Het_AC") {
    
    
    Custom <- TRUE
    
    Alpha[(
      Alpha$A_M == A$Allele[1] & Alpha$A_F == c("T") &
        Alpha$B_M == "T" & Alpha$B_F == c("T") &
        Alpha$C_M == C$Allele[1] & Alpha$C_F == c("T") &
        Alpha$grna_dep_M == levels(Alpha$grna_dep_M)[1] &
        Alpha$grna_dep_F == levels(Alpha$grna_dep_F)[1] &  
        Alpha$Dep_M == levels(Alpha$Dep_M)[1] &
        Alpha$Dep_F == 0 &
        Alpha$SD == 0
    ) # Note that deposition is set to be only paternal
    , length(Alpha[1,])] <- value
    
    running_sum <-  value
    
  } else if (type == "Standing_r1") {
    
    Custom <- TRUE
    
    Alpha[(
      Alpha$A_M == "T" & Alpha$A_F == "r1" &
        Alpha$B_M == "T" & Alpha$B_F == "T" &
        Alpha$C_M == "T" & Alpha$C_F == "T" &
        Alpha$grna_dep_M == levels(Alpha$grna_dep_M)[1] &
        Alpha$grna_dep_F == levels(Alpha$grna_dep_F)[1] &  
        Alpha$Dep_M == levels(Alpha$Dep_M)[1] &
        Alpha$Dep_F == 0 &
        Alpha$SD == 0
    ) # Note that deposition is set to be only paternal
    , length(Alpha[1,])] <- value * (1-value) * 2
    
    Alpha[(
      Alpha$A_M == "r1" & Alpha$A_F == "r1" &
        Alpha$B_M == "T" & Alpha$B_F == "T" &
        Alpha$C_M == "T" & Alpha$C_F == "T" &
        Alpha$grna_dep_M == levels(Alpha$grna_dep_M)[1] &
        Alpha$grna_dep_F == levels(Alpha$grna_dep_F)[1] &  
        Alpha$Dep_M == levels(Alpha$Dep_M)[1] &
        Alpha$Dep_F == 0 &
        Alpha$SD == 0
    ) # Note that deposition is set to be only paternal
    , length(Alpha[1,])] <- value^2
    
    running_sum <-  1-(1-value)^2
    
  } else if (type == "Hom_Sep") {
    Homozygotes_Freq <- value
    Hetrozygotes_Freq <- 0
    
    Sep <- TRUE
    
  } else if (type == "Het_Sep") {
    
    Homozygotes_Freq <- 0
    Hetrozygotes_Freq <- value
    
    Sep <- TRUE
    
  } else if (type == "WT") {
    
    Homozygotes_Freq <- value
    Hetrozygotes_Freq <- 0
    
    Custom <- TRUE
    
    Alpha[(
      Alpha$A_M == "T" & Alpha$A_F == "T" &
        Alpha$B_M == "T" & Alpha$B_F == "T" &
        Alpha$C_M == "T" & Alpha$C_F == "T" &
        Alpha$grna_dep_M == levels(Alpha$grna_dep_M)[1] &
        Alpha$grna_dep_F == levels(Alpha$grna_dep_F)[1] &  
        Alpha$Dep_M == levels(Alpha$Dep_M)[1] &
        Alpha$Dep_F == 0 &
        Alpha$SD == 0
    ) # Note that deposition is set to be only paternal
    , length(Alpha[1,])] <- value
    
  } else {
    stop("No valid genotype preset selected")
  }
  
  
  if (Sep == FALSE & Custom == FALSE) {
    # Here we change the homozygote drive frequency
    Alpha[(
      Alpha$A_M == A$Allele[1] & Alpha$A_F == A$Allele[1] &
        Alpha$B_M == B$Allele[1] & Alpha$B_F == B$Allele[1] &
        Alpha$C_M == C$Allele[1] & Alpha$C_F == C$Allele[1] &
        Alpha$grna_dep_M == levels(Alpha$grna_dep_M)[1] &
        Alpha$grna_dep_F == levels(Alpha$grna_dep_F)[1] &  
        Alpha$Dep_M == levels(Alpha$Dep_M)[1] &
        Alpha$Dep_F == levels(Alpha$Dep_F)[1] &
        Alpha$SD == 0
    ) # Note that deposition is set to be from both parents (but does not matter, as they are homo and have no 'T' that can be cut)
    , length(Alpha[1,])] <- Homozygotes_Freq
    
    # Here we change the heterozygote drive frequency
    Alpha[(
      Alpha$A_M == A$Allele[1] & Alpha$A_F == c("T") &
        Alpha$B_M == B$Allele[1] & Alpha$B_F == c("T") &
        Alpha$C_M == C$Allele[1] & Alpha$C_F == c("T") &
        Alpha$grna_dep_M == levels(Alpha$grna_dep_M)[1] &
        Alpha$grna_dep_F == levels(Alpha$grna_dep_F)[1] &  
        Alpha$Dep_M == levels(Alpha$Dep_M)[1] &
        Alpha$Dep_F == 0 & 
        Alpha$SD == 0
    ) # Note that deposition is set to be only paternal
    , length(Alpha[1,])] <- Hetrozygotes_Freq
    
    running_sum <- running_sum + Homozygotes_Freq + Hetrozygotes_Freq
    
  } else if (Custom == FALSE) {
    
    if (A$Allele[1] != "T") {
      
      # Here we change the homozygote drive frequency
      Alpha[(
        Alpha[, 1] == A$Allele[1] & Alpha$A_F == A$Allele[1] &
          Alpha$B_M == "T" & Alpha$B_F == "T" &
          Alpha$C_M == "T" & Alpha$C_F == "T" &
          Alpha$grna_dep_M == levels(Alpha$grna_dep_M)[1] &
          Alpha$grna_dep_F == levels(Alpha$grna_dep_F)[1] &  
          Alpha$Dep_M == levels(Alpha$Dep_M)[1] &
          Alpha$Dep_F == levels(Alpha$Dep_F)[1] &
          Alpha$SD == 0
      ) # Note that deposition is set to be from both parents (but does not matter, as they are homo and have no 'T' that can be cut)
      , length(Alpha[1,])] <- Homozygotes_Freq
      
      # Here we change the heterozygote drive frequency
      Alpha[(
        Alpha[, 1] == A$Allele[1] & Alpha$A_F == "T" &
          Alpha$B_M == "T" & Alpha$B_F == "T" &
          Alpha$C_M == "T" & Alpha$C_F == "T" &
          Alpha$grna_dep_M == levels(Alpha$grna_dep_M)[1] &
          Alpha$grna_dep_F == levels(Alpha$grna_dep_F)[1] &  
          Alpha$Dep_M == levels(Alpha$Dep_M)[1] &
          Alpha$Dep_F == levels(Alpha$Dep_F)[1] &
          Alpha$SD == 0
      ) # Note that deposition is set to be from both parents (but does not matter, as they are homo and have no 'T' that can be cut)
      , length(Alpha[1,])] <- Hetrozygotes_Freq
      
      running_sum <- running_sum + Homozygotes_Freq + Hetrozygotes_Freq
      
    }
    
    if (B$Allele[1] != "T") {
      
      Alpha[(
        Alpha[, 1] == "T" & Alpha$A_F == "T" &
          Alpha$B_M == B$Allele[1] & Alpha$B_F == B$Allele[1] &
          Alpha$C_M == "T" & Alpha$C_F == "T" &
          Alpha$grna_dep_M == levels(Alpha$grna_dep_M)[1] &
          Alpha$grna_dep_F == levels(Alpha$grna_dep_F)[1] &        
          Alpha$Dep_M == levels(Alpha$Dep_M)[1] &
          Alpha$Dep_F == levels(Alpha$Dep_F)[1] &
          Alpha$SD == 0
      ) # Note that deposition is set to be from both parents (but does not matter, as they are homo and have no 'T' that can be cut)
      , length(Alpha[1,])] <- Homozygotes_Freq
      
      
      Alpha[(
        Alpha[, 1] == "T" & Alpha$A_F == "T" &
          Alpha$B_M == B$Allele[1] & Alpha$B_F == "T" &
          Alpha$C_M == "T" & Alpha$C_F == "T" &
          Alpha$grna_dep_M == levels(Alpha$grna_dep_M)[1] &
          Alpha$grna_dep_F == levels(Alpha$grna_dep_F)[1] &        
          Alpha$Dep_M == levels(Alpha$Dep_M)[1] &
          Alpha$Dep_F == levels(Alpha$Dep_F)[1] &
          Alpha$SD == 0
      ) # Note that deposition is set to be from both parents (but does not matter, as they are homo and have no 'T' that can be cut)
      , length(Alpha[1,])] <- Hetrozygotes_Freq
      
      running_sum <- running_sum + Homozygotes_Freq + Hetrozygotes_Freq
      
    }
    
    if (C$Allele[1] != "T") {
      
      Alpha[(
        Alpha[, 1] == "T" & Alpha$A_F == "T" &
          Alpha$B_M == "T" & Alpha$B_F == "T" &
          Alpha$C_M == C$Allele[1] & Alpha$C_F == C$Allele[1] &
          Alpha$grna_dep_M == levels(Alpha$grna_dep_M)[1] &
          Alpha$grna_dep_F == levels(Alpha$grna_dep_F)[1] &  
          Alpha$Dep_M == levels(Alpha$Dep_M)[1] &
          Alpha$Dep_F == levels(Alpha$Dep_F)[1] &
          Alpha$SD == 0
      ) # Note that deposition is set to be from both parents (but does not matter, as they are homo and have no 'T' that can be cut)
      , length(Alpha[1,])] <- Homozygotes_Freq
      
      Alpha[(
        Alpha[, 1] == "T" & Alpha$A_F == "T" &
          Alpha$B_M == "T" & Alpha$B_F == "T" &
          Alpha$C_M == C$Allele[1] & Alpha$C_F == C$Allele[1] &
          Alpha$grna_dep_M == levels(Alpha$grna_dep_M)[1] &
          Alpha$grna_dep_F == levels(Alpha$grna_dep_F)[1] &  
          Alpha$Dep_M == levels(Alpha$Dep_M)[1] &
          Alpha$Dep_F == levels(Alpha$Dep_F)[1] &
          Alpha$SD == 0
      ) # Note that deposition is set to be from both parents (but does not matter, as they are homo and have no 'T' that can be cut)
      , length(Alpha[1,])] <- Hetrozygotes_Freq
      
      
      running_sum <- running_sum + Homozygotes_Freq + Hetrozygotes_Freq
      
    }
  }
  
  Alpha$Freq <- Alpha$Freq + (Freq * (1-running_sum))

  return(Alpha$Freq/sum(Alpha$Freq))
}

##  ............................................................................
##  Conversion matrix set-up                                                ####

# Here we set up the:
# nuclease expression (C) and deposition (D) conversion matrices
# one C for each sex, only one D for both
Conversion_D_and_G <-
  function(Alpha, nuc_parameters_list, geno_list) {
    # By default all our simulations have three loci
    # two alleles for each locus
    numb_alleles <- 6
    
    # we use letter to refer to the locus, and numbers to the allele
    allele_names <- c("A_M", "A_F",
                      "B_M", "B_F",
                      "C_M", "C_F")
    
    
    
    # parameters that change genotype to genotype
    # we set some defaults that can change in the
    # genotype loop
    # we use nonsense values here so that it is clear when
    # they are not correctly changed in the loop
    df_g <- data.frame(
      # genotype specific
      t_t_allele_hdr = as.integer(1),
      grna_present_exp = rep(as.integer(1), numb_alleles), # changed in loop
      cas9_present_exp = as.integer(1), # changed in loop
      cut_male = as.double(-0.5), # changed in loop
      cut_female = as.double(-0.5), # changed in loop
      cut_soma_male = as.double(-0.5), # changed in loop
      cut_soma_female = as.double(-0.5), # changed in loop
      grna_present_dep_male = as.integer(1), # changed in loop
      grna_present_dep_female = as.integer(1), # changed in loop
      cas9_present_dep_male = as.integer(1),
      cas9_present_dep_female = as.integer(1),
      cas9grna_complex_cut_male_dep = as.double(-0.5),
      cas9grna_complex_cut_female_dep = as.double(-0.5),
      cas9grna_complex_cut_both_dep = as.double(-0.5),
      grna_exp_but_not_dep_male = as.integer(1), # changed in loop
      grna_exp_but_not_dep_female = as.integer(1), # changed in loop
      shadowdrive_cut_male_dep = as.double(-0.5),
      shadowdrive_cut_female_dep = as.double(-0.5),
      shadowdrive_cut_both_dep = as.double(-0.5)
      # with shadow drive dep hdr can depend on the genotype
    )
    # here we name the 6 alleles
    row.names(df_g) <- allele_names
    
    # here we create a vector with the
    # type 2 to type 1 mutation ratio
    # in our simulations this only matters
    # for the A locus. B is neutral, C is never cut
    # we use the same ratio for the germline and embryo, except
    # the hdr is much lower in the embryo, meaning more resistant
    # alleles are created with the given r2:r1 ratio
    ratio_type2_to_type1 <- c(
      geno_list$A$nonFunctional_NHEJ_Fraction_Locus[1],
      geno_list$A$nonFunctional_NHEJ_Fraction_Locus[1],
      geno_list$B$nonFunctional_NHEJ_Fraction_Locus[1],
      geno_list$B$nonFunctional_NHEJ_Fraction_Locus[1],
      geno_list$C$nonFunctional_NHEJ_Fraction_Locus[1],
      geno_list$C$nonFunctional_NHEJ_Fraction_Locus[1]
    )

    # locus-specific HDR penalty (value per allele, but identical within a locus)
    hdr_penalty_by_locus <- c(
      geno_list$A$HDR_Penalty_Locus[1],
      geno_list$A$HDR_Penalty_Locus[1],
      geno_list$B$HDR_Penalty_Locus[1],
      geno_list$B$HDR_Penalty_Locus[1],
      geno_list$C$HDR_Penalty_Locus[1],
      geno_list$C$HDR_Penalty_Locus[1]
    )
    
    # parameters that don't change by genotype
    df_o <- data.frame(
      # applies to all genotypes
      cut_expressed_male = nuc_parameters_list$Germline_Cutrate_Male,
      cut_expressed_soma_male = nuc_parameters_list$Somatic_Cutrate_Male,
      hdr_germline_male = nuc_parameters_list$Germline_HDR_Male * hdr_penalty_by_locus,
      cut_expressed_female = nuc_parameters_list$Germline_Cutrate_Female,
      cut_expressed_soma_female = nuc_parameters_list$Somatic_Cutrate_Female,
      hdr_germline_female = nuc_parameters_list$Germline_HDR_Female * hdr_penalty_by_locus,
      r2r1_ratio_germline = ratio_type2_to_type1,# varies by allele
      cut_dep_male = nuc_parameters_list$Deposition_Cutrate_From_Male_Parent,
      cut_dep_female = nuc_parameters_list$Deposition_Cutrate_From_Female_Parent,
      cut_dep_both = as.double(-5), # calculated below
      hdr_dep = nuc_parameters_list$Deposition_HDR * hdr_penalty_by_locus,
      r2r1_ratio_embryo = ratio_type2_to_type1,# varies by allele
      Shadow_Drive_Cas9_Deposition_Only = nuc_parameters_list$Shadow_Drive_Cas9_Deposition_Only,
      shadow_drive_hdr = nuc_parameters_list$Deposition_HDR * hdr_penalty_by_locus,
      T_T_hdr_off = nuc_parameters_list$Homozygous_Wildtype_No_HDR
    )
    
    # If shadow drive is 'TRUE', we set its
    # hdr to be that of the germline
    # if it is 'FALSE', it still occurs,
    # but it is functionally identical to complex deposition
    if (nuc_parameters_list$Shadow_Drive_Cas9_Deposition_Only) {
      df_o$shadow_drive_hdr <- nuc_parameters_list$Germline_HDR_Male * hdr_penalty_by_locus
      
      # here we check if the germline hdr rates
      # are the same for males and females
      # it is ok if they are different,
      # so long as shadow drive is not also used.
      # This check can be narrowed to only throw an error if
      # they are different when shadow drive is set to occur
      # however, in our simulations they should not be different anyway
      # narrow this check if needed, or adapt the shadow drive process
      if (nuc_parameters_list$Germline_HDR_Male !=
          nuc_parameters_list$Germline_HDR_Female)
      {
        # shadow drive is not correctly set-up if the
        # hdr rates are different between the sexes
        stop("male and female HDR rates not the same with shadow drive")
      }
    }
    
    # here we name the 6 alleles
    row.names(df_o) <- allele_names
    
    
    # here we calculate the cut rate if both
    # maternal and paternal deposition occur
    df_o$cut_dep_both <- (1 -
                            ((1 - nuc_parameters_list$Deposition_Cutrate_From_Male_Parent) *
                               (1 - nuc_parameters_list$Deposition_Cutrate_From_Female_Parent)
                            ))
    
    # in our simulations no paternal deposition occurs
    # so if the combined deposition is not equal to that of
    # the maternal rate, something has gone wrong
    # adapt this is you do want paternal deposition
    if (round(df_o$cut_dep_both[1], 5) !=
        round(nuc_parameters_list$Deposition_Cutrate_From_Female_Parent[1],5)) {
      stop("combined deposition is not the same as maternal")
      # this should not occur in our specific simulations
    }
    
    ##### Conversion matrices set up #######
    
    # here we pre initialize the three matrices this function creates
    
    # the first is the deposition based cutting conversion matrix
    # deposition conversion is independent of the sex experiencing the deposition
    # here we set up a matrix of the intended size
    dep_conv_matrix <- sparseMatrix(i = length(Alpha[, 1]), j = length(Alpha[, 1]), x = 0, dims = c(length(Alpha[, 1]), length(Alpha[, 1])))
    
    # the second matrix is the male germline expression matrix
    # expression conversion is not sex independent
    # so one matrix for each sex
    exp_conv_male_matrix <- sparseMatrix(i = length(Alpha[, 1]), j = length(Alpha[, 1]), x = 0, dims = c(length(Alpha[, 1]), length(Alpha[, 1])))
    
    # the third matrix is the female germline expression matrix
    exp_conv_female_matrix <- sparseMatrix(i = length(Alpha[, 1]), j = length(Alpha[, 1]), x = 0, dims = c(length(Alpha[, 1]), length(Alpha[, 1])))

    exp_soma_conv_male_matrix <- sparseMatrix(i = length(Alpha[, 1]), j = length(Alpha[, 1]), x = 0, dims = c(length(Alpha[, 1]), length(Alpha[, 1])))
    
    # the third matrix is the female germline expression matrix
    exp_soma_conv_female_matrix <- sparseMatrix(i = length(Alpha[, 1]), j = length(Alpha[, 1]), x = 0, dims = c(length(Alpha[, 1]), length(Alpha[, 1])))
    
    ############ populating the matrices
    # This only needs to be done once per run,
    # as it is based on all the possible genotypes and parameters
    # (not dependent on the current state of the population)
    
    
    # Looping once over every possible diploid genotype
    # this includes different deposition states
    # we the figure out what other genotypes it
    # can be converted into by the drive
    for (i in seq_len(length(Alpha[, 1]))) {
      
      # Here we set the presence of the gRNA to 0
      # We will then set it to 1 if the gRNA is present
      df_g$grna_present_exp <- 0
      df_g$grna_present_dep_male <- 0
      df_g$grna_present_dep_female <- 0
      
      # We loop over loci_names (q)
      # and check if each of the gRNAs is expressed
      # in the current genotype (i)
      for (q in seq_len(numb_alleles)) {
        
        
        # Here we check if there at 'T/T' pairs
        # In this case, we assume HDR does not occur
        # We assume repeated cutting occurs until a non-HDR event
        if (Alpha[i, q] == "T" &
            Alpha[i, q + c(1, -1, 1, -1, 1, -1)[q]] == "T" &
            df_o$T_T_hdr_off[q]) {
          df_g$t_t_allele_hdr[q] <- 0} else {
            df_g$t_t_allele_hdr[q] <- 1}
        
        # This if statement checks if the combined string of
        # the genotype contains the gRNA for this particular loci (q)
        # If the gRNA is present we set the gRNA status to 1, else it stays 0
        # here we check if the gRNA has been expressed for the particular allele
        if (grepl(c("A", "A", "B", "B", "C", "C")[q], paste(sapply(Alpha[i, seq_len(numb_alleles)], as.character), collapse = ""), fixed = TRUE)) 
        {df_g$grna_present_exp[q] <- 1} else
        {df_g$grna_present_exp[q] <- 0}
        
        # here we check if the gRNA has been deposited for the particular allele
        # paternal deposition
        if(grepl(c("A", "A", "B", "B", "C", "C")[q],
                 Alpha[i, "grna_dep_M"],
                 fixed = TRUE)) {
          df_g$grna_present_dep_male[q] <- 1} else {
            df_g$grna_present_dep_male[q] <- 0}
        
        # here we check if the gRNA has been deposited for the particular allele
        # maternal deposition
        if (grepl(c("A", "A", "B", "B", "C", "C")[q],
                  Alpha[i, "grna_dep_F"],
                  fixed = TRUE)) {
          df_g$grna_present_dep_female[q] <- 1} else {
            df_g$grna_present_dep_female[q] <- 0}
      }
      
      # here we check for gRNA expression without deposition
      # this would allow for shadow drive
      df_g$grna_exp_but_not_dep_male <-  as.integer(as.logical((df_g$grna_present_exp == 1 & df_g$grna_present_dep_male == 0)))
      
      df_g$grna_exp_but_not_dep_female <-  as.integer(as.logical((df_g$grna_present_exp == 1 & df_g$grna_present_dep_female == 0)))
      
      df_g$cas9_present_exp <- 0
      # Check if there is a 'N' allele in this genotype,
      # if so change the nuclease exp presence to 1
      if (grepl("N",
                paste(sapply(Alpha[i, seq_len(numb_alleles)], as.character),
                      collapse = "")))
      {df_g$cas9_present_exp <- 1} else 
      {df_g$cas9_present_exp <- 0}
      
      # cut rate for the particular allele
      df_g$cut_male <- df_g$cas9_present_exp *
        df_g$grna_present_exp *
        df_o$cut_expressed_male
      
      # cut rate for the particular allele
      df_g$cut_female <- df_g$cas9_present_exp *
        df_g$grna_present_exp *
        df_o$cut_expressed_female

              # cut rate for the particular allele
      df_g$cut_soma_male <- df_g$cas9_present_exp *
        df_g$grna_present_exp *
        df_o$cut_expressed_soma_male
      
      # cut rate for the particular allele
      df_g$cut_soma_female <- df_g$cas9_present_exp *
        df_g$grna_present_exp *
        df_o$cut_expressed_soma_female
      
      
      df_g$cas9_present_dep_male <- 0
      df_g$cas9_present_dep_female <- 0
      
      # Here we check if this genotype has deposition
      # if so, we set deposition to "1"
      # first for paternal dep
      if ((as.numeric(levels(Alpha[i, "Dep_M"]))[Alpha[i, "Dep_M"]] > 0))
      {df_g$cas9_present_dep_male <- 1} else
      {df_g$cas9_present_dep_male <- 0}
      
      # first for maternal dep
      if ((as.numeric(levels(Alpha[i, "Dep_F"]))[Alpha[i, "Dep_F"]] > 0))
      {df_g$cas9_present_dep_female <- 1} else
      {df_g$cas9_present_dep_female <- 0}
      
      # cut rate by deposition for the particular allele
      df_g$cas9grna_complex_cut_male_dep <-
        df_g$cas9_present_dep_male *
        df_g$grna_present_dep_male *
        df_o$cut_dep_male
      
      # shadow drive!
      # expressed gRNA, with deposited Cas9
      # cut rate by deposition for the particular allele
      df_g$shadowdrive_cut_male_dep <- 
        df_g$cas9_present_dep_male *
        df_g$grna_exp_but_not_dep_male *
        df_o$cut_dep_male
      
      # cut rate by deposition for the particular allele
      df_g$cas9grna_complex_cut_female_dep <-
        df_g$cas9_present_dep_female *
        df_g$grna_present_dep_female *
        df_o$cut_dep_female
      
      # shadow drive!
      # expressed gRNA, with deposited Cas9
      # cut rate by deposition for the particular allele
      df_g$shadowdrive_cut_female_dep <-
        (df_g$cas9_present_dep_female *
           df_g$grna_exp_but_not_dep_female *
           df_o$cut_dep_female)
      
      
      # Here we calculate the male and female deposition based cut-rate
      # In our simulations, only maternal deposition occurs
      # it is calculated by multiplying the 1-cutrate values
      # for each sex
      
      df_g$cas9grna_complex_cut_both_dep <-
        1 - (
          (1 - df_g$cas9grna_complex_cut_male_dep) *
            (1 - df_g$cas9grna_complex_cut_female_dep)
        )
      
      df_g$shadowdrive_cut_both_dep <-
        1 - (
          (1 - df_g$shadowdrive_cut_male_dep) *
            (1 - df_g$shadowdrive_cut_female_dep)
        )
      
      # Here we check if shadow drive and complex deposition
      # occur simultaneously this should not happen in our simulations
      # In general it may be possible for shadow dep to come from the male
      # and complex dep from the female, but we have no paternal dep
      # adapt if needed for your uses
      if (sum((
        df_g$cas9grna_complex_cut_both_dep > 1) &
        (df_g$shadowdrive_cut_both_dep > 1)) > 1) {
        stop("shadow drive somehow occurring with cas9:gRNA complex deposition")
      }
      
      # we have now set up all the parameters for this conversion
      # we will now define the conversion outcomes
      
      # We will create a row with just the
      # conversion outcomes for the particular genotype
      # we are considering in this loop
      
      # We start with no change and set that to 1
      genotype_conversion_to <- Alpha[i, ]
      
      # We will create the conversion genotype allele by allele
      # We first reset all the alleles
      genotype_conversion_to[seq_len(numb_alleles)] <- ""
      
      # We will store the conversion probabilities
      # for each genotype as a column in the matrix
      genotype_conversion_to$dep_conversion_prob <- 1
      genotype_conversion_to$male_germline_conv_prob <- 1
      genotype_conversion_to$female_germline_conv_prob <- 1
      genotype_conversion_to$male_soma_conv_prob <- 1
      genotype_conversion_to$female_soma_conv_prob <- 1
      
      genotype_conversion_to$valid <- FALSE
      
      conversion_outcomes <- genotype_conversion_to[rep(1,(4^numb_alleles)),]
      conversion_outcomes[1,]$valid <- TRUE
      
      #genotype_conversion_to_before <- conversion_outcomes
      #temp_outcomes_1 <- conversion_outcomes
      
      # this is a table that we use to switch the parameters
      # per conversion outcome
      # this is used in the main loop later
      # the structure is:
      # (a(b-x))^c
      # we then use different values for a, b, c
      # (1(1-x))^1 = 1-x
      # (-1(0-x))^1 = x
      # (0(0-x))^0 = 1
      df_a <- data.frame(
        allele_outcome = c("T", "r1", "r2", "Homologous allele"),
        cut_a = c(1, -1, -1, -1), # 1-cut, cut, cut, cut
        cut_b = c(1, 0, 0, 0),    # 1-cut, cut, cut, cut
        cut_c = c(1, 1, 1, 1),    # 1-cut, cut, cut, cut
        hdr_a = c(0, 1, 1, -1),   # 1, 1-hdr, 1-hdr, hdr
        hdr_b = c(0, 1, 1, 0),    # 1, 1-hdr, 1-hdr, hdr
        hdr_c = c(0, 1, 1, 1),    # 1, 1-hdr, 1-hdr, hdr
        r21_a = c(0, 1, -1, 0),   # 1, 1-r2, r1, 1
        r21_b = c(0, 1, 0, 0),    # 1, 1-r2, r1, 1
        r21_c = c(0, 1, 1, 0)     # 1, 1-r2, r1, 1
      )
      
      # io keeps track of the amount of outcomes we have
      io <- 1
      
      
      # We loop over each allele (z)
      # of the current genotype (i)
      # and build up the probability of each possible
      # conversion outcome, allele by allele
      # For the first loop of z the
      # possible conversions is only the starting genotype
      # and this stays so until an allele is a T (and can be cut)
      for (z in allele_names) {
        
        z_index <-  which(allele_names == z)
        
        # We check if the current allele is a T
        # The only allele that can change
        if (Alpha[i, z] == "T" &
            (sum(df_g[z, c("cas9grna_complex_cut_both_dep", "shadowdrive_cut_both_dep", "cut_male", "cut_female")]) > 0)) {
          
          df_a$allele_outcome[4] <- as.character(Alpha[i, z_index + c(1, -1, 1, -1, 1, -1)[z_index]])
          
          
          # If the current allele is a "T",
          # we create a copy of the previous possible conversion
          # outcomes (genotype_conversion_to)
          # all the previous possible alleles replicated but
          # now with each possible conversion of the T (no-cut, f, r, h)
          # As such, each new T allele increases the possible alleles by 4
          # and this can occur for every one of the alleles (z)
          
          # genotype_conversion_to_before stores the possible outcomes
          # before taking account the new alleles changes
          #genotype_conversion_to_before <- conversion_outcomes
          
          
          conversion_outcomes[1:(io*4),] <- conversion_outcomes[1:io,]
          
          for (k in seq_len(length(df_a$allele_outcome))) {
            
            
            # temp_outcomes_1 will store the possible outcomes
            # including the new allele, it will be reset for each of the
            # 4 possible conversions
            #temp_outcomes_1 <- genotype_conversion_to_before
            
            
            # Here we define what this genotype is being converted to
            #temp_outcomes_1[, z] <- df_a$allele_outcome[k]
            conversion_outcomes[(1 + (io * (k - 1))):(io + (io * (k - 1))), z_index] <- df_a$allele_outcome[k]
            
            # The probability of the new allele is multiplied
            # by the probability given the previous alleles (temp_Outcomes_1)
            
            # Deposition (Cas9:gRNA complex + shadow) based cutting
            conversion_outcomes[(1 + (io * (k - 1))):(io + (io * (k - 1))),]$dep_conversion_prob <-
              (((((df_a$cut_a[k] * (df_a$cut_b[k] - df_g[z, "cas9grna_complex_cut_both_dep"]))^df_a$cut_c[k]) *
                   ((df_a$hdr_a[k] * (df_a$hdr_b[k] - (df_o[z, "hdr_dep"] * df_g[z, "t_t_allele_hdr"])))^df_a$hdr_c[k]) *
                   ((df_a$r21_a[k] * (df_a$r21_b[k] - df_o[z, "r2r1_ratio_embryo"]))^df_a$r21_c[k])))*(1 - df_a$cut_b[k]) + # PLUS! what the shadow drive does
                 (1 - df_g[z, "cas9grna_complex_cut_both_dep"]) *
                 ((df_a$cut_a[k] * (df_a$cut_b[k] - df_g[z, "shadowdrive_cut_both_dep"]))^df_a$cut_c[k]) *
                 ((df_a$hdr_a[k] * (df_a$hdr_b[k] - (df_o[z, "shadow_drive_hdr"] * df_g[z, "t_t_allele_hdr"])))^df_a$hdr_c[k]) *
                 ((df_a$r21_a[k] * (df_a$r21_b[k] - df_o[z, "r2r1_ratio_embryo"]))^df_a$r21_c[k])) * 
              conversion_outcomes[(1 + (io * (k - 1))):(io + (io * (k - 1))),]$dep_conversion_prob
            
            # Male germline
            conversion_outcomes[(1 + (io * (k - 1))):(io + (io * (k - 1))),]$male_germline_conv_prob <-
              ((df_a$cut_a[k] * (df_a$cut_b[k] - df_g[z, "cut_male"]))^df_a$cut_c[k]) *
              ((df_a$hdr_a[k] * (df_a$hdr_b[k] - (df_o[z, "hdr_germline_male"] * df_g[z, "t_t_allele_hdr"])))^df_a$hdr_c[k]) *
              ((df_a$r21_a[k] * (df_a$r21_b[k] - df_o[z, "r2r1_ratio_germline"]))^df_a$r21_c[k]) *
              conversion_outcomes[(1 + (io * (k - 1))):(io + (io * (k - 1))),]$male_germline_conv_prob
            
            # Female germline
            conversion_outcomes[(1 + (io * (k - 1))):(io + (io * (k - 1))),]$female_germline_conv_prob <-
              ((df_a$cut_a[k] * (df_a$cut_b[k] - df_g[z, "cut_female"]))^df_a$cut_c[k]) *
              ((df_a$hdr_a[k] * (df_a$hdr_b[k] - (df_o[z, "hdr_germline_female"] * df_g[z, "t_t_allele_hdr"])))^df_a$hdr_c[k]) *
              ((df_a$r21_a[k] * (df_a$r21_b[k] - df_o[z, "r2r1_ratio_germline"]))^df_a$r21_c[k]) *
              conversion_outcomes[(1 + (io * (k - 1))):(io + (io * (k - 1))),]$female_germline_conv_prob


            # Male soma
            conversion_outcomes[(1 + (io * (k - 1))):(io + (io * (k - 1))),]$male_soma_conv_prob <-
              ((df_a$cut_a[k] * (df_a$cut_b[k] - df_g[z, "cut_soma_male"]))^df_a$cut_c[k]) *
              ((df_a$hdr_a[k] * (df_a$hdr_b[k] - (df_o[z, "hdr_dep"] * df_g[z, "t_t_allele_hdr"])))^df_a$hdr_c[k]) *
              ((df_a$r21_a[k] * (df_a$r21_b[k] - df_o[z, "r2r1_ratio_embryo"]))^df_a$r21_c[k]) *
              conversion_outcomes[(1 + (io * (k - 1))):(io + (io * (k - 1))),]$male_soma_conv_prob
            
            # Female soma
            conversion_outcomes[(1 + (io * (k - 1))):(io + (io * (k - 1))),]$female_soma_conv_prob <-
              ((df_a$cut_a[k] * (df_a$cut_b[k] - df_g[z, "cut_soma_female"]))^df_a$cut_c[k]) *
              ((df_a$hdr_a[k] * (df_a$hdr_b[k] - (df_o[z, "hdr_dep"] * df_g[z, "t_t_allele_hdr"])))^df_a$hdr_c[k]) *
              ((df_a$r21_a[k] * (df_a$r21_b[k] - df_o[z, "r2r1_ratio_embryo"]))^df_a$r21_c[k]) *
              conversion_outcomes[(1 + (io * (k - 1))):(io + (io * (k - 1))),]$female_soma_conv_prob
            
            # The new allele is appended to the conversion outcomes
            #conversion_outcomes <-
            #rbind(conversion_outcomes, temp_outcomes_1)
          }
          
          io <- io*4
          
          # If the current allele is not a T the overall conversion probabilities should stay the same
        } else if (Alpha[i, z] != "T" |
                   (sum(df_g[z,c("cas9grna_complex_cut_both_dep", "shadowdrive_cut_both_dep", "cut_male", "cut_female")]) == 0)) {
          
          
          # The original unchanged allele is
          # added to the list of outcome but with no change in probabilities
          # no additional conversion genotypes are created
          # the new allele is only added to the existing conversion outcomes
          #conversion_outcomes[, z] <- Alpha[i, z]
          
          conversion_outcomes[1:io, z_index] <- as.character(Alpha[i, z])
          
          # Error checking
        } else {
          stop("The 'if else' is not set up correctly")
        }
      }
      
      # Here we remove the outcomes that have a 0 probability
      conversion_outcomes_valid <- conversion_outcomes[(rowSums(conversion_outcomes[,c(
        "dep_conversion_prob",
         "male_germline_conv_prob",
          "female_germline_conv_prob",
          "male_soma_conv_prob",
          "female_soma_conv_prob"
          )]) > 0 & conversion_outcomes$valid),]
      
      # At the end of the process of adding each allele the total conversion probability should sum to 1
      if (abs((sum(conversion_outcomes_valid[,c("dep_conversion_prob", "male_germline_conv_prob", "female_germline_conv_prob", "male_soma_conv_prob", "female_soma_conv_prob")]) / 5) -
              1) > 1.5e-14) {
        print("Sum of possible outcomes is not 1")
        print(sum(conversion_outcomes_valid[,c("dep_conversion_prob", "male_germline_conv_prob", "female_germline_conv_prob", "male_soma_conv_prob", "female_soma_conv_prob")]) / 5)
        print(i)
        stop("Sum of possible outcomes is not 1")
      }
      
      # We should now have a list of all the possible conversion outcomes and probabilities for the current genotype
      # We now Update D and G with those conversions
      # We go through the possible conversion outcomes 1 by 1
      for (k in seq_len(length(conversion_outcomes_valid[, 1]))) {
      
        # For the current conversion outcome, we exclude all the indexes of the genotypes that are not the same by comparing them allele by allele (and deposition state)
        index_of_conv <- 
          which((
            Alpha[, 1] == conversion_outcomes_valid[k, 1] &
              Alpha[, 2] == conversion_outcomes_valid[k, 2] &
              Alpha[, 3] == conversion_outcomes_valid[k, 3] &
              Alpha[, 4] == conversion_outcomes_valid[k, 4] &
              Alpha[, 5] == conversion_outcomes_valid[k, 5] &
              Alpha[, 6] == conversion_outcomes_valid[k, 6] &
              Alpha[, 7] == conversion_outcomes_valid[k, 7] &
              Alpha[, 8] == conversion_outcomes_valid[k, 8] &
              Alpha[, 9] == conversion_outcomes_valid[k, 9] &
              Alpha[, 10] == conversion_outcomes_valid[k, 10] &
              Alpha[, 11] == conversion_outcomes_valid[k, 11]
            ))
        
        if (length(index_of_conv) > 1) {
          stop("Multiple outcomes match conversion result, this should be unique")
        }
        
        # We should now be left with a single outcome that matches our current possible outcome
        
        # We now change the conversion probability of the outcome genotype in D and H to match that of the probability of our possible outcome
        # We do addition to the value of D/H, because the HDR outcome can often be the same as a type 1 or 2 conversion. So the list of possible outcome is not necessarily unique meaning different processes have the same outcomes
        dep_conv_matrix[i, index_of_conv] <-
          conversion_outcomes_valid[k, "dep_conversion_prob"] + dep_conv_matrix[i, index_of_conv]
        
        exp_conv_male_matrix[i, index_of_conv] <-
          conversion_outcomes_valid[k, "male_germline_conv_prob"] + exp_conv_male_matrix[i, index_of_conv]
        
        exp_conv_female_matrix[i, index_of_conv] <-
          conversion_outcomes_valid[k, "female_germline_conv_prob"] + exp_conv_female_matrix[i, index_of_conv]

        exp_soma_conv_male_matrix[i, index_of_conv] <-
          conversion_outcomes_valid[k, "male_soma_conv_prob"] + exp_soma_conv_male_matrix[i, index_of_conv]
        
        exp_soma_conv_female_matrix[i, index_of_conv] <-
          conversion_outcomes_valid[k, "female_soma_conv_prob"] + exp_soma_conv_female_matrix[i, index_of_conv]
      }
    }
    
    if (round(((
      sum(dep_conv_matrix) + sum(exp_conv_male_matrix) + sum(exp_conv_female_matrix)
    ) / length(Alpha[, 1])),5) != 3) {
      # Can occur when a gRNA targets a locus without non-T alleles
      stop("D or G freq != 1")
    }

    # Per-matrix per-row validation: each row of each conversion matrix must sum to 1
    conv_matrices <- list(
      dep = dep_conv_matrix,
      exp_male = exp_conv_male_matrix,
      exp_female = exp_conv_female_matrix,
      soma_male = exp_soma_conv_male_matrix,
      soma_female = exp_soma_conv_female_matrix
    )
    for (mat_name in names(conv_matrices)) {
      row_sums <- rowSums(conv_matrices[[mat_name]])
      bad_rows <- which(abs(row_sums - 1) > 1e-10)
      if (length(bad_rows) > 0) {
        stop("Conversion matrix '", mat_name, "' has rows not summing to 1. ",
             "First bad row: ", bad_rows[1], " (sum = ", row_sums[bad_rows[1]], ")")
      }
    }

    D_G_List <-
      list(dep_conv_matrix = dep_conv_matrix,
           exp_conv_male_matrix = exp_conv_male_matrix,
           exp_conv_female_matrix = exp_conv_female_matrix,
           exp_soma_conv_male_matrix = exp_soma_conv_male_matrix,
           exp_soma_conv_female_matrix = exp_soma_conv_female_matrix
           )

    return(D_G_List)
  }

Conversion_H <- function(Alpha, geno_list, Delta, Dep) {
  
  # We consider only the non-sex converted case
  Delta <- Delta[Delta$SD == 0,]
  rownames(Delta) <- NULL
  Alpha <- Alpha[Alpha$SD == 0,]
  rownames(Alpha) <- NULL
  
  Alleles <- 6
  
  H <- Matrix(matrix(0, nrow = length(Alpha[, 1]), ncol = length(Delta[, 1])), sparse = TRUE)
  
  colnames(H) <-
    apply(Delta[,-6], 1, paste, collapse = ".")
  
  rownames(H) <-
    apply(Alpha[,-11], 1, paste, collapse = ".")
  
  for (k in 1:length(Alpha[, 1])) {
    # We make a copy of Alpha just so we can look-up the index of each conversion outcome genotype
    Lookup_Outcomes <- Delta
    
    # Check if there is a 'N' allele in this genotype, if so change the nuclease cut rate from 0
    if ((grepl("N", paste(sapply(Alpha[k, 1:Alleles], as.character), collapse = "")))) 
    {Expres <- Dep[1]} else 
    {Expres <- 0}
    
    if (length(levels(Delta$gRNA_Deposition)) > 1) {
      
      # Here we check if any of the gRNAs are expressed
      if (any(grepl("A", unlist(Alpha[k, 1:Alleles]), fixed = TRUE)))
      {a_grna <- "A"} else 
      {a_grna <- "_"}
      
      if (any(grepl("B", unlist(Alpha[k, 1:Alleles]), fixed = TRUE)))
      {b_grna <- "B"} else 
      {b_grna <- "_"}
      
      if (any(grepl("C", unlist(Alpha[k, 1:Alleles]), fixed = TRUE)))
      {c_grna <- "C"} else 
      {c_grna <- "_"}
      
      all_grnas <- paste(paste(a_grna, b_grna, sep = ""), c_grna, sep = "")
      
    } else if (length(levels(Delta$gRNA_Deposition)) == 1) {
      # grna status seemingly not being recorded
      # This will be done to speed up the calculation if deposition 
      # rate is 0 anyway
      all_grnas <- "___"
    } else {
      stop("gRNA check in H not set up correctly")
    }
    
    if(length(Lookup_Outcomes[(as.character(Lookup_Outcomes$Dep) == as.character(Expres) &
                               as.character(Lookup_Outcomes$gRNA_Deposition) == all_grnas),]$A) == 0) {
      warning("the deposition state is seemingly impossible")
      # needs to be update for multiple alleles
    }
    
    Lookup_Outcomes <-
      Lookup_Outcomes[(as.character(Lookup_Outcomes$Dep) == as.character(Expres) &
                         as.character(Lookup_Outcomes$gRNA_Deposition) == all_grnas),]
    
    # For the current conversion outcome, we exclude all the indexes 
    # of the genotypes that are not the same by comparing them allele by allele (and deposition state)
    # we have already taken care of deposition
    for (l in seq_len(length(Delta[1,]) - 4)) {
      # Set the new Lookup_Outcomes value to only include the alleles that have the same allele as the possible outcome
      Lookup_Outcomes <-
        Lookup_Outcomes[(
          as.character(Lookup_Outcomes[, l]) == as.character(Alpha[k, (l * 2) - 1]) |
            as.character(Lookup_Outcomes[, l]) == as.character(Alpha[k, (l * 2)])),]
    }
    
    if(length(Lookup_Outcomes[,1]) == 0) {
      stop("No haploid genotypes that match diploid source")
    }
    
    weight <- 1 / length(Lookup_Outcomes[, 1])
    
    for (e in 1:length(Lookup_Outcomes[, 1])) {
      # Set the new Lookup_Outcomes value to only include the alleles that have the same allele as the possible outcome
      H[k, as.double(rownames(Lookup_Outcomes))] <- weight
    }
  }
  
  if (sum(rowSums(H) != 1) > 0) {
    stop("H conversion matrix not set correctly")
  }
  
  return(H)
}

Fitness_Costs <- function(Alpha, Fit_Values, Sex) {
  
  # Checks 'Sex' to select the right sex-specific fitness values
  # Set s to 0 if male, or 1 if female
  if (Sex == "M") {
    s <- 0
  } else if (Sex == "F") {
    s <- 1
  } else {
    stop("No valid sex selected")
  }
  
  # We get the number of alleles
  Alleles <- length(Alpha[1,]) - 3
  
  # Fit_Sex will store the fitness costs per genotype for the sex this function is called for
  Fit_Sex <- rep(1, length(Alpha[, Alleles + 1]))
  
  Sex_Bias <- Fit_Sex
  
  merged_table <- rbind(Fit_Values[[1]],Fit_Values[[2]],Fit_Values[[3]])
  
  vector_complex_fitness <- 1
  
  for(i in 2:10) {
    
    if((any(merged_table[, grep(paste0("Homo_Phenotype",i,"_Carrier_Fitness"), colnames(merged_table))]<1) |
    any(merged_table[, grep(paste0("Homo_Phenotype",i,"_Homozygous_Fitness"), colnames(merged_table))]<1) |
    any(merged_table[, grep(paste0("Homo_Phenotype",i,"_Sex_Ratio"), colnames(merged_table))]<1))) {
      
      vector_complex_fitness <- cbind(vector_complex_fitness, i)
      
    }
  }
  
  # We select the fitness values fed into the function
  A <- Fit_Values[[1]]$Allele
  Ar <-
    Fit_Values[[1]][, grep("Multiplicative_Fitness_Males", colnames(Fit_Values[[1]])) + s]
  Ar_add <-
    Fit_Values[[1]][, grep("Carrier_Fitness_Males", colnames(Fit_Values[[1]])) + s]
  
  Asr_add <-
    Fit_Values[[1]][, grep("Carrier_Sex_Ratio_Males", colnames(Fit_Values[[1]])) + s]
  
  # Initialize lists to store phenotype data for A
  A_Homo_Fit_Type <- list()
  A_Homo_Fit_Cost <- list()
  A_Carrier_Fit_Cost <- list()
  A_Homo_SexRatio <- list()
  
  # Get data for all 10 phenotypes for A
  for(i in vector_complex_fitness) {
    A_Homo_Fit_Type[[i]] <- 
      Fit_Values[[1]][, grep(paste0("Homo_Phenotype",i,"_Male"), colnames(Fit_Values[[1]])) + s]
    A_Carrier_Fit_Cost[[i]] <-
      Fit_Values[[1]][, grep(paste0("Homo_Phenotype",i,"_Carrier_Fitness"), colnames(Fit_Values[[1]]))]
    A_Homo_Fit_Cost[[i]] <-
      Fit_Values[[1]][, grep(paste0("Homo_Phenotype",i,"_Homozygous_Fitness"), colnames(Fit_Values[[1]]))]
    A_Homo_SexRatio[[i]] <-
      Fit_Values[[1]][, grep(paste0("Homo_Phenotype",i,"_Sex_Ratio"), colnames(Fit_Values[[1]]))]
  }
  
  B <- Fit_Values[[2]]$Allele
  Br <-
    Fit_Values[[2]][, grep("Multiplicative_Fitness_Males", colnames(Fit_Values[[2]])) + s]
  Br_add <-
    Fit_Values[[2]][, grep("Carrier_Fitness_Males", colnames(Fit_Values[[2]])) + s]
  
  Bsr_add <-
    Fit_Values[[2]][, grep("Carrier_Sex_Ratio_Males", colnames(Fit_Values[[2]])) + s]
  
  # Initialize lists to store phenotype data for B  
  B_Homo_Fit_Type <- list()
  B_Homo_Fit_Cost <- list()
  B_Carrier_Fit_Cost <- list()
  B_Homo_SexRatio <- list()
  
  # Get data for all 10 phenotypes for B
  for(i in vector_complex_fitness) {
    B_Homo_Fit_Type[[i]] <-
      Fit_Values[[2]][, grep(paste0("Homo_Phenotype",i,"_Male"), colnames(Fit_Values[[2]])) + s]
    B_Carrier_Fit_Cost[[i]] <-
      Fit_Values[[2]][, grep(paste0("Homo_Phenotype",i,"_Carrier_Fitness"), colnames(Fit_Values[[2]]))]
    B_Homo_Fit_Cost[[i]] <-
      Fit_Values[[2]][, grep(paste0("Homo_Phenotype",i,"_Homozygous_Fitness"), colnames(Fit_Values[[2]]))]
    B_Homo_SexRatio[[i]] <-
      Fit_Values[[2]][, grep(paste0("Homo_Phenotype",i,"_Sex_Ratio"), colnames(Fit_Values[[2]]))]
  }
  
  C <- Fit_Values[[3]]$Allele
  Cr <-
    Fit_Values[[3]][, grep("Multiplicative_Fitness_Males", colnames(Fit_Values[[3]])) + s]
  Cr_add <-
    Fit_Values[[3]][, grep("Carrier_Fitness_Males", colnames(Fit_Values[[3]])) + s]
  Csr_add <-
    Fit_Values[[3]][, grep("Carrier_Sex_Ratio_Males", colnames(Fit_Values[[3]])) + s]
  
  # Initialize lists to store phenotype data for C
  C_Homo_Fit_Type <- list()
  C_Homo_Fit_Cost <- list()
  C_Carrier_Fit_Cost <- list()
  C_Homo_SexRatio <- list()
  
  # Get data for all 10 phenotypes for C
  for(i in vector_complex_fitness) {
    C_Homo_Fit_Type[[i]] <-
      Fit_Values[[3]][, grep(paste0("Homo_Phenotype",i,"_Male"), colnames(Fit_Values[[3]])) + s]
    C_Carrier_Fit_Cost[[i]] <-
      Fit_Values[[3]][, grep(paste0("Homo_Phenotype",i,"_Carrier_Fitness"), colnames(Fit_Values[[3]]))]
    C_Homo_Fit_Cost[[i]] <-
      Fit_Values[[3]][, grep(paste0("Homo_Phenotype",i,"_Homozygous_Fitness"), colnames(Fit_Values[[3]]))]
    C_Homo_SexRatio[[i]] <-
      Fit_Values[[3]][, grep(paste0("Homo_Phenotype",i,"_Sex_Ratio"), colnames(Fit_Values[[3]]))]
  }
  
  # We combine the allele names of A, B and C
  All_Alleles_Name <- c(A, B, C)
  
  # We combine the multaplicative fitness costs of all the alleles
  All_ArBrCr <- c(Ar, Br, Cr)
  
  # We combine the carrier fitness costs of all the alleles
  All_ArBrCr_add <- c(Ar_add, Br_add, Cr_add)
  
  # We combine the carrier sex ratio of all the alleles
  All_AsrBsrCsr_add <- c(Asr_add, Bsr_add, Csr_add)
  
  # Initialize lists for combined phenotype data
  All_Homo_Fit_Type <- list()
  All_Homo_Fit_Cost <- list()
  All_Carrier_Fit_Cost <- list()
  All_Homo_SexRatio <- list()
  
  # Combine phenotype data for all alleles
  for(i in vector_complex_fitness) {
    All_Homo_Fit_Type[[i]] <- c(A_Homo_Fit_Type[[i]], B_Homo_Fit_Type[[i]], C_Homo_Fit_Type[[i]])
    All_Homo_Fit_Cost[[i]] <- c(A_Homo_Fit_Cost[[i]], B_Homo_Fit_Cost[[i]], C_Homo_Fit_Cost[[i]])
    All_Carrier_Fit_Cost[[i]] <- c(A_Carrier_Fit_Cost[[i]], B_Carrier_Fit_Cost[[i]], C_Carrier_Fit_Cost[[i]])
    All_Homo_SexRatio[[i]] <- c(A_Homo_SexRatio[[i]], B_Homo_SexRatio[[i]], C_Homo_SexRatio[[i]])
    
    if (
      length(unique(A_Homo_Fit_Cost[[i]])) > 1 |
      length(unique(B_Homo_Fit_Cost[[i]])) > 1 |
      length(unique(C_Homo_Fit_Cost[[i]])) > 1 |
      length(unique(A_Carrier_Fit_Cost[[i]])) > 1 |
      length(unique(B_Carrier_Fit_Cost[[i]])) > 1 |
      length(unique(C_Carrier_Fit_Cost[[i]])) > 1 |
      length(unique(A_Homo_SexRatio[[i]])) > 1 |
      length(unique(B_Homo_SexRatio[[i]])) > 1 |
      length(unique(C_Homo_SexRatio[[i]])) > 1
        ) {
      stop(
        paste0(
          "Homozygous KO fitness costs must be the same for all alleles in phenotype ",
          i,
          ", otherwise combinations cannot be extrapolated"
        )
      )
    }
  }
  
  # This matrix will hold the individual reduced fitness 'contributions' of each allele
  # This will be used to assign 'blame'
  Individual_Fitness_Values <-
    matrix(data = 1,
           nrow = length(Alpha[, 1]),
           ncol = 6)
  
  Individual_SexRatio_Values <- Individual_Fitness_Values
  
  # Hap_Stuff will note (with a change to '1') if either allele is haplosufficient and is initiated at 0
  # Only the columns corresponding to 'r' (1,3,5) will be used
  Hap_Suff <-
    matrix(0, ncol = 6, nrow = length(Alpha[, Alleles + 1]))
  
  Hap_Suff <- list(
    Hap_Suff,
    Hap_Suff,
    Hap_Suff,
    Hap_Suff,
    Hap_Suff,
    Hap_Suff,
    Hap_Suff,
    Hap_Suff,
    Hap_Suff,
    Hap_Suff
  )
  
  
  # 'r' is the index for the diploid alleles locus (A=1, B=3, C=5)
  r <- 1
  
  # Here we loop over every allele at every locus, starting with the first allele of A
  for (q in 1:length(All_Alleles_Name)) {
    # If both match, we apply the individual fitness cost to the first allele twice
    Individual_Fitness_Values[(Alpha[, r + 0] == All_Alleles_Name[q] &
                                 Alpha[, r + 1] == All_Alleles_Name[q]), r] <-
      Individual_Fitness_Values[(Alpha[, r + 0] == All_Alleles_Name[q] &
                                   Alpha[, r + 1] == All_Alleles_Name[q]), r] * All_ArBrCr_add[q] * All_ArBrCr[q] * All_ArBrCr[q]
    
    # If both alleles match, we apply the sex-ratio carrier effect to the first allele
    Individual_SexRatio_Values[(Alpha[, r + 0] == All_Alleles_Name[q] &
                                 Alpha[, r + 1] == All_Alleles_Name[q]), r] <-
      Individual_SexRatio_Values[(Alpha[, r + 0] == All_Alleles_Name[q] &
                                   Alpha[, r + 1] == All_Alleles_Name[q]), r] * All_AsrBsrCsr_add[q]
    
    # If only the first matches, it gets all the blame (r)
    Individual_Fitness_Values[(Alpha[, r + 0] == All_Alleles_Name[q] &
                                 Alpha[, r + 1] != All_Alleles_Name[q]), r] <-
      Individual_Fitness_Values[(Alpha[, r + 0] == All_Alleles_Name[q] &
                                   Alpha[, r + 1] != All_Alleles_Name[q]), r] * All_ArBrCr_add[q] * All_ArBrCr[q]
    
    # If only the second matches, it gets all the blame (r+1)
    Individual_Fitness_Values[(Alpha[, r + 0] != All_Alleles_Name[q] &
                                 Alpha[, r + 1] == All_Alleles_Name[q]), r + 1] <-
      Individual_Fitness_Values[(Alpha[, r + 0] != All_Alleles_Name[q] &
                                   Alpha[, r + 1] == All_Alleles_Name[q]), r + 1] * All_ArBrCr_add[q] * All_ArBrCr[q]
    
    # If only the first matches, it gets all the blame (r)
    Individual_SexRatio_Values[(Alpha[, r + 0] == All_Alleles_Name[q] &
                                 Alpha[, r + 1] != All_Alleles_Name[q]), r] <-
      Individual_SexRatio_Values[(Alpha[, r + 0] == All_Alleles_Name[q] &
                                   Alpha[, r + 1] != All_Alleles_Name[q]), r] * All_AsrBsrCsr_add[q] 
    
    # If only the second matches, it gets all the blame (r+1)
    Individual_SexRatio_Values[(Alpha[, r + 0] != All_Alleles_Name[q] &
                                 Alpha[, r + 1] == All_Alleles_Name[q]), r + 1] <-
      Individual_SexRatio_Values[(Alpha[, r + 0] != All_Alleles_Name[q] &
                                   Alpha[, r + 1] == All_Alleles_Name[q]), r + 1] * All_AsrBsrCsr_add[q] 
    
    # For each phenotype, update Hap_Suff
    for(i in vector_complex_fitness) { 
      Hap_Suff[[i]][Alpha[, r + 0] == All_Alleles_Name[q], r + 0] <-
        Hap_Suff[[i]][Alpha[, r + 0] == All_Alleles_Name[q], r + 0] + All_Homo_Fit_Type[[i]][q]
      
      Hap_Suff[[i]][Alpha[, r + 1] == All_Alleles_Name[q], r + 0] <-
        Hap_Suff[[i]][Alpha[, r + 1] == All_Alleles_Name[q], r + 0] + All_Homo_Fit_Type[[i]][q]
    }
    
    if (q == length(c(A, B))) {
      r <- 5
    } else if (q == length(A)) {
      r <- 3
    }
  }
  
  
  # for (r in c(1, 3, 5)) {
    # For each phenotype
    for(i in vector_complex_fitness) {
      
      # For the genotypes where no haplosufficent allele is present we apply the recessive fitness costs
      Individual_SexRatio_Values[(rowSums(Hap_Suff[[i]]) < 1),colSums(Hap_Suff[[i]])>0] <-
        Individual_SexRatio_Values[(rowSums(Hap_Suff[[i]]) < 1),colSums(Hap_Suff[[i]])>0]* 
        All_Homo_SexRatio[[i]][1]^(1/sum(colSums(Hap_Suff[[i]])>0))
      
      Individual_Fitness_Values[(rowSums(Hap_Suff[[i]]) < 1),colSums(Hap_Suff[[i]])>0] <-
        Individual_Fitness_Values[(rowSums(Hap_Suff[[i]]) < 1),colSums(Hap_Suff[[i]])>0]* 
        All_Homo_Fit_Cost[[i]][1]^(1/sum(colSums(Hap_Suff[[i]])>0))
      
      Individual_Fitness_Values[(rowSums(Hap_Suff[[i]]) < 2),colSums(Hap_Suff[[i]])>0] <-
        Individual_Fitness_Values[(rowSums(Hap_Suff[[i]]) < 2),colSums(Hap_Suff[[i]])>0]* 
        All_Carrier_Fit_Cost[[i]][1]^(1/sum(colSums(Hap_Suff[[i]])>0))
    }
    
    # For all fitness costs, not just the haplosufficent ones:
    # If individuals are homozygous, we move all the fitness cost to one allele, and set the other to 1
    Individual_Fitness_Values[(Alpha[, r + 0] == Alpha[, r + 1]), r + 0] <-
      Individual_Fitness_Values[(Alpha[, r + 0] == Alpha[, r + 1]), r + 0] *
      Individual_Fitness_Values[(Alpha[, r + 0] == Alpha[, r + 1]), r + 1]
    
    Individual_Fitness_Values[(Alpha[, r + 0] == Alpha[, r + 1]), r + 1] <- 1
    
    Individual_SexRatio_Values[(Alpha[, r + 0] == Alpha[, r + 1]), r + 0] <-
      Individual_SexRatio_Values[(Alpha[, r + 0] == Alpha[, r + 1]), r + 0] *
      Individual_SexRatio_Values[(Alpha[, r + 0] == Alpha[, r + 1]), r + 1]
    
    Individual_SexRatio_Values[(Alpha[, r + 0] == Alpha[, r + 1]), r + 1] <- 1
  
  return(list(Fitness_Costs = Individual_Fitness_Values,
         Sex_Ratios = Individual_SexRatio_Values))
}
##  ............................................................................
##  Simulations                                                             ####
# The generation by generation genotype frequency simulation
Main_Model <- function(Model_Matrices, replicate) {
  
  
  Genotypes <- 1:length(Model_Matrices$Alpha_M$Freq)

  # Validate migration source dimensions if present
  if (!is.null(Model_Matrices$Migration_Source_Data_M)) {
    if (nrow(Model_Matrices$Migration_Source_Data_M) != length(Genotypes)) {
      stop("Migration source has ", nrow(Model_Matrices$Migration_Source_Data_M),
           " genotypes but model has ", length(Genotypes))
    }
  }

  # Pre-initializing the data-frames we will use
  Germline_M <- Genotypes
  Germline_F <- Genotypes
  
  Embry_DepOnly <- Genotypes
  Embryo <- Genotypes
  
  Zyg_Fit_M <- Genotypes
  Zyg_Fit_F <- Genotypes
  
  Zyg_DepFit_M <- Genotypes
  Zyg_DepFit_F <- Genotypes
  
  Haploid_M <- Genotypes
  Haploid_F <- Genotypes
  
  Fit_Male <- apply(Model_Matrices$Fit_Male, 1, prod)
  Fit_Female <- apply(Model_Matrices$Fit_Female, 1, prod)
  
  SRatio_Male <- apply(Model_Matrices$SRatio_Male, 1, prod)
  SRatio_Female <- apply(Model_Matrices$SRatio_Female, 1, prod)
  
  starting_pop_M <- as.numeric(Model_Matrices$Alpha_M$Freq)
  starting_pop_F <- as.numeric(Model_Matrices$Alpha_F$Freq)
  
  suppressMessages(# this is to suppress the 'New names: warning spam'
    Outcome <- bind_cols(Model_Matrices$Alpha_M[, -13],
                         as.data.frame(
                           matrix(
                             data = NA,
                             ncol = Model_Matrices$Generations,
                             nrow = length(Model_Matrices$Alpha_M[, 1])
                           )
                         )))
  
  
  # These are the lists we use to record the outcomes
  Outcomes_M_PostFit <- Outcome
  Outcomes_F_PostFit <- Outcome
  Outcomes_MF_PostFit <- Outcome
  
  Outcomes_M_PreFit <- Outcome
  Outcomes_F_PreFit <- Outcome
  Outcomes_MF_PreFit <- Outcome
  
  Single_Outcomes <- as.data.frame(matrix(
    data = NA,
    ncol = Model_Matrices$Generations,
    nrow = 5)
  )
  
  rownames(Single_Outcomes) <- 
    c("Sex_Ratio_PreFit",
      "Sex_Ratio_PostFit",
      "Genetic_Load_M",
      "Genetic_Load_F",
      "Genetic_Load_MF"
      )
  
  Single_Outcomes$Sex <- c("Male&Female","Male&Female","Male","Female","Male&Female")
  Single_Outcomes$Pre_or_Post <- c("Pre","Post","Post","Post","Post")
  Single_Outcomes$Locus <- "General"
  Single_Outcomes$Outcome <- c("Sex_Ratio","Sex_Ratio","Genetic_Load","Genetic_Load","Genetic_Load")
  
  # cbind(Model_Matrices$Alpha_M,Fit_Male, (rowSums(t(t(Model_Matrices$D) * Fit_Male))), Fit_Female,(rowSums(t(t(Model_Matrices$D) * Fit_Female))))
  # cbind(Model_Matrices$Alpha_M,as.data.frame(as.matrix(Model_Matrices$D)))
  
  # Here we combine the deposition/somatic matrix and fitness cost vector
  D_hat_M <-
   rowSums(t(t(Model_Matrices$D %*% Model_Matrices$exp_soma_conv_male) * Fit_Male))
  # cbind(Model_Matrices$Alpha_M,as.numeric(D_hat_M))
  D_hat_F <-
    rowSums(t(t(Model_Matrices$D %*% Model_Matrices$exp_soma_conv_female) * Fit_Female))
  # cbind(Model_Matrices$Alpha_F,as.numeric(D_hat_F))
  
  D_no_Fit <- Model_Matrices$D
  
  # this is the offset for where the freq data should start
  # the first columns describe the genotype
  off_set_df <- 12
  
  Release_File <- Model_Matrices$Release_File
  
  Haploid_M <- as.numeric(as.numeric(Model_Matrices$Alpha_M[Model_Matrices$Alpha_M$SD == 0,]$Freq) %*% Model_Matrices$H_M)
  Haploid_F <- as.numeric(as.numeric(Model_Matrices$Alpha_F[Model_Matrices$Alpha_M$SD == 0,]$Freq) %*% Model_Matrices$H_F)

  # Guard: verify haploid gamete frequencies sum to 1
  if (abs(sum(Haploid_M) - 1) > 1e-8) warning("Haploid_M frequencies do not sum to 1: ", sum(Haploid_M))
  if (abs(sum(Haploid_F) - 1) > 1e-8) warning("Haploid_F frequencies do not sum to 1: ", sum(Haploid_F))

  Haploid_M_SD <- rep(0, length(Haploid_M))
  Haploid_F_SD <- rep(0, length(Haploid_F))
  
  debug <- FALSE
  
  for (t in 1:(Model_Matrices$Generations)) {
    
    
    if (debug) {
      print(cbind(colnames(Model_Matrices$H_M),Haploid_M,Haploid_M_SD,Haploid_F,Haploid_F_SD))
      print(paste0("Haploid_M:",sum(Haploid_M),"Haploid_M_SD:",sum(Haploid_M_SD),"Haploid_F:",sum(Haploid_F),"Haploid_F_SD:",sum(Haploid_F_SD)))
    }
    
    # Here we set the new value for alpha, which applies to both sexes
    Embryo_M <- as.vector(as.numeric(Haploid_M %*% t(as.numeric(Haploid_F)))) # Xy - XX -> M_50/F_50
    Embryo_F <- as.vector(as.numeric(Haploid_M %*% t(as.numeric(Haploid_F)))) # Xy - XX -> M_50/F_50
    
    if (sum(Model_Matrices$Alpha_M$SD) > 0) {
      
      # Sex converted
      Embryo_M <- Embryo_M + 0.5 * as.vector(as.numeric(Haploid_M %*% t(as.numeric(Haploid_F_SD)))) + # Xy - Xy -> M_50/F_25/Dead_25
      0.5 * as.vector(as.numeric(Haploid_M_SD %*% t(as.numeric(Haploid_F_SD))))  # XX - Xy -> M_50/F_50
      
      Embryo_F <- Embryo_F + 0.25 * as.vector(as.numeric(Haploid_M %*% t(as.numeric(Haploid_F_SD)))) + # Xy - Xy -> M_50/F_25/Dead_25
        as.vector(as.numeric(Haploid_M_SD %*% t(as.numeric(Haploid_F)))) + # XX - XX - > F_100 
        0.5 * as.vector(as.numeric(Haploid_M_SD %*% t(as.numeric(Haploid_F_SD))))  # XX - Xy -> M_50/F_50
      
      Embryo_M <- c(Embryo_M, rep(0, length(Embryo_M)))
      Embryo_F <- c(Embryo_F, rep(0, length(Embryo_F)))
    }
    
    if (sum(Embryo_M) > 0) {Embryo_M <- Embryo_M/sum(Embryo_M)}
    if (sum(Embryo_F) > 0) {Embryo_F <- Embryo_F/sum(Embryo_F)}
    
    # Here we are doing code to do releases at the adult stage after fitness costs
    # In most cases this will reflect what people do
    if (t %in% as.numeric(Release_File$Generation[Release_File$Release_PreorPost_Fitness_Costs == "Pre"])) {
      
      for (r in 1:length(Release_File$Release_Type[Release_File$Release_PreorPost_Fitness_Costs == "Pre"])) {
        
        if (Release_File$Release_PreorPost_Fitness_Costs[Release_File$Release_PreorPost_Fitness_Costs == "Pre"][r] == "Pre" &
         t == Release_File$Generation[Release_File$Release_PreorPost_Fitness_Costs == "Pre"][r]) {
          
          
        Embryo_M <-
          Default_Geno_Freq(
            Freq = Embryo_M,
            Alpha = Model_Matrices$Alpha_M,
            A = Model_Matrices$geno_list$A,
            B = Model_Matrices$geno_list$B,
            C = Model_Matrices$geno_list$C,
            type = Release_File$Release_Type[r],
            value = as.numeric(Release_File$Male_Release_Fraction[r])
          )
        
        Embryo_F <-
          Default_Geno_Freq(
            Freq = Embryo_F,
            Alpha = Model_Matrices$Alpha_F,
            A = Model_Matrices$geno_list$A,
            B = Model_Matrices$geno_list$B,
            C = Model_Matrices$geno_list$C,
            type = Release_File$Release_Type[r],
            value = as.numeric(Release_File$Female_Release_Fraction[r])
          )
        }
      }
    }
    
    
    if (sum(Model_Matrices$Alpha_M$SD) > 0) {
      
      
      if (sum(Embryo_M[as.logical(Model_Matrices$Alpha_M$SD)] + Embryo_F[as.logical(Model_Matrices$Alpha_M$SD)]) > 0) {warning("sexconverts appear before step")}
      
      
      # To account for sex conversion, we have to keep track of genetically male and female individuals
      # We process the phenotypic females, and then keep track of the ratio of 
      Embryo_M_tmp <- c(Embryo_M[!as.logical(Model_Matrices$Alpha_M$SD)] * SRatio_Male[!as.logical(Model_Matrices$Alpha_M$SD)],
                    Embryo_M[as.logical(Model_Matrices$Alpha_M$SD)] + Embryo_F[!as.logical(Model_Matrices$Alpha_M$SD)] * (1-SRatio_Female[!as.logical(Model_Matrices$Alpha_M$SD)]))
      
      Embryo_F_tmp <- c(Embryo_F[!as.logical(Model_Matrices$Alpha_M$SD)] * SRatio_Female[!as.logical(Model_Matrices$Alpha_M$SD)],
                    Embryo_F[as.logical(Model_Matrices$Alpha_M$SD)] + Embryo_M[!as.logical(Model_Matrices$Alpha_M$SD)] * (1-SRatio_Male[!as.logical(Model_Matrices$Alpha_M$SD)]))
    
      Embryo_M <- Embryo_M_tmp
      Embryo_F <- Embryo_F_tmp
      
      if (debug) {
      print(cbind(Outcomes_F_PreFit[,1:11],Embryo_M_tmp)[Embryo_M_tmp>0,])
      print(cbind(Outcomes_F_PreFit[,1:11],Embryo_F_tmp)[Embryo_F_tmp>0,])
      }
      }
    
    
    Embryo_M_Geno_Phenotype_Ratio <- Embryo_M
    Embryo_F_Geno_Phenotype_Ratio <- Embryo_F

    # Guard: check for population collapse before dividing by embryo sums
    if (sum(Embryo_M) == 0 | sum(Embryo_F) == 0) {
      warning("Population collapsed at generation ", t, " (embryo frequencies sum to zero)")
      break
    }

    Outcomes_F_PreFit[, off_set_df + t] <- (Embryo_F / sum(Embryo_F))
    Outcomes_M_PreFit[, off_set_df + t] <- (Embryo_M / sum(Embryo_M))
    Outcomes_MF_PreFit[, off_set_df + t] <- (Embryo_F + Embryo_M) / (sum(Embryo_F) + sum(Embryo_M))

    # This is the embryonic sex ratio. This is 0.5 unless sex conversion occurs
    Single_Outcomes[c("Sex_Ratio_PreFit"), t] <- sum(Embryo_M) / sum(Embryo_F + Embryo_M)
    Zyg_Fit_M <- as.numeric(Embryo_M * D_hat_M) #Embryo_M * (rowSums(t(t(Model_Matrices$D) * Fit_Male)))
    Zyg_Fit_F <- as.numeric(Embryo_F * D_hat_F) #Embryo_F * (rowSums(t(t(Model_Matrices$D) * Fit_Female)))


    Single_Outcomes[c("Genetic_Load_M"), t] <- sum(Zyg_Fit_M) / sum(Embryo_M)
    Single_Outcomes[c("Genetic_Load_F"), t] <- sum(Zyg_Fit_F) / sum(Embryo_F)
    Single_Outcomes[c("Genetic_Load_MF"), t] <- sum(Zyg_Fit_M + Zyg_Fit_F) / sum(Embryo_M + Embryo_F)

    if (sum(Zyg_Fit_M) == 0 | sum(Zyg_Fit_F) == 0) {warning("pop died out")}
    if (sum(Zyg_Fit_M) == 0 | sum(Zyg_Fit_F) == 0) break
    
    # If this is replicate 0, we do a deterministic simulation
    # If it is not the 0th simulation, we bottleneck the frequencies by the population size
    if (replicate != 0) {
      male_pop <- rbinom(1, as.numeric(Release_File$Adult_Pop_Size[1]), (sum(Zyg_Fit_M) /
                                                                           (sum(Zyg_Fit_M) + sum(Zyg_Fit_F))))
      # Here we make a random selection of the different genotypes based on the populations size
      Zyg_Fit_M <- rmultinom(n = 1, size = male_pop, prob = Zyg_Fit_M)[, 1]
      Zyg_Fit_F <- rmultinom(n = 1, size = as.numeric(Release_File$Adult_Pop_Size[1]) - male_pop, prob = Zyg_Fit_F
      )[, 1]
      
    }
    
    # if (sum(Zyg_Fit_M) == 0 | sum(Zyg_Fit_F) == 0) break
    
    
    # Here we are doing code to do releases at the adult stage after fitness costs
    # In most cases this will reflect what people do
    if (t %in% as.numeric(Release_File$Generation[Release_File$Release_PreorPost_Fitness_Costs == "Post"])) {
      
      for (r in 1:length(Release_File$Release_Type)) {
        
        if (t == Release_File$Generation[r] &
            "Post"== Release_File$Release_PreorPost_Fitness_Costs[r]) {
        
        Zyg_Fit_M <-
          Default_Geno_Freq(
            Freq = Zyg_Fit_M,
            Alpha = Model_Matrices$Alpha_M,
            A = Model_Matrices$geno_list$A,
            B = Model_Matrices$geno_list$B,
            C = Model_Matrices$geno_list$C,
            type = Release_File$Release_Type[r],
            value = as.numeric(Release_File$Male_Release_Fraction[r])
          )
        
        Zyg_Fit_F <-
          Default_Geno_Freq(
            Freq = Zyg_Fit_F,
            Alpha = Model_Matrices$Alpha_F,
            A = Model_Matrices$geno_list$A,
            B = Model_Matrices$geno_list$B,
            C = Model_Matrices$geno_list$C,
            type = Release_File$Release_Type[r],
            value = as.numeric(Release_File$Female_Release_Fraction[r])
          )
        }
      }
      
      if (replicate != 0) {
        
        # # Here we add an extra bottleneck to account for the variability specific to a release
        # Zyg_Fit_M <- rmultinom(
        #   n = 1,
        #   size = as.numeric(Release_File$Release_Pop_Effective_Size[r]),
        #   prob = Zyg_Fit_M
        # ) / as.numeric(Release_File$Release_Pop_Effective_Size[r])
        # 
        # Zyg_Fit_F <- rmultinom(
        #   n = 1,
        #   size = as.numeric(Release_File$Release_Pop_Effective_Size[r]),
        #   prob = Zyg_Fit_F
        # ) / as.numeric(Release_File$Release_Pop_Effective_Size[r])
      }
    }
    
    ###################
    
    # This is the adult sex ratio. This is includes sex-conversion and fitness costs
    Single_Outcomes[c("Sex_Ratio_PostFit"), t] <- sum(Zyg_Fit_M) / sum(Zyg_Fit_F + Zyg_Fit_M)

    # One-way migration: blend in source population post-fitness frequencies
    if (!is.null(Model_Matrices$Migration_Source_Data_M) &&
        Model_Matrices$Migration_Rate > 0) {

      if (t <= ncol(Model_Matrices$Migration_Source_Data_M)) {
        m <- Model_Matrices$Migration_Rate
        source_freq_M <- Model_Matrices$Migration_Source_Data_M[, t]
        source_freq_F <- Model_Matrices$Migration_Source_Data_F[, t]

        # Blend: (1-m) of current post-fitness adults + m of source post-fitness adults
        Zyg_Fit_M <- (1 - m) * Zyg_Fit_M + m * source_freq_M * sum(Zyg_Fit_M)
        Zyg_Fit_F <- (1 - m) * Zyg_Fit_F + m * source_freq_F * sum(Zyg_Fit_F)
      } else {
        message("Migration source exhausted at generation ", t, " - no migration applied")
      }
    }

    ###################
    # Here we save the outcomes
    Outcomes_F_PostFit[, off_set_df + t] <- (Zyg_Fit_F / sum(Zyg_Fit_F))
    Outcomes_M_PostFit[, off_set_df + t] <- Zyg_Fit_M / sum(Zyg_Fit_M)
    Outcomes_MF_PostFit[, off_set_df + t] <- (Zyg_Fit_F + Zyg_Fit_M) / sum(Zyg_Fit_F + Zyg_Fit_M)
    
    # Here we normalize the losses. We only track the relative frequency, not absolute
    Zyg_DepFit_M <- as.numeric(t(Zyg_Fit_M) %*% D_no_Fit)
    Zyg_DepFit_F <- as.numeric(t(Zyg_Fit_F) %*% D_no_Fit)
    
    # Applying the expression based conversion
    Germline_M <-
      as.numeric(t(Zyg_DepFit_M) %*% Model_Matrices$exp_conv_male)
    Germline_F <-
      as.numeric(t(Zyg_DepFit_F) %*% Model_Matrices$exp_conv_female)

    
    if (debug) {
      print(cbind(Outcomes_F_PreFit[,1:11],Germline_M)[Germline_M>0,])
      print(sum(Germline_M))
      print(cbind(Outcomes_F_PreFit[,1:11],Germline_F)[Germline_F>0,])
      print(sum(Germline_F))
    }
    
    
    #####################################################
    # Changes Gametes
    
    # Non-sex converted
    Haploid_M <- as.numeric(Germline_M[!as.logical(Model_Matrices$Alpha_M$SD)] %*% Model_Matrices$H_M)
    Haploid_F <- as.numeric(Germline_F[!as.logical(Model_Matrices$Alpha_M$SD)] %*% Model_Matrices$H_F)
    
    if (sum(Model_Matrices$Alpha_M$SD) > 0) {
      # Sex converted
      Haploid_M_SD <- as.numeric(Germline_M[as.logical(Model_Matrices$Alpha_M$SD)] %*% Model_Matrices$H_M)
      Haploid_F_SD <- as.numeric(Germline_F[as.logical(Model_Matrices$Alpha_M$SD)] %*% Model_Matrices$H_F) 
    }
  }
  
  return(
    list(
      Outcomes_F_PreFit = Outcomes_F_PreFit,
      Outcomes_M_PreFit = Outcomes_M_PreFit,
      Outcomes_MF_PreFit = Outcomes_MF_PreFit,
      Outcomes_F_PostFit = Outcomes_F_PostFit,
      Outcomes_M_PostFit = Outcomes_M_PostFit,
      Outcomes_MF_PostFit = Outcomes_MF_PostFit,
      # Fit_Loss = Fit_Loss,
      Single_Outcomes = Single_Outcomes
      # Fit_Deposition = Fit_Deposition,
      # Conversion_Change = Conversion_Change
    )
  )
}

#   ____________________________________________________________________________
#   Connected populations (chain) with back-and-forth migration            ####
#
# This replaces the earlier sequential, one-way two-population approach. All
# populations are advanced one generation in lockstep and migrants are then
# exchanged at the gamete (haploid) stage, exactly as in the two-population
# model of Verkuijl et al. 2022 (PLoS Genetics; OSF project 6a2eg), generalised
# to n populations arranged in a chain.
#
# Topology (chain of n populations), migration rate m each generation:
#   Pop 1 (release end):  H1 <- (1-m) H1 + m H2
#   Pop p (1 < p < n):    Hp <- (1-m) Hp + (m/2) H(p-1) + (m/2) H(p+1)
#   Pop n (far end):      Hn <- (1-m) Hn + (m/2) H(n-1) + (m/2) WT
# The far end is supplemented with wild-type gametes (the "general population"),
# preserving the wild-type supplementation of the 2022 model. Each population
# receives a total immigrant fraction m, so gamete frequencies remain
# normalised (every source vector itself sums to 1).
#
# Because all populations share the same drive, the conversion / fitness
# matrices are built once (Model_Matrices). Populations differ only in their
# Release_File: the release population (Pop 1) carries the release, the others
# are wild-type with zero release.
#
# Arguments
#   Model_Matrices : output of All_Main_Calculations() (shared by all pops).
#   release_files  : list of per-population Release_File data frames; its length
#                    sets the number of populations.
#   migration_rate : m above (0 reduces this to n independent populations).
#   wt_supplementation : if TRUE the far end receives m/2 wild-type gametes; if
#                    FALSE the far end receives its full immigrant fraction m
#                    from its single neighbour.
#   replicate      : 0 = deterministic; != 0 = per-population stochastic
#                    bottleneck (matching Main_Model).
#
# Returns a list (pop_1, pop_2, ...) where each element has the same shape as
# a Main_Model() result, so Outcomes_by_Allele() can consume each directly.
Main_Model_connected <- function(Model_Matrices,
                                 release_files,
                                 migration_rate,
                                 wt_supplementation = TRUE,
                                 replicate = 0) {

  n_pop       <- length(release_files)
  Genotypes   <- 1:length(Model_Matrices$Alpha_M$Freq)
  Generations <- Model_Matrices$Generations
  off_set_df  <- 12

  # ---- Shared derived quantities (identical to Main_Model setup) ----
  Fit_Male      <- apply(Model_Matrices$Fit_Male, 1, prod)
  Fit_Female    <- apply(Model_Matrices$Fit_Female, 1, prod)
  SRatio_Male   <- apply(Model_Matrices$SRatio_Male, 1, prod)
  SRatio_Female <- apply(Model_Matrices$SRatio_Female, 1, prod)

  D_hat_M <- rowSums(t(t(Model_Matrices$D %*% Model_Matrices$exp_soma_conv_male)   * Fit_Male))
  D_hat_F <- rowSums(t(t(Model_Matrices$D %*% Model_Matrices$exp_soma_conv_female) * Fit_Female))
  D_no_Fit <- Model_Matrices$D

  SD_present <- sum(Model_Matrices$Alpha_M$SD) > 0
  SD_logical <- as.logical(Model_Matrices$Alpha_M$SD)

  # ---- Wild-type gamete pools (used as the migration supplementation source) --
  Haploid_M_WT <- as.numeric(as.numeric(Model_Matrices$Alpha_M[Model_Matrices$Alpha_M$SD == 0, ]$Freq) %*% Model_Matrices$H_M)
  Haploid_F_WT <- as.numeric(as.numeric(Model_Matrices$Alpha_F[Model_Matrices$Alpha_M$SD == 0, ]$Freq) %*% Model_Matrices$H_F)
  Haploid_SD_WT <- rep(0, length(Haploid_M_WT))

  # ---- Per-population recorder templates ----
  suppressMessages(
    Outcome <- bind_cols(Model_Matrices$Alpha_M[, -13],
                         as.data.frame(matrix(
                           data = NA,
                           ncol = Generations,
                           nrow = length(Model_Matrices$Alpha_M[, 1])
                         )))
  )
  Single_template <- as.data.frame(matrix(data = NA, ncol = Generations, nrow = 5))
  rownames(Single_template) <- c("Sex_Ratio_PreFit", "Sex_Ratio_PostFit",
                                 "Genetic_Load_M", "Genetic_Load_F", "Genetic_Load_MF")
  Single_template$Sex         <- c("Male&Female", "Male&Female", "Male", "Female", "Male&Female")
  Single_template$Pre_or_Post <- c("Pre", "Post", "Post", "Post", "Post")
  Single_template$Locus       <- "General"
  Single_template$Outcome     <- c("Sex_Ratio", "Sex_Ratio", "Genetic_Load", "Genetic_Load", "Genetic_Load")

  # ---- Per-population state + recorders ----
  Hm <- Hf <- Hm_sd <- Hf_sd <- vector("list", n_pop)
  rec <- vector("list", n_pop)
  for (p in 1:n_pop) {
    Hm[[p]]    <- Haploid_M_WT
    Hf[[p]]    <- Haploid_F_WT
    Hm_sd[[p]] <- Haploid_SD_WT
    Hf_sd[[p]] <- Haploid_SD_WT
    rec[[p]] <- list(
      F_PreFit  = Outcome, M_PreFit  = Outcome, MF_PreFit  = Outcome,
      F_PostFit = Outcome, M_PostFit = Outcome, MF_PostFit = Outcome,
      Single = Single_template,
      dead = FALSE
    )
  }

  # ---- One generation for a single population ----------------------------------
  # Faithful port of the Main_Model generation body, minus the old one-way
  # migration block. Returns the new (pre-migration) gametes plus the columns to
  # record. `dead = TRUE` signals collapse so the caller can freeze that pop.
  step_pop <- function(Haploid_M, Haploid_F, Haploid_M_SD, Haploid_F_SD, Release_File, t) {

    # New embryo frequencies from the gamete pools (applies to both sexes)
    Embryo_M <- as.vector(as.numeric(Haploid_M %*% t(as.numeric(Haploid_F))))
    Embryo_F <- as.vector(as.numeric(Haploid_M %*% t(as.numeric(Haploid_F))))

    if (SD_present) {
      Embryo_M <- Embryo_M + 0.5 * as.vector(as.numeric(Haploid_M %*% t(as.numeric(Haploid_F_SD)))) +
        0.5 * as.vector(as.numeric(Haploid_M_SD %*% t(as.numeric(Haploid_F_SD))))

      Embryo_F <- Embryo_F + 0.25 * as.vector(as.numeric(Haploid_M %*% t(as.numeric(Haploid_F_SD)))) +
        as.vector(as.numeric(Haploid_M_SD %*% t(as.numeric(Haploid_F)))) +
        0.5 * as.vector(as.numeric(Haploid_M_SD %*% t(as.numeric(Haploid_F_SD))))

      Embryo_M <- c(Embryo_M, rep(0, length(Embryo_M)))
      Embryo_F <- c(Embryo_F, rep(0, length(Embryo_F)))
    }

    if (sum(Embryo_M) > 0) {Embryo_M <- Embryo_M / sum(Embryo_M)}
    if (sum(Embryo_F) > 0) {Embryo_F <- Embryo_F / sum(Embryo_F)}

    # Releases applied before fitness costs
    if (t %in% as.numeric(Release_File$Generation[Release_File$Release_PreorPost_Fitness_Costs == "Pre"])) {
      for (r in 1:length(Release_File$Release_Type[Release_File$Release_PreorPost_Fitness_Costs == "Pre"])) {
        if (Release_File$Release_PreorPost_Fitness_Costs[Release_File$Release_PreorPost_Fitness_Costs == "Pre"][r] == "Pre" &
            t == Release_File$Generation[Release_File$Release_PreorPost_Fitness_Costs == "Pre"][r]) {
          Embryo_M <- Default_Geno_Freq(
            Freq = Embryo_M, Alpha = Model_Matrices$Alpha_M,
            A = Model_Matrices$geno_list$A, B = Model_Matrices$geno_list$B, C = Model_Matrices$geno_list$C,
            type = Release_File$Release_Type[r], value = as.numeric(Release_File$Male_Release_Fraction[r]))
          Embryo_F <- Default_Geno_Freq(
            Freq = Embryo_F, Alpha = Model_Matrices$Alpha_F,
            A = Model_Matrices$geno_list$A, B = Model_Matrices$geno_list$B, C = Model_Matrices$geno_list$C,
            type = Release_File$Release_Type[r], value = as.numeric(Release_File$Female_Release_Fraction[r]))
        }
      }
    }

    # Sex conversion bookkeeping (genetic vs phenotypic sex)
    if (SD_present) {
      if (sum(Embryo_M[SD_logical] + Embryo_F[SD_logical]) > 0) {warning("sexconverts appear before step")}

      Embryo_M_tmp <- c(Embryo_M[!SD_logical] * SRatio_Male[!SD_logical],
                        Embryo_M[SD_logical] + Embryo_F[!SD_logical] * (1 - SRatio_Female[!SD_logical]))
      Embryo_F_tmp <- c(Embryo_F[!SD_logical] * SRatio_Female[!SD_logical],
                        Embryo_F[SD_logical] + Embryo_M[!SD_logical] * (1 - SRatio_Male[!SD_logical]))
      Embryo_M <- Embryo_M_tmp
      Embryo_F <- Embryo_F_tmp
    }

    # Guard: population collapse before dividing by embryo sums
    if (sum(Embryo_M) == 0 | sum(Embryo_F) == 0) {
      warning("Population collapsed at generation ", t, " (embryo frequencies sum to zero)")
      return(list(dead = TRUE))
    }

    F_PreFit  <- Embryo_F / sum(Embryo_F)
    M_PreFit  <- Embryo_M / sum(Embryo_M)
    MF_PreFit <- (Embryo_F + Embryo_M) / (sum(Embryo_F) + sum(Embryo_M))

    sr_pre <- sum(Embryo_M) / sum(Embryo_F + Embryo_M)

    Zyg_Fit_M <- as.numeric(Embryo_M * D_hat_M)
    Zyg_Fit_F <- as.numeric(Embryo_F * D_hat_F)

    gl_m  <- sum(Zyg_Fit_M) / sum(Embryo_M)
    gl_f  <- sum(Zyg_Fit_F) / sum(Embryo_F)
    gl_mf <- sum(Zyg_Fit_M + Zyg_Fit_F) / sum(Embryo_M + Embryo_F)

    if (sum(Zyg_Fit_M) == 0 | sum(Zyg_Fit_F) == 0) {warning("pop died out")}
    if (sum(Zyg_Fit_M) == 0 | sum(Zyg_Fit_F) == 0) return(list(dead = TRUE))

    # Stochastic bottleneck by population size (skipped for replicate 0)
    if (replicate != 0) {
      male_pop <- rbinom(1, as.numeric(Release_File$Adult_Pop_Size[1]),
                         (sum(Zyg_Fit_M) / (sum(Zyg_Fit_M) + sum(Zyg_Fit_F))))
      Zyg_Fit_M <- rmultinom(n = 1, size = male_pop, prob = Zyg_Fit_M)[, 1]
      Zyg_Fit_F <- rmultinom(n = 1, size = as.numeric(Release_File$Adult_Pop_Size[1]) - male_pop, prob = Zyg_Fit_F)[, 1]
    }

    # Releases applied after fitness costs
    if (t %in% as.numeric(Release_File$Generation[Release_File$Release_PreorPost_Fitness_Costs == "Post"])) {
      for (r in 1:length(Release_File$Release_Type)) {
        if (t == Release_File$Generation[r] & "Post" == Release_File$Release_PreorPost_Fitness_Costs[r]) {
          Zyg_Fit_M <- Default_Geno_Freq(
            Freq = Zyg_Fit_M, Alpha = Model_Matrices$Alpha_M,
            A = Model_Matrices$geno_list$A, B = Model_Matrices$geno_list$B, C = Model_Matrices$geno_list$C,
            type = Release_File$Release_Type[r], value = as.numeric(Release_File$Male_Release_Fraction[r]))
          Zyg_Fit_F <- Default_Geno_Freq(
            Freq = Zyg_Fit_F, Alpha = Model_Matrices$Alpha_F,
            A = Model_Matrices$geno_list$A, B = Model_Matrices$geno_list$B, C = Model_Matrices$geno_list$C,
            type = Release_File$Release_Type[r], value = as.numeric(Release_File$Female_Release_Fraction[r]))
        }
      }
    }

    sr_post <- sum(Zyg_Fit_M) / sum(Zyg_Fit_F + Zyg_Fit_M)

    # NOTE: the old one-way migration block lived here; migration is now handled
    # between populations in the caller, at the gamete stage below.

    F_PostFit  <- Zyg_Fit_F / sum(Zyg_Fit_F)
    M_PostFit  <- Zyg_Fit_M / sum(Zyg_Fit_M)
    MF_PostFit <- (Zyg_Fit_F + Zyg_Fit_M) / sum(Zyg_Fit_F + Zyg_Fit_M)

    # Deposition, then expression-based germline conversion
    Zyg_DepFit_M <- as.numeric(t(Zyg_Fit_M) %*% D_no_Fit)
    Zyg_DepFit_F <- as.numeric(t(Zyg_Fit_F) %*% D_no_Fit)

    Germline_M <- as.numeric(t(Zyg_DepFit_M) %*% Model_Matrices$exp_conv_male)
    Germline_F <- as.numeric(t(Zyg_DepFit_F) %*% Model_Matrices$exp_conv_female)

    # New gamete pools (carried to the next generation, after migration)
    new_Haploid_M <- as.numeric(Germline_M[!SD_logical] %*% Model_Matrices$H_M)
    new_Haploid_F <- as.numeric(Germline_F[!SD_logical] %*% Model_Matrices$H_F)
    new_Haploid_M_SD <- Haploid_M_SD
    new_Haploid_F_SD <- Haploid_F_SD
    if (SD_present) {
      new_Haploid_M_SD <- as.numeric(Germline_M[SD_logical] %*% Model_Matrices$H_M)
      new_Haploid_F_SD <- as.numeric(Germline_F[SD_logical] %*% Model_Matrices$H_F)
    }

    list(
      dead = FALSE,
      Haploid_M = new_Haploid_M, Haploid_F = new_Haploid_F,
      Haploid_M_SD = new_Haploid_M_SD, Haploid_F_SD = new_Haploid_F_SD,
      F_PreFit = F_PreFit, M_PreFit = M_PreFit, MF_PreFit = MF_PreFit,
      F_PostFit = F_PostFit, M_PostFit = M_PostFit, MF_PostFit = MF_PostFit,
      sr_pre = sr_pre, sr_post = sr_post, gl_m = gl_m, gl_f = gl_f, gl_mf = gl_mf
    )
  }

  # ---- Chain migration at the gamete stage (snapshot, then exchange) -----------
  mix_chain <- function(H, wt) {
    if (n_pop == 1 || migration_rate <= 0) return(H)
    m <- migration_rate
    out <- H
    for (p in 1:n_pop) {
      if (p == 1) {
        out[[p]] <- (1 - m) * H[[1]] + m * H[[2]]
      } else if (p == n_pop) {
        if (wt_supplementation) {
          out[[p]] <- (1 - m) * H[[p]] + 0.5 * m * H[[p - 1]] + 0.5 * m * wt
        } else {
          out[[p]] <- (1 - m) * H[[p]] + m * H[[p - 1]]
        }
      } else {
        out[[p]] <- (1 - m) * H[[p]] + 0.5 * m * H[[p - 1]] + 0.5 * m * H[[p + 1]]
      }
    }
    out
  }

  # ---- Lockstep generation loop -----------------------------------------------
  for (t in 1:Generations) {
    new_Hm <- new_Hf <- new_Hm_sd <- new_Hf_sd <- vector("list", n_pop)

    for (p in 1:n_pop) {
      if (rec[[p]]$dead) {
        new_Hm[[p]] <- Hm[[p]]; new_Hf[[p]] <- Hf[[p]]
        new_Hm_sd[[p]] <- Hm_sd[[p]]; new_Hf_sd[[p]] <- Hf_sd[[p]]
        next
      }

      s <- step_pop(Hm[[p]], Hf[[p]], Hm_sd[[p]], Hf_sd[[p]], release_files[[p]], t)

      if (isTRUE(s$dead)) {
        rec[[p]]$dead <- TRUE
        new_Hm[[p]] <- Hm[[p]]; new_Hf[[p]] <- Hf[[p]]
        new_Hm_sd[[p]] <- Hm_sd[[p]]; new_Hf_sd[[p]] <- Hf_sd[[p]]
        next
      }

      rec[[p]]$F_PreFit[,  off_set_df + t] <- s$F_PreFit
      rec[[p]]$M_PreFit[,  off_set_df + t] <- s$M_PreFit
      rec[[p]]$MF_PreFit[, off_set_df + t] <- s$MF_PreFit
      rec[[p]]$F_PostFit[,  off_set_df + t] <- s$F_PostFit
      rec[[p]]$M_PostFit[,  off_set_df + t] <- s$M_PostFit
      rec[[p]]$MF_PostFit[, off_set_df + t] <- s$MF_PostFit
      rec[[p]]$Single["Sex_Ratio_PreFit",  t] <- s$sr_pre
      rec[[p]]$Single["Sex_Ratio_PostFit", t] <- s$sr_post
      rec[[p]]$Single["Genetic_Load_M",  t] <- s$gl_m
      rec[[p]]$Single["Genetic_Load_F",  t] <- s$gl_f
      rec[[p]]$Single["Genetic_Load_MF", t] <- s$gl_mf

      new_Hm[[p]] <- s$Haploid_M; new_Hf[[p]] <- s$Haploid_F
      new_Hm_sd[[p]] <- s$Haploid_M_SD; new_Hf_sd[[p]] <- s$Haploid_F_SD
    }

    # Exchange migrants between populations (gamete stage)
    Hm    <- mix_chain(new_Hm,    Haploid_M_WT)
    Hf    <- mix_chain(new_Hf,    Haploid_F_WT)
    Hm_sd <- mix_chain(new_Hm_sd, Haploid_SD_WT)
    Hf_sd <- mix_chain(new_Hf_sd, Haploid_SD_WT)
  }

  # ---- Package each population like a Main_Model() result ----------------------
  results <- vector("list", n_pop)
  for (p in 1:n_pop) {
    results[[p]] <- list(
      Outcomes_F_PreFit  = rec[[p]]$F_PreFit,
      Outcomes_M_PreFit  = rec[[p]]$M_PreFit,
      Outcomes_MF_PreFit = rec[[p]]$MF_PreFit,
      Outcomes_F_PostFit  = rec[[p]]$F_PostFit,
      Outcomes_M_PostFit  = rec[[p]]$M_PostFit,
      Outcomes_MF_PostFit = rec[[p]]$MF_PostFit,
      Single_Outcomes = rec[[p]]$Single
    )
  }
  names(results) <- paste0("pop_", 1:n_pop)
  results
}

# Takes in the individual parameters and calls the conversion matrix functions
# Creates a 'Model_Matrices' list that has all the information to run a simulation
All_Main_Calculations <- function(Parameters) {
  
  Fit_Sex_Male <-
    Fitness_Costs(Parameters$Alpha, Parameters$geno_list, "M")
  Fit_Sex_Female <-
    Fitness_Costs(Parameters$Alpha, Parameters$geno_list, "F")
  
  Fit_Male <-
    Fit_Sex_Male$Fitness_Costs
  Fit_Female <-
    Fit_Sex_Female$Fitness_Costs
  
  SRatio_Male <-
    Fit_Sex_Male$Sex_Ratios
  SRatio_Female <-
    Fit_Sex_Female$Sex_Ratios
  
  D_and_Gs <-
    Conversion_D_and_G(Parameters$Alpha,
                       Parameters$nuc_rates_list,
                       Parameters$geno_list)

  
  D_conv <- unlist(D_and_Gs$dep_conv_matrix)
  
  exp_conv_male <- unlist(D_and_Gs$exp_conv_male_matrix)

  exp_conv_female <- unlist(D_and_Gs$exp_conv_female_matrix)

  exp_soma_conv_male <- unlist(D_and_Gs$exp_soma_conv_male_matrix)

  exp_soma_conv_female <- unlist(D_and_Gs$exp_soma_conv_female_matrix)
  
  H_F <-
    Conversion_H(Parameters$Alpha,
                 Parameters$geno_list,
                 Parameters$Delta_F,
                 Parameters$nuc_rates_list$F_Dep)
  H_M <-
    Conversion_H(Parameters$Alpha,
                 Parameters$geno_list,
                 Parameters$Delta_M,
                 Parameters$nuc_rates_list$M_Dep)
  
  return(
    list(
      Alpha_M = Parameters$Alpha_M,
      Alpha_F = Parameters$Alpha_F,
      Delta_M = Parameters$Delta_M,
      Delta_F = Parameters$Delta_F,
      D = D_conv,
      exp_conv_male = exp_conv_male,
      exp_conv_female = exp_conv_female,
      exp_soma_conv_male = exp_soma_conv_male,
      exp_soma_conv_female = exp_soma_conv_female,
      H_M = H_M,
      H_F = H_F,
      Generations = Parameters$Generations,
      Fit_Male = Fit_Male,
      Fit_Female = Fit_Female,
      SRatio_Male = SRatio_Male,
      SRatio_Female = SRatio_Female,
      nuc_rates_list = Parameters$nuc_rates_list,
      Release_File = Parameters$Release_File,
      geno_list = Parameters$geno_list,
      Migration_Rate = Parameters$Release_File$Migration_Rate[1],
      Migration_Source_Data_M = Parameters$Migration_Source_Data_M,
      Migration_Source_Data_F = Parameters$Migration_Source_Data_F

    )
  )
}

# Takes the high-level input selections (e.g, 'Replacement')
# Sets up the right matrices
all_default_parameters <-
  function(Allele_Fit_File_Name,
           Nuc_File_Name,
           Release_File_Name,
           project_root
           ) {
    # 1.
    ### Default Nuc rate set-up####
    # Imports the pre-set default values for cut and repair rates
    nuc_rates_list <- Default_Nuc_Rates(Nuc_File_Name,project_root)

    ### Changes to the default nuc rates go here ###
    
    # 2.
    ### Default Genotypes set-up ####
    # Imports the pre-set default values for genotype and fitness
    geno_list <- Default_Geno(Allele_Fit_File_Name,project_root)
    ### Changes to the default genotypes and fitness rates go here ###
    
    ### ###
    # Here we create a single vector that contains all the allele names
    all_alleles <- c(geno_list$A$Allele,
                     geno_list$B$Allele,
                     geno_list$C$Allele)
    # Here we check for elements that express gRNAs against each loci
    # two options: the gRNA has been deposited, or it has not
    if (any(grepl('_A', all_alleles, fixed = TRUE))) {a_grna <- c("A", "_")} else {a_grna <- "_"}
    
    if (any(grepl('_B', all_alleles, fixed = TRUE))) {b_grna <- c("B", "_")} else {b_grna <- "_"}
    
    if (any(grepl('_C', all_alleles, fixed = TRUE))) {c_grna <- c("C", "_")} else {c_grna <- "_"}
    
    # Here we generate a vector of the possible gRNAs that can be deposited
    grna_dep_options <- expand.grid(A = as.factor(a_grna[!is.na(a_grna)]),
                                    B = as.factor(b_grna[!is.na(b_grna)]),
                                    C = as.factor(c_grna[!is.na(c_grna)]))
    # Here we create a string of each row of the grna_dep_options matrix
    gRNA_Deposition <- apply(grna_dep_options, 1, paste, collapse = "")
    
    # If there is no Cas9 dep, gRNA dep does not matter either
    if ((length(nuc_rates_list$M_Dep) > 1) & nuc_rates_list$gRNA_Deposition) {grna_dep_m <- gRNA_Deposition} else {grna_dep_m <- "___"}
    if ((length(nuc_rates_list$F_Dep) > 1) & nuc_rates_list$gRNA_Deposition) {grna_dep_f <- gRNA_Deposition} else {grna_dep_f <- "___"}
    
    # Create all the possible haploid combinations of the provided alleles
    # We add a final column to both Delta and Alpha that will hold the genotype Freqs
    
    Male_no_SD <- TRUE
    Female_no_SD <- TRUE
    
    for (allele in c("A", "B", "C")) {
      
      Male_no_SD <- all(
        all(geno_list[[allele]]$Carrier_Sex_Ratio_Males == 1),
        !(!all(geno_list[[allele]]$Homo_Phenotype1_Male) & !all(geno_list[[allele]]$Homo_Phenotype1_Sex_Ratio == 1)),
        !(!all(geno_list[[allele]]$Homo_Phenotype2_Male) & !all(geno_list[[allele]]$Homo_Phenotype2_Sex_Ratio == 1)),
        !(!all(geno_list[[allele]]$Homo_Phenotype3_Male) & !all(geno_list[[allele]]$Homo_Phenotype3_Sex_Ratio == 1)),
        !(!all(geno_list[[allele]]$Homo_Phenotype4_Male) & !all(geno_list[[allele]]$Homo_Phenotype4_Sex_Ratio == 1)),
        !(!all(geno_list[[allele]]$Homo_Phenotype5_Male) & !all(geno_list[[allele]]$Homo_Phenotype5_Sex_Ratio == 1)),
        !(!all(geno_list[[allele]]$Homo_Phenotype6_Male) & !all(geno_list[[allele]]$Homo_Phenotype6_Sex_Ratio == 1)),
        !(!all(geno_list[[allele]]$Homo_Phenotype7_Male) & !all(geno_list[[allele]]$Homo_Phenotype7_Sex_Ratio == 1)),
        !(!all(geno_list[[allele]]$Homo_Phenotype8_Male) & !all(geno_list[[allele]]$Homo_Phenotype8_Sex_Ratio == 1)),
        !(!all(geno_list[[allele]]$Homo_Phenotype9_Male) & !all(geno_list[[allele]]$Homo_Phenotype9_Sex_Ratio == 1)),
        !(!all(geno_list[[allele]]$Homo_Phenotype10_Male) & !all(geno_list[[allele]]$Homo_Phenotype10_Sex_Ratio == 1)),
        Male_no_SD)
      
      Female_no_SD <- all(
        all(geno_list[[allele]]$Carrier_Sex_Ratio_Females == 1),
        !(!all(geno_list[[allele]]$Homo_Phenotype1_Female) & !all(geno_list[[allele]]$Homo_Phenotype1_Sex_Ratio == 1)),
        !(!all(geno_list[[allele]]$Homo_Phenotype2_Female) & !all(geno_list[[allele]]$Homo_Phenotype2_Sex_Ratio == 1)),
        !(!all(geno_list[[allele]]$Homo_Phenotype3_Female) & !all(geno_list[[allele]]$Homo_Phenotype3_Sex_Ratio == 1)),
        !(!all(geno_list[[allele]]$Homo_Phenotype4_Female) & !all(geno_list[[allele]]$Homo_Phenotype4_Sex_Ratio == 1)),
        !(!all(geno_list[[allele]]$Homo_Phenotype5_Female) & !all(geno_list[[allele]]$Homo_Phenotype5_Sex_Ratio == 1)),
        !(!all(geno_list[[allele]]$Homo_Phenotype6_Female) & !all(geno_list[[allele]]$Homo_Phenotype6_Sex_Ratio == 1)),
        !(!all(geno_list[[allele]]$Homo_Phenotype7_Female) & !all(geno_list[[allele]]$Homo_Phenotype7_Sex_Ratio == 1)),
        !(!all(geno_list[[allele]]$Homo_Phenotype8_Female) & !all(geno_list[[allele]]$Homo_Phenotype8_Sex_Ratio == 1)),
        !(!all(geno_list[[allele]]$Homo_Phenotype9_Female) & !all(geno_list[[allele]]$Homo_Phenotype9_Sex_Ratio == 1)),
        !(!all(geno_list[[allele]]$Homo_Phenotype10_Female) & !all(geno_list[[allele]]$Homo_Phenotype10_Sex_Ratio == 1)),
        Female_no_SD)
    }
    
    Male_no_SD <- all(Male_no_SD, Female_no_SD)
    Female_no_SD <- all(Male_no_SD, Female_no_SD)
    
    
    Delta_M <- cbind(
      expand.grid(
        SD = if (Male_no_SD) {0} else {c(0, 1)},
        Dep = as.factor(nuc_rates_list$M_Dep),
        gRNA_Deposition = as.factor(grna_dep_m),
        C = as.factor(geno_list$C$Allele),
        B = as.factor(geno_list$B$Allele),
        A = as.factor(geno_list$A$Allele)
      ),
      Freq = as.double(0)
    )
    
    # Here we reorder the columns
    Delta_M <-
      Delta_M[, c("A", "B", "C", "gRNA_Deposition", "Dep", "SD", "Freq")]
    
    Delta_F <- cbind(
      expand.grid(
        SD = if (Female_no_SD) {0} else {c(0, 1)},
        Dep = as.factor(nuc_rates_list$F_Dep),
        gRNA_Deposition = as.factor(grna_dep_f),
        C = as.factor(geno_list$C$Allele),
        B = as.factor(geno_list$B$Allele),
        A = as.factor(geno_list$A$Allele)
      ),
      Freq = as.double(0)
    )
    
    # Here we reorder the columns
    Delta_F <-
      Delta_F[, c("A", "B", "C", "gRNA_Deposition", "Dep", "SD", "Freq")]
    
    # Here we create a vector of the possible diploid genotypes
    Alpha <- expand.grid(
      Dep_M = as.factor(nuc_rates_list$M_Dep),
      grna_dep_M = as.factor(grna_dep_m),
      C_M = as.factor(geno_list$C$Allele),
      B_M = as.factor(geno_list$B$Allele),
      A_M = as.factor(geno_list$A$Allele),
      Dep_F = as.factor(nuc_rates_list$F_Dep),
      grna_dep_F = as.factor(grna_dep_f),
      C_F = as.factor(geno_list$C$Allele),
      B_F = as.factor(geno_list$B$Allele),
      A_F = as.factor(geno_list$A$Allele),
      SD = if (Male_no_SD & Female_no_SD) {0} else {c(0, 1)}
    )
    Alpha$Freq <- as.double(0)
    Alpha <-
      Alpha[, c(
        "A_M",
        "A_F",
        "B_M",
        "B_F",
        "C_M",
        "C_F",
        "grna_dep_M",
        "grna_dep_F",
        "Dep_M",
        "Dep_F",
        "SD",
        "Freq"
      )]
    
    
    Alpha[(
      Alpha[, 1] == c("T") & Alpha$A_F == c("T") &
        Alpha$B_M == c("T") & Alpha$B_F == c("T") &
        Alpha$C_M == c("T") & Alpha$C_F == c("T") &
        Alpha$grna_dep_M == levels(Alpha$grna_dep_M)[1] &
        Alpha$grna_dep_F == levels(Alpha$grna_dep_F)[1] &  
        Alpha$Dep_M == c(0) & Alpha$Dep_F == c(0) &
        Alpha$SD == 0
    ),]$Freq <- 1
    
    
    default_file_location <- if (file.exists("Parameter Files/Default_Release.csv")) {"Parameter Files/Default_Release.csv"} else {file.path(project_root,"Default_Release.csv")}
    
    ## --- 1.  Read defaults ----
    defaults <- read_csv(
      default_file_location,
      col_names = c("Parameter", "Default"),
      col_types = cols(.default = col_character()),
      na = c("", "NA")
    )
    
    ## --- 2.  Read parameter file (any width) ----
    
    if (file.exists(file.path("Parameter Files", Release_File_Name))) {
      params <- read_csv(
        file.path("Parameter Files", Release_File_Name),
        col_names = FALSE,
        na = c("", "NA"),
        col_types = cols(.default = col_character())
      )
    } else {
      message(paste0("Using: ",default_file_location))
      params <- defaults
    }
    
    names(params)[1] <- "Parameter"              # first col = row names
    
    ## --- 3.  Abort early if param file contains unknown rows ----
    unknown <- setdiff(params$Parameter, defaults$Parameter)
    
    if (length(unknown)>0)
      warning("Unknown parameter(s) in parameter file: ",
           paste0(Release_File_Name,":  ", unknown, collapse = ", "))
    
    
    patched <- defaults %>%
      left_join(params, by = "Parameter") %>% # Suffixes .x and .y automatically on conflict
      {
        # df_after_join is the data passed into this block (the result of left_join)
        df_after_join <- .
        
        # Determine the actual name of the 'Default' column from the 'defaults' dataframe.
        # It's "Default.x" if 'params' also had a 'Default' column, otherwise it's "Default".
        defaults_source_column_name <- if ("Default.x" %in% names(df_after_join)) "Default.x" else "Default"
        
        # Identify the columns that need to be patched.
        # These are all columns EXCEPT "Parameter" and the defaults_source_column_name.
        # These target columns might have come from 'params' (and potentially suffixed with .y)
        # or were other data columns in 'defaults' (not its 'Default' column).
        columns_to_patch <- setdiff(names(df_after_join), c("Parameter", defaults_source_column_name))
        
        # If there are columns to patch, apply the default values
        if (length(columns_to_patch) > 0) {
          df_mutated <- mutate(
            df_after_join,
            across(
              all_of(columns_to_patch),
              # For each value (.x) in the current column being processed by across:
              # if it's NA or empty, replace it with the value from defaults_source_column_name
              # .data[[defaults_source_column_name]] robustly accesses the column
              ~ if_else(is.na(.) | . == "", .data[[defaults_source_column_name]], .)
            )
          )
        } else {
          # No columns to patch (e.g., if params only had 'Parameter' or was empty)
          df_mutated <- df_after_join
        }
        df_mutated
      } %>%
      {
        # df_after_mutation is the data passed into this block
        df_after_mutation <- .
        
        # Again, determine the name of the default column from 'defaults' to remove it
        defaults_source_column_name_to_remove <- if ("Default.x" %in% names(df_after_mutation)) "Default.x" else "Default"
        
        df_selected <- select(df_after_mutation, -all_of(defaults_source_column_name_to_remove))
        
        # Rename any columns that ended with ".y" (came from params and had a name conflict)
        # back to their original names. e.g., "Default.y" becomes "Default".
        columns_suffixed_y <- grep("\\.y$", names(df_selected), value = TRUE)
        if (length(columns_suffixed_y) > 0) {
          new_names_for_y_columns <- sub("\\.y$", "", columns_suffixed_y)
          df_renamed <- rename_with(df_selected,
                                    ~ new_names_for_y_columns[match(., columns_suffixed_y)],
                                    .cols = all_of(columns_suffixed_y))
        } else {
          df_renamed <- df_selected
        }
        df_renamed
      } %>%
      arrange(match(Parameter, defaults$Parameter)) # Ensure original row order based on defaults
    
    ## --- 5.  Force correct types (everything read as character so far) ----
    numeric_variables <- c(
      "Simulation_Length",
      "Replicate_Simulations",
      "Generation",
      "Male_Release_Fraction",
      "Female_Release_Fraction",
      "Adult_Pop_Size",
      "Release_Pop_Effective_Size",
      "Migration_Rate"
    )
    
    logical_variables <- c(
      "Score_Seperate_From_Seed",
      "Sex_Mached"
    )
    
    ## t(patched)[2:ncol(patched), ] is what you already have
    mat <- t(patched)[-1, ]        # same thing, just a bit clearer

    if (length((patched[1,]))>2) {
      ## turn every column of that matrix into one list element
      Release_File <- asplit(mat, 2)     # MARGIN = 2  →  split by columns
    } else {
      Release_File <- mat
    }

    ## make each element a plain vector (no dim/row names)
    Release_File <- lapply(Release_File, unname)
    
    ## (optional) restore the column names as list names
    names(Release_File) <- patched$Parameter
    
    ## (optional) convert "20", "150", … to numeric/logic where possible
    Release_File <- lapply(Release_File, type.convert, as.is = TRUE)
    
    Release_File[numeric_variables] <- lapply(Release_File[numeric_variables], as.numeric)
    Release_File[logical_variables] <- lapply(Release_File[logical_variables], as.logical)
    
    string_variables <- c("Release_Type", "Post", "Release_PreorPost_Fitness_Costs", "Migration_Source")
    
    ## ---------------- 2.  FLAG UNDECLARED COLUMNS -------------------------------
    expected_cols   <- c(numeric_variables, logical_variables, string_variables)
    unexpected_cols <- setdiff(names(Release_File), expected_cols)
    
    if (length(unexpected_cols) > 0) {
      message(
        "The following columns in Release_File are NOT in numeric_variables ",
        "or logical_variables: ",
        paste(unexpected_cols, collapse = ", ")
      )
      stop()
    }
    
    Generations <- as.numeric(Release_File$Simulation_Length[1])
    Adult_Pop_Size <- as.numeric(Release_File$Adult_Pop_Size[1])

    # This is an optional check
    if (FALSE) {
      # Here we check the haploid and diploid vectors are set-up correctly
      # We fill in the 2 haploid and 1 diploid frequency vectors according to their alleles
      # We then multiply the haploid vectors and check that we get the same outcome as when we individually filled in the diploid vector (slower)
      
      A_M <- as.list(2 * (1:length(levels(Delta_M$A))))
      names(A_M) <- levels(Delta_M$A)
      
      B_M <- as.list(3 * (1:length(levels(Delta_M$B))))
      names(B_M) <- levels(Delta_M$B)
      
      C_M <- as.list(5 * (1:length(levels(Delta_M$C))))
      names(C_M) <- levels(Delta_M$C)
      
      Dep_M <- as.list(7 * (1:length(levels(Delta_M$Dep))))
      names(Dep_M) <- levels(Delta_M$Dep)
      
      grna_dep_M <- as.list(11 * (1:length(levels(Delta_M$gRNA_Deposition))))
      names(grna_dep_M) <- levels(Delta_M$gRNA_Deposition)
      
      
      A_F <- as.list(13 * (1:length(levels(Delta_F$A))))
      names(A_F) <- levels(Delta_F$A)
      
      B_F <- as.list(17 * (1:length(levels(Delta_F$B))))
      names(B_F) <- levels(Delta_F$B)
      
      C_F <- as.list(19 * (1:length(levels(Delta_F$C))))
      names(C_F) <- levels(Delta_F$C)
      
      Dep_F <- as.list(23 * (1:length(levels(Delta_F$Dep))))
      names(Dep_F) <- levels(Delta_F$Dep)
      
      grna_dep_F <- as.list(29 * (1:length(levels(Delta_F$gRNA_Deposition))))
      names(grna_dep_F) <- levels(Delta_F$gRNA_Deposition)
      
      Delta_M$Freq <-
        as.double(A_M[as.character(Delta_M$A)]) * as.double(B_M[as.character(Delta_M$B)]) * as.double(C_M[as.character(Delta_M$C)]) * as.double(Dep_M[as.character(Delta_M$Dep)]) * as.double(grna_dep_M[as.character(Delta_M$gRNA_Deposition)])
      
      Delta_F$Freq <-
        as.double(A_F[as.character(Delta_F$A)]) * as.double(B_F[as.character(Delta_F$B)]) * as.double(C_F[as.character(Delta_F$C)]) * as.double(Dep_F[as.character(Delta_F$Dep)]) * as.double(grna_dep_F[as.character(Delta_F$gRNA_Deposition)])
      
      
      Alpha$Freq <-
        as.double(A_M[as.character(Alpha$A_M)]) * as.double(A_F[as.character(Alpha$A_F)]) *
        as.double(B_M[as.character(Alpha$B_M)]) * as.double(B_F[as.character(Alpha$B_F)]) *
        as.double(C_M[as.character(Alpha$C_M)]) * as.double(C_F[as.character(Alpha$C_F)]) *
        as.double(Dep_M[as.character(Alpha$Dep_M)]) * as.double(Dep_F[as.character(Alpha$Dep_F)]) *
        as.double(grna_dep_M[as.character(Alpha$grna_dep_M)]) * as.double(grna_dep_F[as.character(Alpha$grna_dep_F)])
      
      # Here we multiply the haploid vectors and check that we get the same outcome as when we individually filled in the diploid vector
      # There is some margin of error for rounding errors
      if (abs(mean(Alpha$Freq - as.vector(as.numeric(
        Delta_M$Freq %*% t(as.numeric(Delta_F$Freq))
      )))) > 1e-10) {
        stop("Haploid and diploid vectors not set up correctly")
      }
      
      # Here we reset the frequencies to 0
      Delta_M$Freq <- 0
      Delta_F$Freq <- 0
      Alpha$Freq <- 0
    }
    
    # 3.
    ### Default Genotypes frequencies ####
    # Imports the default Genotypes frequencies
    
    Alpha_M <- Alpha
    Alpha_M$Freq <- Alpha_M$Freq*0
    Alpha_F <- Alpha_M
    
    Alpha_M[(
      Alpha_M[, 1] == c("T") & Alpha_M$A_F == c("T") &
        Alpha_M$B_M == c("T") & Alpha_M$B_F == c("T") &
        Alpha_M$C_M == c("T") & Alpha_M$C_F == c("T") &
        Alpha_M$grna_dep_M == levels(Alpha_M$grna_dep_M)[1] &
        Alpha_M$grna_dep_F == levels(Alpha_M$grna_dep_F)[1] &  
        Alpha_M$Dep_M == c(0) & Alpha_M$Dep_F == c(0) &
        Alpha_M$SD == 0
    ),]$Freq <- 1
    
    Alpha_F[(
      Alpha_F[, 1] == c("T") & Alpha_F$A_F == c("T") &
        Alpha_F$B_M == c("T") & Alpha_F$B_F == c("T") &
        Alpha_F$C_M == c("T") & Alpha_F$C_F == c("T") &
        Alpha_F$grna_dep_M == levels(Alpha_F$grna_dep_M)[1] &
        Alpha_F$grna_dep_F == levels(Alpha_F$grna_dep_F)[1] &  
        Alpha_F$Dep_M == c(0) & Alpha_F$Dep_F == c(0) &
        Alpha_F$SD == 0
    ),]$Freq <- 1
  
    
    # Validate migration rate
    if (!is.null(Release_File$Migration_Rate) &&
        (Release_File$Migration_Rate[1] < 0 || Release_File$Migration_Rate[1] > 1)) {
      stop("Migration_Rate must be between 0 and 1, got: ", Release_File$Migration_Rate[1])
    }

    # Load migration source data if specified
    Migration_Source_Data_M <- NULL
    Migration_Source_Data_F <- NULL
    if (!is.null(Release_File$Migration_Source) &&
        Release_File$Migration_Source[1] != "none" &&
        Release_File$Migration_Rate[1] > 0) {

      source_prefix <- Release_File$Migration_Source[1]

      # Try Outcomes/ first (where run_simulations saves), then project_root
      for (base_dir in c("Outcomes", file.path(project_root, "Outcomes"), project_root)) {
        m_file <- file.path(base_dir, paste0(source_prefix, "_M.csv"))
        f_file <- file.path(base_dir, paste0(source_prefix, "_F.csv"))
        if (file.exists(m_file) && file.exists(f_file)) {
          Migration_Source_Data_M <- as.matrix(read.csv(m_file, header = TRUE))
          Migration_Source_Data_F <- as.matrix(read.csv(f_file, header = TRUE))
          break
        }
      }

      if (is.null(Migration_Source_Data_M)) {
        warning("Migration source files not found for prefix: ", source_prefix)
      }
    }

    Parameters <- list(
      Alpha = Alpha,
      Alpha_M = Alpha_M,
      Alpha_F = Alpha_F,
      Delta_M = Delta_M,
      Delta_F = Delta_F,
      Generations = Generations,
      geno_list = geno_list,
      nuc_rates_list = nuc_rates_list,
      Adult_Pop_Size = Adult_Pop_Size,
      Release_File = Release_File,
      Migration_Source_Data_M = Migration_Source_Data_M,
      Migration_Source_Data_F = Migration_Source_Data_F
    )

    return(Parameters)
  }

##  ............................................................................
### Post simulation processing                                              ####

# Calculates the individual allele by allele frequency from the combined genotype outcomes matrix
Outcomes_by_Allele <- function(Main_Outcomes, Parameters, replicate) {
  

  # Here we grab the allele legends of each row
  Allele_Lookup <- Main_Outcomes$Outcomes_MF_PreFit[, 1:11]
  
  # Here we store the possible alleles of A, B, and C
  A <- Parameters$geno_list$A$Allele
  B <- Parameters$geno_list$B$Allele
  C <- Parameters$geno_list$C$Allele
  
  # Here we combine all three allele possibilities into one vector
  Allele_Names <- c(A, B, C)
  # We create a vector of the same length that correspond to the locus of each allele
  Locus_Names <-
    c(rep('A', length(A)), rep('B', length(B)), rep('C', length(C)))
  # We create a vector of the same length that corresponds to the column position of the locus in the data frame
  Locus_Allele_Pos <-
    c(rep(1, length(A)), rep(3, length(B)), rep(5, length(C)))
  
  Numb_Alleles <- length(Allele_Names)
  
  # We preinitialise the data frames that will hold the allele by allele outcomes
  # We create a spot for each time point, for each allele
  # we ignore deposition states
  Outcomes_Comb <-
    data.frame(
      t = rep(1:Parameters$Generations, 6*Numb_Alleles),
      Allele_Freq = 0.0,
      Carrier_Freq = 0.0,
      Het_Freq = 0.0,
      Hom_Freq = 0.0,
      Outcome = "",
      Locus = "",
      Sex = "",
      Pre_or_Post_Fit = "",
      Replicate = replicate
    )
  
  Fit_Comb <- data.frame(
    t = rep(1:Parameters$Generations, Numb_Alleles),
    Allele_Freq = 0.0,
    Outcome = "",
    Locus = "",
    Replicate = replicate
  )
  
  Conv_MF <- data.frame(
    t = rep(1:Parameters$Generations, Numb_Alleles),
    Allele_Freq = 0.0,
    Outcome = "",
    Locus = "",
    Replicate = replicate
  )
  
  Tot_Death <- data.frame(
    t = rep(1:Parameters$Generations, Numb_Alleles),
    Allele_Freq = 0.0,
    Outcome = "",
    Locus = "",
    Replicate = replicate
  )
  
  safeColSums <- function(x) {
    if (is.null(nrow(x))) {
      
      if (is.list(x)) {
        return(0)
      }
      return(x)  # Keeps the single row as a matrix
    } else {
      return(colSums(x))
    }
  }
  
  
  # Here we save the total amount of generations + the descriptive columns
  tot_length <- length(Main_Outcomes$Outcomes_MF_PreFit[1,])
  descript_cols_end <- 13 # the number of descriptive columns
  
  # Here we set the threshold frequency above which the drive is considered effective
  # Threshold <- 0.9
  
  data_tabel <- list(
    "1" = Main_Outcomes$Outcomes_MF_PreFit,
    "2" = Main_Outcomes$Outcomes_F_PreFit,
    "3" = Main_Outcomes$Outcomes_M_PreFit,
    "4" = Main_Outcomes$Outcomes_MF_PostFit,
    "5" = Main_Outcomes$Outcomes_F_PostFit,
    "6" = Main_Outcomes$Outcomes_M_PostFit
  )
  
  Allele_rows <- (1:Parameters$Generations)
  
  sex <- c("Male&Female", "Female", "Male", "Male&Female", "Female", "Male")
  pre_post <- c("Pre", "Pre", "Pre", "Post", "Post", "Post")
  
  for (y in 1:length(sex)) {
  # Here we loop over every allele and locus combination
  for (x in 1:Numb_Alleles) {
    Outcomes_Comb[Allele_rows,]$Sex <- sex[y]
    Outcomes_Comb[Allele_rows,]$Pre_or_Post_Fit <- pre_post[y]
    Outcomes_Comb[Allele_rows,]$t <- 1:Parameters$Generations
    Outcomes_Comb[Allele_rows,]$Outcome <- Allele_Names[x]
    Outcomes_Comb[Allele_rows,]$Locus <- Locus_Names[x]
    
    Outcomes_Comb[Allele_rows,]$Allele_Freq <-
      (safeColSums(sapply(
        data_tabel[[y]][Allele_Lookup[, Locus_Allele_Pos[x]] == Allele_Names[x], descript_cols_end:tot_length], as.double
      )) + safeColSums(sapply(
        data_tabel[[y]][Allele_Lookup[, Locus_Allele_Pos[x] + 1] == Allele_Names[x], descript_cols_end:tot_length], as.double
      ))) / 2
    
    
    Outcomes_Comb[Allele_rows,]$Carrier_Freq <- 
      safeColSums(sapply(data_tabel[[y]][ (Allele_Lookup[, Locus_Allele_Pos[x]] == Allele_Names[x]) |
                                                   (Allele_Lookup[, Locus_Allele_Pos[x] + 1] == Allele_Names[x]), descript_cols_end:tot_length], as.double))
    
    Outcomes_Comb[Allele_rows, ]$Het_Freq <- safeColSums(sapply(data_tabel[[y]][xor(Allele_Lookup[, Locus_Allele_Pos[x]] == Allele_Names[x],
                                                                                    Allele_Lookup[, Locus_Allele_Pos[x] + 1] == Allele_Names[x]), descript_cols_end:tot_length], as.double))
    
    Outcomes_Comb[Allele_rows, ]$Hom_Freq <- safeColSums(sapply(data_tabel[[y]][((Allele_Lookup[, Locus_Allele_Pos[x]] == Allele_Names[x]) &
                                                                                   (Allele_Lookup[, Locus_Allele_Pos[x] + 1] == Allele_Names[x])), descript_cols_end:tot_length], as.double))
    
    
    Allele_rows <- Allele_rows + (Parameters$Generations)
    
  }
  }
  return(
    list(
      Outcomes_MF = Outcomes_Comb,
      Fit_MF = Fit_Comb,
      Conv_MF = Conv_MF,
      Tot_Death = Tot_Death
    )
  )
}

# This function adjusts the parameters for Range and SA figs
Parameter_Change <-
  function(Parameters, Value, Factor, Graphs_Current) {
    if (Factor == "no-change") {
      
    } else {
      stop("Unknown variable specified as factor in range graph")
    }
    
    return(Parameters)
  }

# Latex commands cant have numbers in them, so here we create a letter system for numbering
alphabet_numbering <-
  paste(
    expand.grid(letters, letters, stringsAsFactors = TRUE)[, 2],
    expand.grid(letters, letters, stringsAsFactors = TRUE)[, 1],
    sep = ""
  )

#   ____________________________________________________________________________
#   Generating the simulation instruction file                              ####

# The below code generates an simulation instruction file
# If that has already been made, this can be skipped (set to FALSE) and simply read in after this section
create_figure_instruc <- function(Figures_to_make) {
  
  # Here we give each figure a unique letter (latex commands don't allow numbers)
  Figures_to_make$Numbering <-
    alphabet_numbering[1:length(Figures_to_make[, 1])]

  # Here we create a name that is unique for the particular simulation
  # It includes the figure numbering
  Figures_to_make$Fig_long_name <-
    paste(paste(Figures_to_make$Fig_name,
                Figures_to_make$Allele_Fit_File, sep = "_"),
          sep = "_")
  
  # Here we expand the figure instructions for each parameter sweep step
  Graphs <-
    as.data.frame(
      lapply(
        Figures_to_make,
        rep,
        Figures_to_make$Steps_F1 * Figures_to_make$Steps_F2
      )
    )
  
  # We add a column that contains the exact file name that will be used
  Graphs$File_Name <- NA
  
  # We add additional columns that help describe the particular simulation in a sweep
  Graphs$i <- NA
  Graphs$j <- NA
  Graphs$Factor_F1_value <- NA
  Graphs$Factor_F2_value <- NA
  
  # Here we loop over every figure, and create the instructions for the parameter sweeps
  # For example, 'Figures_to_make' has the instruction for every figure. But some figures are composed of SA_size*SA_size simulations
  # In the below loop we fill in those specific parameter sweep instructions
  for (fig_idex in seq_len(length(Figures_to_make$Run_Fig))) {
    
    
    # Here we create a list of the possible parameter steps we will use
    # Here we set the names of the steps list
    Steps_list <- list('i' = 0, 'j' = 0)
    
    # Here we populate the steps list with the steps we will take (if any)
    for (s in 1:2) {
      
      range_type <- c(Figures_to_make$Range_Type_F1[fig_idex], Figures_to_make$Range_Type_F2[fig_idex])[s]
      step_start <- c(Figures_to_make$Start_F1[fig_idex], Figures_to_make$Start_F2[fig_idex])[s]
      step_end <- c(Figures_to_make$Stop_F1[fig_idex], Figures_to_make$Stop_F2[fig_idex])[s]
      step_number <- c(Figures_to_make$Steps_F1[fig_idex], Figures_to_make$Steps_F2[fig_idex])[s]
      step_exponent <- c(Figures_to_make$Exponent_F1[fig_idex], Figures_to_make$Exponent_F2[fig_idex])[s]
      
      if (range_type == "Lin") {
        Steps_list[[s]] <- 
          seq(step_start, step_end, length.out = step_number)
        
      } else if (range_type == "Exp") {
        Steps_list[[s]] <-
          step_exponent ^ seq(step_start, step_end, length.out = step_number)
        
      } else if (range_type == "log") {
        Steps_list[[s]] <- log(seq(step_start, step_end, length.out = step_number))/log(step_end)
        
      } else if (range_type == "none") {
        Steps_list[[s]] <- 1
        
      } else {
        stop("Unknown F1 step factor (e.i, not Linear/Exponent/log)")
      }
    }
    
    print(paste('Fig index:', fig_idex))
    
    unique_combos <- expand.grid(i = seq_len(length(Steps_list[[1]])),
                                 j = seq_len(length(Steps_list[[2]])))
    
    value_i <- Steps_list[[1]][unique_combos$i]
    value_j <- Steps_list[[2]][unique_combos$j]
    
    Factor1_name <-
      paste(Figures_to_make$Factor_F1[fig_idex],
            round(value_i, 3),
            sep = "_")
    Factor2_name <-
      paste(Figures_to_make$Factor_F2[fig_idex],
            round(value_j, 3),
            sep = "_")
    
    dir.create("Outcomes",showWarnings = FALSE)
    for (i in paste0("Outcomes/", Graphs[Graphs$Numbering == Figures_to_make$Numbering[fig_idex],]$Project_Folder)) {
      dir.create(i,showWarnings = FALSE)
    }
    Graphs[Graphs$Numbering == Figures_to_make$Numbering[fig_idex],]$File_Name <-
      paste0(Graphs[Graphs$Numbering == Figures_to_make$Numbering[fig_idex],]$Project_Folder,"/",
            paste(
              Graphs[Graphs$Numbering == Figures_to_make$Numbering[fig_idex],]$Fig_long_name,
              apply(data.frame(Factor1_name, Factor2_name), 1, paste, collapse = "_"),
              sep = ""),sep = "")
    
    Graphs[Graphs$Numbering == Figures_to_make$Numbering[fig_idex],]$i <- unique_combos$i
    Graphs[Graphs$Numbering == Figures_to_make$Numbering[fig_idex],]$j <- unique_combos$j
    
    Graphs[Graphs$Numbering == Figures_to_make$Numbering[fig_idex],]$Factor_F1_value <- value_i
    Graphs[Graphs$Numbering == Figures_to_make$Numbering[fig_idex],]$Factor_F2_value <- value_j
    
  }
  
  # Graphs <- Graphs[order(Graphs$Figure_type),]
  # write.csv(apply(Graphs, 2, as.character), "Graphs.csv")
  
  return(Graphs)
}

##  ............................................................................
##  Conversion matrix simplification                                        ####

# Function to sort alleles within a genotype to make order irrelevant
sort_genotype <- function(genotype) {

  # Split genotype into alleles
  alleles <- strsplit(genotype, "/")[[1]]
  
  # Sort alleles alphabetically
  sorted_alleles <- sort(alleles)
  
  # Recombine into genotype string
  paste(sorted_alleles, collapse = "/")
}

# Function to identify if a genotype is in the identity matrix
is_identity_match <- function(genotype1, genotype2) {
  # Sort both genotypes
  sorted1 <- sort_genotype(genotype1)
  sorted2 <- sort_genotype(genotype2)
  
  # Check if they match after sorting
  sorted1 == sorted2
}

# Function to collapse similar genotypes using wildcards
collapse_genotype <- function(genotype) {
  # Split genotype into alleles
  
  
  alleles <- strsplit(genotype, "/")[[1]]
  
  # Replace non-T alleles with wildcard
  collapsed_alleles <- sapply(alleles, function(allele) {
    if (allele == "T") {
      return("T")
    } else {
      return("(non-T)")
    }
  })
  
  # Recombine into genotype string
  paste(collapsed_alleles, collapse = "/")
}

# Function to simplify conversion matrix
simplify_conversion_matrix <- function(conversion_matrix, genotypes) {
  
  
  genotypes$identity <- FALSE
  
  # Loop through each row of conversion matrix
  for(i in 1:nrow(conversion_matrix)) {

    matrix_row <- conversion_matrix[,i]
    matrix_row[i] <- matrix_row[i] - 1
    
    # Add genotype to appropriate vector
    if(all(matrix_row == 0)) {
      genotypes$identity[i] <- TRUE
    } 
  }

  genotypes <- data.frame(lapply(genotypes, as.character), stringsAsFactors=FALSE)
  genotypes$identity <- as.logical(genotypes$identity)
  genotypes_new <- genotypes
  
  for (allele in colnames(genotypes[1:10])) {
    for (i in unique(genotypes[, allele])[!unique(genotypes[genotypes$identity, ][, allele]) %in% unique(genotypes[!genotypes$identity, ][, allele])]) {
      genotypes_new[genotypes_new[, allele] == i, allele] <- gsub(
        x = genotypes[genotypes[, allele] == i, allele],
        pattern = i,
        replacement = paste0(
          "(", paste0(
            unique(genotypes[, allele])[!unique(genotypes[genotypes$identity, ][, allele]) %in% 
                                        unique(genotypes[!genotypes$identity, ][, allele])],
                      collapse = "/"), ")")
      )
    }
  }
  
  return(genotypes_new)
}


# ---- Summary Table & Meta-Table Helper Functions ----

# Converts a parameter set to a canonical string for identity comparison
hash_parameter_set <- function(param_set) {
  if (is.data.frame(param_set)) {
    # For data.frames (e.g., geno_list$A): sort by row, serialize all values
    df <- param_set[order(rownames(param_set)), , drop = FALSE]
    paste(apply(df, 1, function(row) paste(sort(names(row)), row[sort(names(row))], sep = "=", collapse = ",")), collapse = "|")
  } else if (is.list(param_set)) {
    # For named lists (nuc_rates_list, Release_File): sort by name, serialize
    nms <- sort(names(param_set))
    paste(sapply(nms, function(nm) {
      val <- param_set[[nm]]
      paste0(nm, "=", paste(val, collapse = ";"))
    }), collapse = "|")
  } else {
    as.character(param_set)
  }
}

# Assigns a profile ID for a given hash, or returns existing one
assign_profile_id <- function(hash_value, profile_registry, data = NULL) {
  idx <- match(hash_value, profile_registry$hashes)
  if (!is.na(idx)) {
    return(list(id = profile_registry$ids[idx], registry = profile_registry))
  }
  new_id <- paste0("Profile_", length(profile_registry$hashes) + 1)
  profile_registry$hashes <- c(profile_registry$hashes, hash_value)
  profile_registry$ids <- c(profile_registry$ids, new_id)
  profile_registry$data[[length(profile_registry$data) + 1]] <- data
  return(list(id = new_id, registry = profile_registry))
}

# Flattens all parameters into a single named list of scalars
flatten_parameters_to_row <- function(Parameters, Graphs_row) {
  flat <- list()

  # Nuc parameters (skip computed fields)
  skip_nuc <- c("F_Dep", "M_Dep")
  for (nm in names(Parameters$nuc_rates_list)) {
    if (nm %in% skip_nuc) next
    val <- Parameters$nuc_rates_list[[nm]]
    if (length(val) == 1) {
      flat[[paste0("Nuc_", nm)]] <- as.character(val)
    } else {
      for (k in seq_along(val)) {
        flat[[paste0("Nuc_", nm, "_", k)]] <- as.character(val[k])
      }
    }
  }

  # Drive design parameters (allele names and structure per locus)
  for (locus in c("A", "B", "C")) {
    geno_df <- Parameters$geno_list[[locus]]
    if (is.null(geno_df) || nrow(geno_df) == 0) next
    allele_names <- if ("Alleles" %in% names(geno_df)) geno_df$Alleles else rownames(geno_df)
    flat[[paste0("Design_", locus, "_Num_Alleles")]] <- as.character(length(allele_names))
    for (ai in seq_along(allele_names)) {
      flat[[paste0("Design_", locus, "_Pos", ai)]] <- as.character(allele_names[ai])
    }
  }

  # Fitness parameters (non-phenotype)
  # Per-locus params (same for all alleles at a locus) vs per-allele params
  phenotype_cols <- grep("^Homo_Phenotype", names(Parameters$geno_list$A), value = TRUE)
  per_locus_params <- c("nonFunctional_NHEJ_Fraction_Locus", "HDR_Penalty_Locus",
                        "Somatic_Cutting_Sensitivity_Locus", "Fitness_Applied_At_Stage_Locus",
                        "Linkage_Disequilibrium")
  skip_cols <- c("Locus", "Alleles")

  for (locus in c("A", "B", "C")) {
    geno_df <- Parameters$geno_list[[locus]]
    if (is.null(geno_df) || nrow(geno_df) == 0) next

    # Per-locus params: extract once (from first allele, since they're the same)
    for (col in intersect(per_locus_params, names(geno_df))) {
      flat[[paste0("Fitness_", locus, "_", col)]] <- as.character(geno_df[1, col])
    }

    # Per-allele params: positional indexing
    for (ai in seq_len(nrow(geno_df))) {
      for (col in names(geno_df)) {
        if (col %in% skip_cols || col %in% phenotype_cols || col %in% per_locus_params) next
        flat[[paste0("Fitness_", locus, "_Pos", ai, "_", col)]] <- as.character(geno_df[ai, col])
      }
    }
  }

  # Phenotype fitness parameters (global cost + per-allele rescue flags)
  # Fitness values are enforced identical across all alleles, so extract once
  for (pheno_i in 1:10) {
    cost_col <- paste0("Homo_Phenotype", pheno_i, "_Homozygous_Fitness")
    if (!(cost_col %in% names(Parameters$geno_list$A))) next
    cost_val <- Parameters$geno_list$A[[cost_col]][1]
    carrier_col <- paste0("Homo_Phenotype", pheno_i, "_Carrier_Fitness")
    carrier_val <- if (carrier_col %in% names(Parameters$geno_list$A)) Parameters$geno_list$A[[carrier_col]][1] else 1
    sr_col <- paste0("Homo_Phenotype", pheno_i, "_Sex_Ratio")
    sr_val <- if (sr_col %in% names(Parameters$geno_list$A)) Parameters$geno_list$A[[sr_col]][1] else 1

    # Skip inactive phenotypes (all values = 1 means no effect)
    if (isTRUE(cost_val == 1 & carrier_val == 1 & sr_val == 1)) next

    flat[[paste0("Phenotype", pheno_i, "_Homozygous_Fitness")]] <- as.character(cost_val)
    flat[[paste0("Phenotype", pheno_i, "_Carrier_Fitness")]] <- as.character(carrier_val)
    flat[[paste0("Phenotype", pheno_i, "_Sex_Ratio")]] <- as.character(sr_val)

    # Build rescue allele lists from Male/Female flags
    male_col <- paste0("Homo_Phenotype", pheno_i, "_Male")
    female_col <- paste0("Homo_Phenotype", pheno_i, "_Female")

    for (sex in c("Male", "Female")) {
      flag_col <- if (sex == "Male") male_col else female_col
      rescue_parts <- c()
      for (locus in c("A", "B", "C")) {
        geno_df <- Parameters$geno_list[[locus]]
        if (is.null(geno_df) || nrow(geno_df) == 0) next
        if (!(flag_col %in% names(geno_df))) next
        allele_names <- if ("Alleles" %in% names(geno_df)) geno_df$Alleles else rownames(geno_df)
        rescue_alleles <- allele_names[as.logical(geno_df[[flag_col]])]
        if (length(rescue_alleles) > 0) {
          rescue_parts <- c(rescue_parts, paste0(locus, ":", paste(rescue_alleles, collapse = ",")))
        }
      }
      flat[[paste0("Phenotype", pheno_i, "_", sex, "_Rescue")]] <- paste(rescue_parts, collapse = "; ")
    }
  }

  # Release parameters
  for (nm in names(Parameters$Release_File)) {
    val <- Parameters$Release_File[[nm]]
    if (length(val) == 1) {
      flat[[paste0("Release_", nm)]] <- as.character(val)
    } else {
      for (k in seq_along(val)) {
        flat[[paste0("Release_", nm, "_", k)]] <- as.character(val[k])
      }
    }
  }

  return(flat)
}

# Generates a human-readable summary table (tall format)
generate_summary_table <- function(Parameters, Graphs_row) {
  rows <- list()
  add_row <- function(section, param, value) {
    rows[[length(rows) + 1]] <<- data.frame(
      Section = section, Parameter = param, Value = as.character(value),
      stringsAsFactors = FALSE
    )
  }

  # Drive Design
  add_row("Drive Design", "Project_Folder", Graphs_row$Project_Folder)
  add_row("Drive Design", "Fig_name", Graphs_row$Fig_name)
  add_row("Drive Design", "Allele_Fit_File", Graphs_row$Allele_Fit_File)
  add_row("Drive Design", "Nuc_File", Graphs_row$Nuc_File)

  for (locus in c("A", "B", "C")) {
    geno_df <- Parameters$geno_list[[locus]]
    if (!is.null(geno_df) && nrow(geno_df) > 0) {
      allele_names <- if ("Alleles" %in% names(geno_df)) geno_df$Alleles else rownames(geno_df)
      add_row("Drive Design", paste0("Locus_", locus, "_Alleles"),
              paste(allele_names, collapse = ", "))
    }
  }

  # Nuclease Rates
  skip_nuc <- c("F_Dep", "M_Dep")
  for (nm in names(Parameters$nuc_rates_list)) {
    if (nm %in% skip_nuc) next
    val <- Parameters$nuc_rates_list[[nm]]
    add_row("Nuclease Rates", nm, paste(val, collapse = ", "))
  }

  # Fitness parameters (non-phenotype)
  phenotype_cols <- grep("^Homo_Phenotype", names(Parameters$geno_list$A), value = TRUE)
  per_locus_params <- c("nonFunctional_NHEJ_Fraction_Locus", "HDR_Penalty_Locus",
                        "Somatic_Cutting_Sensitivity_Locus", "Fitness_Applied_At_Stage_Locus",
                        "Linkage_Disequilibrium")
  for (locus in c("A", "B", "C")) {
    geno_df <- Parameters$geno_list[[locus]]
    if (is.null(geno_df) || nrow(geno_df) == 0) next
    section <- paste0("Fitness - ", locus)

    # Per-locus parameters (once)
    for (col in intersect(per_locus_params, names(geno_df))) {
      add_row(section, col, as.character(geno_df[1, col]))
    }

    # Per-allele parameters
    allele_names <- if ("Alleles" %in% names(geno_df)) geno_df$Alleles else rownames(geno_df)
    for (ai in seq_along(allele_names)) {
      allele <- allele_names[ai]
      for (col in names(geno_df)) {
        if (col %in% c("Locus", "Alleles") || col %in% phenotype_cols || col %in% per_locus_params) next
        add_row(section, paste0(allele, ": ", col), as.character(geno_df[ai, col]))
      }
    }
  }

  # Phenotype fitness (global cost + rescue alleles)
  for (pheno_i in 1:10) {
    cost_col <- paste0("Homo_Phenotype", pheno_i, "_Homozygous_Fitness")
    if (!(cost_col %in% names(Parameters$geno_list$A))) next
    cost_val <- Parameters$geno_list$A[[cost_col]][1]
    carrier_col <- paste0("Homo_Phenotype", pheno_i, "_Carrier_Fitness")
    carrier_val <- if (carrier_col %in% names(Parameters$geno_list$A)) Parameters$geno_list$A[[carrier_col]][1] else 1
    sr_col <- paste0("Homo_Phenotype", pheno_i, "_Sex_Ratio")
    sr_val <- if (sr_col %in% names(Parameters$geno_list$A)) Parameters$geno_list$A[[sr_col]][1] else 1

    if (isTRUE(cost_val == 1 & carrier_val == 1 & sr_val == 1)) next

    section <- paste0("Phenotype ", pheno_i)
    add_row(section, "Homozygous_Fitness", cost_val)
    add_row(section, "Carrier_Fitness", carrier_val)
    add_row(section, "Sex_Ratio", sr_val)

    male_col <- paste0("Homo_Phenotype", pheno_i, "_Male")
    female_col <- paste0("Homo_Phenotype", pheno_i, "_Female")
    for (sex in c("Male", "Female")) {
      flag_col <- if (sex == "Male") male_col else female_col
      rescue_parts <- c()
      for (locus in c("A", "B", "C")) {
        geno_df <- Parameters$geno_list[[locus]]
        if (is.null(geno_df) || nrow(geno_df) == 0 || !(flag_col %in% names(geno_df))) next
        allele_names <- if ("Alleles" %in% names(geno_df)) geno_df$Alleles else rownames(geno_df)
        rescue_alleles <- allele_names[as.logical(geno_df[[flag_col]])]
        if (length(rescue_alleles) > 0) {
          rescue_parts <- c(rescue_parts, paste0(locus, ":", paste(rescue_alleles, collapse = ",")))
        }
      }
      add_row(section, paste0(sex, "_Rescue"), paste(rescue_parts, collapse = "; "))
    }
  }

  # Release Parameters
  for (nm in names(Parameters$Release_File)) {
    val <- Parameters$Release_File[[nm]]
    add_row("Release", nm, paste(val, collapse = ", "))
  }

  # Sweep Configuration
  add_row("Sweep", "Factor_F1", as.character(Graphs_row$Factor_F1))
  add_row("Sweep", "Factor_F1_Value", as.character(Graphs_row$Factor_F1_value))
  add_row("Sweep", "Factor_F2", as.character(Graphs_row$Factor_F2))
  add_row("Sweep", "Factor_F2_Value", as.character(Graphs_row$Factor_F2_value))

  do.call(rbind, rows)
}

# Generates two meta-tables from accumulated rows:
# 1. Profiles-only table (compact overview)
# 2. Detailed table (profiles + individual varying parameters)
generate_meta_tables <- function(meta_accumulator) {
  if (length(meta_accumulator$rows) == 0) return(list(profiles = data.frame(), detailed = data.frame()))

  # Combine all rows (bind_rows handles different column sets with NA fill)
  meta_df <- bind_rows(lapply(meta_accumulator$rows, function(row) {
    as_tibble(lapply(row, function(x) as.character(x)))
  }))

  # Identity columns always kept (but only if they vary or are identifiers)
  identity_cols <- c("Numbering", "Project_Folder", "Fig_name")
  profile_cols <- c("Profile_Design", "Profile_Fitness", "Profile_Nuc", "Profile_Release")

  # Find columns that vary across simulations
  varying_cols <- sapply(names(meta_df), function(nm) {
    col <- meta_df[[nm]]
    # For Phenotype columns, treat NA vs non-NA as variation
    # (phenotype active in some sims but not others is meaningful)
    if (startsWith(nm, "Phenotype")) {
      length(unique(col)) > 1
    } else {
      # For other columns, only compare non-NA values
      vals <- na.omit(col)
      length(unique(vals)) > 1
    }
  })
  varying_names <- names(which(varying_cols))

  # Profiles-only table: identity + profiles + any varying sweep/config columns
  profile_keep <- unique(c(identity_cols, profile_cols,
                           intersect(varying_names, names(meta_df))))
  # For profiles table, only include non-individual-parameter varying columns
  param_prefix <- c("Design_", "Fitness_", "Phenotype", "Nuc_", "Release_")
  is_individual_param <- sapply(varying_names, function(nm) any(startsWith(nm, param_prefix)))
  profile_varying <- varying_names[!is_individual_param]
  profile_ordered <- c(identity_cols, profile_cols, sort(profile_varying))
  profile_ordered <- intersect(profile_ordered, names(meta_df))

  # Detailed table: identity + profiles + all varying columns
  detailed_varying <- sort(varying_names[!varying_names %in% c(identity_cols, profile_cols)])
  detailed_ordered <- c(identity_cols, profile_cols, detailed_varying)
  detailed_ordered <- intersect(detailed_ordered, names(meta_df))

  list(
    profiles = meta_df[, profile_ordered, drop = FALSE],
    detailed = meta_df[, detailed_ordered, drop = FALSE]
  )
}

# Generates profile lookup tables in tall format
generate_profile_lookups <- function(profile_registry, type) {
  if (length(profile_registry$data) == 0) return(data.frame())

  rows <- list()

  for (i in seq_along(profile_registry$ids)) {
    pid <- profile_registry$ids[i]
    dat <- profile_registry$data[[i]]

    if (type == "design" && is.list(dat)) {
      # dat is a list of character vectors (allele names per locus)
      for (locus in names(dat)) {
        allele_names <- dat[[locus]]
        for (ai in seq_along(allele_names)) {
          rows[[length(rows) + 1]] <- data.frame(
            Profile_ID = pid, Locus = locus, Position = ai,
            Allele = allele_names[ai],
            stringsAsFactors = FALSE
          )
        }
      }
    } else if (type == "fitness" && is.list(dat)) {
      # dat is a list of data.frames (fitness values per locus, no allele names)
      for (locus in names(dat)) {
        geno_df <- dat[[locus]]
        if (is.null(geno_df) || nrow(geno_df) == 0) next
        for (ai in seq_len(nrow(geno_df))) {
          for (col in names(geno_df)) {
            rows[[length(rows) + 1]] <- data.frame(
              Profile_ID = pid, Locus = locus, Position = ai,
              Parameter = col, Value = as.character(geno_df[ai, col]),
              stringsAsFactors = FALSE
            )
          }
        }
      }
    } else if (type == "release" && is.list(dat)) {
      for (nm in names(dat)) {
        val <- dat[[nm]]
        if (length(val) == 1) {
          rows[[length(rows) + 1]] <- data.frame(
            Profile_ID = pid, Release_Event = 1, Parameter = nm,
            Value = as.character(val), stringsAsFactors = FALSE
          )
        } else {
          for (k in seq_along(val)) {
            rows[[length(rows) + 1]] <- data.frame(
              Profile_ID = pid, Release_Event = k, Parameter = nm,
              Value = as.character(val[k]), stringsAsFactors = FALSE
            )
          }
        }
      }
    } else if (is.list(dat)) {
      # nuc type
      for (nm in names(dat)) {
        if (nm %in% c("F_Dep", "M_Dep")) next
        val <- dat[[nm]]
        rows[[length(rows) + 1]] <- data.frame(
          Profile_ID = pid, Parameter = nm,
          Value = paste(as.character(val), collapse = ", "),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(rows) == 0) return(data.frame())
  do.call(rbind, rows)
}


run_simulations <- function(Graphs, project_root) {

# cl <- makeCluster(detectCores() - 2)
# registerDoParallel(cl)
# on.exit(stopCluster(cl))

# Meta-table accumulator
meta_accumulator <- list(
  rows = list(),
  profile_design  = list(hashes = character(), ids = character(), data = list()),
  profile_fitness = list(hashes = character(), ids = character(), data = list()),
  profile_nuc     = list(hashes = character(), ids = character(), data = list()),
  profile_release = list(hashes = character(), ids = character(), data = list())
)
seen_groups <- character()

# Figure creation loop ####
# We iterate over every figure pre-set, and store the simulation data as a csv file
# foreach (fig_idex = ((1:length(Graphs$Run_Fig))[Graphs$Run_Fig]), .packages = c("Matrix", "dplyr", "tidyr")) %dopar% {
for (fig_idex in ((1:length(Graphs$Run_Fig))[Graphs$Run_Fig])) {
  
  print(Graphs$File_Name[fig_idex])
  print(paste(
    Graphs$i[fig_idex] + ((Graphs$j[fig_idex] - 1) * Graphs$Steps_F1[fig_idex]),
    (Graphs$Steps_F1[fig_idex] * Graphs$Steps_F2[fig_idex]),
    sep = "/"
  ))
  
  # if (!file.exists(paste(Graphs$File_Name[fig_idex], ".csv", sep = "_Alleles"))) {
  if(TRUE) {
    # Here we load the basic parameters
    # These may be changed if we are doing 'range' or 'SA' fig
    
    Project_Folder <- paste0(Graphs$Project_Folder[fig_idex],"/")
    
    Allele_Fit_File = paste0(Project_Folder, Graphs$Allele_Fit_File[fig_idex], ".csv")
    Nuc_File = paste0(Project_Folder, Graphs$Nuc_File[fig_idex], ".csv")
    Release_File = paste0(Project_Folder, Graphs$Release_File[fig_idex], ".csv")
    
    Parameters <-
      all_default_parameters(
        Allele_Fit_File = Allele_Fit_File,
        Nuc_File = Nuc_File,
        Release_File = Release_File, 
        project_root = project_root
      )
    
    # Accumulate meta-table row (base parameters, before sweep adjustment)
    if (!(Graphs$Numbering[fig_idex] %in% seen_groups)) {
      seen_groups <- c(seen_groups, Graphs$Numbering[fig_idex])

      # Build design-only data (allele names/structure) for design profile
      design_data <- list()
      for (locus in c("A", "B", "C")) {
        geno_df <- Parameters$geno_list[[locus]]
        if (!is.null(geno_df) && nrow(geno_df) > 0) {
          allele_names <- if ("Alleles" %in% names(geno_df)) geno_df$Alleles else rownames(geno_df)
          design_data[[locus]] <- allele_names
        }
      }

      # Build fitness-only data (values without allele names) for fitness profile
      fitness_data <- list()
      for (locus in c("A", "B", "C")) {
        geno_df <- Parameters$geno_list[[locus]]
        if (!is.null(geno_df) && nrow(geno_df) > 0) {
          fitness_data[[locus]] <- geno_df[, !names(geno_df) %in% c("Locus", "Alleles"), drop = FALSE]
        }
      }

      # Hash and assign profiles for each parameter class
      design_hash <- hash_parameter_set(design_data)
      result_d <- assign_profile_id(design_hash, meta_accumulator$profile_design, design_data)
      meta_accumulator$profile_design <- result_d$registry

      fitness_hash <- hash_parameter_set(fitness_data)
      result_f <- assign_profile_id(fitness_hash, meta_accumulator$profile_fitness, fitness_data)
      meta_accumulator$profile_fitness <- result_f$registry

      nuc_hash <- hash_parameter_set(Parameters$nuc_rates_list)
      result_n <- assign_profile_id(nuc_hash, meta_accumulator$profile_nuc, Parameters$nuc_rates_list)
      meta_accumulator$profile_nuc <- result_n$registry

      release_hash <- hash_parameter_set(Parameters$Release_File)
      result_r <- assign_profile_id(release_hash, meta_accumulator$profile_release, Parameters$Release_File)
      meta_accumulator$profile_release <- result_r$registry

      # Flatten base parameters and build meta row
      flat_params <- flatten_parameters_to_row(Parameters, Graphs[fig_idex, ])
      meta_row <- c(
        list(
          Numbering = Graphs$Numbering[fig_idex],
          Project_Folder = Graphs$Project_Folder[fig_idex],
          Fig_name = Graphs$Fig_name[fig_idex],
          Profile_Design = result_d$id,
          Profile_Fitness = result_f$id,
          Profile_Nuc = result_n$id,
          Profile_Release = result_r$id,
          Factor_F1 = as.character(Graphs$Factor_F1[fig_idex]),
          Factor_F1_Start = as.character(Graphs$Start_F1[fig_idex]),
          Factor_F1_Stop = as.character(Graphs$Stop_F1[fig_idex]),
          Factor_F1_Steps = as.character(Graphs$Steps_F1[fig_idex]),
          Factor_F2 = as.character(Graphs$Factor_F2[fig_idex]),
          Factor_F2_Start = as.character(Graphs$Start_F2[fig_idex]),
          Factor_F2_Stop = as.character(Graphs$Stop_F2[fig_idex]),
          Factor_F2_Steps = as.character(Graphs$Steps_F2[fig_idex])
        ),
        flat_params
      )
      meta_accumulator$rows[[length(meta_accumulator$rows) + 1]] <- meta_row
    }

    # Here we loop over the possible parameters steps
    # For the most basic allele figure this is run once
    # For the range fig i varies between some values
    # For SA figs combinations of i and j are tested
    
    i <- Graphs$i[fig_idex]
    j <- Graphs$j[fig_idex]
    
    
    # These names get appended to the file names
    # We reset them here so they can change based on the factors varied
    Factor1_name <- ""
    Factor2_name <- ""
    
    # Here check if factors are set for this figure, if so the parameters are adjusted accordingly for range and SA figs
    if (Graphs$Factor_F1[fig_idex] != "none") {
      Parameters <-
        Parameter_Change(Parameters,
                         Graphs$Factor_F1_value[fig_idex],
                         Graphs$Factor_F1[fig_idex],
                         Graphs[fig_idex,])
      
    }
    if (Graphs$Factor_F2[fig_idex] != "none") {
      Parameters <-
        Parameter_Change(Parameters,
                         Graphs$Factor_F2_value[fig_idex],
                         Graphs$Factor_F2[fig_idex],
                         Graphs[fig_idex,])
      
    }
    
    
    
    # Write per-simulation summary table (after sweep adjustment)
    Outcome_Folder <- paste0("Outcomes/", Graphs[fig_idex, ]$Project_Folder[1], "/")
    dir.create(Outcome_Folder, recursive = TRUE, showWarnings = FALSE)
    summary_df <- generate_summary_table(Parameters, Graphs[fig_idex, ])
    write.csv(summary_df,
              paste0(Outcome_Folder, "Summary_", Graphs$Fig_long_name[fig_idex], ".csv"),
              row.names = FALSE)

    # Here we use the default parameters to create the conversion matrices
    Model_Matrices <- All_Main_Calculations(Parameters)
    
    for (r in 0:Parameters$Release_File$Replicate_Simulations[1]) {
      
      message(paste0(r,"/",Parameters$Release_File$Replicate_Simulations[1]))
      
    Main_Outcomes <-
      Main_Model(Model_Matrices, r) # Single run of the model

        Outcomes_by_Allele_List <-
      Outcomes_by_Allele(Main_Outcomes, Parameters, r) # Here we reformat the data to be per allele (instead of a vector of all separate genotype)
    
    
    ## 1. Identify the time‑point columns (works no matter how many you have)
    tp_cols <- grep("^V\\d+$", names(Main_Outcomes$Single_Outcomes), value = TRUE)
    
    ## 2. Reshape, duplicate into freq columns, and give missing metadata defaults
    single_long <- Main_Outcomes$Single_Outcomes %>% 
      pivot_longer(cols = all_of(tp_cols),
                   names_to  = "t",
                   values_to = "value",
                   values_drop_na = FALSE) %>% 
      mutate(
        t = parse_number(t),                    # "V1" → 1, keeps numeric order
        Allele_Freq  = value,                   # copy into the four freq fields
        Carrier_Freq = value,
        Het_Freq     = value,
        Hom_Freq     = value,
        Pre_or_Post_Fit = Pre_or_Post,          # rename to match MF table
        Replicate   = r,                       # or whatever ID you use
        Factor_F1   = 0,
        Factor_F2   = 1
      ) %>% 
      select(                                    # put columns in the same order
        t, Allele_Freq, Carrier_Freq, Het_Freq, Hom_Freq,
        Outcome, Locus, Sex, Pre_or_Post_Fit,
        Replicate, Factor_F1, Factor_F2
      )
    
    ## 3. Append to the master table
    Outcomes_by_Allele_List$Outcomes_MF <- bind_rows(
      Outcomes_by_Allele_List$Outcomes_MF,
      single_long
    )
    
    Outcomes_by_Allele_List$Outcomes_MF$Factor_F1 <-
      Graphs$Factor_F1_value[fig_idex]
    Outcomes_by_Allele_List$Outcomes_MF$Factor_F2 <-
      Graphs$Factor_F2_value[fig_idex]
    
    # Here we store the outcomes per genotype data
    
    Project_Folder <- paste0("Outcomes/", Graphs[fig_idex,]$Project_Folder[1],"/")
    
    if (produce_conversion_table_files) {
      
      # # Simplify the matrix
      # simplified_matrix <- simplify_conversion_matrix(Model_Matrices$D, Model_Matrices$Alpha_M)
      
      Alpha_Allele_Names <- do.call(paste, as.data.frame(Model_Matrices$Alpha_M[1:11], stringsAsFactors=FALSE))
      
      Matrix_row_column_names <- function(Matrix,Alpha_Allele_Names) {
      Matrix <- cbind(Alpha_Allele_Names, as.data.frame(as.matrix(Matrix)))
      names(Matrix) <- c("Alleles",Alpha_Allele_Names)
      row.names(Matrix) <- Alpha_Allele_Names
      return(Matrix)
      }
      
      ### Changes to the default nuc rates go here ###

      write.csv(
        t(rbind(Model_Matrices$geno_list$A,Model_Matrices$geno_listB,Model_Matrices$geno_list$C)),
        paste0(Project_Folder, "Allele_", Graphs$Fig_long_name[fig_idex], ".csv")
      )
      
      write.csv(
        Model_Matrices$nuc_rates_list,
        paste0(Project_Folder, "Nuc_", Graphs$Fig_long_name[fig_idex], ".csv")
      )
      
      write.csv(
        Model_Matrices$Release_File,
        paste0(Project_Folder, "Release_", Graphs$Fig_long_name[fig_idex], ".csv")
      )
      
      
      Experimental_Data = paste0("Parameter Files/", Graphs$Project_Folder[fig_idex],"/", 'Experimental_Data', ".csv")

      if (file.exists(Experimental_Data)) {
        write.csv(
          read.csv(Experimental_Data, header = TRUE),
          paste0(Project_Folder, "Experimental_", Graphs$Fig_long_name[fig_idex], ".csv")
        )
      }
      
      # Here we store the conversion matrices
      write.csv(
        Matrix_row_column_names(Model_Matrices$D, Alpha_Allele_Names),
        paste0(Project_Folder, "Dep_", Graphs$Fig_long_name[fig_idex], ".csv")
      )

      write.csv(
        Matrix_row_column_names(Model_Matrices$exp_conv_male, Alpha_Allele_Names),
        paste0(Project_Folder, "ExpM_", Graphs$Fig_long_name[fig_idex], ".csv")
      )
      write.csv(
        Matrix_row_column_names(Model_Matrices$exp_conv_female, Alpha_Allele_Names),
        paste0(Project_Folder, "ExpF_", Graphs$Fig_long_name[fig_idex], ".csv")
      )
      write.csv(
        cbind(Model_Matrices$Alpha_M, as.matrix(Model_Matrices$H_M)),
        paste0(Project_Folder, "HaploidM_", Graphs$Fig_long_name[fig_idex], ".csv")
      )
      write.csv(
        cbind(Model_Matrices$Alpha_F, as.matrix(Model_Matrices$H_F)),
        paste0(Project_Folder, "HaploidF_", Graphs$Fig_long_name[fig_idex], ".csv")
      )
      
      # apply(Model_Matrices$Fit_Male, 1, prod)
      write.csv(
        cbind(Model_Matrices$Alpha_M, as.matrix(Model_Matrices$Fit_Male)),
        paste0(Project_Folder, "FitnessM_", Graphs$Fig_long_name[fig_idex], ".csv")
      )
      write.csv(
        cbind(Model_Matrices$Alpha_F, as.matrix(Model_Matrices$Fit_Female)),
        paste0(Project_Folder, "FitnessF_", Graphs$Fig_long_name[fig_idex], ".csv")
      )
    }
    
    if (r == 0) {
      # Here we store the outcomes by allele data
      write.csv(
        Outcomes_by_Allele_List$Outcomes_MF,
        paste0(Project_Folder, "Outcomes_", Graphs$Fig_long_name[fig_idex], ".csv")
        )

      # Save sex-specific raw genotype frequencies from deterministic run for migration source use
      off <- 12  # genotype descriptor columns
      for (sex_label in c("M", "F")) {
        outcome_key <- paste0("Outcomes_", sex_label, "_PostFit")
        geno_freq <- as.matrix(Main_Outcomes[[outcome_key]][, (off + 1):ncol(Main_Outcomes[[outcome_key]])])
        geno_freq <- geno_freq[, !apply(is.na(geno_freq), 2, all), drop = FALSE]
        colnames(geno_freq) <- paste0("V", 1:ncol(geno_freq))

        write.csv(
          geno_freq,
          paste0(Project_Folder, "Genotype_Freq_", Graphs$Fig_long_name[fig_idex], "_", sex_label, ".csv"),
          row.names = FALSE
        )
      }
    } else {
      
      write.table(Outcomes_by_Allele_List$Outcomes_MF,
        paste0(Project_Folder, "Outcomes_", Graphs$Fig_long_name[fig_idex], ".csv"),
        sep = ",", append = T, col.names = FALSE)
    }
  
    }
  } else {
    print("Allele outcome file exists, skipping")
  }
}

# Generate and write meta-tables and profile lookups
meta_tables <- generate_meta_tables(meta_accumulator)
if (nrow(meta_tables$profiles) > 0) {
  write.csv(meta_tables$profiles, paste0("Outcomes/Meta_Table.csv"), row.names = FALSE)
  write.csv(meta_tables$detailed, paste0("Outcomes/Meta_Table_Detailed.csv"), row.names = FALSE)
  print(paste("Meta tables written with", nrow(meta_tables$profiles), "simulation groups"))

  for (type in c("design", "fitness", "nuc", "release")) {
    registry <- meta_accumulator[[paste0("profile_", type)]]
    lookup_df <- generate_profile_lookups(registry, type)
    if (nrow(lookup_df) > 0) {
      write.csv(lookup_df,
                paste0("Outcomes/Profile_Lookup_", tools::toTitleCase(type), ".csv"),
                row.names = FALSE)
    }
  }
  print("Profile lookup tables written")
}

print("Done with simulations.")
}
