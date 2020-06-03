#############################################
#  Functions for LipidView_analysis script  #
#############################################

lipid_plot <- function(lipid_var, profile_var) {
  profile_var <- sym(profile_var)

  data %>%
    filter(class == lipid_var) %>%
    mutate(
      parsed_species = str_remove(species, lipid_var),
      parsed_species = str_remove(parsed_species, "O-"),
      parsed_species = str_split(parsed_species, ":|;"),
      chain_length = map_chr(parsed_species, 1),
      dbonds = map_chr(parsed_species, 2),
      oh = ifelse(str_detect(species, ";"), map_chr(parsed_species, 3), NA)
    ) %>%
    group_by(samples) %>%
    mutate(uM_sum = sum(uM)) %>%
    mutate(mol_percent = uM * 100 / uM_sum) %>%
    group_by(!!profile_var, samples) %>%
    summarise(mol_percent = sum(mol_percent)) %>%
    ggplot(aes_string(
      x = profile_var,
      y = "mol_percent",
      fill = "samples"
    )) +
    geom_col(position = "dodge")
}

####

lipid_barplot <- function(data, lipid_var, profile_var) {
  
  profile_var <- sym(profile_var)
  
  blacklist <- c("class", "CA", "Chol", "SAA", "Ergo")
  
  plot_data <- data %>%
    filter(if (lipid_var == "class") {class == class} else {class == lipid_var}) %>%
    mutate(species = if (!(lipid_var %in% blacklist)) {str_remove(species, lipid_var)} else {species},
           parsed_species = if (!(lipid_var %in% blacklist)) {str_remove(species, "O-")} else {NA},
           parsed_species = if (!(lipid_var %in% blacklist)) {str_split(parsed_species, ":|;")} else {NA},
           chain_length = if (!(lipid_var %in% blacklist)) {map_chr(parsed_species, 1)} else {NA},
           dbonds = if (!(lipid_var %in% blacklist)) {map_chr(parsed_species, 2)} else {NA},
           oh = if (!(lipid_var %in% blacklist)) {if (any(str_detect(species, ";"))) {map_chr(parsed_species, 3)} else {NA}} else {NA}) %>%
    group_by(samples) %>%
    mutate(uM_sum = sum(uM)) %>%
    mutate(mol_percent = uM * 100 / uM_sum)
  
  
  if (profile_var == sym("species") & nrow(unique(plot_data[,profile_var])) > 1) {
    species_order <- fct_inorder(unique(plot_data$species), ordered = NA)
  }
  
  plot_data <- group_by(plot_data, !!profile_var, samples) %>%
    summarise(mol_percent = sum(mol_percent), .groups = "drop")
  
  ggplot(plot_data, aes_string(x = profile_var, y = "mol_percent", fill = "samples")) +
    geom_col(position = "dodge") +
    ggtitle(paste(lipid_var, profile_var, "profile")) +
    theme(axis.text.x = element_text(angle = 45),
          plot.title = element_text(hjust = 0.5)) +
    if (profile_var == sym("species") & nrow(unique(plot_data[,profile_var])) > 1) {
      scale_x_discrete(limits = species_order)
    }
}

####

write_data_to_excelsheet <- function(df1, df2, df3) {
  
  writeData(wb, f, as.character(files[f]), startCol = 1, startRow = 1)
  writeData(wb, f, df1, startCol = 2, startRow = 2, rowNames = FALSE)
  if (any(paste0("e", f, " [uM]") == colnames(df2))) {
  writeData(wb, f, df2[, -1], startCol = (2 + ncol(df1) + 1), startRow = 2, rowNames = FALSE)
  }
  writeData(wb, f, as.character(paste0(f, " mol%")), startCol = 1, startRow = (nrow(df1) + 5))
  writeData(wb, f, df3, startCol = 1, startRow = (nrow(df1) + 6))
  insertPlot(wb, f, width = 15, height = 5, fileType = "png", units = "in",startRow = ((2*nrow(df1))+9))
}

####