library(data.table)
library(qqman)
library(htmltools)

manhattan_plot <- function(args) {
  
  inputs <- args[1:(length(args)-2)]
  output1 <- args[length(args)-1]
  output2 <- args[length(args)]
  
  inputs <- args[1:(length(args) - 2)]
  output_merged <- args[length(args) - 1]
  output_report <- args[length(args)]

  all_tables <- lapply(inputs, function(file) {
    basename_file <- basename(file)  # quita el path, solo el nombre
    prefix <- sub("_association.*\\.assoc$", "", basename_file)
    
    dt <- fread(file, select = c("SNP", "BP", "CHR", "P"))
    dt[, SNP_clean := sub(paste0("^", prefix, "_"), "", SNP)]
    cols_to_rename <- setdiff(names(dt), c("SNP", "SNP_clean", "BP", "CHR"))
    setnames(dt, cols_to_rename, paste0(cols_to_rename, "_", prefix))
    dt[, SNP := SNP_clean][, SNP_clean := NULL]
    return(dt)
  })

  merged_dt <- Reduce(function(x, y) merge(x, y, by = c("SNP", "BP", "CHR"), all = TRUE), all_tables)

  pvalue_cols <- grep("^P_", names(merged_dt), value = TRUE)
  logP_dt <- merged_dt[, lapply(.SD, function(x) -log10(as.numeric(x))), .SDcols = pvalue_cols]
  merged_dt[, mean_logP := rowMeans(logP_dt, na.rm = TRUE)]
  merged_dt[, combined_mean_P := 10^(-mean_logP)]
  merged_dt[, c("mean_logP") := NULL]

  fwrite(merged_dt, output_merged, sep = "\t")

  top_hits <- merged_dt[order(combined_mean_P)][1:30, c("CHR", "BP", "combined_mean_P", pvalue_cols), with=FALSE]

  top_hits_transformed <- copy(top_hits)
  for (col in c("combined_mean_P", pvalue_cols)) {
    new_col <- paste0("logP_", col)
    top_hits_transformed[, (new_col) := fifelse(is.na(get(col)), NA_real_, -log10(as.numeric(get(col))))]
  }
  cols_to_show <- c("CHR", "BP", paste0("logP_", c("combined_mean_P", pvalue_cols)))
  final_table <- top_hits_transformed[, ..cols_to_show]
  setorder(final_table, CHR, BP)

  html_table <- tags$table(
    style = "border-collapse: collapse; width: 90%; margin: auto;",
    tags$caption(
      style = "font-weight: bold; font-size: 1.2em; margin-bottom: 10px;",
      "Top 30 SNPs (-log10 P-values)"
    ),
    tags$thead(
      tags$tr(
        lapply(names(final_table), function(colname) {
          tags$th(style = "border: 1px solid #999; padding: 8px; background-color: #f2f2f2;", colname)
        })
      )
    ),
    tags$tbody(
      lapply(1:nrow(final_table), function(i) {
        bg_color <- if (i %% 2 == 0) "#ffffff" else "#f9f9f9"
        tags$tr(
          style = paste0("background-color: ", bg_color, ";"),
          lapply(final_table[i], function(val) {
            value <- ifelse(is.na(val), "-", format(val, digits=4))
            tags$td(style = "border: 1px solid #999; padding: 8px; text-align: center;", value)
          })
        )
      })
    )
  )

  tmp_plot <- tempfile(fileext = ".png")
  png(tmp_plot, width = 1400, height = 800)
  manhattan(merged_dt[, .(CHR, BP, SNP, P = combined_mean_P)],
            suggestiveline = -log10(1e-5),
            genomewideline = -log10(5e-8),
            col = c("grey30", "skyblue3"),
            cex = 0.6,
            main = "Manhattan Plot (Combined)")
  dev.off()

  intro_section <- tags$div(
    style = "width:80%; margin:auto; font-size:1.1em; line-height:1.5;",
    tags$p("This report summarizes the combined results of multiple genome-wide association studies (GWAS) performed across several datasets."),
    tags$p("P-values from individual analyses were integrated by computing the mean of their -log10-transformed values to obtain an aggregated measure of association."),
    tags$p("The following sections present a Manhattan plot of the aggregated results and a table listing the most significant variants identified.")
  )

  manhattan_caption <- tags$p(
    style = "width:80%; margin:auto; font-size:1em; line-height:1.5;",
    "The Manhattan plot displays the -log10 transformed combined p-values along the genome. Each point represents a single nucleotide polymorphism (SNP), plotted according to its chromosomal position. Alternating colors distinguish different chromosomes."
  )

  table_caption <- tags$p(
    style = "width:80%; margin:auto; font-size:1em; line-height:1.5; margin-top:20px;",
    "The table below lists the top 30 SNPs based on the aggregated combined p-values. For each SNP, its chromosomal position is indicated, along with the -log10 transformed p-values from each individual dataset as well as the combined value. Missing values are represented by '-'."
  )

  report <- tags$html(
    tags$head(
      tags$title("GWAS Mini Report"),
      tags$style("body { font-family: Arial, sans-serif; }")
    ),
    tags$body(
      tags$h1("GWAS Combined Analysis", style = "text-align: center;"),
      intro_section,
      tags$h2("Manhattan Plot", style = "text-align: center;"),
      manhattan_caption,
      tags$div(
        style = "text-align: center;",
        tags$img(src=basename(tmp_plot), style="max-width:90%; height:auto;")
      ),
      tags$h2("Top 30 SNPs", style = "text-align: center; margin-top: 40px;"),
      table_caption,
      html_table
    )
  )

  report_dir <- dirname(output_report)
  file.copy(tmp_plot, file.path(report_dir, basename(tmp_plot)), overwrite = TRUE)

  save_html(report, file = output_report)
  
}

#run as Rscript
args <- commandArgs(trailingOnly = TRUE)
manhattan_plot(args)