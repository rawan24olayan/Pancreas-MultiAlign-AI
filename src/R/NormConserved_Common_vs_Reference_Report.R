# Load Libraries (Install if missing: install.packages(c("ggplot2", "dplyr", "gridExtra")))
library(ggplot2)
library(dplyr)
library(gridExtra)

# 1. Load Agent Findings
report_path <- "results/final_research_report.csv"
if (!file.exists(report_path)) {
  stop("❌ Error: final_research_report.csv not found. Run the Agent first!")
}
df <- read.csv(report_path)

# 2. Visual: Disease Association Bar Plot
p1 <- ggplot(df, aes(x = reorder(Gene, Confidence), y = Confidence, fill = Confidence)) +
  geom_bar(stat = "identity", width = 0.7) +
  coord_flip() +
  scale_fill_gradient(low = "#3498db", high = "#e74c3c") +
  labs(title = "BENGAL Phase 3: Clinical Validation",
       subtitle = "UK Biobank / Open Targets Association Scores",
       x = "Conserved Orthologs (Human/Mouse)",
       y = "Association Confidence Score") +
  theme_minimal()

# 3. Visual: Data Table
table_theme <- ttheme_minimal(
  core = list(bg_params = list(fill = c("white", "#f9f9f9"), col = NA)),
  colhead = list(fg_params = list(col = "white"), bg_params = list(fill = "#2c3e50"))
)
p2 <- tableGrob(df, rows = NULL, theme = table_theme)

# 4. Save to PDF (Standard Bio-Journal Format)
pdf("results/BENGAL_Discovery_Report.pdf", width = 8.5, height = 11)
grid.arrange(p1, p2, ncol = 1, heights = c(1, 1))
dev.off()

cat("✨ Success: BENGAL_Discovery_Report.pdf generated in results/\n")
