# Converts the RetiGene gene atlas (Table S1 of Am J Hum Genet 2025;112(10):2253-2265,
# doi:10.1016/j.ajhg.2025.08.017, https://retigene.erdc.info) into PriorR panels:
#   Rscript Program_files/scripts/import_retigene.R <TableS1.xlsx>     (run from the app/ folder)
# Writes Program_files/Panels/retigene.tsv: retained genes as "green" (confidence 3) and
# candidate genes as "amber" (confidence 2); discarded genes are left out. Needs readxl.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: import_retigene.R <TableS1.xlsx>")

t <- readxl::read_excel(args[1], skip = 2, col_types = "text", .name_repair = "minimal")
t <- t[!is.na(t$Gene) & t$Category != "Discarded", ]
retained <- t$Retained == "Yes"
confidence <- ifelse(retained, "3", "2")
category_cols <- grep("^Category [0-9]$", names(t), value = TRUE)

panels <- list()
add <- function(id, name, rows) {
  panels[[id]] <<- data.frame(panel_id = id, panel_name = name, version = "2025",
                              gene = t$Gene[rows], confidence = confidence[rows],
                              mode_of_inheritance = t$Inheritance[rows], stringsAsFactors = FALSE)
}

add("all", "RetiGene - all inherited retinal disease genes", rep(TRUE, nrow(t)))
add("non_syndromic", "RetiGene - non-syndromic", retained & t$Category %in% c("Non-syndromic", "Both"))
add("syndromic", "RetiGene - syndromic", retained & t$Category %in% c("Syndromic", "Both"))
subsets <- sort(unique(na.omit(unlist(t[retained, category_cols]))))
for (s in subsets) {
  rows <- retained & apply(t[, category_cols], 1, function(r) s %in% r)
  add(gsub("[^A-Za-z0-9]+", "_", tolower(s)), paste("RetiGene -", s), rows)
}

out <- do.call(rbind, unname(panels))
data.table::fwrite(out, file.path("Program_files", "Panels", "retigene.tsv"), sep = "\t")
message(length(panels), " RetiGene panels: ", sum(retained), " retained genes, ", sum(!retained), " candidates")
