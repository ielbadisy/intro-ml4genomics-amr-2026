# Build the clean, isolate-level dataset used by lab-tutorial/ from the raw
# BV-BRC downloads in raw/. One row per isolate, phenotype for 8 antibiotics
# (incl. ciprofloxacin, the future ML target) as columns.
#
# By default this re-derives the table from the cached raw files (raw/*.csv),
# which is fast and fully reproducible (fixed seed). Set REFRESH_DATA <- TRUE
# to re-query the BV-BRC API instead (requires network access).
#
# Run from data/: Rscript build_clean_data.R

library(httr2)
library(jsonlite)
library(basetable)

compact <- function(x) Filter(Negate(is.null), x)

REFRESH_DATA <- FALSE

bvbrc_base <- "https://www.bv-brc.org/api"

antibiotic_panel <- c(
  "ciprofloxacin", "ceftazidime", "piperacillin/tazobactam", "cefotaxime",
  "gentamicin", "meropenem", "ampicillin", "trimethoprim/sulfamethoxazole"
)
ab_slug <- c(
  "ciprofloxacin"                 = "ciprofloxacin",
  "ceftazidime"                   = "ceftazidime",
  "piperacillin/tazobactam"       = "pip_taz",
  "cefotaxime"                    = "cefotaxime",
  "gentamicin"                    = "gentamicin",
  "meropenem"                     = "meropenem",
  "ampicillin"                    = "ampicillin",
  "trimethoprim/sulfamethoxazole" = "sxt"
)

amr_fields <- c(
  "genome_id", "genome_name", "antibiotic", "resistant_phenotype",
  "measurement", "measurement_value", "measurement_unit", "measurement_sign",
  "laboratory_typing_method", "testing_standard", "testing_standard_year"
)

## ---- Step 1: AMR phenotypes (long format, one row per isolate x antibiotic)

make_amr_request <- function(ab, n_max = 3000) {
  q <- paste(
    "eq(taxon_id,562)",
    "eq(evidence,Laboratory%20Method)",
    sprintf("eq(antibiotic,%s)", utils::URLencode(ab, reserved = TRUE)),
    sprintf("select(%s)", paste(amr_fields, collapse = ",")),
    sprintf("limit(%d)", n_max),
    sep = "&"
  )
  url <- paste0(bvbrc_base, "/genome_amr/?", q)
  request(url) |> req_headers(accept = "application/json") |> req_retry(max_tries = 3)
}

parse_amr_response <- function(resp) {
  txt <- resp_body_string(resp)
  if (nchar(trimws(txt)) == 0 || txt == "[]") return(NULL)
  fromJSON(txt, flatten = TRUE) |> as.data.frame(stringsAsFactors = FALSE)
}

raw_path <- "raw/ecoli_amr_raw.csv"

if (!REFRESH_DATA && file.exists(raw_path)) {
  ecoli_amr_raw <- read.csv(raw_path, stringsAsFactors = FALSE, colClasses = c(genome_id = "character"))
} else {
  amr_requests   <- map(antibiotic_panel, make_amr_request, n_max = 3000)
  amr_responses  <- req_perform_parallel(amr_requests, on_error = "continue", progress = FALSE)
  ecoli_amr_raw  <- amr_responses |> map(parse_amr_response) |> compact() |> rbindfill()
  dir.create("raw", showWarnings = FALSE)
  write.csv(ecoli_amr_raw, raw_path, row.names = FALSE)
}

## ---- Step 2: clean, dedupe, select the isolate sample

amr_clean <- ecoli_amr_raw |>
  transform(antibiotic_slug = ab_slug[antibiotic]) |>
  subset(resistant_phenotype %in% c("Resistant", "Intermediate", "Susceptible")) |>
  orderrows(by = c("genome_id", "antibiotic_slug", "testing_standard_year"), decreasing = c(FALSE, FALSE, TRUE)) |>
  uniquerows(cols = c("genome_id", "antibiotic_slug"), .keep_all = TRUE)

coverage <- count(amr_clean, by = "genome_id", name = "n_antibiotics_tested")

has_cip <- amr_clean |>
  subset(antibiotic_slug == "ciprofloxacin") |>
  (\(d) d$genome_id)()

eligible_ids <- coverage |>
  subset(n_antibiotics_tested >= 3 & genome_id %in% has_cip) |>
  (\(d) d$genome_id)()

# Sorted explicitly: `count()`'s output order (and, before basetable 1.4.1,
# the tie-break among genome_ids with an equal n_antibiotics_tested) is an
# implementation detail of the grouping engine, not a documented contract.
# Fixing the order here means the sample drawn below is reproducible given
# the seed, regardless of which basetable version produced `eligible_ids`.
eligible_ids <- sort(eligible_ids)

set.seed(42)
target_n   <- min(700, length(eligible_ids))
sample_ids <- sort(sample(eligible_ids, target_n))

## ---- Step 3: genomic/collection metadata for the sampled isolates

meta_fields <- c(
  "genome_id", "genome_name", "collection_year", "isolation_country",
  "isolation_source", "host_name", "host_group", "genome_length",
  "gc_content", "mlst", "contigs", "genome_quality"
)

make_metadata_request <- function(ids) {
  ids_str <- paste(ids, collapse = ",")
  q <- paste(
    sprintf("in(genome_id,(%s))", ids_str),
    sprintf("select(%s)", paste(meta_fields, collapse = ",")),
    sprintf("limit(%d)", length(ids) + 10),
    sep = "&"
  )
  url <- paste0(bvbrc_base, "/genome/?", q)
  request(url) |> req_headers(accept = "application/json") |> req_retry(max_tries = 3)
}

meta_path <- "raw/ecoli_metadata_raw.csv"

if (!REFRESH_DATA && file.exists(meta_path)) {
  ecoli_meta_raw <- read.csv(meta_path, stringsAsFactors = FALSE, colClasses = c(genome_id = "character"))
} else {
  id_batches     <- base::split(sample_ids, ceiling(seq_along(sample_ids) / 150))
  meta_requests  <- map(id_batches, make_metadata_request)
  meta_responses <- req_perform_parallel(meta_requests, on_error = "continue", progress = FALSE)
  ecoli_meta_raw <- meta_responses |> map(parse_amr_response) |> compact() |> rbindfill()
  dir.create("raw", showWarnings = FALSE)
  write.csv(ecoli_meta_raw, meta_path, row.names = FALSE)
}

## ---- Step 4: long -> wide (one row per isolate), join, save

amr_sub <- subset(amr_clean, genome_id %in% sample_ids)

amr_wide <- towide(
  amr_sub, names = "antibiotic_slug", values = "resistant_phenotype",
  idcols = "genome_id", fun = function(x) x[1]
)

ecoli_amr_isolate <- ecoli_meta_raw |>
  uniquerows(cols = "genome_id", .keep_all = TRUE) |>
  merge(amr_wide, by = "genome_id", all.x = TRUE)

out_path <- "ecoli_amr_isolate_v2.csv"
write.csv(ecoli_amr_isolate, out_path, row.names = FALSE)

cat(sprintf("Wrote %s: %d isolates x %d columns\n", out_path, nrow(ecoli_amr_isolate), ncol(ecoli_amr_isolate)))

if ("package:basetable" %in% search()) {
  detach("package:basetable", unload = TRUE, force = TRUE)
}
