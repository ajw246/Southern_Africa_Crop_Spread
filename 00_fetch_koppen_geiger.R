# Download and verify the climate grid archive

ZIP_URL  <- "https://ndownloader.figshare.com/files/61012822"
ZIP_MD5  <- "7fc2f5a15d4f5fe0ce59c9a9b502aa09"
ZIP_SIZE <- 130618411
DEST    <- file.path("data", "koppen_geiger")
PERIODS <- c("1901_1930", "1931_1960", "1961_1990", "1991_2020")
GRID    <- "koppen_geiger_0p00833333.tif"
EXPECTED_SHA256 <- c(
  "1901_1930" = "bcce8f83c6c850a34705ff49dbbeb19baae2a409ca2492e18e30ace1786893ef",
  "1931_1960" = "cc71784dcb61635a4913cb42111bdf98ab7fba3157b5f9f3e297da633f3e33d3",
  "1961_1990" = "ccdc6c70b20a62eb09b9dac715c136c48a7d3721e2fccaf19b52ad282916fefd",
  "1991_2020" = "2130f0071dfb2904947d8ec3a0d807fac71004df76e769262004f1602e4d6a13",
  "legend"    = "2ede2ad270a036cc11c31705a2c1dbf0314a8cf011fc972cd4a9665e3339e5e5"
)

# Check that the repository data folder exists

if (!dir.exists("data")) {
  stop("No data/ folder here.\n",
       "Set R's working directory to the repository root - the folder that\n",
       "contains this script, data/ and figures/.", call. = FALSE)
}

# Skip the download when all required files are present

wanted <- c(file.path(DEST, PERIODS, GRID), file.path(DEST, "legend.txt"))
if (all(file.exists(wanted))) {
  message("The Koppen-Geiger grids are already in ", DEST, " - nothing to do.")
  message("To force a fresh download, delete that folder and run this again.")
} else {

# Download, verify and unpack the required climate grids

  old_timeout <- getOption("timeout")
  on.exit(options(timeout = old_timeout), add = TRUE)
  options(timeout = 1800)
  zip_path <- file.path(tempdir(), "koppen_geiger_tif.zip")
  on.exit(unlink(zip_path), add = TRUE)
  message("Downloading the Koppen-Geiger archive (about ",
          round(ZIP_SIZE / 1024^2), " MB).")
  message("This is a one-off; it usually takes a few minutes.")
  ok <- try(
    utils::download.file(ZIP_URL, zip_path, mode = "wb", quiet = FALSE),
    silent = TRUE
  )
  if (inherits(ok, "try-error") || !file.exists(zip_path)) {
    stop("The download failed.\n",
         "Check that this computer is online and that figshare.com is\n",
         "reachable, then run this script again.  The file can also be\n",
         "fetched by hand from https://doi.org/10.6084/m9.figshare.21789074",
         call. = FALSE)
  }
  got_md5 <- unname(tools::md5sum(zip_path))
  if (!identical(got_md5, ZIP_MD5)) {
    stop("The downloaded archive does not match its published checksum.\n",
         "  expected MD5: ", ZIP_MD5, "\n",
         "  got MD5:      ", got_md5, "\n",
         "The download may have been truncated - try again.  If it keeps\n",
         "failing, Figshare may have published a new version; see\n",
         "DATA_SOURCES.md for the version this code was built against.",
         call. = FALSE)
  }
  message("Checksum matches the published archive.")
  unzip_dir <- file.path(tempdir(), "kg_unzip")
  unlink(unzip_dir, recursive = TRUE)
  dir.create(unzip_dir, showWarnings = FALSE)
  on.exit(unlink(unzip_dir, recursive = TRUE), add = TRUE)
  contents <- utils::unzip(zip_path, list = TRUE)$Name
  find_one <- function(pattern, what) {
    hits <- grep(pattern, contents, value = TRUE)
    if (length(hits) != 1) {
      stop("Expected exactly one ", what, " in the archive, found ",
           length(hits), ".\n",
           "The archive layout may have changed.  Unpack it by hand and copy\n",
           "the four historical period folders of the 1 km (0p00833333)\n",
           "GeoTIFF set, plus legend.txt, into ", DEST, "/.",
           call. = FALSE)
    }
    hits
  }
  targets <- vapply(
    PERIODS,
    function(p) find_one(paste0("(^|/)", p, "/", GRID, "$"),
                         paste("grid for", p)),
    character(1)
  )
  legend_in_zip <- find_one("(^|/)legend\\.txt$", "legend.txt")
  message("Unpacking ", length(targets) + 1, " files.")
  utils::unzip(zip_path, files = c(targets, legend_in_zip), exdir = unzip_dir)
  for (p in PERIODS) {
    dir.create(file.path(DEST, p), recursive = TRUE, showWarnings = FALSE)
    file.copy(file.path(unzip_dir, targets[[p]]),
              file.path(DEST, p, GRID), overwrite = TRUE)
  }
  file.copy(file.path(unzip_dir, legend_in_zip),
            file.path(DEST, "legend.txt"), overwrite = TRUE)
  missing <- wanted[!file.exists(wanted)]
  if (length(missing)) {
    stop("These files did not make it out of the archive:\n  ",
         paste(missing, collapse = "\n  "), call. = FALSE)
  }
  if (requireNamespace("digest", quietly = TRUE)) {
    got <- c(
      vapply(PERIODS, function(p)
        digest::digest(file.path(DEST, p, GRID), algo = "sha256", file = TRUE),
        character(1)),
      legend = digest::digest(file.path(DEST, "legend.txt"),
                              algo = "sha256", file = TRUE)
    )
    bad <- names(got)[got != EXPECTED_SHA256[names(got)]]
    if (length(bad)) {
      warning("These files do not match the SHA-256 recorded in ",
              "DATA_SOURCES.md: ", paste(bad, collapse = ", "),
              ". Figshare may have published a new version.", call. = FALSE)
    } else {
      message("All five files match the SHA-256 recorded in DATA_SOURCES.md.")
    }
  }
  written <- sum(file.size(wanted))
  message("Done. ", round(written / 1024^2), " MB written to ", DEST, "/")
  message("01_data_wrangling.R and Figures 5, S1 and S3 can now be run.")
}
