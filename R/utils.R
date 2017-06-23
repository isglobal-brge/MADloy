doc_date <- function() {
    format(Sys.Date(), "%d %B %Y")
}

pkg_ver <- function(name) {
    paste(name, utils::packageVersion(name))
}

# ! Return the samples names in a MADloy object
samples <- function(x) {
    x$par$files
}
