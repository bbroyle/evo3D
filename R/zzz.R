.onAttach <- function(libname, pkgname) {
  if (!requireNamespace("msa", quietly = TRUE)) {
    packageStartupMessage(
      "Optional package 'msa' is not installed.\n",
      "For MSA alignment features, install it with:\n",
      "  BiocManager::install('msa')"
    )
  }
}