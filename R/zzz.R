.onAttach <- function(libname, pkgname) {
  if (!("IMPRINTS.CETSA" %in% installed.packages())) {
    packageStartupMessage(
      "Please install the last version of `IMPRINTS.CETSA` with",
      " `devtools::install_github('nkdailingyun/IMPRINTS.CETSA')`",
      " or by procuring you the tar.gz file."
    )
  }
  if (!("IMPRINTS.CETSA.app" %in% installed.packages())) {
    packageStartupMessage(
      "Please install the last version of `IMPRINTS.CETSA.app` with",
      " `devtools::install_github('mgerault/IMPRINTS.CETSA.app')`",
      " or by procuring you the tar.gz file."
    )
  }

  packageStartupMessage(
    "\n",
    "Welcome to IMPRINTS.PhosphoQP package! \n",
    "To access the documentation type browseVignettes(package = 'IMPRINTS.PhosphoQP') or \n",
    "with vignette('IMPRINTS.PhosphoQP_doc', package = 'IMPRINTS.PhosphoQP')"
    )
}
