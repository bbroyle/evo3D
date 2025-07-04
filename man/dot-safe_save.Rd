% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/WRAPPER_MODULE.R
\name{.safe_save}
\alias{.safe_save}
\title{Safely Resolve a Non-Clashing Output Path}
\usage{
.safe_save(path, is_dir = FALSE)
}
\arguments{
\item{path}{Character. Desired path.}

\item{is_dir}{Logical. Is it a directory?}

\item{max_tries}{Integer. How many numbered attempts before fallback to POSIXct timestamp}

\item{tag}{Character. Optional tag to append to the path.}
}
\value{
A list with:
  \item{path}{Safe file path}
  \item{systime}{Used timestamp (in fallback or NULL if not used)}
}
\description{
Avoids overwriting files/directories by appending numeric suffixes.
If no safe name is found, uses a fallback with a consistent timestamp.
}
\keyword{internal}
