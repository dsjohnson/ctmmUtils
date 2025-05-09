% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/as_telem.R
\name{as_telem}
\alias{as_telem}
\title{Convert Wildlife Computers data imported with
`\link[ctmmUtils]{read_wc_dirs}` to a `telemetry` object from the
`ctmm` package.}
\usage{
as_telem(x, ...)
}
\arguments{
\item{x}{An sf data frame output by the function `\link[ctmmUtils]{read_wc_dirs}`.}

\item{...}{Additional arguments to be passed to `\link[ctmm]{as.telemetry}`}
}
\description{
Convert Wildlife Computers data imported with
`\link[ctmmUtils]{read_wc_dirs}` to a `telemetry` object from the
`ctmm` package.
}
\author{
Josh M. London, Devin S. Johnson
}
