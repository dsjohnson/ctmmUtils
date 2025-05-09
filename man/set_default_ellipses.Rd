% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/set_default_ellipses.R
\name{set_default_ellipses}
\alias{set_default_ellipses}
\title{Set Default Values for MIssing ARGOS Diagnostic Data}
\usage{
set_default_ellipses(x)
}
\arguments{
\item{x}{An sf data frame output by the function `\link[ctmmUtils]{read_wc_dirs}`.}
}
\description{
This function will fill in missing ARGOS diagnostic ellipse data for either type of location, `Argos` or `FastGPS`. This is necessary for converting the
data to `telemetry` format for the `ctmm` package. The default values used for non-KF Argos data (old location classes) are those in Vincent et al. (2002; Table 1 nonfilterd). For
`FastGPS` data, the values in Dujon et al. (2014; Table 1, 95\%-tile unfiltered.).

Because Dujon et al. (2014) only used error magnitude, we assumed normally distributed error positions, \eqn{N(0,\sigma^2)}. Thus, the error magitudes are Rayleigh distributed. So
we used the CDF of the Rayleigh distribution and the 95\%-tile given in the table to solve for the appropriate \eqn{\sigma}.
}
