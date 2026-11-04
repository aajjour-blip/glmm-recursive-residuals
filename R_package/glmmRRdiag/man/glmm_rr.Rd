
\name{glmm_rr}
\alias{glmm_rr}
\title{Posterior-predictive recursive residual diagnostics for a GLMM}
\usage{
glmm_rr(fit, data = NULL, group = NULL, order = NULL, index_vars = NULL,
        B = 0, n_quad = 25, seed = NULL, alpha = 0.05,
        progress = interactive())
}
\arguments{
\item{fit}{A fitted \code{lme4} GLMM object inheriting from \code{merMod}.}
\item{data}{Optional data frame used to fit the model.}
\item{group}{Character name of the clustering variable.}
\item{order}{Optional character name of the within-cluster ordering variable.}
\item{index_vars}{Optional character vector of variables used for indexed cumulative residual processes.}
\item{B}{Number of parametric refit-calibration replicates. Use \code{B = 0} for observed residuals only.}
\item{n_quad}{Number of normal quadrature nodes for random-intercept posterior integration.}
\item{seed}{Optional random seed.}
\item{alpha}{Nominal decision level for the detected flag.}
\item{progress}{Logical; print progress during calibration.}
}
\description{
Computes posterior-predictive recursive residuals for supported random-intercept generalized linear mixed models and optionally performs parametric refit calibration.
}
\value{An object of class \code{glmm_rr}.}
\details{
The manuscript prototype supports \code{lme4::glmer}-style random-intercept Bernoulli/binomial and Poisson models, with experimental Gamma support. The returned object contains residuals, observed component statistics, optional refit-calibrated component p-values, the raw minimum component p-value, and the calibrated minimum-p omnibus statistic \code{p_combo}.
}
\examples{
\dontrun{
library(lme4)
data(cbpp, package = "lme4")
cbpp$y <- with(cbpp, incidence / size) # Bernoulli-count responses need expansion for this prototype.
}
}
