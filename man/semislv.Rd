% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/solve.R
\name{semislv}
\alias{semislv}
\title{Solve Semi-parametric estimation by implicit profiling}
\usage{
semislv(
  theta,
  lambda,
  Phi_fn,
  Psi_fn,
  jac = list(),
  intermediates = list(),
  method = "implicit",
  diy = FALSE,
  control = list(max_iter = 100, tol = 0.001),
  save = list(time = TRUE, path = FALSE),
  ...
)
}
\arguments{
\item{theta}{the initial value of parametric part}

\item{lambda}{the initial value of non-parametric part}

\item{Phi_fn}{the equation function highly relevant to the parametric part}

\item{Psi_fn}{the equation function highly relevant to the non-parametric part}

\item{jac}{a list containing some of deterivate info of Phi_der_theta_fn, Psi_der_theta_fn, Phi_der_lambda_fn, Psi_der_lambda_fn,}

\item{intermediates}{a list containing the important variables for diy mode}

\item{method}{"implicit" or "iterative"}

\item{diy}{a bool value to decide to parse user designed function}

\item{control}{a list like list(max_iter = 100, tol = 1e-3) to control the early stop}

\item{save}{a list like list(time = FALSE, path = FALSE) to control saving setting}

\item{...}{static parameter for Phi_fn, Psi_fn. Diy execution function.}
}
\value{
A save space containing final iteration result and iteration path
}
\description{
Solve Semi-parametric estimation by implicit profiling
}
\examples{
Phi_fn <- function(theta, lambda, alpha) 2 * theta + alpha * lambda
Psi_fn <- function(theta, lambda, alpha) 2 * lambda + alpha * theta
# build quasi jacobiean by package NumDeriv
res <- semislv(1, 1, Phi_fn, Psi_fn, alpha = 1)
res <- semislv(1, 1, Phi_fn, Psi_fn, method = "iterative", alpha = 1)
# parsing all mathematical Jacobian function by user
res <- semislv(1, 1, Phi_fn, Psi_fn, jac = list(
        Phi_der_theta_fn = function(theta, lambda, alpha) 2,
        Phi_der_lambda_fn = function(theta, lambda, alpha) alpha,
        Psi_der_theta_fn = function(theta, lambda, alpha) alpha,
        Psi_der_lambda_fn = function(theta, lambda, alpha) 2
), method = "implicit", alpha = 1)
res <- semislv(1, 1, Phi_fn, Psi_fn, jac = list(
        Phi_der_theta_fn = function(theta, lambda, alpha) 2,
        Psi_der_lambda_fn = function(theta, lambda, alpha) 2
), method = "iterative", alpha = 1)
# parsing partial mathemetical user-provided Jacobian, the rest will be generated by the NumDeriv
res <- semislv(1, 1, Phi_fn, Psi_fn,
        jac = list(Phi_der_theta_fn = function(theta, lambda, alpha) 2),
        method = "implicit", alpha = 1
)
res <- semislv(1, 1, Phi_fn, Psi_fn,
        jac = list(Phi_der_theta_fn = function(theta, lambda, alpha) 2),
        method = "iterative", alpha = 1
)
# use some package or solve the updating totally by the user
# Cases: (1) use thirty party package (2) save the intermediates
# use diy = True, then the package will be just a wrapper for your personalise code
# diy is an advanced mode for researchers, see more examples in our vigettee documents
}
