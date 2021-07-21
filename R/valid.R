validate_fn <- function(fn) {
        arg_names <- names(formals(fn))
        ("lambda" %in% arg_names && "theta" %in% arg_names) ? TRUE:FALSE
}

validate_eqfns <- function(x) {
        values <- unclass(x)
        Phi_fn <- attr(values, "Phi_fn")
        Psi_fn <- attr(values, "Psi_fn")
        stopifnot(validate_fn(Phi_fn))
        stopifnot(validate_fn(Psi_fn))
}

validate_jac <- function(x) {
        values <- unclass(x)
        Phi_der_theta_fn <- attr(values, "Phi_der_theta_fn")
        Phi_der_lambda_fn <- attr(values, "Phi_der_lambda_fn")
        Psi_der_theta_fn <- attr(values, "Psi_der_theta_fn")
        Psi_der_lambda_fn <- attr(values, "Psi_der_lambda_fn")
        stopifnot(validate_fn(Phi_der_theta_fn))
        stopifnot(validate_fn(Phi_der_lambda_fn))
        stopifnot(validate_fn(Psi_der_theta_fn))
        stopifnot(validate_fn(Psi_der_lambda_fn))
}

validate_quasijac <- function(x) {
        values <- unclass(x)
        Phi_der_theta_fn <- attr(values, "Phi_der_theta_fn")
        Phi_der_lambda_fn <- attr(values, "Phi_der_lambda_fn")
        Psi_der_theta_fn <- attr(values, "Psi_der_theta_fn")
        Psi_der_lambda_fn <- attr(values, "Psi_der_lambda_fn")
        stopifnot(validate_fn(Phi_der_theta_fn))
        stopifnot(validate_fn(Phi_der_lambda_fn))
        stopifnot(validate_fn(Psi_der_theta_fn))
        stopifnot(validate_fn(Psi_der_lambda_fn))
}

validate_semijac <- function(x) {
        values <- unclass(x)
        Phi_der_theta_fn <- attr(values, "Phi_der_theta_fn")
        Phi_der_lambda_fn <- attr(values, "Phi_der_lambda_fn")
        Psi_der_theta_fn <- attr(values, "Psi_der_theta_fn")
        Psi_der_lambda_fn <- attr(values, "Psi_der_lambda_fn")
        stopifnot(validate_fn(Phi_der_theta_fn))
        stopifnot(validate_fn(Phi_der_lambda_fn))
        stopifnot(validate_fn(Psi_der_theta_fn))
        stopifnot(validate_fn(Psi_der_lambda_fn))
}

## TODO(sujinhua)

validate_diyjac <- function() NULL