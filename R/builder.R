


new_eqfns <- function(Phi_fn = function(theta, lambda, ...) NULL,
                      Psi_fn = function(theta, lambda, ...) NULL) {
        stopifnot(is.function(Phi_fn))
        stopifnot(is.function(Psi_fn))
        data <- list(Phi_fn = Phi_fn, Psi_fn = Psi_fn)
        structure(data, class = "eqfns")
}

new_jac <- function(Phi_der_theta_fn = function(theta, lambda, ...) NULL,
                    Phi_der_lambda_fn = function(theta, lambda, ...) NULL,
                    Psi_der_theta_fn = function(theta, lambda, ...) NULL,
                    Psi_der_lambda_fn = function(theta, lambda, ...) NULL) {
        stopifnot(is.function(Phi_der_theta_fn))
        stopifnot(is.function(Phi_der_lambda_fn))
        stopifnot(is.function(Psi_der_theta_fn))
        stopifnot(is.function(Psi_der_lambda_fn))
        print(Phi_der_theta_fn)
        data <- list(
                Phi_der_theta_fn = Phi_der_theta_fn,
                Phi_der_lambda_fn = Phi_der_lambda_fn,
                Psi_der_theta_fn = Psi_der_theta_fn,
                Psi_der_lambda_fn = Psi_der_lambda_fn
        )
        x <- structure(data,
                class = "jac"
        )
        print(class(x))
        x
}


new_quasijac <- function(Phi_fn = function(theta, lambda, ...) NULL,
                         Psi_fn = function(theta, lambda, ...) NULL) {
        stopifnot(is.function(Phi_fn))
        stopifnot(is.function(Psi_fn))
        Phi_der_theta_fn <- function(theta, lambda, ...) numDeriv::jacobian(func = function(x) Phi_fn(theta = x, lambda = lambda, ...), x = theta)
        Phi_der_lambda_fn <- function(theta, lambda, ...) numDeriv::jacobian(func = function(x) Phi_fn(theta = theta, lambda = x, ...), x = lambda)
        Psi_der_theta_fn <- function(theta, lambda, ...) numDeriv::jacobian(func = function(x) Psi_fn(theta = x, lambda = lambda, ...), x = theta)
        Psi_der_lambda_fn <- function(theta, lambda, ...) numDeriv::jacobian(func = function(x) Psi_fn(theta = theta, lambda = x, ...), x = lambda)
        data <- list(
                Phi_der_theta_fn = Phi_der_theta_fn,
                Phi_der_lambda_fn = Phi_der_lambda_fn,
                Psi_der_theta_fn = Psi_der_theta_fn,
                Psi_der_lambda_fn = Psi_der_lambda_fn
        )
        structure(
                data,
                class = "quasijac"
        )
}

new_semijac <- function(Phi_fn = function(theta, lambda, ...) NULL, Psi_fn = function(theta, lambda, ...) NULL, ...) {
        stopifnot(is.function(Phi_fn))
        stopifnot(is.function(Psi_fn))
        quasijac <- new_quasijac(Phi_fn = Phi_fn, Psi_fn = Psi_fn)
        der_name_list <- c("Phi_der_theta_fn", "Phi_der_lambda_fn", "Psi_der_theta_fn", "Psi_der_lambda_fn")
        args <- list(...)
        quasi_names <- list()
        purrr::map(der_name_list, function(x) if (x %in% names(args)) eval(parse(text = paste("args$", x, "<<-", "quasijac$", x, ";quasi_names<<-", x))))
        data <- list(
                Phi_der_theta_fn = args$Phi_der_theta_fn,
                Phi_der_lambda_fn = args$Phi_der_lambda_fn,
                Psi_der_theta_fn = args$Psi_der_theta_fn,
                Psi_der_lambda_fn = args$Psi_der_lambda_fn
        )
        if (length(quasi_names) == 0) {
                quasijac
        } else {
                structure(
                        data,
                        quasi_names = quasi_names, class = "semijac"
                )
        }
}

# use key map to wrap it
new_diyjac <- function(...) {
        args <- list(...)
        ordered_fn <- list()
        fn_args <- list()
        map_fn <- function(x) {
                if (is.function(x)) {
                        ordered_fn[x] <<- args[x]
                } else {
                        fn_args[x] <<- args[x]
                }
        }
        purrr::walk(names(args), map_fn)
        fn_args_name <- names(fn_args)
        intermedials <- purrr::map(fn_args_name, function(x) eval(x))
        return_fn <- ifelse(("method" %in% fn_args_name) && method == "implicit", function() list(Phi_der_theta = Phi_der_theta, Phi_der_lambda = Phi_der_lambda, Psi_der_theta = Psi_der_theta, Psi_der_lambda = Psi_der_lambda), function() list(Phi_der_theta = Phi_der_theta, Psi_der_lambda = Psi_der_lambda))
        data <- list(ordered_fn = ordered_fn, intermedials = intermedials, return_fn = return_fn)
        structure(
                data,
                class = "diyjac"
        )
}

get_subclass <- function(method, jac_class) {
        first <- ifelse(method == "iterative", "IT", "IP")
        second <- ifelse(class(jac_class) == "diyjac", "HM", "AT")
        paste(first, second, sep = "")
}

new_iterspace <- function(initials = list(), Phi_fn, Psi_fn, jac_like, control) {
        stopifnot("theta" %in% names(initials))
        stopifnot("lambda" %in% names(initials))
        stopifnot("method" %in% names(initials))
        sub_class_name <- get_subclass(initials$method, class(jac_like))
        eqfns <- list(Phi_fn = Phi_fn, Psi_fn = Psi_fn)
        data <- list(initials = initials, eqfns = eqfns, jac_like = jac_like, iter_step = 0, update_delta = list(theta = NULL, lambda = NULL), parameters = list(lambda = initials$lambda, theta = initials$theta), control = control)
        structure(data, class = paste("iterspace.", sub_class_name, sep = ""))
}

new_savespace <- function(save.path = FALSE, ...) {
        data <- list(save.path = save.path, save.time = FALSE, save.der = FALSE, cache = list(...), res_path = list())
        # save.time = save.time, save.der = save.der,
        structure(data, class = "savespace")
}