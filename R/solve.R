#' Solve Semi-parametric estimation by implicit profiling
#'
#' @param theta the initial value of parametric part
#' @param lambda the initial value of non-parametric part
#' @param Phi_fn the equation function highly relevant to the parametric part
#' @param Psi_fn the equation function highly relevant to the non-parametric part
#' @param jac a list containing some of deterivate info of Phi_der_theta_fn, Psi_der_theta_fn, Phi_der_lambda_fn, Psi_der_lambda_fn,
#' @param method "implicit" or "iterative"
#' @param diy a bool value to decide to parse user designed function
#' @param control a list like list(max_iter = 100, tol = 1e-3) to control the early stop
#' @param ...  static parameter for Phi_fn, Psi_fn. Diy execution function.
#' @return A save space containing final iteration result and iteration path
#' @examples
#' Phi_fn <- function(theta, lambda, alpha) 2 * theta + alpha * lambda
#' Psi_fn <- function(theta, lambda, alpha) 2 * lambda + alpha * theta
#' res <- semislv(1, 1, Phi_fn, Psi_fn, method = "implicit", alpha = 1)
#' @export
semislv <- function(theta, lambda, Phi_fn, Psi_fn, jac = list(), intermediates = list(), method = c("iterative", "implicit"), diy = FALSE, control = list(max_iter = 100, tol = 1e-3), ...) {
        if (diy) {
                args <- rlang::dots_list(theta = theta, lambda = lambda, method = method, intermediates = intermediates, !!!list(...), .homonyms = "first")
                jac_like <- do.call(new_diyjac, args)
        } else if (length(jac) == 4) {
                jac_like <- do.call(new_jac, jac)
        } else {
                args <- rlang::dots_list(Phi_fn = Phi_fn, Psi_fn = Psi_fn, !!!jac, !!!list(...), .homonyms = "first")
                jac_like <- do.call(new_semijac, args)
        }
        ## jac_like is only a function with class called jac
        initials <- list(...)
        initials$theta <- theta
        initials$lambda <- lambda
        initials$method <- method

        iterspace <- new_iterspace(initials = initials, Phi_fn = Phi_fn, Psi_fn = Psi_fn, jac_like = jac_like, control = control)
        iterspace$tol <- control$tol
        savespace <- new_savespace(save.path = TRUE)
        t0 <- Sys.time()
        for (i in 1:control$max_iter) {
                iterspace <- update(iterspace)
                if (iterspace$iter_over) break
                savespace <- savestats(iterspace = iterspace, savespace = savespace, step = i)
                # TODO(sujinhua): to unify the results
        }
        savespace$run.time <- Sys.time() - t0
        return(savespace)
}