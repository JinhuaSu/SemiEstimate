#' @export
semislv <- function(theta, lambda, Phi_fn, Psi_fn, jac = list(), intermediates = list(), method = c("iterative", "implicit"), diy = FALSE, control = list(max_iter = 100, tol = 1e-3), ...) {
        if (diy) {
                args <- rlang::dots_list(theta = theta, lambda = lambda, method = method, intermediates = intermediates, !!!list(...), .homonyms = "first")
                jac_like <- do.call(new_diyjac, args)
        } else if (length(jac) == 4) {
                jac_like <- do.call(new_jac, jac)
        } else {
                jac$Phi_fn < Phi_fn
                jac$Psi_fn < Psi_fn
                jac_like <- do.call(new_semijac, jac)
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
                savespace <- savestats(iterspace = iterspace, savespace = savespace, step = i)
                if (iterspace$iter_over) break
        }
        savespace$run.time <- Sys.time() - t0
        return(savespace)
}