semislv <- function(theta, lambda, Phi_fn, Psi_fn, jac = list(), method = c("iterative", "implicit"), jacobian = FALSE, control = list(max_iter = 100, tol = 1e-3), ...) {
        stopifnot(is.list(jac))
        # Phi_der_theta_fn = jac$Phi_der_theta_fn, Phi_der_lambda_fn = jac$Phi_der_lambda_fn, Psi_der_theta_fn = jac$Psi_der_theta_fn, Psi_der_lambda_fn = jac$Psi_der_lambda_fn
        print(length(jac))
        jac_like <- ifelse(length(jac) == 4, do.call(new_jac, jac), new_semijac(Phi_fn = Phi_fn, Psi_fn = Psi_fn, !!!jac))
        print(jac_like)
        print(class(jac_like))
        initials <- list(...)
        initials["theta"] <- theta
        initials["lambda"] <- lambda
        initials["method"] <- method
        iterspace <- new_iterspace(initials = initials, Phi_fn = Phi_fn, Psi_fn = Psi_fn, jac_like = jac_like, control = control)
        iterspace$tol <- control$tol
        savespace <- new_savespace()
        for (i in 1:control$max_iter) {
                print(class(iterspace))
                flag <- update(iterspace)
                savestats(iterspace = iterspace, savespace = savespace)
                if (flag) break
        }
        return(savespace)
}