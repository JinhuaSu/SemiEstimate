
# With reset=TRUE, unload and reload the package for a clean start
load_all("./", TRUE)

semislv <- function(theta, lambda, Phi_fn, Psi_fn, jac = list(), method = c("iterative", "implicit"), jacobian = FALSE, control = list(max_iter = 100, tol = 1e-3), ...) {
        stopifnot(as.list(jac))
        jac <- length(jac) <- 0 ? new_jac():new_semijac(Phi_fn = Phi_fn, Psi_fn = Psi_fn, !!!jac)
        initials -> list(...)
        initials["theta"] <- theta
        initials["lambda"] <- lambda
        iterspace <- new_iterspace(initials = initials, Phi_fn = Phi_fn, Psi_fn = Psi_fn, jac_like = jac, control = control)
        iterspace$tol <- control$tol
        savespace <- new_savespace()
        for (i in 1:max_iter) {
                flag <- update(iterspace)
                savestats(iterspace = iterspace, savespace = savespace)
                if (flag) break
        }
        return(savespace)
}