
update <- function(iterspace) {
        UseMethod(update)
}

small_update <- function(vec, tol){
	return(all(abs(vec) < tol))
}

update.iterspace.ITAT <- function(iterspace) {

	Phi <- iterspace$eqfns$Phi_fn(theta = iterspace$parameters$theta, lambda = iterspace$parameters$lambda, !!!parameters)
	Phi_der_theta <- iterspace$jac_like$Phi_der_theta_fn(theta = iterspace$parameters$theta, lambda = iterspace$parameters$lambda, !!!parameters)
	iterapce$update_delta$theta <- -1 * solve(Phi_der_theta,Phi)

	iterspace$parameters$theta <<- theta + iterapce$update_delta$theta
	Psi <- iterspace$eqfns$Psi_fn(theta = iterspace$parameters$theta, lambda = iterspace$parameters$lambda, !!!parameters)
	Psi_der_lambda <- iterspace$jac_like$Phi_der_theta_fn(theta = iterspace$parameters$theta, lambda = iterspace$parameters$lambda, !!!parameters)
	iterapce$update_delta$lambda <- -1 * solve(Psi_der_lambda,Psi)
	iterspace$parameters$lambda <<- lambda + iterapce$update_delta$lambda
	return(small_update(iterapce$update_delta$theta, iterspace$tol) && small_update(iterapce$update_delta$lambda, iterspace$tol))
}

update.iterspace.ITHM <- function(iterspace) {
	Phi <- iterspace$eqfns$Phi_fn(theta = iterspace$parameters$theta, lambda = iterspace$parameters$lambda, !!!parameters)
	purrr::walk(ordered_fn, function(x) x(!!!iterspace$initials, !!!iterspace$parameters))
	der_list <- iterspace$jac_like$return_fn()
	Phi_der_theta <- der_list$Phi_der_theta
	iterapce$update_delta$theta <- -1 * solve(Phi_der_theta,Phi)

	Psi <- iterspace$eqfns$Psi_fn(theta = iterspace$parameters$theta, lambda = iterspace$parameters$lambda, !!!parameters)
	Psi_der_lambda <- der_list$Psi_der_lambda
	iterapce$update_delta$lambda <- -1 * solve(Psi_der_lambda,Psi)
	iterspace$parameters$theta <<- theta + iterapce$update_delta$theta
	iterspace$parameters$lambda <<- lambda + iterapce$update_delta$lambda
	return((small_update(iterapce$update_delta$theta, iterspace$tol) && small_update(iterapce$update_delta$lambda, iterspace$tol))

}

update.iterspace.IPAT <- function(iterspace) {
	Phi <- iterspace$eqfns$Phi_fn(theta = iterspace$parameters$theta, lambda = iterspace$parameters$lambda, !!!parameters)
	Phi_der_theta <- iterspace$jac_like$Phi_der_theta_fn(theta = iterspace$parameters$theta, lambda = iterspace$parameters$lambda, !!!parameters)
	Phi_der_lambda <- iterspace$jac_like$Phi_der_lambda_fn(theta = iterspace$parameters$theta, lambda = iterspace$parameters$lambda, !!!parameters)
	Psi_der_theta <- iterspace$jac_like$Psi_der_theta_fn(theta = iterspace$parameters$theta, lambda = iterspace$parameters$lambda, !!!parameters)
	iterapce$update_delta$theta <- -1 * solve(Phi_der_theta + Phi_der_lambda * solve(Psi_der_lambda,Psi_der_theta),Phi)

	Psi <- iterspace$eqfns$Psi_fn(theta = iterspace$parameters$theta, lambda = iterspace$parameters$lambda, !!!parameters)
	Psi_der_lambda <- iterspace$jac_like$Phi_der_theta_fn(theta = iterspace$parameters$theta, lambda = iterspace$lambda, !!!parameters)
	iterapce$update_delta$lambda <- -1 * solve(Psi_der_lambda,Psi)
	iterspace$parameters$theta <<- theta + iterapce$update_delta$theta
	iterspace$parameters$lambda <<- lambda + iterapce$update_delta$lambda
	return((small_update(iterapce$update_delta$theta, iterspace$tol) && small_update(iterapce$update_delta$lambda, iterspace$tol))
}

update.iterspace.IPHM <- function(iterspace) {
	Phi <- iterspace$eqfns$Phi_fn(theta = iterspace$parameters$theta, lambda = iterspace$parameters$lambda, !!!parameters)
	Phi_der_theta <- iterspace$jac_like$Phi_der_theta_fn(theta = iterspace$parameters$theta, lambda = iterspace$parameters$lambda, !!!parameters)
	Phi_der_lambda <- iterspace$jac_like$Phi_der_lambda_fn(theta = iterspace$parameters$theta, lambda = iterspace$parameters$lambda, !!!parameters)
	Psi_der_theta <- iterspace$jac_like$Psi_der_theta_fn(theta = iterspace$parameters$theta, lambda = iterspace$parameters$lambda, !!!parameters)
	iterapce$update_delta$theta <- -1 * solve(Phi_der_theta + Phi_der_lambda * solve(Psi_der_lambda,Psi_der_theta),Phi)
	Psi <- iterspace$eqfns$Psi_fn(theta = iterspace$parameters$theta, lambda = iterspace$parameters$lambda, !!!parameters)
	Psi_der_lambda <- iterspace$jac_like$Phi_der_theta_fn(theta = iterspace$parameters$theta, lambda = iterspace$parameters$lambda, !!!parameters)
	iterapce$update_delta$lambda <- -1 * solve(Psi_der_lambda,Psi)
	iterspace$parameters$theta <<- theta + iterapce$update_delta$theta
	iterspace$parameters$lambda <<- lambda + iterapce$update_delta$lambda
	return((small_update(iterapce$update_delta$theta, iterspace$tol) && small_update(iterapce$update_delta$lambda,iterspace$tol))
}

savestats <- function(savespace, iterspace, step) {
	savespace$iterspace <<- iterspace
	if(save.path){
		savespace$res_path[step] <<- iterspace
	}
}