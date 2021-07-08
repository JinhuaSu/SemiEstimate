
Newton <- function(beta0, alpha, tol = 1e-7) {
        H_GS <- matrix(c(2, alpha, alpha, 2), nrow = 2, ncol = 2)
        beta_GS <- beta0
        step_GS <- 0
        series <- NULL
        t0 <- Sys.time()
        f <- function(x, y) {
                return(c(2 * x + alpha * y, 2 * y + alpha * x))
        }

        while (TRUE) {
                bscore <- f(beta_GS[1], beta_GS[2])
                if (all(abs(bscore) < tol)) break
                beta_GS <- beta_GS - solve(H_GS, bscore)
                step_GS <- step_GS + 1
                series <- rbind(series, beta_GS)
        }

        run_time_GS <- Sys.time() - t0
        return(list(
                beta = beta_GS, step = step_GS, run_time = run_time_GS,
                series = series
        ))
}

It <- function(beta0, alpha, tol = 1e-7) {
        beta_IPS <- beta0
        step_IPS <- 0
        series_IPS <- NULL
        direction <- NULL
        t0 <- Sys.time()
        while (TRUE) {
                x <- beta_IPS[1]
                y <- beta_IPS[2]

                yscore <- 2 * y + alpha * x
                y <- y - yscore / 2
                xscore <- 2 * x + alpha * y
                x <- x - xscore / 2
                if (all(abs(c(xscore, yscore)) < tol)) break
                beta_IPS <- c(x, y)
                series_IPS <- rbind(series_IPS, beta_IPS)
                direction <- rbind(direction, c(xscore, yscore))
                step_IPS <- step_IPS + 1
        }

        run_time_IPS <- Sys.time() - t0
        return(list(
                beta = beta_IPS, step = step_IPS, run_time = run_time_IPS,
                series = series_IPS, direction = direction
        ))
}

build_ip_der <- function(Phi, Psi, Phi_der_theta, Psi_der_theta, Phi_der_lambda, Psi_der_lambda) {

}

build_it_der <- function(Phi, Psi, Phi_der_theta, Psi_der_lambda) {

}

# 1. Psi
# 2. Psi_der_lambda
# 3. Psi_der_theta
# 4. Phi
# 5. Phi_der_theta
# 6. Phi_der_lambda

der_class <- methods::setRefClass("der_class",
        fields = list(
                theta = NULL,
                lambda = NULL,
                data =
                        list(
                                Z = Z,
                                delta = NULL,
                                KC = KC,
                                n = n,
                                KCd = KCd
                        ),
                step = 0,
                intermediates =
                        list(
                                hC = lambda,
                                beta = theta,
                                dif = NULL,
                                hCd = NULL,
                                lp = NULL,
                                gij = NULL,
                                wZbar = NULL,
                                tmp = NULL,
                                hscore = NULL,
                                hHess = NULL,
                                Zbar_up = NULL,
                                Z_bar = NULL,
                                gi = NULL,
                                bscore = NULL,
                                part1 = NULL,
                                part2 = NULL,
                                theta_Hess = NULL
                        )
        ),
        methods =
                list(
                        Psi = function() {
                                expit <- function(d) {
                                        return(1 / (1 + exp(-d)))
                                }
                                intermediates.lp <<- drop(data.Z %*% intermediates.beta)
                                intermediates.gij <<- expit(outer(intermediates.hC, intermediates.lp, "+"))
                                intermediates.tmp <<- data.KC * intermediates.gij
                                intermediates.wZbar <<- intermediates.tmp * (1 - intermediates.gij)
                                intermediates.hscore <<- apply(intermediates.tmp, 1, sum) - data.KCd
                                return(intermediates.hscore)
                        },
                        Psi_der_lambda = function() {
                                intermediates.hHess <<- apply(intermediates.wZbar, 1, sum)
                                return(intermediates.hHess)
                        },
                        Psi_der_theta = function() {
                                intermediates.Zbar_up <<- intermediates.wZbar %*% data.Z
                                intermediates.Zbar <<- intermediates.Zbar_up / hHess
                                return(intermediates.Zbar_up)
                        },
                        Phi = function() {
                                intermediates.gi <<- expit(intermediates.hC + intermediates.lp)
                                intermediates.bscore <<- drop(t(Z) %*% (delta - gi))
                                return(intermediates.bscore)
                        },
                        Phi_der_theta = function() {
                                intermediates.part1 <<- t(intermediates.gi * (1 - intermediates.gi) * Z)
                        },
                        Phi_der_lambda = function() {
                                intermediates.part2 <<- intermediates.Zbar - data.Z
                        },
                        theta_Hess = function() {
                                intermediates.theta_Hess <<- intermediates.part1 %*% intermediates.part2
                                step <<- step + 1
                                return(intermediates.thet)
                        }
                )
)

der_sa <- der_class(
        theta = theta0,
        lambda = lambda0,
        data =
                list(
                        Z = Z,
                        delta = NULL,
                        KC = KC,
                        n = n,
                        KCd = KCd
                ),
        step = 0,
        intermediates = list(
                hC = lambda,
                beta = theta,
                dif = NULL,
                hCd = NULL,
                lp = NULL,
                gij = NULL,
                wZbar = NULL,
                tmp = NULL,
                hscore = NULL,
                hHess = NULL,
                Zbar_up = NULL,
                Z_bar = NULL,
                gi = NULL,
                bscore = NULL,
                part1 = NULL,
                part2 = NULL,
                theta_Hess = NULL
        )
)


iter_class <- methods::setRefClass("ip_iter", fields = list(der_class_inst = der_sa, max_iter = max_iter, result = list(theta = c(der_class_inst.theta), lambda = c(der_class_inst.lambda))), methods = list(update = function() {
        psi <- der_class_inst.Psi()
        der_class_inst.Psi_der_lambda()
        lambda_Hess <- der_class_inst.Psi_der_theta()
        lambda_update <- -1 * solve(lambda_Hess) * psi
        result.lambda <<- c(result.lambda, result.lambda[length(result.lambda)] + lambda_update)
        phi <- der_class_inst.Phi()
        der_class_inst.Phi_der_theta()
        der_class_inst.Phi_der_lambda()
        theta_Hess <- der_class_inst.theta_Hess()
        theta_update <- -1 * solve(theta_Hess) * phi
        result.theta <<- c(result.theta, result.theta[length(result.theta)] + theta_update)
}, iter_over = function() {
        return(der_class_inst.step >= max_iter)
}))

iter_sa <- iter_class()

build_class <- methods::setRefClass("build_ip", fields = list(iter_class_inst = NULL, ), methods = list(function() {
        return(iter_class_inst.result)
}))

build_sa <- build_class(iter_class_inst = iter_sa)


semi_solve <- function(theta0, lambda0, data, Phi, Psi, Phi_der_theta = NULL, Psi_der_theta = NULL, Phi_der_lambda = NULL, Psi_der_lambda = NULL, der_class = NULL, method = "IP", lambda_flag = FALSE, lambda_loop = FALSE, save_intermediate = FALSE, record_time = FALSE, max_iter = 100, tol = 1e-7) {
        # build a der_class(reference class) from expressions or approximation
        if (der_class == NULL) {
                if (method == "IP") {
                        der_class <- build_ip_der(Phi, Psi, Phi_der_theta, Psi_der_theta, Phi_der_lambda, Psi_der_lambda)
                } else {
                        der_class <- build_it_der(Phi, Psi, Phi_der_theta, Psi_der_lambda)
                }
        }
        # build a iter_class(reference class) from der_class(S3)
        if (method == "IP") {
                iter_class <- build_ip_iter(theta0, lambda0, lambda_flag, lambda_loop, save_intermediate, tol)
        } else {
                iter_class <- build_it_iter(theta0, lambda0, lambda_flag, lambda_loop, save_intermediate, tol)
        }
        if (record_time) {
                t0 <- Sys.time()
        }
        for (iter in 1:max_iter) {
                iter_class.update()
                if (iter_class.iter_over()) break
        }
        # get res_class(S4) from iter_class
        if (record_time) {
                t <- Sys.time() - t0
                res_class <- build_res(iter_class, t)
        } else {
                res_class <- build_res(iter_class)
        }
        # rewrite the print for res_class
        return(res_class)
}