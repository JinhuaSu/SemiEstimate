## semi_solve()

# =====================================
# ======= z = x^2 + y^2 + alpha xy=====
# =====================================
# test simple jac
Phi_fn <- function(theta, lambda, alpha) 2 * theta + alpha * lambda

Psi_fn <- function(theta, lambda, alpha) 2 * lambda + alpha * theta
res <- semislv(1, 1, Phi_fn, Psi_fn, jac = list(Phi_der_theta_fn = function(theta, lambda, alpha) 2, Phi_der_lambda_fn = function(theta, lambda, alpha) alpha, Psi_der_theta_fn = function(theta, lambda, alpha) alpha, Psi_der_lambda_fn = function(theta, lambda, alpha) 2), method = "implicit", jacobian = TRUE,alpha=1)
res <- semislv(1, 1, Phi_fn, Psi_fn, jac = list(Phi_der_theta_fn = function(theta, lambda, alpha) 2, Phi_der_lambda_fn = function(theta, lambda, alpha) alpha, Psi_der_theta_fn = function(theta, lambda, alpha) alpha, Psi_der_lambda_fn = function(theta, lambda, alpha) 2), method = "iterative", jacobian = TRUE,alpha=1)
Newton <- function(beta0, alpha) {
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
                if (all(abs(bscore) < 1e-7)) break
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

It <- function(beta0, alpha) {
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
                if (all(abs(c(xscore, yscore)) < 1e-7)) break
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

Ip <- function(beta0, alpha) {
        beta_DPS <- beta0
        step_DPS <- 0
        series_DPS <- NULL
        t0 <- Sys.time()
        while (TRUE) {
                x <- beta_DPS[1]
                y <- beta_DPS[2]
                yscore <- 2 * y + alpha * x
                y <- y - yscore / 2
                xscore <- 2 * x + alpha * y
                x <- x - 2 * xscore / (4 - alpha^2)
                if (all(abs(c(xscore, yscore)) < 1e-7)) break
                beta_DPS <- c(x, y)
                series_DPS <- rbind(series_DPS, beta_DPS)
                step_DPS <- step_DPS + 1
        }

        run_time_DPS <- Sys.time() - t0
        return(list(
                beta = beta_DPS, step = step_DPS, run_time = run_time_DPS,
                series = series_DPS
        ))
}