## semi_solve()

# =====================================
# ======= z = x^2 + y^2 + alpha xy=====
# =====================================
# test simple jac


# provided jocabiean and mathematical
Phi_fn <- function(theta, lambda, alpha) 2 * theta + alpha * lambda
Psi_fn <- function(theta, lambda, alpha) 2 * lambda + alpha * theta
res <- semislv(1, 1, Phi_fn, Psi_fn,method = "implicit", alpha = 1)
res <- semislv(1, 1, Phi_fn, Psi_fn, jac = list(Phi_der_theta_fn = function(theta, lambda, alpha) 2, Phi_der_lambda_fn = function(theta, lambda, alpha) alpha, Psi_der_theta_fn = function(theta, lambda, alpha) alpha, Psi_der_lambda_fn = function(theta, lambda, alpha) 2), method = "implicit", alpha = 1)
res <- semislv(1, 1, Phi_fn, Psi_fn, jac = list(Phi_der_theta_fn = function(theta, lambda, alpha) 2, Phi_der_lambda_fn = function(theta, lambda, alpha) alpha, Psi_der_theta_fn = function(theta, lambda, alpha) alpha, Psi_der_lambda_fn = function(theta, lambda, alpha) 2), method = "iterative", alpha = 1)

Newton <- function(beta0, alpha) {
        H_GS <- matrix(c(2, alpha, alpha, 2), nrow = 2, ncol = 2) # Hessian Matrix
        beta_GS <- beta0
        step_GS <- 0
        series <- NULL
        t0 <- Sys.time()
        f <- function(x, y) {
                return(c(2 * x + alpha * y, 2 * y + alpha * x))
        } # 要求解的两个式子

        while (TRUE) {
                bscore <- f(beta_GS[1], beta_GS[2])
                if (all(abs(bscore) < 1e-7)) break
                beta_GS <- beta_GS - solve(H_GS, bscore) # 牛顿法迭代式
                step_GS <- step_GS + 1
                series <- rbind(series, beta_GS)
        }

        run_time_GS <- Sys.time() - t0
        return(list(
                beta = beta_GS, step = step_GS, run_time = run_time_GS,
                series = series
        ))
}

get_fit_from_raw -> function(raw_fit) {
        fit <- list()
        fit$beta <- c(raw_fit$parameters$theta, raw_fit$parameters$lambda)
        fit$step <- raw_fit$step
        fit$run_time <- raw_fit$run.time
        series <- c(res$iterspace$initials$theta, res$iterspace$initials$lambda)
        purrr::walk(raw_fit$res_path, function(x) series <<- rbind(series, c(x$parameters$theta, x$parameters$lambda)))
        fit$series <- series
        fit
}


run_Ip <- function(intermediates, theta, lambda, alpha) {
        x <- theta
        y <- lambda
        yscore <- 2 * y + alpha * x
        intermediates$y_delta <- -yscore / 2 # IP lambda迭代式
        xscore <- 2 * x + alpha * y
        intermediates$x_delta <- -2 * xscore / (4 - alpha^2) # IP theta迭代式
        intermediates
}

theta_delta_Ip <- function(intermediates) {
        intermediates$theta_delta <- intermediates$x_delta
        intermediates
}

lambda_delta_Ip <- function(intermediates) {
        intermediates$lambda_delta <- intermediates$y_delta
        intermediates
}

run_It <- function(intermediates, theta, lambda, alpha) {
        x <- theta
        y <- lambda
        yscore <- 2 * y + alpha * x
        intermediates$y_delta <- -yscore / 2 # IP lambda迭代式
        xscore <- 2 * x + alpha * y
        intermediates$x_delta <- -2 * xscore / (4 - alpha^2) # IP theta迭代式
        intermediates
}

theta_delta_It <- function(intermediates) {
        intermediates$theta_delta <- intermediates$x_delta
        intermediates
}

lambda_delta_It <- function(intermediates) {
        intermediates$lambda_delta <- intermediates$y_delta
        intermediates
}
# =================================
# =========== 均匀撒点测试 ========
# =================================

j <- 1
step_all <- list()
series_all <- list()
direction_all <- list()
for (k in c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8)) {
        step <- list()
        series <- list()
        direction <- list()
        length <- 10
        theta <- seq(0, 2 * base::pi, length.out = length)
        alpha <- k
        for (i in 1:10) {
                C <- i^2
                x <- (sqrt(C * (1 - alpha / 2)) * cos(theta) + sqrt(C * (1 + alpha / 2)) * sin(theta)) / sqrt(2 - alpha^2 / 2)
                y <- (sqrt(C * (1 - alpha / 2)) * cos(theta) - sqrt(C * (1 + alpha / 2)) * sin(theta)) / sqrt(2 - alpha^2 / 2)
                sub_step <- matrix(nrow = 3, ncol = length)
                sub_series <- list()
                k1 <- list()
                k2 <- list()
                k3 <- list()
                sub_direction <- list()
                for (ii in 1:length) {
                        beta0 <- c(x[ii], y[ii])
                        Newton_fit <- Newton(beta0, alpha)
                        ## It_fit = It(beta0, alpha)
                        ## Ip_fit = Ip(beta0, alpha)
                        ## Ip_raw_fit <- semislv(theta = beta0[1], lambda = beta0[2], Phi_fn, Psi_fn, jac = list(Phi_der_theta_fn = function(theta, lambda, alpha) 2, Phi_der_lambda_fn = function(theta, lambda, alpha) alpha, Psi_der_theta_fn = function(theta, lambda, alpha) alpha, Psi_der_lambda_fn = function(theta, lambda, alpha) 2), method = "implicit", jacobian = TRUE, alpha = alpha)
                        Ip_raw_fit <- semislv(theta = beta0[1], lambda = beta0[2], Phi_fn, Psi_fn, method = "implicit", diy = TRUE, run_Ip = run_Ip, theta_delta = theta_delta, lambda_delta = lambda_delta, alpha = alpha)
                        Ip_fit <- get_fit_from_raw(Ip_raw_fit)
                        It_raw_fit <- semislv(theta = beta0[1], lambda = beta0[2], Phi_fn, Psi_fn, jac = list(Phi_der_theta_fn = function(theta, lambda, alpha) 2, Phi_der_lambda_fn = function(theta, lambda, alpha) alpha, Psi_der_theta_fn = function(theta, lambda, alpha) alpha, Psi_der_lambda_fn = function(theta, lambda, alpha) 2), method = "iterative", jacobian = TRUE, alpha = alpha)
                        It_fit <- get_fit_from_raw(It_raw_fit)
                        sub_step[, ii] <- c(Newton_fit$step, It_fit$step, Ip_fit$step)
                        k1[[ii]] <- Newton_fit$series
                        k2[[ii]] <- It_fit$series
                        k3[[ii]] <- Ip_fit$series
                        sub_direction[[ii]] <- It_fit$direction
                }
                step[[i]] <- sub_step
                sub_series[["Newton"]] <- k1
                sub_series[["It"]] <- k2
                sub_series[["Ip"]] <- k3
                series[[i]] <- sub_series
                direction[[i]] <- sub_direction
        }
        step_all[[j]] <- step
        series_all[[j]] <- series
        direction_all[[j]] <- direction
        j <- j + 1
}

for (i in 1:10) {
        k <- step_all[[i]]
        s <- NULL
        for (j in 1:length(k)) {
                s <- cbind(s, k[[j]])
        }
        print(apply(s, 1, mean))
}

# [1] 1 2 2
# [1] 1.00 3.78 3.00
# [1] 1.00 4.38 3.00
# [1] 1.00 5.22 3.00
# [1] 1.00 6.25 3.00
# [1] 1.00 7.48 3.00
# [1] 1.00 9.38 3.00
# [1]  1.00 12.29  3.00
# [1]  1.00 17.88  3.00
# [1]  1.0 34.6  3.0

k <- step_all[[9]]
k <- lapply(k, apply, 1, mean)
s <- NULL
for (i in 1:10) {
        s <- rbind(s, k[[i]])
}
print(s)

#       [,1] [,2] [,3]
#  [1,]    1 14.4    3
#  [2,]    1 16.2    3
#  [3,]    1 17.1    3
#  [4,]    1 17.5    3
#  [5,]    1 18.2    3
#  [6,]    1 18.4    3
#  [7,]    1 18.9    3
#  [8,]    1 19.2    3
#  [9,]    1 19.4    3
# [10,]    1 19.5    3