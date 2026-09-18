# Derivatives and constrained fitting for the abundance profile test.

.sacat_abundance_effect_derivatives <- function(beta, X_b, gh) {
    n <- nrow(X_b)
    p <- ncol(X_b) + 1L
    X0 <- cbind(X_b[, 1L], Group = 0, X_b[, -1L, drop = FALSE])
    X1 <- X0
    X1[, 2L] <- 1
    displacement <- sqrt(2) * exp(beta[p + 1L]) * gh$node
    weight <- gh$weight / sum(gh$weight)
    eta0 <- outer(drop(X0 %*% beta[seq_len(p)]), displacement, `+`)
    eta1 <- eta0 + beta[2L]
    first0 <- plogis(-eta0)
    first1 <- plogis(-eta1)
    second0 <- -plogis(eta0) * first0
    second1 <- -plogis(eta1) * first1
    value <- mean(drop((plogis(eta1, log.p = TRUE) -
        plogis(eta0, log.p = TRUE)) %*% weight))
    location_gradient <- (crossprod(X1, drop(first1 %*% weight)) -
        crossprod(X0, drop(first0 %*% weight))) / n
    scale_gradient <- mean(drop((first1 - first0) %*%
        (weight * displacement)))
    gradient <- c(location_gradient, scale_gradient)
    hessian <- matrix(0, p + 1L, p + 1L)
    hessian[seq_len(p), seq_len(p)] <-
        (crossprod(X1, X1 * drop(second1 %*% weight)) -
         crossprod(X0, X0 * drop(second0 %*% weight))) / n
    cross <- (crossprod(X1, drop(second1 %*% (weight * displacement))) -
        crossprod(X0, drop(second0 %*% (weight * displacement)))) / n
    hessian[seq_len(p), p + 1L] <-
        hessian[p + 1L, seq_len(p)] <- cross
    hessian[p + 1L, p + 1L] <- scale_gradient +
        mean(drop((second1 - second0) %*% (weight * displacement^2)))
    list(value = value, gradient = gradient, hessian = hessian)
}

.sacat_abundance_profile_kernel <- function(y, N, eta, sigma, X, gh) {
    nodes <- count_log_hy_adaptive_ref(
        y, N, eta, sigma, gh, return_nodes = TRUE
    )
    weight <- exp(nodes$log_terms - nodes$log_term_normalizer)
    displacement <- nodes$x_node - eta
    probability <- plogis(nodes$x_node)
    first <- y - N * probability
    second <- -N * probability * (1 - probability)
    scale_first <- first * displacement
    mean_first <- rowSums(weight * first)
    mean_scale_first <- rowSums(weight * scale_first)
    centered_first <- first - mean_first
    centered_scale_first <- scale_first - mean_scale_first
    eta_eta <- rowSums(weight * second) +
        rowSums(weight * centered_first^2)
    eta_scale <- rowSums(weight * second * displacement) +
        rowSums(weight * centered_first * centered_scale_first)
    scale_scale <- rowSums(weight *
        (second * displacement^2 + scale_first)) +
        rowSums(weight * centered_scale_first^2)
    score <- cbind(X * mean_first, mean_scale_first)
    q <- ncol(score)
    location <- seq_len(q - 1L)
    hessian <- array(0, c(length(y), q, q))
    for (i in seq_along(y)) {
        hessian[i, location, location] <- eta_eta[i] * tcrossprod(X[i, ])
        hessian[i, location, q] <- hessian[i, q, location] <-
            eta_scale[i] * X[i, ]
        hessian[i, q, q] <- scale_scale[i]
    }
    list(log_h = nodes$log_hy, score = score, hessian = hessian)
}

.sacat_abundance_profile_derivatives <- function(beta, y, N, X, gh) {
    positive <- which(y > 0)
    X <- X[positive, , drop = FALSE]
    N <- N[positive]
    eta <- drop(X %*% beta[-length(beta)])
    sigma <- exp(beta[length(beta)])
    observed <- .sacat_abundance_profile_kernel(
        y[positive], N, eta, sigma, X, gh
    )
    zero <- .sacat_abundance_profile_kernel(
        rep(0, length(positive)), N, eta, sigma, X, gh
    )
    log_detection <- zt_log1mexp(zero$log_h)
    ratio <- exp(zero$log_h - log_detection)
    score <- observed$score + zero$score * ratio
    information <- matrix(0, length(beta), length(beta))
    for (i in seq_along(positive)) {
        information <- information - observed$hessian[i, , ] -
            ratio[i] * zero$hessian[i, , ] -
            ratio[i] * (1 + ratio[i]) * tcrossprod(zero$score[i, ])
    }
    list(score = colSums(score), information = information,
         nll = -sum(observed$log_h - log_detection))
}

.sacat_abundance_profile_decode <- function(lambda, target, X_b, gh) {
    p <- length(lambda) + 1L
    nuisance <- setdiff(seq_len(p), 2L)
    beta <- numeric(p)
    beta[nuisance] <- lambda
    beta[2L] <- target
    for (iteration in seq_len(12L)) {
        effect <- .sacat_abundance_effect_derivatives(beta, X_b, gh)
        error <- effect$value - target
        if (abs(error) < 1e-11) return(list(beta = beta, effect = effect))
        slope <- effect$gradient[2L]
        if (!is.finite(slope) || slope < 1e-10) break
        beta[2L] <- beta[2L] - error / slope
        if (!is.finite(beta[2L]) || abs(beta[2L]) > 40) break
    }
    fn <- function(value) {
        beta[2L] <- value
        .sacat_abundance_effect_derivatives(beta, X_b, gh)$value - target
    }
    beta[2L] <- uniroot(fn, c(-40, 40), tol = 1e-11)$root
    list(beta = beta,
         effect = .sacat_abundance_effect_derivatives(beta, X_b, gh))
}

.sacat_abundance_profile_numeric_fit <- function(fit, y, N, X_b, target, gh) {
    X <- fit$solver_diagnostics$X_eta
    p <- length(fit$theta)
    nuisance <- setdiff(seq_len(p), 2L)
    positive_count <- sum(y > 0)
    lower <- c(-80, rep(-40, p - 3L), log(1e-6))
    upper <- c(40, rep(40, p - 3L), log(128))
    objective <- function(lambda) {
        decoded <- .sacat_abundance_profile_decode(lambda, target, X_b, gh)
        -sum(zt_beta_loglik_by_sample_inference(decoded$beta, y, N, X, gh))
    }
    derivative <- function(lambda, hessian = FALSE) {
        steps <- zt_inference_steps(lambda, length(y), lower, upper)
        value <- objective(lambda)
        gradient <- as.numeric(zt_central_derivative_matrix_fixed(
            objective, lambda, steps$step, value, steps$scheme
        ))
        information <- if (hessian) {
            zt_central_hessian_fixed(
                objective, lambda, steps$step, value, steps$scheme
            )
        } else NULL
        list(nll = value, gradient = gradient, information = information)
    }
    fit_start <- function(start) {
        optimization <- optim(
            start[nuisance],
            fn = function(value) objective(value) / positive_count,
            gr = function(value) derivative(value)$gradient / positive_count,
            method = "L-BFGS-B", lower = lower, upper = upper,
            control = list(maxit = 200L, factr = 1e4, pgtol = 1e-8)
        )
        lambda <- optimization$par
        result <- derivative(lambda, TRUE)
        iterations <- 0L
        for (iteration in seq_len(15L)) {
            if (max(abs(result$gradient)) / positive_count <= 1e-7) break
            chol(result$information)
            direction <- -drop(solve(result$information, result$gradient))
            accepted <- FALSE
            for (step in 2^-(0:12)) {
                trial <- lambda + step * direction
                if (any(trial <= lower) || any(trial >= upper)) next
                candidate <- tryCatch(
                    derivative(trial, TRUE), error = function(e) NULL
                )
                if (!is.null(candidate) && is.finite(candidate$nll) &&
                    candidate$nll <= result$nll + 1e-9 &&
                    (candidate$nll < result$nll ||
                     max(abs(candidate$gradient)) < max(abs(result$gradient)))) {
                    lambda <- trial
                    result <- candidate
                    accepted <- TRUE
                    iterations <- iterations + 1L
                    break
                }
            }
            if (!accepted) break
        }
        decoded <- .sacat_abundance_profile_decode(lambda, target, X_b, gh)
        result$theta <- decoded$beta
        result$effect <- decoded$effect$value
        result$gradient <- -result$gradient
        result$score_residue <- max(abs(result$gradient)) / positive_count
        result$constraint_error <- abs(result$effect - target)
        result$iterations <- iterations
        result$fallback <- TRUE
        result$numerical_derivatives <- TRUE
        if (!is.finite(result$score_residue) || result$score_residue > 1e-6) {
            stop("constrained nuisance fit did not converge")
        }
        if (!is.finite(result$constraint_error) || result$constraint_error > 1e-8) {
            stop("abundance effect constraint did not converge")
        }
        chol(result$information)
        result
    }
    candidate <- tryCatch(fit_start(fit$theta), error = function(e) NULL)
    if (!is.null(candidate)) return(candidate)
    candidates <- lapply(c(0.5, 2), function(sigma) {
        start <- fit$theta
        start[p] <- log(sigma)
        tryCatch(fit_start(start), error = function(e) NULL)
    })
    candidates <- Filter(Negate(is.null), candidates)
    if (!length(candidates)) stop("constrained abundance maximum unavailable")
    candidates[[which.min(vapply(candidates, `[[`, numeric(1), "nll"))]]
}

.sacat_abundance_profile_fit <- function(fit, y, N, X_b, target, gh) {
    X <- fit$solver_diagnostics$X_eta
    p <- length(fit$theta)
    nuisance <- setdiff(seq_len(p), 2L)
    positive_count <- sum(y > 0)
    lower <- c(-80, rep(-40, p - 3L), log(1e-6))
    upper <- c(40, rep(40, p - 3L), log(128))
    cache <- new.env(parent = emptyenv())
    cache$lambda <- NULL
    evaluate <- function(lambda) {
        if (identical(lambda, cache$lambda)) return(cache$result)
        decoded <- .sacat_abundance_profile_decode(lambda, target, X_b, gh)
        effect <- decoded$effect
        derivatives <- .sacat_abundance_profile_derivatives(
            decoded$beta, y, N, X, gh
        )
        J <- matrix(0, p, p - 1L)
        J[cbind(nuisance, seq_len(p - 1L))] <- 1
        J[2L, ] <- -effect$gradient[nuisance] / effect$gradient[2L]
        information <- crossprod(J, (derivatives$information +
            derivatives$score[2L] / effect$gradient[2L] *
                effect$hessian) %*% J)
        result <- list(
            theta = decoded$beta, effect = effect$value,
            information = (information + t(information)) / 2,
            gradient = drop(crossprod(J, derivatives$score)),
            nll = derivatives$nll
        )
        cache$lambda <- lambda
        cache$result <- result
        result
    }
    fit_start <- function(start) {
        lambda <- start[nuisance]
        result <- evaluate(lambda)
        iterations <- 0L
        fallback <- FALSE
        for (iteration in seq_len(15L)) {
            iterations <- iteration
            if (max(abs(result$gradient)) / positive_count < 1e-7) break
            direction <- tryCatch({
                chol(result$information)
                drop(solve(result$information, result$gradient))
            }, error = function(e) NULL)
            if (is.null(direction)) break
            accepted <- FALSE
            for (step in 2^-(0:12)) {
                trial <- lambda + step * direction
                if (any(trial <= lower) || any(trial >= upper)) next
                candidate <- tryCatch(evaluate(trial), error = function(e) NULL)
                if (!is.null(candidate) && is.finite(candidate$nll) &&
                    candidate$nll <= result$nll + 1e-9 &&
                    (candidate$nll < result$nll ||
                     max(abs(candidate$gradient)) < max(abs(result$gradient)))) {
                    lambda <- trial
                    result <- candidate
                    accepted <- TRUE
                    break
                }
            }
            if (!accepted) break
        }
        if (max(abs(result$gradient)) / positive_count >= 1e-7) {
            fallback <- TRUE
            optimization <- optim(
                lambda,
                fn = function(value) evaluate(value)$nll / positive_count,
                gr = function(value) -evaluate(value)$gradient / positive_count,
                method = "L-BFGS-B", lower = lower, upper = upper,
                control = list(maxit = 200L, factr = 1e4, pgtol = 1e-8)
            )
            result <- evaluate(optimization$par)
        }
        result$score_residue <- max(abs(result$gradient)) / positive_count
        result$constraint_error <- abs(result$effect - target)
        result$iterations <- iterations
        result$fallback <- fallback
        if (!is.finite(result$score_residue) || result$score_residue > 1e-6) {
            stop("constrained nuisance fit did not converge")
        }
        if (!is.finite(result$constraint_error) || result$constraint_error > 1e-8) {
            stop("abundance effect constraint did not converge")
        }
        chol(result$information)
        result
    }
    candidate <- tryCatch(fit_start(fit$theta), error = function(e) NULL)
    if (!is.null(candidate)) return(candidate)
    candidates <- lapply(c(0.5, 2), function(sigma) {
        start <- fit$theta
        start[p] <- log(sigma)
        tryCatch(fit_start(start), error = function(e) NULL)
    })
    candidates <- Filter(Negate(is.null), candidates)
    if (!length(candidates)) {
        return(.sacat_abundance_profile_numeric_fit(
            fit, y, N, X_b, target, gh
        ))
    }
    candidates[[which.min(vapply(candidates, `[[`, numeric(1), "nll"))]]
}

.sacat_abundance_profile_test <- function(fit, y, N, X_b, background,
                                          variance, gh) {
    restricted <- .sacat_abundance_profile_fit(fit, y, N, X_b, background, gh)
    full_nll <- -sum(zt_beta_loglik_by_sample_inference(
        fit$theta, y, N, fit$solver_diagnostics$X_eta, gh
    ))
    likelihood_ratio <- 2 * (restricted$nll - full_nll)
    if (!is.finite(likelihood_ratio) || likelihood_ratio < -1e-7) {
        stop("constrained likelihood exceeds unrestricted likelihood")
    }
    likelihood_ratio <- max(0, likelihood_ratio)
    if (fit$raw_delta == background) likelihood_ratio <- 0
    gradient <- fit$solver_diagnostics$effect_gradient
    curvature <- .sacat_abundance_equilibrated_solve(
        fit$solver_diagnostics$information, matrix(gradient, ncol = 1L),
        condition_limit = 1e12, backward_tolerance = 1e-8
    )
    if (!isTRUE(curvature$ok)) stop("profile curvature solve failed")
    curvature_variance <- sum(gradient * curvature$solution)
    if (!is.finite(curvature_variance) || curvature_variance <= 0 ||
        !is.finite(variance) || variance <= 0) {
        stop("profile test variance is nonpositive")
    }
    statistic <- likelihood_ratio * curvature_variance / variance
    list(
        p = pchisq(statistic, df = 1, lower.tail = FALSE),
        z = sign(fit$raw_delta - background) * sqrt(statistic),
        statistic = statistic, likelihood_ratio = likelihood_ratio,
        curvature_variance = curvature_variance, fit = restricted
    )
}

.sacat_abundance_profile_tests <- function(result, fits, Y, N, X_b, gh,
                                           cluster = NULL) {
    index <- which(result$formed)
    test_one <- function(j) {
        tryCatch(
            .sacat_abundance_profile_test(
                fits[[j]], as.numeric(Y[, j]), N, X_b,
                result$background[j], result$se[j]^2, gh
            ),
            error = function(e) list(error = conditionMessage(e))
        )
    }
    tests <- .sacat_taxon_lapply(index, test_one, cluster)
    profile_fits <- vector("list", ncol(Y))
    likelihood_ratio <- curvature_variance <- rep(NA_real_, ncol(Y))
    for (k in seq_along(index)) {
        j <- index[k]
        test <- tests[[k]]
        if (!is.null(test$error)) {
            result$formed[j] <- FALSE
            result$reason[j] <- paste0("abundance_profile_failed: ", test$error)
            next
        }
        result$p[j] <- test$p
        result$z[j] <- test$z
        likelihood_ratio[j] <- test$likelihood_ratio
        curvature_variance[j] <- test$curvature_variance
        profile_fits[[j]] <- test$fit
    }
    if (!is.null(result$diagnostics)) {
        table <- result$diagnostics$taxon
        table$corrected_p_value <- result$p
        table$formed <- result$formed
        table$reason <- result$reason
        table$profile_likelihood_ratio <- likelihood_ratio
        table$profile_curvature_variance <- curvature_variance
        table$profile_statistic <- result$z^2
        result$diagnostics$taxon <- table
        names(profile_fits) <- colnames(Y)
        result$diagnostics$profile_fits <- profile_fits
    }
    result$background <- NULL
    result
}
