test_that("the default structural engine forms a regular result", {
    n <- 80L
    group <- rep(c(0, 1), each = n / 2L)
    depth <- rep(c(4000, 5500, 7000, 8500), length.out = n)
    index <- seq_len(n)
    y <- as.integer(7 + ((index * 7L) %% 13L) + 3L * group)
    y[index %% 5L == 0L | index %% 11L == 0L] <- 0L

    fit <- DASRA:::zt_count_structural_test(
        y = y,
        N = depth,
        g = group,
        keep_fit = TRUE
    )

    expect_true(fit$tested)
    expect_true(fit$regular)
    expect_identical(fit$reason, "ok")
    expect_true(is.finite(fit$p) && fit$p >= 0 && fit$p <= 1)
    expect_true(is.finite(fit$diagnostics$score_z))
    expect_gt(fit$diagnostics$V, 0)
    expect_true(all(is.finite(fit$gamma)))
    expect_true(all(fit$gamma >= 0 & fit$gamma <= 1))
    expect_true(all(is.finite(fit$rho)))
    expect_true(all(fit$rho > 0 & fit$rho < 1))
    expect_identical(fit$diagnostics$fit$quadrature$Q, 1001L)
})

test_that("the Q=1001 structural quadrature rule preserves its invariants", {
    gh <- DASRA:::make_structural_gh_rule(1001L)
    midpoint <- (gh$Q + 1L) / 2L
    largest_log_weight <- max(gh$log_weight)
    log_weight_sum <- largest_log_weight + log(sum(exp(
        gh$log_weight - largest_log_weight
    )))

    expect_identical(gh$Q, 1001L)
    expect_length(gh$node, 1001L)
    expect_length(gh$log_weight, 1001L)
    expect_true(all(is.finite(gh$node)))
    expect_true(all(is.finite(gh$log_weight)))
    expect_true(all(diff(gh$node) > 0))
    expect_equal(gh$node, -rev(gh$node), tolerance = 1e-10)
    expect_equal(
        gh$log_weight, rev(gh$log_weight), tolerance = 1e-8
    )
    expect_equal(gh$node[midpoint], 0, tolerance = 1e-12)
    expect_equal(log_weight_sum, 0, tolerance = 1e-12)
    expect_equal(sum(gh$weight), 1, tolerance = 1e-12)
    expect_equal(gh$weight, rev(gh$weight), tolerance = 1e-12)
    expect_equal(sum(gh$weight * gh$node), 0, tolerance = 1e-13)
    expect_equal(sum(gh$weight * gh$node^2), 0.5, tolerance = 1e-12)
})
