## -----------------------------------------------------------------------------------------------------------------------
#| label: setup
library(MASS)
library(ggplot2)
theme_set(theme_bw())
library(lavaan)
library(mirt)
library(splines)


## -----------------------------------------------------------------------------------------------------------------------
#| label: example-1
set.seed(1360)
# Simulation condition
num_obs <- 50  # small sample per group
# Simulate scalar invariant data with same mean
lambda <- c(0.9, 0.7, 0.5)
nu <- c(0, 0, 0)
theta <- diag(1 - lambda^2)
# Group 1
psi1 <- 1.6
alpha1 <- 0
# Group 2
psi2 <- 0.4
alpha2 <- 0.8
# Function for data generation
gendata <- function(nobs, lambda1, lambda2 = lambda1,
                    nu1, nu2 = nu1, theta1, theta2 = theta1,
                    psi1, psi2, alpha1, alpha2) {
    zero_p <- 0 * lambda  # zero vector of length p (for convenience)
    psi_theta1 <- rbind(
        c(psi1, zero_p),
        cbind(zero_p, theta1)
    )
    eta_eps1 <- MASS::mvrnorm(nobs,
        mu = c(alpha1, 0 * lambda1),
        Sigma = psi_theta1, empirical = TRUE
    )
    eta1 <- eta_eps1[, 1]
    eps1 <- eta_eps1[, -1]
    y1 <- t(nu1 + t(tcrossprod(eta1, lambda1) + eps1))
    psi_theta2 <- rbind(
        c(psi2, zero_p),
        cbind(zero_p, theta2)
    )
    eta_eps2 <- MASS::mvrnorm(nobs,
        mu = c(alpha2, 0 * lambda),
        Sigma = psi_theta2, empirical = TRUE
    )
    eta2 <- eta_eps2[, 1]
    eps2 <- eta_eps2[, -1]
    y2 <- t(nu2 + t(tcrossprod(eta2, lambda2) + eps2))
    out <- rbind(cbind(eta1, y1, group = 1), cbind(eta2, y2, group = 2))
    out <- data.frame(out)
    out$group <- factor(out$group)
    colnames(out) <- c("eta", paste0("y", seq_along(lambda)), "group")
    out
}
dat_y <- gendata(num_obs, lambda1 = lambda, nu1 = nu, theta1 = theta,
                 psi1 = psi1, psi2 = psi2, alpha1 = alpha1, alpha2 = alpha2)


## -----------------------------------------------------------------------------------------------------------------------
#| label: fig-ss-inv
#| fig-cap: Sum scores against true latent factor under full invariance.
# Plot
text_range <- c(1, .9, .8, .7)
yrange <- range(rowSums(dat_y[, 2:4]))
ycoords <- yrange[1] + diff(yrange) * text_range
lab_df <- data.frame(
    x = -3,
    y = ycoords,
    par = rep(c("Slope", "Intercept"), 2),
    group = factor(rep(1:2, each = 2))
)
lab_df$val <- c(2.1, 0, 2.1, 0)
p_ss_inv <- ggplot(
    dat_y,
    aes(x = eta, y = rowSums(dat_y[, 2:4]), color = group)
) +
    geom_smooth(aes(linetype = group), method = "lm", se = FALSE,
                alpha = 0.5, fullrange = TRUE) +
    geom_point(aes(shape = group), size = 1.5) +
    labs(x = expression(eta), y = "sum scores") +
    scale_shape(solid = FALSE) +
    geom_label(
        aes(
            x = x, y = y,
            label = paste0(par, " (Group ", group, ") = ", val)
        ),
        data = lab_df,
        hjust = 0, 
        size = 2.5, 
        show.legend = FALSE
    )


## -----------------------------------------------------------------------------------------------------------------------
#| label: example-1-strict-invariance
strict_fit <- cfa(
    "
    f =~ NA * y1 + y2 + y3
    f ~~ c(1.6, 0.4) * f
    ",
    data = dat_y,
    # std.lv = TRUE,
    group = "group",
    group.equal = c("loadings", "intercepts", "residuals"),
    likelihood = "wishart"
)
# summary(strict_fit)


## -----------------------------------------------------------------------------------------------------------------------
#| include: false
# Function for scoring matrix
compute_a_reg <- function(lambda, theta, psi) {
    covy <- lambda %*% psi %*% t(lambda) + theta
    ginvcovy <- MASS::ginv(covy)
    tlam_invcov <- crossprod(lambda, ginvcovy)
    psi %*% tlam_invcov
}
# scoring matrix for Group 1
a1_reg <- compute_a_reg(matrix(lambda), psi = psi1, theta = theta)
slp1_reg <- a1_reg %*% lambda
int1_reg <- (1 - slp1_reg) * alpha1
# scoring matrix for Group 2
a2_reg <- compute_a_reg(matrix(lambda), psi = psi2, theta = theta)
slp2_reg <- a2_reg %*% lambda
int2_reg <- (1 - slp2_reg) * alpha2


## -----------------------------------------------------------------------------------------------------------------------
#| label: fig-fs-inv
#| fig-cap: Regression factor scores against true latent factor under full invariance.
# Regression factor scores
dat_y$fsyr <- lavPredict(strict_fit, assemble = TRUE)$f
# Plot
yrange <- range(dat_y$fsyr)
ycoords <- yrange[1] + diff(yrange) * text_range
lab_df$y <- ycoords
lab_df$val <- round(c(slp1_reg, int1_reg, slp2_reg, int2_reg), 2)
p_fs_inv <- ggplot(
    dat_y,
    aes(x = eta, y = fsyr, color = group)
) +
    geom_smooth(aes(linetype = group), method = "lm", se = FALSE,
                alpha = 0.5, fullrange = TRUE) +
    geom_point(aes(shape = group), size = 1.5) +
    labs(x = expression(eta), y = "regression factor scores") +
    scale_shape(solid = FALSE) +
    geom_label(
        aes(
            x = x, y = y,
            label = paste0(par, " (Group ", group, ") = ", val)
        ),
        data = lab_df,
        hjust = 0, 
        size = 2.5, 
        show.legend = FALSE
    )


## -----------------------------------------------------------------------------------------------------------------------
#| label: example-2
set.seed(1360)
nu2 <- c(0, 0, 0.5)  # new intercepts
dat_ex2 <- gendata(num_obs,
    lambda1 = lambda, nu1 = nu, nu2 = nu2, theta1 = theta,
    psi1 = psi1, psi2 = psi2, alpha1 = alpha1, alpha2 = alpha2
)


## -----------------------------------------------------------------------------------------------------------------------
#| label: fig-ss-pinv
#| fig-cap: Sum scores against true latent factor under partial invariance.
# Plot
yrange <- range(rowSums(dat_ex2[, 2:4]))
ycoords <- yrange[1] + diff(yrange) * text_range
lab_df$y <- ycoords
lab_df$val <- c(2.1, 0, 2.1, 0.5)
p_ss_pinv <- ggplot(
    dat_ex2,
    aes(x = eta, y = rowSums(dat_ex2[, 2:4]), color = group)
) +
    geom_smooth(aes(linetype = group), method = "lm", se = FALSE,
                alpha = 0.5, fullrange = TRUE) +
    geom_point(aes(shape = group), size = 1.5) +
    labs(x = expression(eta), y = "sum scores") +
    scale_shape(solid = FALSE) +
    geom_label(
        aes(
            x = x, y = y,
            label = paste0(par, " (Group ", group, ") = ", val)
        ),
        data = lab_df,
        hjust = 0, 
        size = 2.5, 
        show.legend = FALSE
    )


## -----------------------------------------------------------------------------------------------------------------------
#| label: example-2-partial-scalar-invariance
pscalar_fit <- cfa(
    "
    f =~ NA * y1 + y2 + y3
    f ~~ c(1.6, NA) * f
    ",
    data = dat_ex2,
    # std.lv = TRUE,
    group = "group",
    group.equal = c("loadings", "intercepts", "residuals"),
    group.partial = c("y3~1", "y3~~y3"),
    likelihood = "wishart")
# summary(pscalar_fit)


## -----------------------------------------------------------------------------------------------------------------------
#| label: fig-fs-pinv
#| fig-cap: Regression factor scores against true latent factor under partial invariance.
# Regression factor scores
dat_ex2$fsyr <- lavPredict(pscalar_fit, assemble = TRUE)$f
# Plot
yrange <- range(dat_ex2$fsyr)
ycoords <- yrange[1] + diff(yrange) * text_range
lab_df$y <- ycoords
lab_df$val <- round(c(slp1_reg, int1_reg, slp2_reg, int2_reg), 2)
p_fs_pinv <- ggplot(
    dat_ex2,
    aes(x = eta, y = fsyr, color = group)
) +
    geom_smooth(aes(linetype = group), method = "lm", se = FALSE,
                alpha = 0.5, fullrange = TRUE) +
    geom_point(aes(shape = group), size = 1.5) +
    labs(x = expression(eta), y = "regression factor scores") +
    scale_shape(solid = FALSE) +
    geom_label(
        aes(
            x = x, y = y,
            label = paste0(par, " (Group ", group, ") = ", val)
        ),
        data = lab_df,
        hjust = 0, 
        size = 2.5, 
        show.legend = FALSE
    )


## -----------------------------------------------------------------------------------------------------------------------
#| label: apafg-ss-vs-fs
#| apa-cap: Observed scores against true latent factor under invariance (a and b) or partial invariance (c and d).
ggpubr::ggarrange(p_ss_inv, p_fs_inv, p_ss_pinv, p_fs_pinv, 
                  common.legend = TRUE, 
                  labels = c("(a)", "(b)", "(c)", "(d)"), 
                  font.label = list(size = 12, color = "black", face = "plain"))
# p_ss_inv + theme(legend.position = "top")
# p_fs_inv + theme(legend.position = "top")
# p_ss_pinv + theme(legend.position = "top")
# p_fs_pinv + theme(legend.position = "top")
# col1 <- ggplot() +
#     annotate("text", x = 1, y = 1, label = "Full Invariance") +
#     theme_void()
# col2 <- ggplot() +
#     annotate("text", x = 1, y = 1, label = "Partial Invariance") +
#     theme_void()
# ggpubr::ggarrange(
#     col1, col2,
#     p_ss_inv, p_fs_inv, p_ss_pinv, p_fs_pinv,
#     ncol = 2, nrow = 3,
#     heights = c(0.3, 1, 1),
#     legend = "right",
#     common.legend = TRUE,
# )


## -----------------------------------------------------------------------------------------------------------------------
#| include: false
# scoring matrix for Group 1
a1_reg <- compute_a_reg(matrix(lambda), psi = 1, theta = theta)
slp1_reg <- a1_reg %*% lambda
int1_reg <- (1 - slp1_reg) * 0.5
# scoring matrix for Group 2
a2_reg <- compute_a_reg(matrix(lambda), psi = 1, theta = theta)
slp2_reg <- a2_reg %*% lambda
int2_reg <- (1 - slp2_reg) * 0.5


## -----------------------------------------------------------------------------------------------------------------------
#| label: fig-fs-inv2
# Function for computing factor scores with common distributions
compute_fscore <- function(y, lambda, nu, theta,
                           psi, alpha) {
    y1c <- t(as.matrix(y))
    meany <- lambda %*% alpha + nu
    y1c <- y1c - as.vector(meany)
    a_mat <- compute_a_reg(lambda, psi = psi, theta = theta)
    t(a_mat %*% y1c + as.vector(alpha))
}
strict_pars <- lavInspect(strict_fit, what = "est")
fscores_reg3 <- list(
    `1` = compute_fscore(
        dat_y[dat_y$group == 1, 2:4],
        lambda = strict_pars[[1]]$lambda,
        theta = strict_pars[[1]]$theta,
        nu = strict_pars[[1]]$nu,
        psi = 1,
        alpha = 0.5
    ),
    `2` = compute_fscore(
        dat_y[dat_y$group == 2, 2:4],
        lambda = strict_pars[[2]]$lambda,
        theta = strict_pars[[2]]$theta,
        nu = strict_pars[[2]]$nu,
        psi = 1,
        alpha = 0.5
    )
)
# Add to data
dat_y$fsyr <- do.call(rbind, fscores_reg3)
# Plot
yrange <- range(dat_y$fsyr)
ycoords <- yrange[1] + diff(yrange) * text_range
lab_df$y <- ycoords
lab_df$val <- round(c(slp1_reg, int1_reg, slp2_reg, int2_reg), 2)
p_fs_inv2 <- ggplot(
    dat_y,
    aes(x = eta, y = fsyr, color = group)
) +
    geom_smooth(aes(linetype = group), method = "lm", se = FALSE,
                alpha = 0.5, fullrange = TRUE) +
    geom_point(aes(shape = group), size = 1.5) +
    labs(x = expression(eta), y = "regression factor scores") +
    scale_shape(solid = FALSE) +
    geom_label(
        aes(
            x = x, y = y,
            label = paste0(par, " (Group ", group, ") = ", val)
        ),
        data = lab_df,
        hjust = 0, 
        size = 2.5, 
        show.legend = FALSE
    )


## -----------------------------------------------------------------------------------------------------------------------
#| include: false
# Use a different common distribution
fscores_reg3b <- list(
    `1` = compute_fscore(
        dat_y[dat_y$group == 1, 2:4],
        lambda = strict_pars[[1]]$lambda,
        theta = strict_pars[[1]]$theta,
        nu = strict_pars[[1]]$nu,
        psi = 100,
        alpha = 50
    ),
    `2` = compute_fscore(
        dat_y[dat_y$group == 2, 2:4],
        lambda = strict_pars[[2]]$lambda,
        theta = strict_pars[[2]]$theta,
        nu = strict_pars[[2]]$nu,
        psi = 100,
        alpha = 50
    )
)
# Add to data
dat_y$fsyr <- do.call(rbind, fscores_reg3b)
# Plot
yrange <- range(dat_y$fsyr)
ycoords <- yrange[1] + diff(yrange) * text_range
lab_df$y <- ycoords
lab_df$val <- round(c(slp1_reg, int1_reg, slp2_reg, int2_reg), 2)
p_fs_inv2b <- ggplot(
    dat_y,
    aes(x = eta, y = fsyr, color = group)
) +
    geom_smooth(aes(linetype = group), method = "lm", se = FALSE,
                alpha = 0.5, fullrange = TRUE) +
    geom_point(aes(shape = group), size = 1.5) +
    labs(x = expression(eta), y = "regression factor scores") +
    scale_shape(solid = FALSE) +
    geom_label(
        aes(
            x = x, y = y,
            label = paste0(par, " (Group ", group, ") = ", val)
        ),
        data = lab_df,
        hjust = 0, 
        size = 2.5, 
        show.legend = FALSE
    )
p_fs_inv2b


## -----------------------------------------------------------------------------------------------------------------------
#| label: example-2-new-loadings
set.seed(1360)
# New loading
lambda2 <- lambda
lambda2[3] <- 10
# New intercepts
nu2 <- c(0, 0, 0.5)
# New uniqueness
theta2 <- theta
theta2[3, 3] <- theta2[3, 3] / 2
dat_ex3 <- gendata(num_obs,
    lambda1 = lambda, lambda2 = lambda2, nu1 = nu, nu2 = nu2,
    theta1 = theta, theta2 = theta2,
    psi1 = psi1, psi2 = psi2, alpha1 = alpha1, alpha2 = alpha2
)


## -----------------------------------------------------------------------------------------------------------------------
#| label: example-2-partial-metric-invariance
pmetric_fit <- cfa(
    "
    f =~ y1 + y2 + y3
    ",
    data = dat_ex3,
    std.lv = TRUE,
    group = "group",
    group.equal = c("loadings", "intercepts", "residuals"),
    group.partial = c("f=~y3", "y3~1", "y3~~y3"))
# summary(pmetric_fit)


## -----------------------------------------------------------------------------------------------------------------------
#| include: false
# scoring matrix for Group 1
a1_reg <- compute_a_reg(matrix(lambda), psi = 1, theta = theta)
slp1_reg <- a1_reg %*% lambda
int1_reg <- (1 - slp1_reg) * 0.5
# scoring matrix for Group 2
a2_reg <- compute_a_reg(matrix(lambda2), psi = 1, theta = theta)
slp2_reg <- a2_reg %*% lambda2
int2_reg <- (1 - slp2_reg) * 0.5


## -----------------------------------------------------------------------------------------------------------------------
#| label: fig-fs-pinv3
pmetric_pars <- lavInspect(pmetric_fit, what = "est")
fscores_reg4 <- list(
    `1` = compute_fscore(
        dat_ex3[dat_ex3$group == 1, 2:4],
        lambda = pmetric_pars[[1]]$lambda,
        theta = pmetric_pars[[1]]$theta,
        nu = pmetric_pars[[1]]$nu,
        psi = 1,
        alpha = 0.5
    ),
    `2` = compute_fscore(
        dat_ex3[dat_ex3$group == 2, 2:4],
        lambda = pmetric_pars[[2]]$lambda,
        theta = pmetric_pars[[2]]$theta,
        nu = pmetric_pars[[2]]$nu,
        psi = 1,
        alpha = 0.5
    )
)
# Add to data
dat_ex3$fsyr <- do.call(rbind, fscores_reg4)
# Plot
yrange <- range(dat_ex3$fsyr)
ycoords <- yrange[1] + diff(yrange) * text_range
lab_df$y <- ycoords
lab_df$val <- round(c(slp1_reg, int1_reg, slp2_reg, int2_reg), 2)
p_fs_pinv2 <- ggplot(
    dat_ex3,
    aes(x = eta, y = fsyr, color = group)
) +
    geom_smooth(aes(linetype = group), method = "lm", se = FALSE,
                alpha = 0.5, fullrange = TRUE) +
    geom_point(aes(shape = group), size = 1.5) +
    labs(x = expression(eta), y = "regression factor scores") +
    scale_shape(solid = FALSE) +
    geom_label(
        aes(
            x = x, y = y,
            label = paste0(par, " (Group ", group, ") = ", val)
        ),
        data = lab_df,
        hjust = 0, 
        size = 2.5, 
        show.legend = FALSE
    )


## -----------------------------------------------------------------------------------------------------------------------
#| include: false
# scoring matrix for Group 1
a1_reg <- R2spa:::compute_a_reg(matrix(lambda), psi = 1, theta = theta)
slp1_reg <- a1_reg %*% lambda
int1_reg <- (1 - slp1_reg) * 0.5
# scoring matrix for Group 2
a2_reg <- R2spa:::compute_a_reg(matrix(lambda2), psi = 1, theta = theta)
slp2_reg <- a2_reg %*% lambda2
int2_reg <- (1 - slp2_reg) * 0.5


## -----------------------------------------------------------------------------------------------------------------------
#| label: apafg-fs-same-prior
#| fig-height: 3.6
#| apa-cap: Regression factor scores against true latent factor under (a) full invariance and (b) partial metric invariance when using the same latent distribution.
ggpubr::ggarrange(p_fs_inv2, p_fs_pinv2, 
                  common.legend = TRUE, 
                  labels = c("(a)", "(b)"), 
                  font.label = list(size = 12, color = "black", face = "plain"))
# p_fs_inv2 + theme(legend.position = "top")
# p_fs_pinv2 + theme(legend.position = "top")




## -----------------------------------------------------------------------------------------------------------------------
#| eval: false
## dat_y$grp <- as.integer(dat_y$group) - 1
## write.table(dat_y[c(1:4, 7)],  # use dummy-coded group
##             file = here::here("analysis", "dat_y.dat"),
##             row.names = FALSE,
##             col.names = FALSE)


## -----------------------------------------------------------------------------------------------------------------------
#| include: false
# Fitting the model in lavaan so that the scoring matrix can be computed
dat_y$grp <- as.integer(dat_y$group) - 1
mimic_fit0 <- cfa(
    "
    f =~ NA * y1 + y2 + y3
    lgroup =~ grp
    f ~ lgroup
    f ~~ 1 * f
    f ~ 0
    grp ~ 0
    lgroup ~ NA * 1
    ",
    data = dat_y,
    meanstructure = TRUE)
# summary(mimic_fit0)


## -----------------------------------------------------------------------------------------------------------------------
#| include: false
# mimic_fs0 <- read.table(here::here("analysis", "mimic_fs0.sav"),
#     col.names = c("y1", "y2", "y3", "group", "f", "f_se")
# )
# dat_y$fsyr <- mimic_fs0$f
dat_y$fsyr <- lavPredict(mimic_fit0)[, 1]
# Confirm lavaan has same factor scores as in Mplus
# all.equal(lavPredict(mimic_fit0)[, 1], mimic_fs0$f,
#           tolerance = 1e-3)
# Verify scoring
pars_mimic0 <- lavInspect(mimic_fit0, what = "est")
a_mimic0 <- with(pars_mimic0, {
    inv_I_minus_beta <- solve(diag(2) - beta)
    sig <- inv_I_minus_beta %*% psi %*% t(inv_I_minus_beta)
    solve(lambda %*% sig %*% t(lambda) + theta, lambda %*% sig)
})
mu_mimic0 <- with(pars_mimic0,
    solve(diag(2) - beta, alpha)
)
# crossprod(a_mimic0,
#     t(lavInspect(mimic_fit0, what = "data")) - 
#     c(pars_mimic0$nu + pars_mimic0$lambda %*% mu_mimic0)
# ) + c(mu_mimic0)  # the same
# So the problem is the lambda matrix!


## -----------------------------------------------------------------------------------------------------------------------
#| label: fig-ss-vs-fs-mimic-scalar
slp_mimic0 <- crossprod(a_mimic0, pars_mimic0$lambda)
int_mimic0 <- (diag(2) - slp_mimic0) %*% mu_mimic0
int1_mimic0 <- int_mimic0 + slp_mimic0 %*% c(0, 0)
int2_mimic0 <- int_mimic0 + slp_mimic0 %*% c(0, 1)
yrange <- range(dat_y$fsyr)
ycoords <- yrange[1] + diff(yrange) * text_range
lab_df$y <- ycoords
lab_df$val <- round(c(slp_mimic0[1], int1_mimic0[1],
                      slp_mimic0[1], int2_mimic0[1]), 2)
# Plot
p_fs_mimic0 <- ggplot(
    dat_y,
    aes(x = eta, y = fsyr, color = group)
) +
    geom_smooth(aes(linetype = group), method = "lm", se = FALSE,
                alpha = 0.5, fullrange = TRUE) +
    geom_point(aes(shape = group), size = 1.5) +
    labs(x = expression(eta), y = "regression factor scores") +
    scale_shape(solid = FALSE) +
    geom_label(
        aes(
            x = x, y = y,
            label = paste0(par, " (Group ", group, ") = ", val)
        ),
        data = lab_df,
        hjust = 0, 
        size = 2.5, 
        show.legend = FALSE
    )


## -----------------------------------------------------------------------------------------------------------------------
#| label: mimic-partial-scalar-invariance
dat_ex2$grp <- as.integer(dat_ex2$group) - 1
mimic_fit <- cfa(
    "
    f =~ NA * y1 + y2 + y3
    lgroup =~ grp
    f ~ lgroup
    y3 ~ lgroup
    f ~~ 1 * f
    f ~ 0
    grp ~ 0
    lgroup ~ NA * 1
    ",
    data = dat_ex2)
# summary(mimic_fit)




## -----------------------------------------------------------------------------------------------------------------------
#| label: fig-ss-vs-fs-mimic-pscalar
# mimic_fs <- read.table(here::here("analysis", "mimic_fs.sav"),
#     col.names = c("y1", "y2", "y3", "group", "f", "f_se")
# )
# dat_ex2$fsyr <- mimic_fs$f
dat_ex2$fsyr <- lavPredict(mimic_fit)[, 1]
# Plot
p_fs_mimic <- ggplot(
    dat_ex2,
    aes(x = eta, y = fsyr, color = group)
) +
    geom_smooth(aes(linetype = group), method = "lm", se = FALSE,
                alpha = 0.5, fullrange = TRUE) +
    geom_point(aes(shape = group), size = 1.5) +
    labs(x = expression(eta), y = "regression factor scores") +
    scale_shape(solid = FALSE) +
    geom_label(
        aes(
            x = x, y = y,
            label = paste0(par, " (Group ", group, ") = ", val)
        ),
        data = lab_df,
        hjust = 0, 
        size = 2.5, 
        show.legend = FALSE
    )








## -----------------------------------------------------------------------------------------------------------------------
#| label: apafg-ss-vs-fs-mimic
#| fig-height: 3.6
#| apa-cap: Regression factor scores against true latent factor under MIMIC with (a) invariance and (b) partial scalar invariance.
ggpubr::ggarrange(p_fs_mimic0, p_fs_mimic, 
                  common.legend = TRUE, 
                  labels = c("(a)", "(b)"), 
                  font.label = list(size = 12, color = "black", face = "plain"))
# p_fs_mimic0 + theme(legend.position = "top")
# p_fs_mimic + theme(legend.position = "top")






## -----------------------------------------------------------------------------------------------------------------------
#| label: example-irt
lambda <- rep(0.8, 9)
nu <- rep(0, 9)
theta <- diag(1 - lambda^2)
dat_y <- gendata(500, lambda1 = lambda, nu1 = nu, theta1 = theta,
                 psi1 = psi1, psi2 = psi2, alpha1 = alpha1, alpha2 = alpha2)
# Cut to binary
names_y <- paste0("y", seq_along(lambda))
names_i <- paste0("i", seq_along(lambda))
y_bin <- t(apply(t(dat_y[, names_y]) > 0, MARGIN = 2, FUN = as.integer))
colnames(y_bin) <- names_i
dat_bin <- cbind(dat_y, y_bin)


## -----------------------------------------------------------------------------------------------------------------------
#| label: fig-eap-invariance-different-distributions
mg_irt <- multipleGroup(dat_bin[names_i], group = dat_bin$group,
                        invariance = c("free_mean", "free_var", "slopes", "intercepts"),
                        technical = list(NCYCLES = 5000),
                        verbose = FALSE)
fs_eap <- fscores(mg_irt, method = "EAP")
dat_bin$fs_eap <- fs_eap
p_eap_inv <- ggplot(
    dat_bin,
    aes(x = eta, y = fs_eap, color = group)
) +
    geom_smooth(aes(linetype = group), se = FALSE, alpha = 0.5,
                method = "lm",
                formula = "y ~ ns(x, df = 3)") +
    geom_point(aes(shape = group), size = 0.5, alpha = 0.5) +
    labs(x = expression(eta), y = "EAP scores") +
    scale_shape(solid = FALSE)


## -----------------------------------------------------------------------------------------------------------------------
#| label: fig-eap-invariance-same-distributions
# Use same priors
fs_eap <- fscores(mg_irt, method = "EAP",
                  mean = c(0, 0), cov = c(1, 1))
dat_bin$fs_eap <- fs_eap
p_eap_inv_sp <- p_eap_inv %+% dat_bin


## -----------------------------------------------------------------------------------------------------------------------
#| label: apafg-eap-inv
#| fig-height: 3.6
#| apa-cap: Expected a posteriori (EAP) scores against true latent factor under full invariance, with (a) group-specific and (b) equal means and variances. The fitted curves are obtained using natural cubic splines.
ggpubr::ggarrange(p_eap_inv, p_eap_inv_sp, 
                  common.legend = TRUE, 
                  labels = c("(a)", "(b)"), 
                  font.label = list(size = 12, color = "black", face = "plain"))
# p_eap_inv
# p_eap_inv_sp


## -----------------------------------------------------------------------------------------------------------------------
#| label: example-irt-partial-invariance
# Noninvariance
nu2 <- nu
nu2[7:9] <- 0.5
lambda2 <- lambda
lambda2[7:9] <- 0.4
dat_ex2 <- gendata(500,
    lambda1 = lambda, nu1 = nu, nu2 = nu2, theta1 = theta,
    psi1 = psi1, psi2 = psi2, alpha1 = alpha1, alpha2 = alpha2
)
# Cut to binary
y_ex2_bin <- t(apply(t(dat_ex2[, names_y]) > 0, MARGIN = 2, FUN = as.integer))
colnames(y_ex2_bin) <- names_i
dat_ex2_bin <- cbind(dat_ex2, y_ex2_bin)


## -----------------------------------------------------------------------------------------------------------------------
#| include: false
#| label: fig-eap-pinv
#| fig-cap: Expected a posteriori (EAP) scores against true latent factor under partial invariance.
mg_irt_ex2 <- multipleGroup(dat_ex2_bin[names_i],
    group = dat_ex2_bin$group,
    invariance = c("free_mean", "free_var", names_i[1:6]),
    technical = list(NCYCLES = 5000),
    verbose = FALSE
)
fs_eap_ex2 <- fscores(mg_irt_ex2, method = "EAP")
dat_ex2_bin$fs_eap <- fs_eap_ex2
p_eap_pinv <- p_eap_inv %+% dat_ex2_bin
p_eap_pinv


## -----------------------------------------------------------------------------------------------------------------------
#| include: false
#| label: fig-eap-pinv-same-distributions
#| fig-cap: Expected a posteriori (EAP) scores against true latent factor under partial invariance, with equal priors across groups.
# Use same priors
fs_eap_ex2 <- fscores(mg_irt_ex2,
    method = "EAP",
    mean = c(0, 0), cov = c(1, 1)
)
dat_ex2_bin$fs_eap <- fs_eap_ex2
p_eap_pinv_sp <- p_eap_pinv %+% dat_ex2_bin
p_eap_pinv_sp


## -----------------------------------------------------------------------------------------------------------------------
#| label: simulate-irt
#| include: false
set.seed(2007)
num_rep <- 5000
num_obs <- 100
etac <- MASS::mvrnorm(num_obs, mu = alpha1, Sigma = psi1, empirical = TRUE)
fs_mat <- matrix(nrow = num_obs * 2, ncol = num_rep)
gendata2 <- function(eta1, eta2, lambda1, lambda2 = lambda1,
                     nu1, nu2 = nu1, theta1, theta2 = theta1) {
    eps1 <- MASS::mvrnorm(length(eta1), mu = 0 * lambda1, Sigma = theta1)
    y1 <- t(nu1 + t(tcrossprod(eta1, lambda1) + eps1))
    eps2 <- MASS::mvrnorm(length(eta2), mu = 0 * lambda2, Sigma = theta2)
    y2 <- t(nu2 + t(tcrossprod(eta2, lambda2) + eps2))
    out <- rbind(cbind(eta1, y1, group = 1), cbind(eta2, y2, group = 2))
    out <- data.frame(out)
    out$group <- factor(out$group)
    names_y <- paste0("y", seq_along(lambda))
    colnames(out) <- c("eta", names_y, "group")
    # Cut to binary
    y_bin <- t(apply(t(out[, names_y]) > 0, MARGIN = 2, FUN = as.integer))
    colnames(y_bin) <- paste0("i", seq_along(lambda))
    cbind(out, y_bin)
}


## -----------------------------------------------------------------------------------------------------------------------
#| label: simulate-irt-loop
#| eval: false
## # Uncomment to run simulations
## # for (r in seq_len(num_rep)) {
## #     dat_y <- gendata2(eta1 = etac, eta2 = etac, lambda1 = lambda,
## #     lambda2 = lambda2, nu1 = nu, nu2 = -nu2, theta1 = theta)
## #     names_i <- paste0("y", 1:9)
## #     mg_irt_sim <- multipleGroup(dat_y[names_i],
## #         group = dat_y$group,
## #         invariance = c("free_mean", "free_var", names_i[1:6]),
## #         technical = list(NCYCLES = 5000),
## #         verbose = FALSE
## #     )
## #     fs <- try(fscores(mg_irt_sim, method = "EAP", mean = c(0, 0), cov = c(1, 1)))
## #     if (inherits(fs, "try-error")) { fs <- NA }
## #     fs_mat[, r] <- fs
## # }
## # saveRDS(fs_mat, "fs_mat.RDS")


## -----------------------------------------------------------------------------------------------------------------------
#| include: false
fs_mat <- readRDS(here::here("fs_mat.RDS"))


## -----------------------------------------------------------------------------------------------------------------------
#| label: apafg-eap-pinv-sim
#| apa-cap: Expected a posteriori (EAP) scores against true latent factor under partial invariance across 5,000 replications. The vertical bars show the 10th and 90th percentiles, and the lines show the median across replications.
data.frame(
    eta = rep(etac, num_rep),
    fs = c(fs_mat),
    group = factor(
        rep(rep(1:2, each = 100), num_rep)
    )
) |>
    ggplot(aes(x = eta, y = fs, color = group)) +
    stat_summary(aes(linetype = group), geom = "linerange", 
                 fun.min = \(x) quantile(x, .1),
                 fun.max = \(x) quantile(x, .9),
                 position = position_dodge2(width = 0.05)) +
    stat_summary(aes(color = group), geom = "line", 
                 fun = median,
                 position = position_dodge2(width = 0.05)) +
    # geom_point(size = 0.1, alpha = 0.1) +
    # geom_smooth(aes(linetype = group), se = FALSE) +
    labs(x = expression(eta), y = "EAP scores")
# matplot(etac, fs_mat[1:100, ])
# matplot(etac, fs_mat[101:200, ])

