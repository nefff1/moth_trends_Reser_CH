# SETUP ------------------------------- ########################################
################################################################################.

# We used R version 4.2.0 (R Core Team 2022)

# ... packages #################################################################
################################################################################.

library(tidyverse); theme_set(theme_classic()) # tidyverse version 2.0.0
library(rstan); options(mc.cores = 4) # rstan version 2.26.22 and Stan version 2.26.1
library(brms) # brms version 2.17.0
library(bayestestR) # bayestestR version 0.13.1
library(mgcv) # mgcv version 1.9-1
library(cowplot) # cowplot version 1.1.1
library(ggh4x) # ggh4x version 0.2.8
library(openxlsx) # openxlsx version 4.2.5.2
library(raster) # raster version 3.6-26
library(giscoR) # giscoR version 0.6.0
library(sf) # sf version 1.0-15
library(ggspatial) # ggspatial version 1.1.9

select <- dplyr::select

# ... set global parameters ####################################################
################################################################################.

min_mass <- min(d_mod_z$mass_tot[d_mod_z$mass_tot > 0])

log_plus_trans <- scales::trans_new(name = "log_plus",
                                    transform = \(x) log(x + min_mass, 10),
                                    inverse = \(x) out <- 10 ^ (x) - min_mass)

v_textsize <- c(axis.title = 8, axis.text = 7, legend.title = 8, additional.text = 6) 

n_iter <- 2000 # number of MCMC iterations

# ... define functions #########################################################
################################################################################.

f_analysis_A_height <- function(formula, data_z, scalings, family, iter, seed, 
                                hours_sel = NULL, run_loo = FALSE) {
  
  data_z <- droplevels(data_z)
  
  set.seed(seed)
  
  out <- list()
  
  response <- as.character(formula[2])
  
  
  l_data <- make_standata(update(formula, ". ~ A * height + ."), 
                          data = data_z)
  
  if (length(hours_sel) == 0){
    hours_sel <- NULL
  }
  
  if (!is.null(hours_sel)){
    l_data_sub <- make_standata(update(formula, ". ~ s(active_hours)"),
                                data = data_z[hours_sel, ])
    
    l_data$N_SEL <- length(hours_sel)
    l_data$SEL <-  hours_sel
    l_data$Xs_add <- l_data_sub$Xs
    l_data$nb_add <- 1
    l_data$knots_add <- l_data$knots_1
    dimnames(l_data$knots_add)[[1]] <- "Zs_add_1"
    l_data$Zs_add_1 <- l_data_sub$Zs_1_1
  }
  
  
  if (family == "zero_inflated_negbinomial"){
    if (!is.null(hours_sel)){
      mod <- stan(file = "Stan_Code/Stan_nb_spline_s1p1_r4.stan",
                  data = l_data,
                  chains = 4, cores = 4,
                  iter = iter)
    } else {
      mod <- stan(file = "Stan_Code/Stan_nb_spline_s1_r4.stan",
                  data = l_data,
                  chains = 4, cores = 4,
                  iter = iter)
    }
  } else if (family == "hurdle_gamma"){
    if (!is.null(hours_sel)){
      mod <- stan(file = "Stan_Code/Stan_hg_spline_s1p1_r4.stan",
                  data = l_data,
                  chains = 4, cores = 4,
                  iter = iter)
    } else {
      mod <- stan(file = "Stan_Code/Stan_hg_spline_s1_r4.stan",
                  data = l_data,
                  chains = 4, cores = 4,
                  iter = iter)
    }
  }
  
  # extract parameter data
  par_names <- names(mod)
  pos <- regexpr("\\[", par_names) 
  pos <- ifelse(pos == -1, 200, pos - 1)
  pars <- unique(substr(par_names, 1, pos))
  pars <- pars[pars %in% c("b", "Intercept", "bs", "bs_add", "zs_1_1",
                           "sds_1_1", "zs_2_1", "sds_2_1", "zs_add_1",
                           "sds_add_1", "shape", "zi", "sd_1", "sd_2",
                           "sd_3", "sd_4", "s_1_1", "s_2_1", "s_add_1")]
  fit <- rstan::extract(mod, pars = pars)
  out$fit <- fit
  if (run_loo) out$loo_lm <- loo(fit)
  
  terms <- brmsterms(update(formula, ". ~ A * height + ."))
  
  # predictions for fixed effects
  out$l_pred_fe <- f_pred_fixed(fit = fit, data = data_z, 
                                terms = terms, scalings = scalings)
  
  
  # predictions for smooths
  out$l_pred_sm <- f_pred_smooths(fit = fit, data = data_z, 
                                  l_data = l_data,
                                  terms = terms, scalings = scalings,
                                  hours_sel = hours_sel)
  
  
  out
}

f_analysis_A_height_splev <- function(formula, data_z, scalings, family, iter, seed, 
                                      hours_sel = NULL) {
  
  data_z <- droplevels(data_z)
  
  set.seed(seed)
  
  out <- list()
  
  response <- as.character(formula[2])
  
  
  l_data <- make_standata(update(formula, ". ~ A * height + ."), 
                          data = data_z)
  
  if (length(hours_sel) == 0){
    hours_sel <- NULL
  }
  
  if (!is.null(hours_sel)){
    l_data_sub <- make_standata(update(formula, ". ~ s(active_hours)"),
                                data = data_z[hours_sel, ])
    
    l_data$N_SEL <- length(hours_sel)
    l_data$SEL <-  hours_sel
    l_data$Xs_add <- l_data_sub$Xs
    l_data$nb_add <- 1
    l_data$knots_add <- l_data$knots_1
    dimnames(l_data$knots_add)[[1]] <- "Zs_add_1"
    l_data$Zs_add_1 <- l_data_sub$Zs_1_1
  }
  
  
  
  if (!is.null(hours_sel)){
    mod <- stan(file = "Stan_Code/Stan_nb_spline_s1p1_r4.stan",
                data = l_data,
                chains = 4, cores = 4,
                iter = iter)
  } else {
    mod <- stan(file = "Stan_Code/Stan_nb_spline_s1_r4.stan",
                data = l_data,
                chains = 4, cores = 4,
                iter = iter)
  }
  
  
  
  # extract parameter data
  par_names <- names(mod)
  pos <- regexpr("\\[", par_names) 
  pos <- ifelse(pos == -1, 200, pos - 1)
  pars <- unique(substr(par_names, 1, pos))
  pars <- pars[pars %in% c("b", "Intercept", "bs", "bs_add", "zs_1_1",
                           "sds_1_1", "zs_2_1", "sds_2_1", "zs_add_1",
                           "sds_add_1", "shape", "zi", "sd_1", "sd_2",
                           "sd_3", "sd_4", "s_1_1", "s_2_1", "s_add_1")]
  fit <- rstan::extract(mod, pars = pars)
  d_coefs <- as.data.frame(fit$b)
  
  terms <- brmsterms(update(formula, ". ~ A * height + ."))
  mm_full <- model.matrix(terms$dpars$mu$fe, data = data_z)[, -1]
  
  names(d_coefs) <- colnames(mm_full)
  d_coefs <- d_coefs |> 
    select(A, height, `A:height`) |> 
    bind_cols(data.frame(Intercept = fit$Intercept))
  out$d_coefs <- d_coefs
  
  
  d_coefs_means <- apply(d_coefs, 2, function(x) ci(x) |> as.data.frame()) |> 
    bind_rows(.id = "var") |>
    bind_cols(mean = apply(d_coefs, 2, mean)) |> 
    select(var, mean, CI_low, CI_high, CI)
  out$d_coefs_means <- d_coefs_means
  
  out
}

f_analysis_A_site <- function(formula, data_z, scalings, family, iter, seed, 
                              hours_sel = NULL) {
  
  set.seed(seed)
  
  out <- list()
  
  response <- as.character(formula[2])
  
  
  l_data <- make_standata(update(formula, ". ~ A + ."), 
                          data = data_z)
  
  if (!is.null(hours_sel)){
    update(formula, . ~ s(active_hours, k = 10))
    
    l_data_sub <- tryCatch({make_standata(update(formula, ". ~ s(active_hours)"),
                                          data = data_z[hours_sel, ])},
                           error = function(e) "error")
    
    if (is.character(l_data_sub)){
      hours_sel <- NULL
      warning(paste("No hours active term for", LOC_i))
    } else {
      l_data$N_SEL <- length(hours_sel)
      l_data$SEL <-  hours_sel
      l_data$Xs_add <- l_data_sub$Xs
      l_data$nb_add <- 1
      l_data$knots_add <- l_data$knots_1
      dimnames(l_data$knots_add)[[1]] <- "Zs_add_1"
      l_data$Zs_add_1 <- l_data_sub$Zs_1_1
    }
  }
  
  if (family == "zero_inflated_negbinomial"){
    if (!is.null(hours_sel)){
      mod_lm <- stan(file = "Stan_Code/Stan_nb_spline_s1p1_r1.stan",
                     data = l_data,
                     chains = 4, cores = 4,
                     iter = iter)
      
    } else {
      
      mod_lm <- stan(file = "Stan_Code/Stan_nb_spline_s1_r1.stan",
                     data = l_data,
                     chains = 4, cores = 4,
                     iter = iter)
      
    }
  } else if (family == "hurdle_gamma"){
    if (!is.null(hours_sel)){
      
      mod_lm <- stan(file = "Stan_Copde/Stan_hg_spline_s1p1_r1.stan",
                     data = l_data,
                     chains = 4, cores = 4,
                     iter = iter)
      
    } else {
      
      mod_lm <- stan(file = "Stan_Copde/Stan_hg_spline_s1_r1.stan",
                     data = l_data,
                     chains = 4, cores = 4,
                     iter = iter)
      
    }
  }
  
  # extract parameter data
  par_names <- names(mod_lm)
  pos <- regexpr("\\[", par_names) 
  pos <- ifelse(pos == -1, 200, pos - 1)
  pars <- unique(substr(par_names, 1, pos))
  pars <- pars[pars %in% c("b", "Intercept", "bs", "bs_add", "zs_1_1",
                           "sds_1_1", "zs_2_1", "sds_2_1", "zs_add_1",
                           "sds_add_1", "shape", "zi", "sd_1", 
                           "s_1_1", "s_2_1", "s_add_1")]
  fit_lm <- rstan::extract(mod_lm, pars = pars)
  out$fit_lm <- fit_lm
  
  out
}

f_pred_fixed <- function(fit, data, terms, scalings) {
  
  mm_full <- model.matrix(terms$dpars$mu$fe, data = data)[, -1]
  
  # effects without interactions -----------------------------------------------.
  
  l_out <- list()
  for (var_i in all.vars(terms$dpars$mu$fe)){
    if (is.numeric(deframe(data[, var_i]))) {
      
      d_pred <- data.frame(expl = seq(min(data[, var_i]), max(data[, var_i]),
                                      length.out = 150)) %>% 
        mutate(expl_c = expl - mean(data[, var_i, drop = T]),
               expl_or = expl * scalings$sd[scalings$var == var_i] +
                 scalings$mean[scalings$var == var_i])
      
      
      m_pred <- matrix(rep(fit$Intercept, each = 150), nrow = 150) + 
        matrix(d_pred$expl_c) %*% fit$b[, which(grepl(var_i, colnames(mm_full)) &
                                                  !grepl(":", colnames(mm_full)))]
      
      d_pred <- cbind(d_pred, m_pred)
      
      d_pred <- d_pred %>% 
        select(-c(expl, expl_c)) %>% 
        pivot_longer(cols = -expl_or,
                     names_to = "iter", values_to = "pred") %>%
        rename(!! sym(var_i) := expl_or) %>% 
        mutate(pred_exp = exp(pred)) %>%
        group_by(!! sym(var_i)) %>% 
        summarise(log_lower = ci(pred, .95)$CI_low,
                  log_upper = ci(pred, .95)$CI_high,
                  log_estimate = mean(pred),
                  lower = ci(pred_exp, .95)$CI_low,
                  upper = ci(pred_exp, .95)$CI_high,
                  estimate = mean(pred_exp),
                  .groups = "drop")
      
    } else {
      
      d_pred <- data.frame(expl = unique(data[, var_i]))
      
      comps <- strsplit(paste(deparse(terms$dpars$mu$fe), collapse = ""), split = " \\+ ")[[1]]
      
      comps <- comps[grepl(var_i, comps)]
      
      mm_data <- model.matrix(as.formula(paste("~", comps)), data = data)[, -1, drop = F]
      mm_pred <- model.matrix(as.formula(paste("~", comps)), data = d_pred)[, -1, drop = F]
      
      for (i in seq_len(ncol(mm_pred))){
        mm_pred[, i] <- mm_pred[, i] - mean(mm_data[, i])
      }
      
      m_pred <- matrix(rep(fit$Intercept, each = nrow(d_pred)), nrow = nrow(d_pred)) + 
        as.matrix(mm_pred) %*% t(as.matrix(fit$b[, which(grepl(var_i, colnames(mm_full)) &
                                                           !grepl(":", colnames(mm_full)))]))
      
      d_pred <- cbind(d_pred, m_pred)
      
      d_pred <- d_pred %>% 
        pivot_longer(cols = - !! enquo(var_i),
                     names_to = "iter", values_to = "pred") %>%
        mutate(pred_exp = exp(pred)) %>%
        group_by(!! sym(var_i)) %>% 
        summarise(log_lower = ci(pred, .95)$CI_low,
                  log_upper = ci(pred, .95)$CI_high,
                  log_estimate = mean(pred),
                  lower = ci(pred_exp, .95)$CI_low,
                  upper = ci(pred_exp, .95)$CI_high,
                  estimate = mean(pred_exp),
                  .groups = "drop")
    }
    
    l_out[[var_i]] <- d_pred
  }
  
  # effects with interactions --------------------------------------------------.
  for (int_i in colnames(mm_full)[grepl(":", colnames(mm_full))]) {
    vars_i <- strsplit(int_i, ":")[[1]]
    
    d_pred <- data.frame()
    for (var_i in vars_i){
      
      if (nrow(d_pred) == 0){
        d_pred <- data.frame(expl = seq(min(data[, var_i]), max(data[, var_i]),
                                        length.out = 150)) %>%
          mutate(expl_c = expl - mean(data[, var_i, drop = T]),
                 expl_or = expl * scalings$sd[scalings$var == var_i] +
                   scalings$mean[scalings$var == var_i]) %>% 
          rename(!! sym(var_i) := expl,
                 !! sym(paste0(var_i, "_C")) := expl_c,
                 !! sym(paste0(var_i, "_OR")) := expl_or)
      } else {
        d_pred <- expand_grid(d_pred, 
                              data.frame(expl = seq(min(data[, var_i]), max(data[, var_i]),
                                                    length.out = 150)) %>% 
                                mutate(expl_c = expl - mean(data[, var_i, drop = T]),
                                       expl_or = expl * scalings$sd[scalings$var == var_i] +
                                         scalings$mean[scalings$var == var_i]) %>% 
                                rename(!! sym(var_i) := expl,
                                       !! sym(paste0(var_i, "_C")) := expl_c,
                                       !! sym(paste0(var_i, "_OR")) := expl_or))
      }
      
    }
    
    
    m_pred <- matrix(rep(fit$Intercept, each = nrow(d_pred)), nrow = nrow(d_pred))
    
    # add main effects
    for (var_i in vars_i){
      m_pred <- m_pred +
        matrix(d_pred[, paste0(var_i, "_C"), drop = T]) %*% fit$b[, which(grepl(var_i, colnames(mm_full)) &
                                                                            !grepl(":", colnames(mm_full)))]
    }
    
    # add interactive effects
    m_pred <- m_pred + 
      matrix(d_pred %>% 
               select(ends_with("_C")) %>% 
               apply(., 1, prod)) %*% fit$b[, which(grepl(int_i, colnames(mm_full)))]
    
    
    d_pred <- cbind(d_pred, m_pred)
    
    d_pred <-
      d_pred %>% 
      select(- c(!!! syms(vars_i)),
             - c(!!! syms(paste0(vars_i, "_C")))) %>% 
      pivot_longer(cols = -ends_with("_OR"),
                   names_to = "iter", values_to = "pred") %>%
      rename_all(~ gsub("_OR$", "", .)) %>% 
      mutate(pred_exp = exp(pred)) %>%
      group_by(!!! syms(vars_i)) %>% 
      summarise(log_lower = ci(pred, .95)$CI_low,
                log_upper = ci(pred, .95)$CI_high,
                log_estimate = mean(pred),
                lower = ci(pred_exp, .95)$CI_low,
                upper = ci(pred_exp, .95)$CI_high,
                estimate = mean(pred_exp),
                .groups = "drop")
    
    l_out[[int_i]] <- d_pred
  }
  
  l_out
}

f_pred_smooths <-  function(fit, data, l_data, terms, scalings, hours_sel){
  l_out <- list()
  
  if (!is.null(hours_sel)){
    vars <- c(all.vars(terms$dpars$mu$sm), "active_hours")
  } else {
    vars <- all.vars(terms$dpars$mu$sm)
  }
  
  
  for (var_i in vars){ 
    
    if (var_i == "active_hours"){
      d_target_unscaled <- data[hours_sel, var_i] %>% 
        rename(var_target = !! sym(var_i))
      
      d_target_scaled <- l_data$Xs_add  %>% 
        as.data.frame()
      
      v_bs <- fit$bs_add[, 1]
      v_s <- fit$s_add_1
      
    } else {
      d_target_unscaled <- data[, var_i] %>% 
        rename(var_target = !! sym(var_i))
      
      d_target_scaled <- l_data$Xs[, which(colnames(l_data$Xs) == paste0("s", var_i, "_1"))]  %>% 
        as.data.frame()
      
      count_i <- which(all.vars(terms$dpars$mu$sm) == var_i)
      
      v_bs <- fit$bs[, count_i]
      v_s <- fit[[paste0("s_", count_i, "_1")]]
      
    }
    
    
    names(d_target_scaled) <- "var_target"
    
    sm_unscaled <- smoothCon(s(var_target),
                             data  = d_target_unscaled,
                             knots = NULL,
                             absorb.cons = TRUE, 
                             modCon = 3,
                             diagonal.penalty = TRUE)
    
    sm_scaled <- smoothCon(s(var_target),
                           data = d_target_scaled,
                           knots = NULL,
                           absorb.cons = TRUE, 
                           modCon = 3,
                           diagonal.penalty = TRUE)
    
    d_newdata_unscaled <- data.frame(var_target = seq(min(d_target_unscaled$var_target), 
                                                      max(d_target_unscaled$var_target), 
                                                      length.out = 365))
    pred_sm <- PredictMat(sm_unscaled[[1]], d_newdata_unscaled)
    
    d_newdata_scaled <- data.frame(var_target = pred_sm[, 9])
    
    pred_sm <- pred_sm[, -9][, c(1,8:2)]
    
    pred <- matrix(rep(fit$Intercept, each = nrow(d_newdata_scaled)),
                   nrow = nrow(d_newdata_scaled)) + 
      as.matrix(d_newdata_scaled$var_target) %*%  v_bs + 
      pred_sm %*% t(v_s)
    
    d_pred <- data.frame(d_newdata_unscaled, pred) %>%
      pivot_longer(-var_target, names_to = "iter", values_to = "pred")
    
    d_pred <- d_pred %>%
      mutate(var_target = var_target * scalings$sd[scalings$var == var_i] +
               scalings$mean[scalings$var == var_i]) %>%
      rename(!! sym(var_i) := var_target) %>% 
      mutate(pred_exp = exp(pred)) %>%
      group_by(!! sym(var_i)) %>% 
      summarise(log_lower = ci(pred, .95)$CI_low,
                log_upper = ci(pred, .95)$CI_high,
                log_estimate = mean(pred),
                lower = ci(pred_exp, .95)$CI_low,
                upper = ci(pred_exp, .95)$CI_high,
                estimate = mean(pred_exp),
                .groups = "drop")
    
    l_out[[var_i]] <- d_pred
  }
  l_out
}

f_A_height_change_numbers <- function(fit, data, terms) {
  
  int_i <- "A:height"
  
  mm_full <- model.matrix(terms$dpars$mu$fe, data = data)[, -1]
  
  # also considers main effects, not just interactive effects!
  
  vars_i <- strsplit(int_i, ":")[[1]]
  
  d_pred <- data.frame()
  for (var_i in vars_i){
    
    if (nrow(d_pred) == 0){
      d_pred <- data.frame(expl = seq(min(data[, var_i]), max(data[, var_i]),
                                      length.out = 150)) %>%
        mutate(expl_c = expl - mean(data[, var_i, drop = T])) %>%
        rename(!! sym(var_i) := expl,
               !! sym(paste0(var_i, "_C")) := expl_c)
    } else {
      d_pred <- expand_grid(d_pred, 
                            data.frame(expl = seq(min(data[, var_i]), max(data[, var_i]),
                                                  length.out = 150)) %>% 
                              mutate(expl_c = expl - mean(data[, var_i, drop = T])) %>%
                              rename(!! sym(var_i) := expl,
                                     !! sym(paste0(var_i, "_C")) := expl_c))
    }
    
  }
  
  
  m_pred <- matrix(rep(fit$Intercept, each = nrow(d_pred)), nrow = nrow(d_pred))
  
  # add main effects
  for (var_i in vars_i){
    m_pred <- m_pred +
      matrix(d_pred[, paste0(var_i, "_C"), drop = T]) %*% fit$b[, which(grepl(var_i, colnames(mm_full)) &
                                                                          !grepl(":", colnames(mm_full)))]
  }
  
  # add interactive effects
  m_pred <- m_pred + 
    matrix(d_pred %>% 
             select(ends_with("_C")) %>% 
             apply(., 1, prod)) %*% fit$b[, which(grepl(int_i, colnames(mm_full)))]
  
  
  d_pred <- cbind(d_pred, m_pred)
  
  
  d_pred %>% 
    filter(A %in% range(A),
           height %in% range(height)) %>% 
    mutate(height = ifelse(height == min(height), "lowest", "highest")) %>% 
    select(- c(!!! syms(paste0(vars_i, "_C")))) %>% 
    arrange(height, A) %>% 
    select(-A) %>% 
    pivot_longer(-height, names_to = "iter", values_to = "pred") %>% 
    group_by(height, iter) %>% 
    summarise(diff_mean = mean(pred[2] - pred[1]),
              .groups = "drop") %>% 
    group_by(height) %>% 
    summarise(diff_log_estimate = mean(diff_mean),
              diff_log_lower = ci(diff_mean)$CI_low,
              diff_log_upper = ci(diff_mean)$CI_high,
              factor_estimate = mean(exp(diff_mean)),
              factor_lower = ci(exp(diff_mean))$CI_low,
              factor_upper = ci(exp(diff_mean))$CI_high,
              .groups = "drop")
  
}

f_apply_Rhat <- function(fit){
  # Rhat function basend on rstan package
  f_Rhat <- \(sims) {
    
    f_tomatrix <- \(obj_draws){
      matrix(as.numeric(obj_draws), ncol = 8, byrow = F) # 4 chains split in two
    }
    
    f_rhat_rfun <- \(sims) {
      chains <- ncol(sims)
      n_samples <- nrow(sims)
      chain_mean <- numeric(chains)
      chain_var <- numeric(chains)
      for (i in seq_len(chains)) {
        chain_mean[i] <- mean(sims[, i])
        chain_var[i] <- var(sims[, i])
      }
      var_between <- n_samples * var(chain_mean)
      var_within <- mean(chain_var)
      sqrt((var_between/var_within + n_samples - 1)/n_samples)
    }
    
    f_z_scale <- \(x){
      S <- length(x)
      r <- rank(x, ties.method = 'average')
      z <- qnorm((r - 1/2)/S)
      z[is.na(x)] <- NA
      if (!is.null(dim(x))) {
        z <- array(z, dim = dim(x), dimnames = dimnames(x))
      }
      z
    }
    
    bulk_rhat <- f_rhat_rfun(f_z_scale(f_tomatrix(sims)))
    sims_folded <- abs(sims - median(sims))
    tail_rhat <- f_rhat_rfun(f_z_scale(f_tomatrix(sims_folded)))
    max(bulk_rhat, tail_rhat)
  }
  
  d_Rhat <- data.frame()
  for (var_i in names(fit)[grepl("Intercept|^b|^bs|^s_", names(fit))]){
    if (length(dim(fit[[var_i]])) > 1){ 
      for (i in seq_len(ncol(fit[[var_i]]))){
        d_Rhat <- d_Rhat |> 
          bind_rows(data.frame(var = paste(var_i, i, sep = "_"),
                               rhat = f_Rhat(fit[[var_i]][, i])))
      }
    } else {
      d_Rhat <- d_Rhat |> 
        bind_rows(data.frame(var = var_i,
                             rhat = f_Rhat(fit[[var_i]])))
    }
  }
  
  d_Rhat
}

f_A_height_plot <- function(pred,
                            data_raw = NULL, mean_yday = NULL, response = NULL,
                            quantiles = c(0, .25, .5, .75, 1)){
  
  # simplified version to take quantiles
  sel <- round(length(unique(pred$height)) * quantiles)
  sel[sel == 0] <- 1
  sel_height <- unique(pred$height)[sel]
  
  pred <- pred %>% 
    filter(height %in% sel_height) %>% 
    rowwise() %>%
    mutate(height_cat = paste0("Q-", quantiles[which(sel_height == height)] * 100, "%")) %>% 
    ungroup() |> 
    mutate(height_cat = factor(height_cat, levels = unique(height_cat),
                               labels = paste0(unique(height_cat), "\n(   ", round(sel_height), "m)")))
  
  
  if (!is.null(mean_yday)){
    pred <- pred %>%
      mutate(Date = as.Date(format(date_decimal(A)), "%Y-%m-%d") + 
               mean_yday) 
  } else {
    pred <- pred %>%
      mutate(Date = A) 
  }
  
  p <- pred %>%
    ggplot(aes(x = Date, col = height_cat))
  
  if (!is.null(data_raw)){
    data_raw <- data_raw %>%
      mutate(Date = Samplingdate) 
    
    p <- p +
      geom_point(data = data_raw, aes_string(y = response), 
                 alpha = .1, size = .4, col = "salmon4")
  }
  
  p <- p +
    geom_ribbon(aes(ymin = lower, ymax = upper), 
                fill = "grey20", alpha = .3, lty = 2) +
    geom_line(aes(y = exp(log_estimate)), size = 1) + # more stable mean estimates
    scale_colour_manual(values = RColorBrewer::brewer.pal(length(quantiles),
                                                          "RdYlBu"))
  
  
  if (!is.null(response)){
    if (response == "mass_tot"){
      p <-  p +
        scale_y_continuous(trans = log_plus_trans)
    } else {
      p <-  p +
        scale_y_continuous(trans = "log1p")
    }
  }
  
  p +
    labs(colour = "Elevation", x = "Year")
}

f_A_height_plot_comb <- function(l_pred,
                                 response = NULL,
                                 quantiles = c(0, .25, .5, .75, 1),
                                 name = NULL){
  
  # simplified version to take quantiles
  sel <- round(length(unique(l_pred[[1]]$height)) * quantiles)
  sel[sel == 0] <- 1
  sel_height <- unique(l_pred[[1]]$height)[sel]
  
  
  d_plot <- lapply(l_pred,
                   function(x) x %>% 
                     filter(height %in% sel_height) %>% 
                     rowwise() %>%
                     mutate(height_cat = paste0("Q-", quantiles[which(sel_height == height)] * 100, "%")) %>% 
                     ungroup() |> 
                     mutate(height_cat = factor(height_cat, levels = unique(height_cat),
                                                labels = paste0(unique(height_cat), 
                                                                "\n(   ", round(sel_height), "m)")))) %>% 
    bind_rows(.id = "Trait") %>% 
    mutate(Trait = factor(Trait, levels = names(l_pred)),
           name = name)
  
  
  p <-
    d_plot %>%
    ggplot(aes(x = A, col = height_cat)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), 
                fill = "grey20", alpha = .3, lty = 2) +
    geom_line(aes(y = exp(log_estimate)), size = 1) +
    scale_colour_manual(values = RColorBrewer::brewer.pal(length(quantiles),
                                                          "RdYlBu"))
  
  if (!is.null(name)) {
    p <- p + facet_nested_wrap(~ name + Trait, scales = "free_y", nrow = 1)
  } else {
    p <- p + facet_wrap(~ Trait, scales = "free_y", nrow = 1)
  }
  
  if (!is.null(response)){
    if (response == "mass_tot"){
      p <-  p +
        scale_y_continuous(trans = log_plus_trans)
    } else {
      p <-  p +
        scale_y_continuous(trans = "log1p")
    }
  }
  
  p +
    labs(colour = "Elevation", x = "Year")
}

f_plot_pred <- function(x, data, response, hours_sel = NULL, line.size = 1){
  var_i <- names(x)[1]
  
  label <- ifelse(var_i %in% names(v_covlabels),
                  v_covlabels[var_i],
                  var_i)
  
  if (var_i == "traptype"){
    x <- x %>% 
      mutate(traptype = factor(traptype, levels = names(v_traplabels),
                               labels = v_traplabels))
    data <- data %>% 
      mutate(traptype = factor(traptype, levels = names(v_traplabels),
                               labels = v_traplabels))
  } else if (var_i == "bulbtype"){
    x <- x %>% 
      mutate(bulbtype = factor(bulbtype, levels = names(v_bulblabels),
                               labels = v_bulblabels))
    data <- data %>% 
      mutate(bulbtype = factor(bulbtype, levels = names(v_bulblabels),
                               labels = v_bulblabels))
  }
  
  
  if (var_i == "active_hours"){
    p <- ggplot(x, aes_string(x = var_i))  +
      geom_point(data = data[hours_sel, ], aes_string(y = response), 
                 alpha = .1, size = .4, col = "salmon4") + 
      geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey20", alpha = .5) +
      geom_line(aes(y = estimate), size = line.size, colour = "darkgreen") 
  } else if (is.numeric(x[, var_i, drop = T])){
    p <- ggplot(x, aes_string(x = var_i)) +
      geom_point(data = data, aes_string(y = response), 
                 alpha = .1, size = .4, col = "salmon4") + 
      geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey20", alpha = .5) +
      geom_line(aes(y = estimate), size = line.size, colour = "darkgreen")
  } else {
    p <- ggplot(x, aes_string(x = var_i)) +
      geom_jitter(data = data, aes_string(y = response), 
                  alpha = .05, size = .4, col = "salmon4", width = .25) + 
      geom_segment(aes_string(y = "lower", yend = "upper", xend = var_i), 
                   col = "grey20", 
                   arrow = arrow(angle = 90, ends = "both", length = unit(.1, "inches"))) +
      geom_point(aes(y = estimate), size = 2, colour = "darkgreen")
  }
  
  p <- p +
    ylab("Estimate") +
    xlab(label) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    coord_cartesian(ylim = c(0, NA))
  
  if (response == "mass_tot"){
    p <-  p +
      scale_y_continuous(trans = log_plus_trans)
  } else {
    p <-  p +
      scale_y_continuous(trans = "log1p")
  }
  
  p
}

f_extract_slopes <- function(fit){
  formula = SCcorr_ric ~
    s(yday) + P_2day + T_2day +
    C(traptype, "contr.sum") +
    C(bulbtype, "contr.sum") +
    n_trap +
    C(sample_previous, "contr.sum") +
    (1 | spattemp_cluster) + (1 | trap_ID) +
    (1 | night_ID) +
    (1 | trap_ID_A)
  terms <- brmsterms(update(formula, ". ~ A * height + ."))
  mm_full <- model.matrix(terms$dpars$mu$fe, data = d_mod_z)[, -1]
  
  d_coefs <- as.data.frame(fit$b)
  
  names(d_coefs) <- colnames(mm_full)
  d_coefs <- d_coefs |> 
    select(A, height, `A:height`) |> 
    bind_cols(data.frame(Intercept = fit$Intercept))
  
  d_coefs
}


f_boot <- function(iter, dat) {
  sel <- sample(1:nrow(dat), nrow(dat), replace = T)
  out <- which(sort(dat$mean[sel]) > 0)[1] - .5
  
  out <- ifelse(is.na(out), nrow(dat) + .5, out)
  
  out
}

# RUN MODELS -------------------- ##############################################
################################################################################.

# ... overall models ###########################################################
################################################################################.

# abundance --------------------------------------------------------------------.
l_abu_A <- f_analysis_A_height(formula = abu_tot ~
                                 s(yday) + P_2day + T_2day +
                                 C(traptype, "contr.sum") +
                                 C(bulbtype, "contr.sum") +
                                 n_trap +
                                 C(sample_previous, "contr.sum") +
                                 (1 | spattemp_cluster) +
                                 (1 | LOC) +
                                 (1 | night_ID) +
                                 (1 | trap_ID_A),
                               data_z = d_mod_z,
                               scalings = filter(d_scalings, data == "full"),
                               family = "zero_inflated_negbinomial",
                               hours_sel = which(d_mod_z$traptype == "p" & d_mod_z$hours_data),
                               iter = n_iter, seed = 126)

# richness ---------------------------------------------------------------------.
l_ric_A <- f_analysis_A_height(formula = SCcorr_ric ~
                                 s(yday) + P_2day + T_2day +
                                 C(traptype, "contr.sum") +
                                 C(bulbtype, "contr.sum") +
                                 n_trap +
                                 C(sample_previous, "contr.sum") +
                                 (1 | spattemp_cluster) + (1 | LOC) +
                                 (1 | night_ID) +
                                 (1 | trap_ID_A),
                               data_z = d_mod_z,
                               scalings = filter(d_scalings, data == "full"),
                               family = "hurdle_gamma",
                               hours_sel = which(d_mod_z$traptype == "p" & d_mod_z$hours_data),
                               iter = n_iter, seed = 8989)

# biomass ----------------------------------------------------------------------.
l_mass_A <- f_analysis_A_height(formula = mass_tot ~
                                  s(yday) + P_2day + T_2day +
                                  C(traptype, "contr.sum") +
                                  C(bulbtype, "contr.sum") +
                                  n_trap +
                                  C(sample_previous, "contr.sum") +
                                  (1 | spattemp_cluster) +
                                  (1 | LOC) +
                                  (1 | night_ID) +
                                  (1 | trap_ID_A),
                                data_z = d_mod_z,
                                scalings = filter(d_scalings, data == "full"),
                                family = "hurdle_gamma",
                                hours_sel = which(d_mod_z$traptype == "p" & d_mod_z$hours_data),
                                iter = n_iter, seed = 811)

# ... per body size group ######################################################
################################################################################.

# abundance --------------------------------------------------------------------.
l_abu_A_small <- f_analysis_A_height(formula = abu_tot ~
                                       s(yday) + P_2day + T_2day +
                                       C(traptype, "contr.sum") +
                                       C(bulbtype, "contr.sum") +
                                       n_trap +
                                       C(sample_previous, "contr.sum") +
                                       (1 | spattemp_cluster) +
                                       (1 | LOC) +
                                       (1 | night_ID) +
                                       (1 | trap_ID_A),
                                     data_z = d_mod_mass_z %>% filter(mass_cat == "small"),
                                     scalings = filter(d_scalings, data == "full"),
                                     family = "zero_inflated_negbinomial",
                                     hours_sel = which(d_mod_mass_z$traptype[d_mod_mass_z$mass_cat == "small"] == "p" &
                                                         d_mod_mass_z$hours_data[d_mod_mass_z$mass_cat == "small"]),
                                     iter = n_iter, seed = 323)
l_abu_A_medium <- f_analysis_A_height(formula = abu_tot ~
                                        s(yday) + P_2day + T_2day +
                                        C(traptype, "contr.sum") +
                                        C(bulbtype, "contr.sum") +
                                        n_trap +
                                        C(sample_previous, "contr.sum") +
                                        (1 | spattemp_cluster) +
                                        (1 | LOC) +
                                        (1 | night_ID) +
                                        (1 | trap_ID_A),
                                      data_z = d_mod_mass_z %>% filter(mass_cat == "medium"),
                                      scalings = filter(d_scalings, data == "full"),
                                      family = "zero_inflated_negbinomial",
                                      hours_sel = which(d_mod_mass_z$traptype[d_mod_mass_z$mass_cat == "medium"] == "p" &
                                                          d_mod_mass_z$hours_data[d_mod_mass_z$mass_cat == "medium"]),
                                      iter = n_iter, seed = 963)
l_abu_A_large <- f_analysis_A_height(formula = abu_tot ~
                                       s(yday) + P_2day + T_2day +
                                       C(traptype, "contr.sum") +
                                       C(bulbtype, "contr.sum") +
                                       n_trap +
                                       C(sample_previous, "contr.sum") +
                                       (1 | spattemp_cluster) +
                                       (1 | LOC) +
                                       (1 | night_ID) +
                                       (1 | trap_ID_A),
                                     data_z = d_mod_mass_z %>% filter(mass_cat == "large"),
                                     scalings = filter(d_scalings, data == "full"),
                                     family = "zero_inflated_negbinomial",
                                     hours_sel = which(d_mod_mass_z$traptype[d_mod_mass_z$mass_cat == "large"] == "p" &
                                                         d_mod_mass_z$hours_data[d_mod_mass_z$mass_cat == "large"]),
                                     iter = n_iter, seed = 97)

# richness ---------------------------------------------------------------------.
l_ric_A_small <- f_analysis_A_height(formula = SCcorr_ric ~
                                       s(yday) + P_2day + T_2day +
                                       C(traptype, "contr.sum") +
                                       C(bulbtype, "contr.sum") +
                                       n_trap +
                                       C(sample_previous, "contr.sum") +
                                       (1 | spattemp_cluster) +
                                       (1 | LOC) +
                                       (1 | night_ID) +
                                       (1 | trap_ID_A),
                                     data_z = d_mod_mass_z %>% filter(mass_cat == "small"),
                                     scalings = filter(d_scalings, data == "full"),
                                     family = "hurdle_gamma",
                                     hours_sel = which(d_mod_mass_z$traptype[d_mod_mass_z$mass_cat == "small"] == "p" &
                                                         d_mod_mass_z$hours_data[d_mod_mass_z$mass_cat == "small"]),
                                     iter = n_iter, seed = 8460)
l_ric_A_medium <- f_analysis_A_height(formula = SCcorr_ric ~
                                        s(yday) + P_2day + T_2day +
                                        C(traptype, "contr.sum") +
                                        C(bulbtype, "contr.sum") +
                                        n_trap +
                                        C(sample_previous, "contr.sum") +
                                        (1 | spattemp_cluster) +
                                        (1 | LOC) +
                                        (1 | night_ID) +
                                        (1 | trap_ID_A),
                                      data_z = d_mod_mass_z %>% filter(mass_cat == "medium"),
                                      scalings = filter(d_scalings, data == "full"),
                                      family = "hurdle_gamma",
                                      hours_sel = which(d_mod_mass_z$traptype[d_mod_mass_z$mass_cat == "medium"] == "p" &
                                                          d_mod_mass_z$hours_data[d_mod_mass_z$mass_cat == "medium"]),
                                      iter = n_iter, seed = 510)
l_ric_A_large <- f_analysis_A_height(formula = SCcorr_ric ~
                                       s(yday) + P_2day + T_2day +
                                       C(traptype, "contr.sum") +
                                       C(bulbtype, "contr.sum") +
                                       n_trap +
                                       C(sample_previous, "contr.sum") +
                                       (1 | spattemp_cluster) +
                                       (1 | LOC) +
                                       (1 | night_ID) +
                                       (1 | trap_ID_A),
                                     data_z = d_mod_mass_z %>% filter(mass_cat == "large"),
                                     scalings = filter(d_scalings, data == "full"),
                                     family = "hurdle_gamma",
                                     hours_sel = which(d_mod_mass_z$traptype[d_mod_mass_z$mass_cat == "large"] == "p" &
                                                         d_mod_mass_z$hours_data[d_mod_mass_z$mass_cat == "large"]),
                                     iter = n_iter, seed = 879)

# biomass ----------------------------------------------------------------------.
l_mass_A_small <- f_analysis_A_height(formula = mass_tot ~
                                        s(yday) + P_2day + T_2day +
                                        C(traptype, "contr.sum") +
                                        C(bulbtype, "contr.sum") +
                                        n_trap +
                                        C(sample_previous, "contr.sum") +
                                        (1 | spattemp_cluster) +
                                        (1 | LOC) +
                                        (1 | night_ID) +
                                        (1 | trap_ID_A),
                                      data_z = d_mod_mass_z %>% filter(mass_cat == "small"),
                                      scalings = filter(d_scalings, data == "full"),
                                      family = "hurdle_gamma",
                                      hours_sel = which(d_mod_mass_z$traptype[d_mod_mass_z$mass_cat == "small"] == "p" &
                                                          d_mod_mass_z$hours_data[d_mod_mass_z$mass_cat == "small"]),
                                      iter = n_iter, seed = 2512)
l_mass_A_medium <- f_analysis_A_height(formula = mass_tot ~
                                         s(yday) + P_2day + T_2day +
                                         C(traptype, "contr.sum") +
                                         C(bulbtype, "contr.sum") +
                                         n_trap +
                                         C(sample_previous, "contr.sum") +
                                         (1 | spattemp_cluster) +
                                         (1 | LOC) +
                                         (1 | night_ID) +
                                         (1 | trap_ID_A),
                                       data_z = d_mod_mass_z %>% filter(mass_cat == "medium"),
                                       scalings = filter(d_scalings, data == "full"),
                                       family = "hurdle_gamma",
                                       hours_sel = which(d_mod_mass_z$traptype[d_mod_mass_z$mass_cat == "medium"] == "p" &
                                                           d_mod_mass_z$hours_data[d_mod_mass_z$mass_cat == "medium"]),
                                       iter = n_iter, seed = 1717)
l_mass_A_large <- f_analysis_A_height(formula = mass_tot ~
                                        s(yday) + P_2day + T_2day +
                                        C(traptype, "contr.sum") +
                                        C(bulbtype, "contr.sum") +
                                        n_trap +
                                        C(sample_previous, "contr.sum") +
                                        (1 | spattemp_cluster) +
                                        (1 | LOC) +
                                        (1 | night_ID) +
                                        (1 | trap_ID_A),
                                      data_z = d_mod_mass_z %>% filter(mass_cat == "large"),
                                      scalings = filter(d_scalings, data == "full"),
                                      family = "hurdle_gamma",
                                      hours_sel = which(d_mod_mass_z$traptype[d_mod_mass_z$mass_cat == "large"] == "p" &
                                                          d_mod_mass_z$hours_data[d_mod_mass_z$mass_cat == "large"]),
                                      iter = n_iter, seed = 112)

# ... per temperature niche group ##############################################
################################################################################.

# abundance --------------------------------------------------------------------.
l_abu_A_cold <- f_analysis_A_height(formula = abu_tot ~
                                      s(yday) + P_2day + T_2day +
                                      C(traptype, "contr.sum") +
                                      C(bulbtype, "contr.sum") +
                                      n_trap +
                                      C(sample_previous, "contr.sum") +
                                      (1 | spattemp_cluster) +
                                      (1 | LOC) +
                                      (1 | night_ID) +
                                      (1 | trap_ID_A),
                                    data_z = d_mod_Tavg_z %>% filter(Tavg_mean_cat == "cold"),
                                    scalings = filter(d_scalings, data == "full"),
                                    family = "zero_inflated_negbinomial",
                                    hours_sel = which(d_mod_Tavg_z$traptype[d_mod_Tavg_z$Tavg_mean_cat == "cold"] == "p" &
                                                        d_mod_Tavg_z$hours_data[d_mod_Tavg_z$Tavg_mean_cat == "cold"]),
                                    iter = n_iter, seed = 1688)
l_abu_A_intermediate <- f_analysis_A_height(formula = abu_tot ~ height +
                                              s(yday) + P_2day + T_2day +
                                              C(traptype, "contr.sum") +
                                              C(bulbtype, "contr.sum") +
                                              n_trap +
                                              C(sample_previous, "contr.sum") +
                                              (1 | spattemp_cluster) +
                                              (1 | LOC) +
                                              (1 | night_ID) +
                                              (1 | trap_ID_A),
                                            data_z = d_mod_Tavg_z %>% filter(Tavg_mean_cat == "intermediate"),
                                            scalings = filter(d_scalings, data == "full"),
                                            family = "zero_inflated_negbinomial",
                                            hours_sel = which(d_mod_Tavg_z$traptype[d_mod_Tavg_z$Tavg_mean_cat == "intermediate"] == "p" &
                                                                d_mod_Tavg_z$hours_data[d_mod_Tavg_z$Tavg_mean_cat == "intermediate"]),
                                            iter = n_iter, seed = 190)
l_abu_A_warm <- f_analysis_A_height(formula = abu_tot ~ height +
                                      s(yday) + P_2day + T_2day +
                                      C(traptype, "contr.sum") +
                                      C(bulbtype, "contr.sum") +
                                      n_trap +
                                      C(sample_previous, "contr.sum") +
                                      (1 | spattemp_cluster) +
                                      (1 | LOC) +
                                      (1 | night_ID) +
                                      (1 | trap_ID_A),
                                    data_z = d_mod_Tavg_z %>% filter(Tavg_mean_cat == "warm"),
                                    scalings = filter(d_scalings, data == "full"),
                                    family = "zero_inflated_negbinomial",
                                    hours_sel = which(d_mod_Tavg_z$traptype[d_mod_Tavg_z$Tavg_mean_cat == "warm"] == "p" &
                                                        d_mod_Tavg_z$hours_data[d_mod_Tavg_z$Tavg_mean_cat == "warm"]),
                                    iter = n_iter, seed = 1881)

# richness ---------------------------------------------------------------------.
l_ric_A_cold <- f_analysis_A_height(formula = SCcorr_ric ~
                                      s(yday) + P_2day + T_2day +
                                      C(traptype, "contr.sum") +
                                      C(bulbtype, "contr.sum") +
                                      n_trap +
                                      C(sample_previous, "contr.sum") +
                                      (1 | spattemp_cluster) +
                                      (1 | LOC) +
                                      (1 | night_ID) +
                                      (1 | trap_ID_A),
                                    data_z = d_mod_Tavg_z %>% filter(Tavg_mean_cat == "cold"),
                                    scalings = filter(d_scalings, data == "full"),
                                    family = "hurdle_gamma",
                                    hours_sel = which(d_mod_Tavg_z$traptype[d_mod_Tavg_z$Tavg_mean_cat == "cold"] == "p" &
                                                        d_mod_Tavg_z$hours_data[d_mod_Tavg_z$Tavg_mean_cat == "cold"]),
                                    iter = n_iter, seed = 871)
l_ric_A_intermediate <- f_analysis_A_height(formula = SCcorr_ric ~
                                              s(yday) + P_2day + T_2day +
                                              C(traptype, "contr.sum") +
                                              C(bulbtype, "contr.sum") +
                                              n_trap +
                                              C(sample_previous, "contr.sum") +
                                              (1 | spattemp_cluster) +
                                              (1 | LOC) +
                                              (1 | night_ID) +
                                              (1 | trap_ID_A),
                                            data_z = d_mod_Tavg_z %>% filter(Tavg_mean_cat == "intermediate"),
                                            scalings = filter(d_scalings, data == "full"),
                                            family = "hurdle_gamma",
                                            hours_sel = which(d_mod_Tavg_z$traptype[d_mod_Tavg_z$Tavg_mean_cat == "intermediate"] == "p" &
                                                                d_mod_Tavg_z$hours_data[d_mod_Tavg_z$Tavg_mean_cat == "intermediate"]),
                                            iter = n_iter, seed = 111)
l_ric_A_warm <- f_analysis_A_height(formula = SCcorr_ric ~
                                      s(yday) + P_2day + T_2day +
                                      C(traptype, "contr.sum") +
                                      C(bulbtype, "contr.sum") +
                                      n_trap +
                                      C(sample_previous, "contr.sum") +
                                      (1 | spattemp_cluster) +
                                      (1 | LOC) +
                                      (1 | night_ID) +
                                      (1 | trap_ID_A),
                                    data_z = d_mod_Tavg_z %>% filter(Tavg_mean_cat == "warm"),
                                    scalings = filter(d_scalings, data == "full"),
                                    family = "hurdle_gamma",
                                    hours_sel = which(d_mod_Tavg_z$traptype[d_mod_Tavg_z$Tavg_mean_cat == "warm"] == "p" &
                                                        d_mod_Tavg_z$hours_data[d_mod_Tavg_z$Tavg_mean_cat == "warm"]),
                                    iter = n_iter, seed = 901)

# biomass ----------------------------------------------------------------------.
l_mass_A_cold <- f_analysis_A_height(formula = mass_tot ~  s(yday) + P_2day + T_2day +
                                       C(traptype, "contr.sum") +
                                       C(bulbtype, "contr.sum") +
                                       n_trap +
                                       C(sample_previous, "contr.sum") +
                                       (1 | spattemp_cluster) +
                                       (1 | LOC) +
                                       (1 | night_ID) +
                                       (1 | trap_ID_A),
                                     data_z = d_mod_Tavg_z %>% filter(Tavg_mean_cat == "cold"),
                                     scalings = filter(d_scalings, data == "full"),
                                     family = "hurdle_gamma",
                                     hours_sel = which(d_mod_Tavg_z$traptype[d_mod_Tavg_z$Tavg_mean_cat == "cold"] == "p" &
                                                         d_mod_Tavg_z$hours_data[d_mod_Tavg_z$Tavg_mean_cat == "cold"]),
                                     iter = n_iter, seed = 964)
l_mass_A_intermediate <- f_analysis_A_height(formula = mass_tot ~  s(yday) + P_2day + T_2day +
                                               C(traptype, "contr.sum") +
                                               C(bulbtype, "contr.sum") +
                                               n_trap +
                                               C(sample_previous, "contr.sum") +
                                               (1 | spattemp_cluster) +
                                               (1 | LOC) +
                                               (1 | night_ID) +
                                               (1 | trap_ID_A),
                                             data_z = d_mod_Tavg_z %>% filter(Tavg_mean_cat == "intermediate"),
                                             scalings = filter(d_scalings, data == "full"),
                                             family = "hurdle_gamma",
                                             hours_sel = which(d_mod_Tavg_z$traptype[d_mod_Tavg_z$Tavg_mean_cat == "intermediate"] == "p" &
                                                                 d_mod_Tavg_z$hours_data[d_mod_Tavg_z$Tavg_mean_cat == "intermediate"]),
                                             iter = n_iter, seed = 111)
l_mass_A_warm <- f_analysis_A_height(formula = mass_tot ~  s(yday) + P_2day + T_2day +
                                       C(traptype, "contr.sum") +
                                       C(bulbtype, "contr.sum") +
                                       n_trap +
                                       C(sample_previous, "contr.sum") +
                                       (1 | spattemp_cluster) +
                                       (1 | LOC) +
                                       (1 | night_ID) +
                                       (1 | trap_ID_A),
                                     data_z = d_mod_Tavg_z %>% filter(Tavg_mean_cat == "warm"),
                                     scalings = filter(d_scalings, data == "full"),
                                     family = "hurdle_gamma",
                                     hours_sel = which(d_mod_Tavg_z$traptype[d_mod_Tavg_z$Tavg_mean_cat == "warm"] == "p" &
                                                         d_mod_Tavg_z$hours_data[d_mod_Tavg_z$Tavg_mean_cat == "warm"]),
                                     iter = n_iter, seed = 255)

# ... per specialisation group #################################################
################################################################################.

# abundance --------------------------------------------------------------------.
l_abu_A_spec_m <- f_analysis_A_height(formula = abu_tot ~
                                        s(yday) + P_2day + T_2day +
                                        C(traptype, "contr.sum") +
                                        C(bulbtype, "contr.sum") +
                                        n_trap +
                                        C(sample_previous, "contr.sum") +
                                        (1 | spattemp_cluster) +
                                        (1 | LOC) +
                                        (1 | night_ID) +
                                        (1 | trap_ID_A),
                                      data_z = d_mod_spec_z %>% filter(Spec == "Monophagous"),
                                      scalings = filter(d_scalings, data == "full"),
                                      family = "zero_inflated_negbinomial",
                                      hours_sel = which(d_mod_spec_z$traptype[d_mod_spec_z$Spec == "Monophagous"] == "p" &
                                                          d_mod_spec_z$hours_data[d_mod_spec_z$Spec == "Monophagous"]),
                                      iter = n_iter, seed = 346)
l_abu_A_spec_o <- f_analysis_A_height(formula = abu_tot ~
                                        s(yday) + P_2day + T_2day +
                                        C(traptype, "contr.sum") +
                                        C(bulbtype, "contr.sum") +
                                        n_trap +
                                        C(sample_previous, "contr.sum") +
                                        (1 | spattemp_cluster) +
                                        (1 | LOC) +
                                        (1 | night_ID) +
                                        (1 | trap_ID_A),
                                      data_z = d_mod_spec_z %>% filter(Spec == "Oligophagous"),
                                      scalings = filter(d_scalings, data == "full"),
                                      family = "zero_inflated_negbinomial",
                                      hours_sel = which(d_mod_spec_z$traptype[d_mod_spec_z$Spec == "Oligophagous"] == "p" &
                                                          d_mod_spec_z$hours_data[d_mod_spec_z$Spec == "Oligophagous"]),
                                      iter = n_iter, seed = 665)
l_abu_A_spec_p <- f_analysis_A_height(formula = abu_tot ~
                                        s(yday) + P_2day + T_2day +
                                        C(traptype, "contr.sum") +
                                        C(bulbtype, "contr.sum") +
                                        n_trap +
                                        C(sample_previous, "contr.sum") +
                                        (1 | spattemp_cluster) +
                                        (1 | LOC) +
                                        (1 | night_ID) +
                                        (1 | trap_ID_A),
                                      data_z = d_mod_spec_z %>% filter(Spec == "Polyphagous"),
                                      scalings = filter(d_scalings, data == "full"),
                                      family = "zero_inflated_negbinomial",
                                      hours_sel = which(d_mod_spec_z$traptype[d_mod_spec_z$Spec == "Polyphagous"] == "p" &
                                                          d_mod_spec_z$hours_data[d_mod_spec_z$Spec == "Polyphagous"]),
                                      iter = n_iter, seed = 863)

# richness ---------------------------------------------------------------------.
l_ric_A_spec_m <- f_analysis_A_height(formula = SCcorr_ric ~
                                        s(yday) + P_2day + T_2day +
                                        C(traptype, "contr.sum") +
                                        C(bulbtype, "contr.sum") +
                                        n_trap +
                                        C(sample_previous, "contr.sum") +
                                        (1 | spattemp_cluster) +
                                        (1 | LOC) +
                                        (1 | night_ID) +
                                        (1 | trap_ID_A),
                                      data_z = d_mod_spec_z %>% filter(Spec == "Monophagous"),
                                      scalings = filter(d_scalings, data == "full"),
                                      family = "hurdle_gamma",
                                      hours_sel = which(d_mod_spec_z$traptype[d_mod_spec_z$Spec == "Monophagous"] == "p" &
                                                          d_mod_spec_z$hours_data[d_mod_spec_z$Spec == "Monophagous"]),
                                      iter = n_iter, seed = 8132)
l_ric_A_spec_o <- f_analysis_A_height(formula = SCcorr_ric ~
                                        s(yday) + P_2day + T_2day +
                                        C(traptype, "contr.sum") +
                                        C(bulbtype, "contr.sum") +
                                        n_trap +
                                        C(sample_previous, "contr.sum") +
                                        (1 | spattemp_cluster) +
                                        (1 | LOC) +
                                        (1 | night_ID) +
                                        (1 | trap_ID_A),
                                      data_z = d_mod_spec_z %>% filter(Spec == "Oligophagous"),
                                      scalings = filter(d_scalings, data == "full"),
                                      family = "hurdle_gamma",
                                      hours_sel = which(d_mod_spec_z$traptype[d_mod_spec_z$Spec == "Oligophagous"] == "p" &
                                                          d_mod_spec_z$hours_data[d_mod_spec_z$Spec == "Oligophagous"]),
                                      iter = n_iter, seed = 551)
l_ric_A_spec_p <- f_analysis_A_height(formula = SCcorr_ric ~
                                        s(yday) + P_2day + T_2day +
                                        C(traptype, "contr.sum") +
                                        C(bulbtype, "contr.sum") +
                                        n_trap +
                                        C(sample_previous, "contr.sum") +
                                        (1 | spattemp_cluster) +
                                        (1 | LOC) +
                                        (1 | night_ID) +
                                        (1 | trap_ID_A),
                                      data_z = d_mod_spec_z %>% filter(Spec == "Polyphagous"),
                                      scalings = filter(d_scalings, data == "full"),
                                      family = "hurdle_gamma",
                                      hours_sel = which(d_mod_spec_z$traptype[d_mod_spec_z$Spec == "Polyphagous"] == "p" &
                                                          d_mod_spec_z$hours_data[d_mod_spec_z$Spec == "Polyphagous"]),
                                      iter = n_iter, seed = 918)

# biomass ----------------------------------------------------------------------.
l_mass_A_spec_m <- f_analysis_A_height(formula = mass_tot ~
                                         s(yday) + P_2day + T_2day +
                                         C(traptype, "contr.sum") +
                                         C(bulbtype, "contr.sum") +
                                         n_trap +
                                         C(sample_previous, "contr.sum") +
                                         (1 | spattemp_cluster) +
                                         (1 | LOC) +
                                         (1 | night_ID) +
                                         (1 | trap_ID_A),
                                       data_z = d_mod_spec_z %>% filter(Spec == "Monophagous"),
                                       scalings = filter(d_scalings, data == "full"),
                                       family = "hurdle_gamma",
                                       hours_sel = which(d_mod_spec_z$traptype[d_mod_spec_z$Spec == "Monophagous"] == "p" &
                                                           d_mod_spec_z$hours_data[d_mod_spec_z$Spec == "Monophagous"]),
                                       iter = n_iter, seed = 810)
l_mass_A_spec_o <- f_analysis_A_height(formula = mass_tot ~
                                         s(yday) + P_2day + T_2day +
                                         C(traptype, "contr.sum") +
                                         C(bulbtype, "contr.sum") +
                                         n_trap +
                                         C(sample_previous, "contr.sum") +
                                         (1 | spattemp_cluster) +
                                         (1 | LOC) +
                                         (1 | night_ID) +
                                         (1 | trap_ID_A),
                                       data_z = d_mod_spec_z %>% filter(Spec == "Oligophagous"),
                                       scalings = filter(d_scalings, data == "full"),
                                       family = "hurdle_gamma",
                                       hours_sel = which(d_mod_spec_z$traptype[d_mod_spec_z$Spec == "Oligophagous"] == "p" &
                                                           d_mod_spec_z$hours_data[d_mod_spec_z$Spec == "Oligophagous"]),
                                       iter = n_iter, seed = 743)
l_mass_A_spec_p <- f_analysis_A_height(formula = mass_tot ~
                                         s(yday) + P_2day + T_2day +
                                         C(traptype, "contr.sum") +
                                         C(bulbtype, "contr.sum") +
                                         n_trap +
                                         C(sample_previous, "contr.sum") +
                                         (1 | spattemp_cluster) +
                                         (1 | LOC) +
                                         (1 | night_ID) +
                                         (1 | trap_ID_A),
                                       data_z = d_mod_spec_z %>% filter(Spec == "Polyphagous"),
                                       scalings = filter(d_scalings, data == "full"),
                                       family = "hurdle_gamma",
                                       hours_sel = which(d_mod_spec_z$traptype[d_mod_spec_z$Spec == "Polyphagous"] == "p" &
                                                           d_mod_spec_z$hours_data[d_mod_spec_z$Spec == "Polyphagous"]),
                                       iter = n_iter, seed = 822)

# ... per overwintering stage ##################################################
################################################################################.

# abundance --------------------------------------------------------------------.
l_abu_A_egg <- f_analysis_A_height(formula = abu_tot ~
                                     s(yday) + P_2day + T_2day +
                                     C(traptype, "contr.sum") +
                                     C(bulbtype, "contr.sum") +
                                     n_trap +
                                     C(sample_previous, "contr.sum") +
                                     (1 | spattemp_cluster) +
                                     (1 | LOC) +
                                     (1 | night_ID) +
                                     (1 | trap_ID_A),
                                   data_z = d_mod_hib_z %>% filter(overwintering_stage == "egg"),
                                   scalings = filter(d_scalings, data == "full"),
                                   family = "zero_inflated_negbinomial",
                                   hours_sel = which(d_mod_hib_z$traptype[d_mod_hib_z$overwintering_stage == "egg"] == "p" &
                                                       d_mod_hib_z$hours_data[d_mod_hib_z$overwintering_stage == "egg"]),
                                   iter = n_iter, seed = 121)
l_abu_A_larva <- f_analysis_A_height(formula = abu_tot ~
                                       s(yday) + P_2day + T_2day +
                                       C(traptype, "contr.sum") +
                                       C(bulbtype, "contr.sum") +
                                       n_trap +
                                       C(sample_previous, "contr.sum") +
                                       (1 | spattemp_cluster) +
                                       (1 | LOC) +
                                       (1 | night_ID) +
                                       (1 | trap_ID_A),
                                     data_z = d_mod_hib_z %>% filter(overwintering_stage == "larva"),
                                     scalings = filter(d_scalings, data == "full"),
                                     family = "zero_inflated_negbinomial",
                                     hours_sel = which(d_mod_hib_z$traptype[d_mod_hib_z$overwintering_stage == "larva"] == "p" &
                                                         d_mod_hib_z$hours_data[d_mod_hib_z$overwintering_stage == "larva"]),
                                     iter = n_iter, seed = 442)
l_abu_A_pupa <- f_analysis_A_height(formula = abu_tot ~
                                      s(yday) + P_2day + T_2day +
                                      C(traptype, "contr.sum") +
                                      C(bulbtype, "contr.sum") +
                                      n_trap +
                                      C(sample_previous, "contr.sum") +
                                      (1 | spattemp_cluster) +
                                      (1 | LOC) +
                                      (1 | night_ID) +
                                      (1 | trap_ID_A),
                                    data_z = d_mod_hib_z %>% filter(overwintering_stage == "pupa"),
                                    scalings = filter(d_scalings, data == "full"),
                                    family = "zero_inflated_negbinomial",
                                    hours_sel = which(d_mod_hib_z$traptype[d_mod_hib_z$overwintering_stage == "pupa"] == "p" &
                                                        d_mod_hib_z$hours_data[d_mod_hib_z$overwintering_stage == "pupa"]),
                                    iter = n_iter, seed = 9877)
l_abu_A_adult <- f_analysis_A_height(formula = abu_tot ~
                                       s(yday) + P_2day + T_2day +
                                       C(traptype, "contr.sum") +
                                       C(bulbtype, "contr.sum") +
                                       n_trap +
                                       C(sample_previous, "contr.sum") +
                                       (1 | spattemp_cluster) +
                                       (1 | LOC) +
                                       (1 | night_ID) +
                                       (1 | trap_ID_A),
                                     data_z = d_mod_hib_z %>% filter(overwintering_stage == "adult"),
                                     scalings = filter(d_scalings, data == "full"),
                                     family = "zero_inflated_negbinomial",
                                     hours_sel = which(d_mod_hib_z$traptype[d_mod_hib_z$overwintering_stage == "adult"] == "p" &
                                                         d_mod_hib_z$hours_data[d_mod_hib_z$overwintering_stage == "adult"]),
                                     iter = n_iter, seed = 2111)

# richness ---------------------------------------------------------------------.
l_ric_A_egg <- f_analysis_A_height(formula = SCcorr_ric ~
                                     s(yday) + P_2day + T_2day +
                                     C(traptype, "contr.sum") +
                                     C(bulbtype, "contr.sum") +
                                     n_trap +
                                     C(sample_previous, "contr.sum") +
                                     (1 | spattemp_cluster) +
                                     (1 | LOC) +
                                     (1 | night_ID) +
                                     (1 | trap_ID_A),
                                   data_z = d_mod_hib_z %>% filter(overwintering_stage == "egg"),
                                   scalings = filter(d_scalings, data == "full"),
                                   family = "hurdle_gamma",
                                   hours_sel = which(d_mod_hib_z$traptype[d_mod_hib_z$overwintering_stage == "egg"] == "p" &
                                                       d_mod_hib_z$hours_data[d_mod_hib_z$overwintering_stage == "egg"]),
                                   iter = n_iter, seed = 12)
l_ric_A_larva <- f_analysis_A_height(formula = SCcorr_ric ~
                                       s(yday) + P_2day + T_2day +
                                       C(traptype, "contr.sum") +
                                       C(bulbtype, "contr.sum") +
                                       n_trap +
                                       C(sample_previous, "contr.sum") +
                                       (1 | spattemp_cluster) +
                                       (1 | LOC) +
                                       (1 | night_ID) +
                                       (1 | trap_ID_A),
                                     data_z = d_mod_hib_z %>% filter(overwintering_stage == "larva"),
                                     scalings = filter(d_scalings, data == "full"),
                                     family = "hurdle_gamma",
                                     hours_sel = which(d_mod_hib_z$traptype[d_mod_hib_z$overwintering_stage == "larva"] == "p" &
                                                         d_mod_hib_z$hours_data[d_mod_hib_z$overwintering_stage == "larva"]),
                                     iter = n_iter, seed = 8797)
l_ric_A_pupa <- f_analysis_A_height(formula = SCcorr_ric ~
                                      s(yday) + P_2day + T_2day +
                                      C(traptype, "contr.sum") +
                                      C(bulbtype, "contr.sum") +
                                      n_trap +
                                      C(sample_previous, "contr.sum") +
                                      (1 | spattemp_cluster) +
                                      (1 | LOC) +
                                      (1 | night_ID) +
                                      (1 | trap_ID_A),
                                    data_z = d_mod_hib_z %>% filter(overwintering_stage == "pupa"),
                                    scalings = filter(d_scalings, data == "full"),
                                    family = "hurdle_gamma",
                                    hours_sel = which(d_mod_hib_z$traptype[d_mod_hib_z$overwintering_stage == "pupa"] == "p" &
                                                        d_mod_hib_z$hours_data[d_mod_hib_z$overwintering_stage == "pupa"]),
                                    iter = n_iter, seed = 194)
l_ric_A_adult <- f_analysis_A_height(formula = SCcorr_ric ~
                                       s(yday) + P_2day + T_2day +
                                       C(traptype, "contr.sum") +
                                       C(bulbtype, "contr.sum") +
                                       n_trap +
                                       C(sample_previous, "contr.sum") +
                                       (1 | spattemp_cluster) +
                                       (1 | LOC) +
                                       (1 | night_ID) +
                                       (1 | trap_ID_A),
                                     data_z = d_mod_hib_z %>% filter(overwintering_stage == "adult"),
                                     scalings = filter(d_scalings, data == "full"),
                                     family = "hurdle_gamma",
                                     hours_sel = which(d_mod_hib_z$traptype[d_mod_hib_z$overwintering_stage == "adult"] == "p" &
                                                         d_mod_hib_z$hours_data[d_mod_hib_z$overwintering_stage == "adult"]),
                                     iter = n_iter, seed = 654)

# biomass ----------------------------------------------------------------------.
l_mass_A_egg <- f_analysis_A_height(formula = mass_tot ~
                                      s(yday) + P_2day + T_2day +
                                      C(traptype, "contr.sum") +
                                      C(bulbtype, "contr.sum") +
                                      n_trap +
                                      C(sample_previous, "contr.sum") +
                                      (1 | spattemp_cluster) +
                                      (1 | LOC) +
                                      (1 | night_ID) +
                                      (1 | trap_ID_A),
                                    data_z = d_mod_hib_z %>% filter(overwintering_stage == "egg"),
                                    scalings = filter(d_scalings, data == "full"),
                                    family = "hurdle_gamma",
                                    hours_sel = which(d_mod_hib_z$traptype[d_mod_hib_z$overwintering_stage == "egg"] == "p" &
                                                        d_mod_hib_z$hours_data[d_mod_hib_z$overwintering_stage == "egg"]),
                                    iter = n_iter, seed = 537)
l_mass_A_larva <- f_analysis_A_height(formula = mass_tot ~
                                        s(yday) + P_2day + T_2day +
                                        C(traptype, "contr.sum") +
                                        C(bulbtype, "contr.sum") +
                                        n_trap +
                                        C(sample_previous, "contr.sum") +
                                        (1 | spattemp_cluster) +
                                        (1 | LOC) +
                                        (1 | night_ID) +
                                        (1 | trap_ID_A),
                                      data_z = d_mod_hib_z %>% filter(overwintering_stage == "larva"),
                                      scalings = filter(d_scalings, data == "full"),
                                      family = "hurdle_gamma",
                                      hours_sel = which(d_mod_hib_z$traptype[d_mod_hib_z$overwintering_stage == "larva"] == "p" &
                                                          d_mod_hib_z$hours_data[d_mod_hib_z$overwintering_stage == "larva"]),
                                      iter = n_iter, seed = 634)
l_mass_A_pupa <- f_analysis_A_height(formula = mass_tot ~
                                       s(yday) + P_2day + T_2day +
                                       C(traptype, "contr.sum") +
                                       C(bulbtype, "contr.sum") +
                                       n_trap +
                                       C(sample_previous, "contr.sum") +
                                       (1 | spattemp_cluster) +
                                       (1 | LOC) +
                                       (1 | night_ID) +
                                       (1 | trap_ID_A),
                                     data_z = d_mod_hib_z %>% filter(overwintering_stage == "pupa"),
                                     scalings = filter(d_scalings, data == "full"),
                                     family = "hurdle_gamma",
                                     hours_sel = which(d_mod_hib_z$traptype[d_mod_hib_z$overwintering_stage == "pupa"] == "p" &
                                                         d_mod_hib_z$hours_data[d_mod_hib_z$overwintering_stage == "pupa"]),
                                     iter = n_iter, seed = 19)
l_mass_A_adult <- f_analysis_A_height(formula = mass_tot ~
                                        s(yday) + P_2day + T_2day +
                                        C(traptype, "contr.sum") +
                                        C(bulbtype, "contr.sum") +
                                        n_trap +
                                        C(sample_previous, "contr.sum") +
                                        (1 | spattemp_cluster) +
                                        (1 | LOC) +
                                        (1 | night_ID) +
                                        (1 | trap_ID_A),
                                      data_z = d_mod_hib_z %>% filter(overwintering_stage == "adult"),
                                      scalings = filter(d_scalings, data == "full"),
                                      family = "hurdle_gamma",
                                      hours_sel = which(d_mod_hib_z$traptype[d_mod_hib_z$overwintering_stage == "adult"] == "p" &
                                                          d_mod_hib_z$hours_data[d_mod_hib_z$overwintering_stage == "adult"]),
                                      iter = n_iter, seed = 944)

# ... sensitivity analyses #####################################################
################################################################################.

# abundance --------------------------------------------------------------------.

# full observation hours data:
l_abu_A_ne <- f_analysis_A_height(formula = abu_tot ~
                                    s(yday) + P_2day + T_2day +
                                    C(traptype, "contr.sum") +
                                    C(bulbtype, "contr.sum") +
                                    n_trap +
                                    C(sample_previous, "contr.sum") +
                                    (1 | spattemp_cluster) +
                                    (1 | LOC) +
                                    (1 | night_ID) +
                                    (1 | trap_ID_A),
                                  data_z = d_mod_ne_z,
                                  scalings = filter(d_scalings, data == "no_estimates"),
                                  family = "zero_inflated_negbinomial",
                                  hours_sel = which(d_mod_ne_z$traptype == "p"),
                                  iter = n_iter, seed = 54987)

# fixed traps only:
l_abu_A_lf <- f_analysis_A_height(formula = abu_tot ~
                                    s(yday) + P_2day + T_2day +
                                    C(traptype, "contr.sum") +
                                    C(bulbtype, "contr.sum") +
                                    C(sample_previous, "contr.sum") +
                                    (1 | spattemp_cluster) +
                                    (1 | LOC) +
                                    (1 | night_ID) +
                                    (1 | trap_ID_A),
                                  data_z = d_mod_lf_z,
                                  scalings = filter(d_scalings, data == "LF"),
                                  family = "zero_inflated_negbinomial",
                                  iter = n_iter, seed = 1213)


# manual traps only:
l_abu_A_p <- f_analysis_A_height(formula = abu_tot ~
                                   s(yday) + P_2day + T_2day +
                                   C(bulbtype, "contr.sum") +
                                   n_trap +
                                   C(sample_previous, "contr.sum") +
                                   (1 | spattemp_cluster) +
                                   (1 | LOC) +
                                   (1 | night_ID) +
                                   (1 | trap_ID_A),
                                 data_z = d_mod_p_z,
                                 scalings = filter(d_scalings, data == "p"),
                                 family = "zero_inflated_negbinomial",
                                 hours_sel = which(d_mod_p_z$hours_data),
                                 iter = n_iter, seed = 6894)

# richness ---------------------------------------------------------------------.

# full observation hours data:
l_ric_A_ne <- f_analysis_A_height(formula = SCcorr_ric ~
                                    s(yday) + P_2day + T_2day +
                                    C(traptype, "contr.sum") +
                                    C(bulbtype, "contr.sum") +
                                    n_trap +
                                    C(sample_previous, "contr.sum") +
                                    (1 | spattemp_cluster) +
                                    (1 | LOC) +
                                    (1 | night_ID) +
                                    (1 | trap_ID_A),
                                  data_z = d_mod_ne_z,
                                  scalings = filter(d_scalings, data == "no_estimates"),
                                  family = "hurdle_gamma",
                                  hours_sel = which(d_mod_ne_z$traptype == "p"),
                                  iter = n_iter, seed = 1211)
# fixed traps only:
l_ric_A_lf <- f_analysis_A_height(formula = SCcorr_ric ~
                                    s(yday) + P_2day + T_2day +
                                    C(traptype, "contr.sum") +
                                    C(bulbtype, "contr.sum") +
                                    C(sample_previous, "contr.sum") +
                                    (1 | spattemp_cluster) + (1 | LOC) +
                                    (1 | night_ID) +
                                    (1 | trap_ID_A),
                                  data_z = d_mod_lf_z,
                                  scalings = filter(d_scalings, data == "LF"),
                                  family = "hurdle_gamma",
                                  iter = n_iter, seed = 45461)
# manual traps only:
l_ric_A_p <- f_analysis_A_height(formula = SCcorr_ric ~
                                   s(yday) + P_2day + T_2day +
                                   C(bulbtype, "contr.sum") +
                                   n_trap +
                                   C(sample_previous, "contr.sum") +
                                   (1 | spattemp_cluster) + (1 | LOC) +
                                   (1 | night_ID) +
                                   (1 | trap_ID_A),
                                 data_z = d_mod_p_z,
                                 scalings = filter(d_scalings, data == "p"),
                                 family = "hurdle_gamma",
                                 hours_sel = which(d_mod_p_z$hours_data),
                                 iter = n_iter, seed = 6894)

# biomass ----------------------------------------------------------------------.

# full observation hours data:
l_mass_A_ne <- f_analysis_A_height(formula = mass_tot ~
                                     s(yday) + P_2day + T_2day +
                                     C(traptype, "contr.sum") +
                                     C(bulbtype, "contr.sum") +
                                     n_trap +
                                     C(sample_previous, "contr.sum") +
                                     (1 | spattemp_cluster) +
                                     (1 | LOC) +
                                     (1 | night_ID) +
                                     (1 | trap_ID_A),
                                   data_z = d_mod_ne_z,
                                   scalings = filter(d_scalings, data == "no_estimates"),
                                   family = "hurdle_gamma",
                                   hours_sel = which(d_mod_ne_z$traptype == "p"),
                                   iter = n_iter, seed = 9101)

# fixed traps only:
l_mass_A_lf <- f_analysis_A_height(formula = mass_tot ~
                                     s(yday) + P_2day + T_2day +
                                     C(traptype, "contr.sum") +
                                     C(bulbtype, "contr.sum") +
                                     C(sample_previous, "contr.sum") +
                                     (1 | spattemp_cluster) +
                                     (1 | LOC) +
                                     (1 | night_ID) +
                                     (1 | trap_ID_A),
                                   data_z = d_mod_lf_z,
                                   scalings = filter(d_scalings, data == "LF"),
                                   family = "hurdle_gamma",
                                   iter = n_iter, seed = 565)

# manual traps only:
l_mass_A_p <- f_analysis_A_height(formula = mass_tot ~
                                    s(yday) + P_2day + T_2day +
                                    C(bulbtype, "contr.sum") +
                                    n_trap +
                                    C(sample_previous, "contr.sum") +
                                    (1 | spattemp_cluster) +
                                    (1 | LOC) +
                                    (1 | night_ID) +
                                    (1 | trap_ID_A),
                                  data_z = d_mod_p_z,
                                  scalings = filter(d_scalings, data == "p"),
                                  family = "hurdle_gamma",
                                  hours_sel = which(d_mod_p_z$hours_data),
                                  iter = n_iter, seed = 6894)

# ... single-site models #######################################################
################################################################################.

d_LOC_sel <- d_mod |> 
  group_by(LOC) |> 
  summarise(n_years = n_distinct(A),
            n_obs = n(),
            .groups = "drop") |> 
  filter(n_years >= 4,
         n_obs >= 10) 


# abundance --------------------------------------------------------------------.

d_slope_A_abu <- data.frame()
for (LOC_i in d_LOC_sel$LOC){
  d_mod_z_i <- d_mod_z |> filter(LOC == LOC_i) |> droplevels()
  
  formula_i <- abu_tot ~
    s(yday) + P_2day + T_2day +
    (1 | A_id)
  
  if (nlevels(d_mod_z_i$traptype) > 1)  {
    formula_i <- update(formula_i, . ~ . + C(traptype, "contr.sum"))
  }
  if (nlevels(d_mod_z_i$bulbtype) > 1)  {
    formula_i <- update(formula_i, . ~ . + C(bulbtype, "contr.sum"))
  }
  if (nlevels(d_mod_z_i$n_trap) > 1)  {
    formula_i <- update(formula_i, . ~ . + n_trap)
  }
  if (nlevels(d_mod_z_i$sample_previous) > 1)  {
    formula_i <- update(formula_i, . ~ . + C(sample_previous, "contr.sum"))
  }
  
  hours_sel_i <- which(d_mod_z_i$traptype == "p" & d_mod_z_i$hours_data)
  
  if (length(hours_sel_i) == 0) hours_sel_i <- NULL
  
  l_abu_A_i <- f_analysis_A_site(formula = formula_i,
                                 data_z = d_mod_z_i,
                                 scalings = filter(d_scalings, data == "full"),
                                 family = "zero_inflated_negbinomial",
                                 hours_sel = hours_sel_i,
                                 iter = n_iter, seed = 126)
  
  
  terms <- brmsterms(update(formula_i, ". ~ A + ."))
  mm_fe <- model.matrix(terms$dpars$mu$fe, data = d_mod_z_i)[, -1]
  
  c_slope_A_abu_i <- l_abu_A_i$fit_lm$b[, which(colnames(mm_fe) == "A"), drop = T]
  
  d_slope_A_abu <- data.frame(LOC = LOC_i, slope = c_slope_A_abu_i, iter = seq_len(length(c_slope_A_abu_i))) |> 
    bind_rows(d_slope_A_abu)
}

# richness ---------------------------------------------------------------------.

d_slope_A_SCcr <- data.frame()
for (LOC_i in d_LOC_sel$LOC){
  d_mod_z_i <- d_mod_z |> filter(LOC == LOC_i) |> droplevels()
  
  formula_i <- SCcorr_ric ~
    s(yday) + P_2day + T_2day +
    (1 | A_id)
  
  if (nlevels(d_mod_z_i$traptype) > 1)  {
    formula_i <- update(formula_i, . ~ . + C(traptype, "contr.sum"))
  }
  if (nlevels(d_mod_z_i$bulbtype) > 1)  {
    formula_i <- update(formula_i, . ~ . + C(bulbtype, "contr.sum"))
  }
  if (nlevels(d_mod_z_i$n_trap) > 1)  {
    formula_i <- update(formula_i, . ~ . + n_trap)
  }
  if (nlevels(d_mod_z_i$sample_previous) > 1)  {
    formula_i <- update(formula_i, . ~ . + C(sample_previous, "contr.sum"))
  }
  
  hours_sel_i <- which(d_mod_z_i$traptype == "p" & d_mod_z_i$hours_dat)
  
  if (length(hours_sel_i) == 0) hours_sel_i <- NULL
  
  l_SCcr_A_i <- f_analysis_A_site(formula = formula_i,
                                  data_z = d_mod_z_i,
                                  scalings = filter(d_scalings, data == "full"),
                                  family = "hurdle_gamma",
                                  hours_sel = hours_sel_i,
                                  iter = n_iter, seed = 126)
  
  
  terms <- brmsterms(update(formula_i, ". ~ A + ."))
  mm_fe <- model.matrix(terms$dpars$mu$fe, data = d_mod_z_i)[, -1]
  
  c_slope_A_SCcr_i <- l_SCcr_A_i$fit_lm$b[, which(colnames(mm_fe) == "A"), drop = T]
  
  d_slope_A_SCcr <- data.frame(LOC = LOC_i, slope = c_slope_A_SCcr_i, iter = seq_len(length(c_slope_A_SCcr_i))) |> 
    bind_rows(d_slope_A_SCcr)
}


# biomass ----------------------------------------------------------------------.

d_slope_A_mass <- data.frame()
for (LOC_i in d_LOC_sel$LOC){
  d_mod_z_i <- d_mod_z |> filter(LOC == LOC_i) |> droplevels()
  
  formula_i <- mass_tot ~
    s(yday) + P_2day + T_2day +
    (1 | A_id)
  
  if (nlevels(d_mod_z_i$traptype) > 1)  {
    formula_i <- update(formula_i, . ~ . + C(traptype, "contr.sum"))
  }
  if (nlevels(d_mod_z_i$bulbtype) > 1)  {
    formula_i <- update(formula_i, . ~ . + C(bulbtype, "contr.sum"))
  }
  if (nlevels(d_mod_z_i$n_trap) > 1)  {
    formula_i <- update(formula_i, . ~ . + n_trap)
  }
  if (nlevels(d_mod_z_i$sample_previous) > 1)  {
    formula_i <- update(formula_i, . ~ . + C(sample_previous, "contr.sum"))
  }
  
  hours_sel_i <- which(d_mod_z_i$traptype == "p" & d_mod_z_i$hours_dat)
  
  if (length(hours_sel_i) == 0) hours_sel_i <- NULL
  
  l_mass_A_i <- f_analysis_A_site(formula = formula_i,
                                  data_z = d_mod_z_i,
                                  scalings = filter(d_scalings, data == "full"),
                                  family = "hurdle_gamma",
                                  hours_sel = hours_sel_i,
                                  iter = n_iter, seed = 126)
  
  
  terms <- brmsterms(update(formula_i, ". ~ A + ."))
  mm_fe <- model.matrix(terms$dpars$mu$fe, data = d_mod_z_i)[, -1]
  
  c_slope_A_mass_i <- l_mass_A_i$fit_lm$b[, which(colnames(mm_fe) == "A"), drop = T]
  
  d_slope_A_mass <- data.frame(LOC = LOC_i, slope = c_slope_A_mass_i, iter = seq_len(length(c_slope_A_mass_i))) |> 
    bind_rows(d_slope_A_mass)
}

# ... single-species models ####################################################
################################################################################.

d_sel <- d_moths |> 
  group_by(Name_std) |> 
  summarise(n_comb = n_distinct(paste(LOC, A)),
            .groups = "drop") |> 
  filter(n_comb >= 100)

for (sp_i in d_sel$Name_std){
  
  d_mod_i <- d_moths %>%
    filter(Name_std == sp_i) |> 
    group_by(LOC, Samplingdate, A_id = A) %>%
    summarise(abu = sum(individualCount),
              .groups = "drop") |> 
    full_join(d_mod_z |> select(-abu_tot),
              by = join_by(LOC, Samplingdate, A_id)) |>
    mutate(abu = ifelse(is.na(abu), 0, abu))
  
  
  out <- f_analysis_A_height_splev(formula = abu ~
                                     s(yday) + P_2day + T_2day +
                                     C(traptype, "contr.sum") +
                                     C(bulbtype, "contr.sum") +
                                     n_trap +
                                     C(sample_previous, "contr.sum") +
                                     (1 | spattemp_cluster) +
                                     (1 | LOC) +
                                     (1 | night_ID) +
                                     (1 | trap_ID_A),
                                   data_z = d_mod_i,
                                   scalings = filter(d_scalings, data == "full"),
                                   family = "zero_inflated_negbinomial",
                                   hours_sel = which(d_mod_i$traptype == "p" & d_mod_i$hours_data),
                                   iter = n_iter,
                                   seed = which(d_sel$Name_std == sp_i))
  
  out$d_coefs <- out$d_coefs |> 
    mutate(Name_std = sp_i)
  out$d_coefs_mean <- out$d_coefs_mean |> 
    mutate(Name_std = sp_i)
  
  # save output temporarily
  saveRDS(out, "Output/singlespecies/splev_", 
          gsub(" |/", "_", sp_i), ".rds")
}


# create combined dataframes ---------------------------------------------------.

files <- list.files("Output/singlespecies", full.names = T)

d_coefs_mean <- lapply(files, function(x) readRDS(x)$d_coefs_mean) |>
  bind_rows()

set.seed(143)
d_coefs_sim <- lapply(files, function(x) readRDS(x)$d_coefs[sample(1:4000, 5000, replace = T), c("A:height", "Name_std")]) |> 
  bind_rows()


# CREATE OUTPUTS -------------- ################################################
################################################################################.

# ... Text / numbers ###########################################################
################################################################################.

f_A_height_change_numbers(l_abu_A$fit,
                          d_mod_z,
                          brmsterms(update(abu_tot ~ s(yday) + P_2day + T_2day + 
                                             C(traptype, "contr.sum") + 
                                             C(bulbtype, "contr.sum") + 
                                             n_trap + C(sample_previous, "contr.sum") + 
                                             (1 | spattemp_cluster) + (1 | LOC) + 
                                             (1 | night_ID) + (1 | trap_ID_A), 
                                           ". ~ A * height + ."))) |> 
  mutate(perc_estimate = (factor_estimate - 1) * 100,
         perc_lower = (factor_lower - 1) * 100,
         perc_upper = (factor_upper - 1) * 100) |> 
  select(-c(diff_log_estimate, diff_log_lower, diff_log_upper)) |> 
  mutate(text_factor = paste0(formatC(factor_estimate, digits = 3, format = "fg", flag = "#"), 
                              " (95%-CI: ", 
                              formatC(factor_lower, digits = 3, format = "fg", flag = "#"), "–", 
                              formatC(factor_upper, digits = 3, format = "fg", flag = "#"), ")"),
         text_perc = paste0(as.character(formatC(perc_estimate, digits = 3, format = "fg", flag = "#")), 
                            "% (95%-CI: ", 
                            formatC(perc_lower, digits = 3, format = "fg", flag = "#"), "% to ", 
                            formatC(perc_upper, digits = 3, format = "fg", flag = "#"), "%)"))

f_A_height_change_numbers(l_ric_A$fit,
                          d_mod_z,
                          brmsterms(update(SCcorr_ric ~ s(yday) + P_2day + T_2day + 
                                             C(traptype, "contr.sum") + 
                                             C(bulbtype, "contr.sum") + 
                                             n_trap + C(sample_previous, "contr.sum") + 
                                             (1 | spattemp_cluster) + (1 | LOC) + 
                                             (1 | night_ID) + (1 | trap_ID_A), 
                                           ". ~ A * height + ."))) |> 
  mutate(perc_estimate = (factor_estimate - 1) * 100,
         perc_lower = (factor_lower - 1) * 100,
         perc_upper = (factor_upper - 1) * 100) |> 
  select(-c(diff_log_estimate, diff_log_lower, diff_log_upper)) |> 
  mutate(text_factor = paste0(formatC(factor_estimate, digits = 3, format = "fg", flag = "#"), 
                              " (95%-CI: ", 
                              formatC(factor_lower, digits = 3, format = "fg", flag = "#"), "–", 
                              formatC(factor_upper, digits = 3, format = "fg", flag = "#"), ")"),
         text_perc = paste0(as.character(formatC(perc_estimate, digits = 3, format = "fg", flag = "#")), 
                            "% (95%-CI: ", 
                            formatC(perc_lower, digits = 3, format = "fg", flag = "#"), "% to ", 
                            formatC(perc_upper, digits = 3, format = "fg", flag = "#"), "%)"))

f_A_height_change_numbers(l_mass_A$fit,
                          d_mod_z,
                          brmsterms(update(mass_tot ~ s(yday) + P_2day + T_2day + 
                                             C(traptype, "contr.sum") + 
                                             C(bulbtype, "contr.sum") + 
                                             n_trap + C(sample_previous, "contr.sum") + 
                                             (1 | spattemp_cluster) + (1 | LOC) + 
                                             (1 | night_ID) + (1 | trap_ID_A), 
                                           ". ~ A * height + ."))) |> 
  mutate(perc_estimate = (factor_estimate - 1) * 100,
         perc_lower = (factor_lower - 1) * 100,
         perc_upper = (factor_upper - 1) * 100) |> 
  select(-c(diff_log_estimate, diff_log_lower, diff_log_upper)) |> 
  mutate(text_factor = paste0(formatC(factor_estimate, digits = 3, format = "fg", flag = "#"), 
                              " (95%-CI: ", 
                              formatC(factor_lower, digits = 3, format = "fg", flag = "#"), "–", 
                              formatC(factor_upper, digits = 3, format = "fg", flag = "#"), ")"),
         text_perc = paste0(as.character(formatC(perc_estimate, digits = 3, format = "fg", flag = "#")), 
                            "% (95%-CI: ", 
                            formatC(perc_lower, digits = 3, format = "fg", flag = "#"), "% to ", 
                            formatC(perc_upper, digits = 3, format = "fg", flag = "#"), "%)"))

f_A_height_change_numbers(l_abu_A_small$fit,
                          d_mod_z,
                          brmsterms(update(abu_tot ~ s(yday) + P_2day + T_2day + 
                                             C(traptype, "contr.sum") + 
                                             C(bulbtype, "contr.sum") + 
                                             n_trap + C(sample_previous, "contr.sum") + 
                                             (1 | spattemp_cluster) + (1 | LOC) + 
                                             (1 | night_ID) + (1 | trap_ID_A), 
                                           ". ~ A * height + ."))) |> 
  mutate(perc_estimate = (factor_estimate - 1) * 100,
         perc_lower = (factor_lower - 1) * 100,
         perc_upper = (factor_upper - 1) * 100) |> 
  select(-c(diff_log_estimate, diff_log_lower, diff_log_upper)) |> 
  mutate(text_factor = paste0(formatC(factor_estimate, digits = 3, format = "fg", flag = "#"), 
                              " (95%-CI: ", 
                              formatC(factor_lower, digits = 3, format = "fg", flag = "#"), "–", 
                              formatC(factor_upper, digits = 3, format = "fg", flag = "#"), ")"),
         text_perc = paste0(as.character(formatC(perc_estimate, digits = 3, format = "fg", flag = "#")), 
                            "% (95%-CI: ", 
                            formatC(perc_lower, digits = 3, format = "fg", flag = "#"), "% to ", 
                            formatC(perc_upper, digits = 3, format = "fg", flag = "#"), "%)"))
f_A_height_change_numbers(l_mass_A_small$fit,
                          d_mod_z,
                          brmsterms(update(mass_tot ~ s(yday) + P_2day + T_2day + 
                                             C(traptype, "contr.sum") + 
                                             C(bulbtype, "contr.sum") + 
                                             n_trap + C(sample_previous, "contr.sum") + 
                                             (1 | spattemp_cluster) + (1 | LOC) + 
                                             (1 | night_ID) + (1 | trap_ID_A), 
                                           ". ~ A * height + ."))) |> 
  mutate(perc_estimate = (factor_estimate - 1) * 100,
         perc_lower = (factor_lower - 1) * 100,
         perc_upper = (factor_upper - 1) * 100) |> 
  select(-c(diff_log_estimate, diff_log_lower, diff_log_upper)) |> 
  mutate(text_factor = paste0(formatC(factor_estimate, digits = 3, format = "fg", flag = "#"), 
                              " (95%-CI: ", 
                              formatC(factor_lower, digits = 3, format = "fg", flag = "#"), "–", 
                              formatC(factor_upper, digits = 3, format = "fg", flag = "#"), ")"),
         text_perc = paste0(as.character(formatC(perc_estimate, digits = 3, format = "fg", flag = "#")), 
                            "% (95%-CI: ", 
                            formatC(perc_lower, digits = 3, format = "fg", flag = "#"), "% to ", 
                            formatC(perc_upper, digits = 3, format = "fg", flag = "#"), "%)"))
f_A_height_change_numbers(l_abu_A_large$fit,
                          d_mod_z,
                          brmsterms(update(abu_tot ~ s(yday) + P_2day + T_2day + 
                                             C(traptype, "contr.sum") + 
                                             C(bulbtype, "contr.sum") + 
                                             n_trap + C(sample_previous, "contr.sum") + 
                                             (1 | spattemp_cluster) + (1 | LOC) + 
                                             (1 | night_ID) + (1 | trap_ID_A), 
                                           ". ~ A * height + ."))) |> 
  mutate(perc_estimate = (factor_estimate - 1) * 100,
         perc_lower = (factor_lower - 1) * 100,
         perc_upper = (factor_upper - 1) * 100) |> 
  select(-c(diff_log_estimate, diff_log_lower, diff_log_upper)) |> 
  mutate(text_factor = paste0(formatC(factor_estimate, digits = 3, format = "fg", flag = "#"), 
                              " (95%-CI: ", 
                              formatC(factor_lower, digits = 3, format = "fg", flag = "#"), "–", 
                              formatC(factor_upper, digits = 3, format = "fg", flag = "#"), ")"),
         text_perc = paste0(as.character(formatC(perc_estimate, digits = 3, format = "fg", flag = "#")), 
                            "% (95%-CI: ", 
                            formatC(perc_lower, digits = 3, format = "fg", flag = "#"), "% to ", 
                            formatC(perc_upper, digits = 3, format = "fg", flag = "#"), "%)"))

f_A_height_change_numbers(l_abu_A_cold$fit,
                          d_mod_z,
                          brmsterms(update(abu_tot ~ s(yday) + P_2day + T_2day + 
                                             C(traptype, "contr.sum") + 
                                             C(bulbtype, "contr.sum") + 
                                             n_trap + C(sample_previous, "contr.sum") + 
                                             (1 | spattemp_cluster) + (1 | LOC) + 
                                             (1 | night_ID) + (1 | trap_ID_A), 
                                           ". ~ A * height + ."))) |> 
  mutate(perc_estimate = (factor_estimate - 1) * 100,
         perc_lower = (factor_lower - 1) * 100,
         perc_upper = (factor_upper - 1) * 100) |> 
  select(-c(diff_log_estimate, diff_log_lower, diff_log_upper)) |> 
  mutate(text_factor = paste0(formatC(factor_estimate, digits = 3, format = "fg", flag = "#"), 
                              " (95%-CI: ", 
                              formatC(factor_lower, digits = 3, format = "fg", flag = "#"), "–", 
                              formatC(factor_upper, digits = 3, format = "fg", flag = "#"), ")"),
         text_perc = paste0(as.character(formatC(perc_estimate, digits = 3, format = "fg", flag = "#")), 
                            "% (95%-CI: ", 
                            formatC(perc_lower, digits = 3, format = "fg", flag = "#"), "% to ", 
                            formatC(perc_upper, digits = 3, format = "fg", flag = "#"), "%)"))
f_A_height_change_numbers(l_ric_A_cold$fit,
                          d_mod_z,
                          brmsterms(update(SCcorr_ric ~ s(yday) + P_2day + T_2day + 
                                             C(traptype, "contr.sum") + 
                                             C(bulbtype, "contr.sum") + 
                                             n_trap + C(sample_previous, "contr.sum") + 
                                             (1 | spattemp_cluster) + (1 | LOC) + 
                                             (1 | night_ID) + (1 | trap_ID_A), 
                                           ". ~ A * height + ."))) |> 
  mutate(perc_estimate = (factor_estimate - 1) * 100,
         perc_lower = (factor_lower - 1) * 100,
         perc_upper = (factor_upper - 1) * 100) |> 
  select(-c(diff_log_estimate, diff_log_lower, diff_log_upper)) |> 
  mutate(text_factor = paste0(formatC(factor_estimate, digits = 3, format = "fg", flag = "#"), 
                              " (95%-CI: ", 
                              formatC(factor_lower, digits = 3, format = "fg", flag = "#"), "–", 
                              formatC(factor_upper, digits = 3, format = "fg", flag = "#"), ")"),
         text_perc = paste0(as.character(formatC(perc_estimate, digits = 3, format = "fg", flag = "#")), 
                            "% (95%-CI: ", 
                            formatC(perc_lower, digits = 3, format = "fg", flag = "#"), "% to ", 
                            formatC(perc_upper, digits = 3, format = "fg", flag = "#"), "%)"))
f_A_height_change_numbers(l_mass_A_cold$fit,
                          d_mod_z,
                          brmsterms(update(mass_tot ~ s(yday) + P_2day + T_2day + 
                                             C(traptype, "contr.sum") + 
                                             C(bulbtype, "contr.sum") + 
                                             n_trap + C(sample_previous, "contr.sum") + 
                                             (1 | spattemp_cluster) + (1 | LOC) + 
                                             (1 | night_ID) + (1 | trap_ID_A), 
                                           ". ~ A * height + ."))) |> 
  mutate(perc_estimate = (factor_estimate - 1) * 100,
         perc_lower = (factor_lower - 1) * 100,
         perc_upper = (factor_upper - 1) * 100) |> 
  select(-c(diff_log_estimate, diff_log_lower, diff_log_upper)) |> 
  mutate(text_factor = paste0(formatC(factor_estimate, digits = 3, format = "fg", flag = "#"), 
                              " (95%-CI: ", 
                              formatC(factor_lower, digits = 3, format = "fg", flag = "#"), "–", 
                              formatC(factor_upper, digits = 3, format = "fg", flag = "#"), ")"),
         text_perc = paste0(as.character(formatC(perc_estimate, digits = 3, format = "fg", flag = "#")), 
                            "% (95%-CI: ", 
                            formatC(perc_lower, digits = 3, format = "fg", flag = "#"), "% to ", 
                            formatC(perc_upper, digits = 3, format = "fg", flag = "#"), "%)"))
f_A_height_change_numbers(l_abu_A_spec_m$fit,
                          d_mod_z,
                          brmsterms(update(abu_tot ~ s(yday) + P_2day + T_2day + 
                                             C(traptype, "contr.sum") + 
                                             C(bulbtype, "contr.sum") + 
                                             n_trap + C(sample_previous, "contr.sum") + 
                                             (1 | spattemp_cluster) + (1 | LOC) + 
                                             (1 | night_ID) + (1 | trap_ID_A), 
                                           ". ~ A * height + ."))) |> 
  mutate(perc_estimate = (factor_estimate - 1) * 100,
         perc_lower = (factor_lower - 1) * 100,
         perc_upper = (factor_upper - 1) * 100) |> 
  select(-c(diff_log_estimate, diff_log_lower, diff_log_upper)) |> 
  mutate(text_factor = paste0(formatC(factor_estimate, digits = 3, format = "fg", flag = "#"), 
                              " (95%-CI: ", 
                              formatC(factor_lower, digits = 3, format = "fg", flag = "#"), "–", 
                              formatC(factor_upper, digits = 3, format = "fg", flag = "#"), ")"),
         text_perc = paste0(as.character(formatC(perc_estimate, digits = 3, format = "fg", flag = "#")), 
                            "% (95%-CI: ", 
                            formatC(perc_lower, digits = 3, format = "fg", flag = "#"), "% to ", 
                            formatC(perc_upper, digits = 3, format = "fg", flag = "#"), "%)"))
f_A_height_change_numbers(l_ric_A_spec_m$fit,
                          d_mod_z,
                          brmsterms(update(SCcorr_ric ~ s(yday) + P_2day + T_2day + 
                                             C(traptype, "contr.sum") + 
                                             C(bulbtype, "contr.sum") + 
                                             n_trap + C(sample_previous, "contr.sum") + 
                                             (1 | spattemp_cluster) + (1 | LOC) + 
                                             (1 | night_ID) + (1 | trap_ID_A), 
                                           ". ~ A * height + ."))) |> 
  mutate(perc_estimate = (factor_estimate - 1) * 100,
         perc_lower = (factor_lower - 1) * 100,
         perc_upper = (factor_upper - 1) * 100) |> 
  select(-c(diff_log_estimate, diff_log_lower, diff_log_upper)) |> 
  mutate(text_factor = paste0(formatC(factor_estimate, digits = 3, format = "fg", flag = "#"), 
                              " (95%-CI: ", 
                              formatC(factor_lower, digits = 3, format = "fg", flag = "#"), "–", 
                              formatC(factor_upper, digits = 3, format = "fg", flag = "#"), ")"),
         text_perc = paste0(as.character(formatC(perc_estimate, digits = 3, format = "fg", flag = "#")), 
                            "% (95%-CI: ", 
                            formatC(perc_lower, digits = 3, format = "fg", flag = "#"), "% to ", 
                            formatC(perc_upper, digits = 3, format = "fg", flag = "#"), "%)"))
f_A_height_change_numbers(l_mass_A_spec_m$fit,
                          d_mod_z,
                          brmsterms(update(mass_tot ~ s(yday) + P_2day + T_2day + 
                                             C(traptype, "contr.sum") + 
                                             C(bulbtype, "contr.sum") + 
                                             n_trap + C(sample_previous, "contr.sum") + 
                                             (1 | spattemp_cluster) + (1 | LOC) + 
                                             (1 | night_ID) + (1 | trap_ID_A), 
                                           ". ~ A * height + ."))) |> 
  mutate(perc_estimate = (factor_estimate - 1) * 100,
         perc_lower = (factor_lower - 1) * 100,
         perc_upper = (factor_upper - 1) * 100) |> 
  select(-c(diff_log_estimate, diff_log_lower, diff_log_upper)) |> 
  mutate(text_factor = paste0(formatC(factor_estimate, digits = 3, format = "fg", flag = "#"), 
                              " (95%-CI: ", 
                              formatC(factor_lower, digits = 3, format = "fg", flag = "#"), "–", 
                              formatC(factor_upper, digits = 3, format = "fg", flag = "#"), ")"),
         text_perc = paste0(as.character(formatC(perc_estimate, digits = 3, format = "fg", flag = "#")), 
                            "% (95%-CI: ", 
                            formatC(perc_lower, digits = 3, format = "fg", flag = "#"), "% to ", 
                            formatC(perc_upper, digits = 3, format = "fg", flag = "#"), "%)"))
f_A_height_change_numbers(l_abu_A_spec_o$fit,
                          d_mod_z,
                          brmsterms(update(abu_tot ~ s(yday) + P_2day + T_2day + 
                                             C(traptype, "contr.sum") + 
                                             C(bulbtype, "contr.sum") + 
                                             n_trap + C(sample_previous, "contr.sum") + 
                                             (1 | spattemp_cluster) + (1 | LOC) + 
                                             (1 | night_ID) + (1 | trap_ID_A), 
                                           ". ~ A * height + ."))) |> 
  mutate(perc_estimate = (factor_estimate - 1) * 100,
         perc_lower = (factor_lower - 1) * 100,
         perc_upper = (factor_upper - 1) * 100) |> 
  select(-c(diff_log_estimate, diff_log_lower, diff_log_upper)) |> 
  mutate(text_factor = paste0(formatC(factor_estimate, digits = 3, format = "fg", flag = "#"), 
                              " (95%-CI: ", 
                              formatC(factor_lower, digits = 3, format = "fg", flag = "#"), "–", 
                              formatC(factor_upper, digits = 3, format = "fg", flag = "#"), ")"),
         text_perc = paste0(as.character(formatC(perc_estimate, digits = 3, format = "fg", flag = "#")), 
                            "% (95%-CI: ", 
                            formatC(perc_lower, digits = 3, format = "fg", flag = "#"), "% to ", 
                            formatC(perc_upper, digits = 3, format = "fg", flag = "#"), "%)"))
f_A_height_change_numbers(l_abu_A_pupa$fit,
                          d_mod_z,
                          brmsterms(update(abu_tot ~ s(yday) + P_2day + T_2day + 
                                             C(traptype, "contr.sum") + 
                                             C(bulbtype, "contr.sum") + 
                                             n_trap + C(sample_previous, "contr.sum") + 
                                             (1 | spattemp_cluster) + (1 | LOC) + 
                                             (1 | night_ID) + (1 | trap_ID_A), 
                                           ". ~ A * height + ."))) |> 
  mutate(perc_estimate = (factor_estimate - 1) * 100,
         perc_lower = (factor_lower - 1) * 100,
         perc_upper = (factor_upper - 1) * 100) |> 
  select(-c(diff_log_estimate, diff_log_lower, diff_log_upper)) |> 
  mutate(text_factor = paste0(formatC(factor_estimate, digits = 3, format = "fg", flag = "#"), 
                              " (95%-CI: ", 
                              formatC(factor_lower, digits = 3, format = "fg", flag = "#"), "–", 
                              formatC(factor_upper, digits = 3, format = "fg", flag = "#"), ")"),
         text_perc = paste0(as.character(formatC(perc_estimate, digits = 3, format = "fg", flag = "#")), 
                            "% (95%-CI: ", 
                            formatC(perc_lower, digits = 3, format = "fg", flag = "#"), "% to ", 
                            formatC(perc_upper, digits = 3, format = "fg", flag = "#"), "%)"))
f_A_height_change_numbers(l_ric_A_pupa$fit,
                          d_mod_z,
                          brmsterms(update(SCcorr_ric ~ s(yday) + P_2day + T_2day + 
                                             C(traptype, "contr.sum") + 
                                             C(bulbtype, "contr.sum") + 
                                             n_trap + C(sample_previous, "contr.sum") + 
                                             (1 | spattemp_cluster) + (1 | LOC) + 
                                             (1 | night_ID) + (1 | trap_ID_A), 
                                           ". ~ A * height + ."))) |> 
  mutate(perc_estimate = (factor_estimate - 1) * 100,
         perc_lower = (factor_lower - 1) * 100,
         perc_upper = (factor_upper - 1) * 100) |> 
  select(-c(diff_log_estimate, diff_log_lower, diff_log_upper)) |> 
  mutate(text_factor = paste0(formatC(factor_estimate, digits = 3, format = "fg", flag = "#"), 
                              " (95%-CI: ", 
                              formatC(factor_lower, digits = 3, format = "fg", flag = "#"), "–", 
                              formatC(factor_upper, digits = 3, format = "fg", flag = "#"), ")"),
         text_perc = paste0(as.character(formatC(perc_estimate, digits = 3, format = "fg", flag = "#")), 
                            "% (95%-CI: ", 
                            formatC(perc_lower, digits = 3, format = "fg", flag = "#"), "% to ", 
                            formatC(perc_upper, digits = 3, format = "fg", flag = "#"), "%)"))
f_A_height_change_numbers(l_mass_A_pupa$fit,
                          d_mod_z,
                          brmsterms(update(mass_tot ~ s(yday) + P_2day + T_2day + 
                                             C(traptype, "contr.sum") + 
                                             C(bulbtype, "contr.sum") + 
                                             n_trap + C(sample_previous, "contr.sum") + 
                                             (1 | spattemp_cluster) + (1 | LOC) + 
                                             (1 | night_ID) + (1 | trap_ID_A), 
                                           ". ~ A * height + ."))) |> 
  mutate(perc_estimate = (factor_estimate - 1) * 100,
         perc_lower = (factor_lower - 1) * 100,
         perc_upper = (factor_upper - 1) * 100) |> 
  select(-c(diff_log_estimate, diff_log_lower, diff_log_upper)) |> 
  mutate(text_factor = paste0(formatC(factor_estimate, digits = 3, format = "fg", flag = "#"), 
                              " (95%-CI: ", 
                              formatC(factor_lower, digits = 3, format = "fg", flag = "#"), "–", 
                              formatC(factor_upper, digits = 3, format = "fg", flag = "#"), ")"),
         text_perc = paste0(as.character(formatC(perc_estimate, digits = 3, format = "fg", flag = "#")), 
                            "% (95%-CI: ", 
                            formatC(perc_lower, digits = 3, format = "fg", flag = "#"), "% to ", 
                            formatC(perc_upper, digits = 3, format = "fg", flag = "#"), "%)"))

# single-species models --------------------------------------------------------.

d_coefs_sim |> 
  group_by(Name_std) |>
  mutate(iter = row_number()) |>
  ungroup() |>
  group_by(iter) |> 
  summarise(perc_pos = mean(`A:height` > 0) * 100,
            .groups = "drop") |> 
  summarise(perc_pos_m = mean(perc_pos),
            perc_pos_lower = ci(perc_pos)$CI_low,
            perc_pos_upper = ci(perc_pos)$CI_high) |> 
  mutate(text = paste0(round(perc_pos_m, 1), " (95%-CI: ", 
                       round(perc_pos_lower, 1), "-",
                       round(perc_pos_upper, 1), ")"))

d_coefs_sim |> 
  group_by(Name_std) |>
  mutate(iter = row_number()) |>
  ungroup() |>
  left_join(d_traits,
            by = "Name_std") |>
  group_by(mass_cat, iter) |> 
  summarise(perc_pos = mean(`A:height` > 0) * 100,
            .groups = "drop") |> 
  group_by(mass_cat) |> 
  summarise(perc_pos_m = mean(perc_pos),
            perc_pos_lower = ci(perc_pos)$CI_low,
            perc_pos_upper = ci(perc_pos)$CI_high,
            .groups = "drop") |> 
  mutate(text = paste0(round(perc_pos_m, 1), " (95%-CI: ", 
                       round(perc_pos_lower, 1), "-",
                       round(perc_pos_upper, 1), ")"))

d_coefs_sim |> 
  group_by(Name_std) |>
  mutate(iter = row_number()) |>
  ungroup() |>
  left_join(d_traits,
            by = "Name_std") |>
  group_by(Tavg_mean_cat, iter) |> 
  summarise(perc_pos = mean(`A:height` > 0) * 100,
            .groups = "drop") |> 
  group_by(Tavg_mean_cat) |> 
  summarise(perc_pos_m = mean(perc_pos),
            perc_pos_lower = ci(perc_pos)$CI_low,
            perc_pos_upper = ci(perc_pos)$CI_high,
            .groups = "drop") |> 
  mutate(text = paste0(round(perc_pos_m, 1), " (95%-CI: ", 
                       round(perc_pos_lower, 1), "-",
                       round(perc_pos_upper, 1), ")"))

d_coefs_sim |> 
  group_by(Name_std) |>
  mutate(iter = row_number()) |>
  ungroup() |>
  left_join(d_traits,
            by = "Name_std") |>
  filter(!is.na(Spec)) |> 
  group_by(Spec, iter) |> 
  summarise(perc_pos = mean(`A:height` > 0) * 100,
            .groups = "drop") |> 
  group_by(Spec) |> 
  summarise(perc_pos_m = mean(perc_pos),
            perc_pos_lower = ci(perc_pos)$CI_low,
            perc_pos_upper = ci(perc_pos)$CI_high,
            .groups = "drop") |> 
  mutate(text = paste0(round(perc_pos_m, 1), " (95%-CI: ", 
                       round(perc_pos_lower, 1), "-",
                       round(perc_pos_upper, 1), ")"))

d_coefs_sim |> 
  group_by(Name_std) |>
  mutate(iter = row_number()) |>
  ungroup() |>
  left_join(d_traits,
            by = "Name_std") |>
  filter(overwintering_stage != "larva/pupa") |> 
  group_by(overwintering_stage, iter) |> 
  summarise(perc_pos = mean(`A:height` > 0) * 100,
            .groups = "drop") |> 
  group_by(overwintering_stage) |> 
  summarise(perc_pos_m = mean(perc_pos),
            perc_pos_lower = ci(perc_pos)$CI_low,
            perc_pos_upper = ci(perc_pos)$CI_high,
            .groups = "drop") |> 
  mutate(text = paste0(round(perc_pos_m, 1), " (95%-CI: ", 
                       round(perc_pos_lower, 1), "-",
                       round(perc_pos_upper, 1), ")"))

# ... Rhat #####################################################################
################################################################################.

d_Rhat_abu <- f_apply_Rhat(l_abu_A$fit) |> 
  mutate(data = "data", trait = "Full") |> 
  bind_rows(f_apply_Rhat(l_abu_A_small$fit) |> 
              mutate(data = "small", trait = "Body size")) |> 
  bind_rows(f_apply_Rhat(l_abu_A_medium$fit) |> 
              mutate(data = "medium", trait = "Body size")) |> 
  bind_rows(f_apply_Rhat(l_abu_A_large$fit) |> 
              mutate(data = "large", trait = "Body size")) |> 
  bind_rows(f_apply_Rhat(l_abu_A_cold$fit) |> 
              mutate(data = "cold", trait = "Temperature niche")) |> 
  bind_rows(f_apply_Rhat(l_abu_A_intermediate$fit) |> 
              mutate(data = "interm.", trait = "Temperature niche")) |> 
  bind_rows(f_apply_Rhat(l_abu_A_warm$fit) |> 
              mutate(data = "warm", trait = "Temperature niche")) |> 
  bind_rows(f_apply_Rhat(l_abu_A_spec_m$fit) |> 
              mutate(data = "mono.", trait = "Specialisation")) |> 
  bind_rows(f_apply_Rhat(l_abu_A_spec_o$fit) |> 
              mutate(data = "oligo.", trait = "Specialisation")) |> 
  bind_rows(f_apply_Rhat(l_abu_A_spec_p$fit) |> 
              mutate(data = "poly.", trait = "Specialisation")) |> 
  bind_rows(f_apply_Rhat(l_abu_A_egg$fit) |> 
              mutate(data = "egg", trait = "Overwintering stage")) |> 
  bind_rows(f_apply_Rhat(l_abu_A_larva$fit) |> 
              mutate(data = "larva", trait = "Overwintering stage")) |> 
  bind_rows(f_apply_Rhat(l_abu_A_pupa$fit) |> 
              mutate(data = "pupa", trait = "Overwintering stage")) |> 
  bind_rows(f_apply_Rhat(l_abu_A_adult$fit) |> 
              mutate(data = "adult", trait = "Overwintering stage")) |> 
  mutate(data = factor(data, levels = unique(data)),
         trait = factor(trait, levels = unique(trait)))

d_Rhat_ric <- f_apply_Rhat(l_ric_A$fit) |> 
  mutate(data = "data", trait = "Full") |> 
  bind_rows(f_apply_Rhat(l_ric_A_small$fit) |> 
              mutate(data = "small", trait = "Body size")) |> 
  bind_rows(f_apply_Rhat(l_ric_A_medium$fit) |> 
              mutate(data = "medium", trait = "Body size")) |> 
  bind_rows(f_apply_Rhat(l_ric_A_large$fit) |> 
              mutate(data = "large", trait = "Body size")) |> 
  bind_rows(f_apply_Rhat(l_ric_A_cold$fit) |> 
              mutate(data = "cold", trait = "Temperature niche")) |> 
  bind_rows(f_apply_Rhat(l_ric_A_intermediate$fit) |> 
              mutate(data = "interm.", trait = "Temperature niche")) |> 
  bind_rows(f_apply_Rhat(l_ric_A_warm$fit) |> 
              mutate(data = "warm", trait = "Temperature niche")) |> 
  bind_rows(f_apply_Rhat(l_ric_A_spec_m$fit) |> 
              mutate(data = "mono.", trait = "Specialisation")) |> 
  bind_rows(f_apply_Rhat(l_ric_A_spec_o$fit) |> 
              mutate(data = "oligo.", trait = "Specialisation")) |> 
  bind_rows(f_apply_Rhat(l_ric_A_spec_p$fit) |> 
              mutate(data = "poly.", trait = "Specialisation")) |> 
  bind_rows(f_apply_Rhat(l_ric_A_egg$fit) |> 
              mutate(data = "egg", trait = "Overwintering stage")) |> 
  bind_rows(f_apply_Rhat(l_ric_A_larva$fit) |> 
              mutate(data = "larva", trait = "Overwintering stage")) |> 
  bind_rows(f_apply_Rhat(l_ric_A_pupa$fit) |> 
              mutate(data = "pupa", trait = "Overwintering stage")) |> 
  bind_rows(f_apply_Rhat(l_ric_A_adult$fit) |> 
              mutate(data = "adult", trait = "Overwintering stage")) |> 
  mutate(data = factor(data, levels = unique(data)),
         trait = factor(trait, levels = unique(trait)))

d_Rhat_mass <- f_apply_Rhat(l_mass_A$fit) |> 
  mutate(data = "data", trait = "Full") |> 
  bind_rows(f_apply_Rhat(l_mass_A_small$fit) |> 
              mutate(data = "small", trait = "Body size")) |> 
  bind_rows(f_apply_Rhat(l_mass_A_medium$fit) |> 
              mutate(data = "medium", trait = "Body size")) |> 
  bind_rows(f_apply_Rhat(l_mass_A_large$fit) |> 
              mutate(data = "large", trait = "Body size")) |> 
  bind_rows(f_apply_Rhat(l_mass_A_cold$fit) |> 
              mutate(data = "cold", trait = "Temperature niche")) |> 
  bind_rows(f_apply_Rhat(l_mass_A_intermediate$fit) |> 
              mutate(data = "interm.", trait = "Temperature niche")) |> 
  bind_rows(f_apply_Rhat(l_mass_A_warm$fit) |> 
              mutate(data = "warm", trait = "Temperature niche")) |> 
  bind_rows(f_apply_Rhat(l_mass_A_spec_m$fit) |> 
              mutate(data = "mono.", trait = "Specialisation")) |> 
  bind_rows(f_apply_Rhat(l_mass_A_spec_o$fit) |> 
              mutate(data = "oligo.", trait = "Specialisation")) |> 
  bind_rows(f_apply_Rhat(l_mass_A_spec_p$fit) |> 
              mutate(data = "poly.", trait = "Specialisation")) |> 
  bind_rows(f_apply_Rhat(l_mass_A_egg$fit) |> 
              mutate(data = "egg", trait = "Overwintering stage")) |> 
  bind_rows(f_apply_Rhat(l_mass_A_larva$fit) |> 
              mutate(data = "larva", trait = "Overwintering stage")) |> 
  bind_rows(f_apply_Rhat(l_mass_A_pupa$fit) |> 
              mutate(data = "pupa", trait = "Overwintering stage")) |> 
  bind_rows(f_apply_Rhat(l_mass_A_adult$fit) |> 
              mutate(data = "adult", trait = "Overwintering stage")) |> 
  mutate(data = factor(data, levels = unique(data)),
         trait = factor(trait, levels = unique(trait)))

d_Rhat_abu |> 
  mutate(response = "Abundance") |> 
  bind_rows(d_Rhat_ric |> 
              mutate(response = "Richness"))|> 
  bind_rows(d_Rhat_mass |> 
              mutate(response = "Biomass")) |> 
  group_by(response, trait, data) |> 
  summarise(prop_thresh = mean(rhat < 1.1)) |> 
  pivot_wider(values_from = prop_thresh, names_from = response)

# ... Table 1 / Table S3 #######################################################
################################################################################.

d_table <- lapply(list(Overall__overall = l_abu_A$fit, 
                       mass_cat__small = l_abu_A_small$fit,
                       mass_cat__medium = l_abu_A_medium$fit,
                       mass_cat__large = l_abu_A_large$fit,
                       Tavg_mean_cat__cold = l_abu_A_cold$fit,
                       Tavg_mean_cat__intermediate = l_abu_A_intermediate$fit,
                       Tavg_mean_cat__warm = l_abu_A_warm$fit,
                       Spec__monophagous = l_abu_A_spec_m$fit,
                       Spec__oligophagous = l_abu_A_spec_o$fit,
                       Spec__polyphagous = l_abu_A_spec_p$fit,
                       overwintering_stage__egg = l_abu_A_egg$fit,
                       overwintering_stage__larva = l_abu_A_larva$fit,
                       overwintering_stage__pupa = l_abu_A_pupa$fit,
                       overwintering_stage__adult = l_abu_A_adult$fit),
                  f_A_height_change_numbers,
                  data = d_mod_z,
                  terms = brmsterms(update(abu_tot ~ s(yday) + P_2day + T_2day + 
                                             C(traptype, "contr.sum") + 
                                             C(bulbtype, "contr.sum") + 
                                             n_trap + C(sample_previous, "contr.sum") + 
                                             (1 | spattemp_cluster) + (1 | trap_ID) + 
                                             (1 | night_ID) + (1 | trap_ID_A), 
                                           ". ~ A * height + ."))) |> 
  bind_rows(.id = "trait__traitvalue") |> 
  separate(trait__traitvalue, into = c("trait", "traitvalue"), sep = "__") |> 
  mutate(var = "Abundance") |> 
  bind_rows(lapply(list(Overall__overall = l_ric_A$fit, 
                        mass_cat__small = l_ric_A_small$fit,
                        mass_cat__medium = l_ric_A_medium$fit,
                        mass_cat__large = l_ric_A_large$fit,
                        Tavg_mean_cat__cold = l_ric_A_cold$fit,
                        Tavg_mean_cat__intermediate = l_ric_A_intermediate$fit,
                        Tavg_mean_cat__warm = l_ric_A_warm$fit,
                        Spec__monophagous = l_ric_A_spec_m$fit,
                        Spec__oligophagous = l_ric_A_spec_o$fit,
                        Spec__polyphagous = l_ric_A_spec_p$fit,
                        overwintering_stage__egg = l_ric_A_egg$fit,
                        overwintering_stage__larva = l_ric_A_larva$fit,
                        overwintering_stage__pupa = l_ric_A_pupa$fit,
                        overwintering_stage__adult = l_ric_A_adult$fit),
                   f_A_height_change_numbers,
                   data = d_mod_z,
                   terms = brmsterms(update(SCcorr_ric ~ s(yday) + P_2day + T_2day + 
                                              C(traptype, "contr.sum") + 
                                              C(bulbtype, "contr.sum") + 
                                              n_trap + C(sample_previous, "contr.sum") + 
                                              (1 | spattemp_cluster) + (1 | trap_ID) + 
                                              (1 | night_ID) + (1 | trap_ID_A), 
                                            ". ~ A * height + ."))) |> 
              bind_rows(.id = "trait__traitvalue") |> 
              separate(trait__traitvalue, into = c("trait", "traitvalue"), sep = "__") |> 
              mutate(var = "Richness")) |> 
  bind_rows(lapply(list(Overall__overall = l_mass_A$fit, 
                        mass_cat__small = l_mass_A_small$fit,
                        mass_cat__medium = l_mass_A_medium$fit,
                        mass_cat__large = l_mass_A_large$fit,
                        Tavg_mean_cat__cold = l_mass_A_cold$fit,
                        Tavg_mean_cat__intermediate = l_mass_A_intermediate$fit,
                        Tavg_mean_cat__warm = l_mass_A_warm$fit,
                        Spec__monophagous = l_mass_A_spec_m$fit,
                        Spec__oligophagous = l_mass_A_spec_o$fit,
                        Spec__polyphagous = l_mass_A_spec_p$fit,
                        overwintering_stage__egg = l_mass_A_egg$fit,
                        overwintering_stage__larva = l_mass_A_larva$fit,
                        overwintering_stage__pupa = l_mass_A_pupa$fit,
                        overwintering_stage__adult = l_mass_A_adult$fit),
                   f_A_height_change_numbers,
                   data = d_mod_z,
                   terms = brmsterms(update(mass_tot ~ s(yday) + P_2day + T_2day + 
                                              C(traptype, "contr.sum") + 
                                              C(bulbtype, "contr.sum") + 
                                              n_trap + C(sample_previous, "contr.sum") + 
                                              (1 | spattemp_cluster) + (1 | trap_ID) + 
                                              (1 | night_ID) + (1 | trap_ID_A), 
                                            ". ~ A * height + ."))) |> 
              bind_rows(.id = "trait__traitvalue") |> 
              separate(trait__traitvalue, into = c("trait", "traitvalue"), sep = "__") |> 
              mutate(var = "Biomass"))

d_scalings_height <- d_scalings |> 
  filter(data == "full",
         var == "height")

d_h_threshold <- lapply(list(Overall__overall = l_abu_A$fit, 
                             mass_cat__small = l_abu_A_small$fit,
                             mass_cat__medium = l_abu_A_medium$fit,
                             mass_cat__large = l_abu_A_large$fit,
                             Tavg_mean_cat__cold = l_abu_A_cold$fit,
                             Tavg_mean_cat__intermediate = l_abu_A_intermediate$fit,
                             Tavg_mean_cat__warm = l_abu_A_warm$fit,
                             Spec__monophagous = l_abu_A_spec_m$fit,
                             Spec__oligophagous = l_abu_A_spec_o$fit,
                             Spec__polyphagous = l_abu_A_spec_p$fit,
                             overwintering_stage__egg = l_abu_A_egg$fit,
                             overwintering_stage__larva = l_abu_A_larva$fit,
                             overwintering_stage__pupa = l_abu_A_pupa$fit,
                             overwintering_stage__adult = l_abu_A_adult$fit),
                        f_extract_slopes) |> 
  bind_rows(.id = "trait__traitvalue") |> 
  separate(trait__traitvalue, into = c("trait", "traitvalue"), sep = "__") |> 
  mutate(var = "Abundance") |> 
  bind_rows(lapply(list(Overall__overall = l_ric_A$fit, 
                        mass_cat__small = l_ric_A_small$fit,
                        mass_cat__medium = l_ric_A_medium$fit,
                        mass_cat__large = l_ric_A_large$fit,
                        Tavg_mean_cat__cold = l_ric_A_cold$fit,
                        Tavg_mean_cat__intermediate = l_ric_A_intermediate$fit,
                        Tavg_mean_cat__warm = l_ric_A_warm$fit,
                        Spec__monophagous = l_ric_A_spec_m$fit,
                        Spec__oligophagous = l_ric_A_spec_o$fit,
                        Spec__polyphagous = l_ric_A_spec_p$fit,
                        overwintering_stage__egg = l_ric_A_egg$fit,
                        overwintering_stage__larva = l_ric_A_larva$fit,
                        overwintering_stage__pupa = l_ric_A_pupa$fit,
                        overwintering_stage__adult = l_ric_A_adult$fit),
                   f_extract_slopes) |> 
              bind_rows(.id = "trait__traitvalue") |> 
              separate(trait__traitvalue, into = c("trait", "traitvalue"), sep = "__") |> 
              mutate(var = "Richness")) |> 
  bind_rows(lapply(list(Overall__overall = l_mass_A$fit, 
                        mass_cat__small = l_mass_A_small$fit,
                        mass_cat__medium = l_mass_A_medium$fit,
                        mass_cat__large = l_mass_A_large$fit,
                        Tavg_mean_cat__cold = l_mass_A_cold$fit,
                        Tavg_mean_cat__intermediate = l_mass_A_intermediate$fit,
                        Tavg_mean_cat__warm = l_mass_A_warm$fit,
                        Spec__monophagous = l_mass_A_spec_m$fit,
                        Spec__oligophagous = l_mass_A_spec_o$fit,
                        Spec__polyphagous = l_mass_A_spec_p$fit,
                        overwintering_stage__egg = l_mass_A_egg$fit,
                        overwintering_stage__larva = l_mass_A_larva$fit,
                        overwintering_stage__pupa = l_mass_A_pupa$fit,
                        overwintering_stage__adult = l_mass_A_adult$fit),
                   f_extract_slopes) |> 
              bind_rows(.id = "trait__traitvalue") |> 
              separate(trait__traitvalue, into = c("trait", "traitvalue"), sep = "__") |> 
              mutate(var = "Biomass")) |> 
  mutate(h_threshold = (-A / `A:height` + mean(d_mod_z$height)) * d_scalings_height$sd + d_scalings_height$mean) |> 
  group_by(var, trait, traitvalue) |> 
  summarise(h_threshold_mean = mean(h_threshold),
            h_threshold_lower = ci(h_threshold)$CI_low,
            h_threshold_upper = ci(h_threshold)$CI_high)

d_h_threshold <- d_h_threshold |> 
  mutate(h_threshold_mean = round(h_threshold_mean),
         h_threshold_mean = case_when(h_threshold_mean < min(d_mod$height) ~ "<min",
                                      h_threshold_mean > max(d_mod$height) ~ ">max",
                                      .default = as.character(h_threshold_mean)),
         h_threshold_lower = round(h_threshold_lower),
         h_threshold_lower = case_when(h_threshold_lower < min(d_mod$height) ~ "<min",
                                       h_threshold_lower > max(d_mod$height) ~ ">max",
                                       .default = as.character(h_threshold_lower)),
         h_threshold_upper = round(h_threshold_upper),
         h_threshold_upper = case_when(h_threshold_upper < min(d_mod$height) ~ "<min",
                                       h_threshold_upper > max(d_mod$height) ~ ">max",
                                       .default = as.character(h_threshold_upper)),
         h_threshold = paste0(h_threshold_mean, "\n(",
                              h_threshold_lower, " – ",
                              h_threshold_upper, ")")) |> 
  select(var, trait, traitvalue, h_threshold)


d_table |> 
  mutate(perc_estimate = (factor_estimate - 1) * 100,
         perc_lower = (factor_lower - 1) * 100,
         perc_upper = (factor_upper - 1) * 100) |> 
  select(-c(diff_log_estimate, diff_log_lower, diff_log_upper)) |> 
  mutate(text_factor = paste0(formatC(factor_estimate, digits = 3, format = "fg", flag = "#"), 
                              " (", 
                              formatC(factor_lower, digits = 3, format = "fg", flag = "#"), "–", 
                              formatC(factor_upper, digits = 3, format = "fg", flag = "#"), ")"),
         text_perc = paste0(case_when(abs(perc_estimate) >= 100 ~ formatC(perc_estimate, format = "f", flag = "+", digits = 0),
                                      abs(perc_estimate) >= 10 ~ formatC(perc_estimate, format = "f", flag = "+", digits = 1),
                                      abs(perc_estimate) >= 1 ~ formatC(perc_estimate, format = "f", flag = "+", digits = 2),
                                      .default = formatC(perc_estimate, digits = 3, format = "fg", flag = "+")), 
                            "% (", 
                            case_when(abs(perc_lower) >= 100 ~ formatC(perc_lower, format = "f", flag = "+", digits = 0),
                                      abs(perc_lower) >= 10 ~ formatC(perc_lower, format = "f", flag = "+", digits = 1),
                                      abs(perc_lower) >= 1 ~ formatC(perc_lower, format = "f", flag = "+", digits = 2),
                                      .default = formatC(perc_lower, digits = 3, format = "fg", flag = "+")), 
                            "%–", 
                            case_when(abs(perc_upper) >= 100 ~ formatC(perc_upper, format = "f", flag = "+", digits = 0),
                                      abs(perc_upper) >= 10 ~ formatC(perc_upper, format = "f", flag = "+", digits = 1),
                                      abs(perc_upper) >= 1 ~ formatC(perc_upper, format = "f", flag = "+", digits = 2),
                                      .default = formatC(perc_upper, digits = 3, format = "fg", flag = "+")), 
                            "%)")) |> 
  select(-c(factor_estimate, factor_lower, factor_upper, perc_estimate, perc_lower, perc_upper)) |> 
  mutate(text = paste(text_factor, text_perc, sep = "\n")) |> 
  select(-c(text_factor, text_perc)) |> 
  pivot_wider(names_from = height, values_from = text) |> 
  left_join(d_h_threshold, by = c("var", "trait", "traitvalue")) |> 
  select(var, everything()) |> 
  write_excel_csv2("Output/Tables/Table1.csv")

# ... Figure 1 #################################################################
################################################################################.

# Europe map -------------------------------------------------------------------.
sf_europe <- gisco_get_countries(region = "Europe")

r_elevEurope <- raster("Data/Elevation_models/dtm_elev.lowestmode_meritdem_m_100m_0..0cm_2000..2018_eumap_epsg3035_v1.01.tif")
r_elevEurope <- aggregate(r_elevEurope, fact = 100, fu = mean)
r_elevEurope <- projectRaster(r_elevEurope, 
                              crs = CRS("+proj=longlat +datum=WGS84 +no_defs"))
r_elevEurope <- mask(r_elevEurope, sf_europe)
# plot(r_elevEurope)

d_elevEurope <- as.data.frame(r_elevEurope, xy = TRUE)
names(d_elevEurope) <- c("x", "y", "elev")
d_elevEurope <- d_elevEurope |> 
  filter(!is.na(elev))

# to cut off overseas areas
cut_pol <-  st_sfc(st_polygon(list(cbind(c(-5, 40, 45, -40, -5),
                                         c(25, 30, 75, 65, 25)))),
                   crs = 4326) %>% 
  st_transform(crs = st_crs(sf_europe)) # harmonize coordinate system

p_europe <-
  sf_europe %>%
  ggplot() +
  geom_raster(data = d_elevEurope, aes(x = x, y = y, fill = elev)) +
  geom_sf(col = "grey5", fill = NA, size = .1) +
  geom_sf(data = sf_europe %>% filter(CNTR_ID == "CH"), fill = NA, 
          col = "red3", size = .5) +
  coord_sf(xlim = c(-9, 20),
           ylim = c(37, 60)) +
  scale_fill_gradient(low = "grey80", high = "grey10") +
  theme(plot.margin = ggplot2::margin(10, 10, 10, 10, "pt"),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.background = element_rect(fill = NA, colour = NA),
        panel.background = element_rect(fill = NA, colour = NA),
        legend.position = "none") 

# maps Switzerland -------------------------------------------------------------.

sf_ch_border <- gisco_get_countries(country = "CH", resolution = "01") |> 
  st_transform(crs = 2056)

# elevation model of Switzerland, available from Federal Office of Topography swisstopo
# https://www.swisstopo.admin.ch/en/height-model-dhm25
r_relief <- raster("Data/Elevation_models/ASCII_GRID_1part/dhm25_grid_raster.asc")
r_relief <- aggregate(r_relief, fact = 10, fu = mean)
crs(r_relief) <- CRS("+init=epsg:21781")
r_relief <- projectRaster(r_relief, crs = CRS("+init=epsg:2056"))
# and again back to data.frame
d_relief <-
  r_relief %>%
  # hide relief outside of Switzerland by masking with country borders
  mask(sf_ch_border) %>%
  # as("SpatialPixelsDataFrame") %>%
  as.data.frame(xy = T) %>%
  rename(value = dhm25_grid_raster) %>%
  filter(!is.na(value))
rm(r_relief)

scale_bar <- data.frame(int5 = 2017,  # matches only this facet
                        geometry = st_sfc(st_linestring(rbind(c(-80, 34), c(-79.5, 34))))) |> 
  st_as_sf(crs = st_crs(sf_ch_border))

sf_samplings <-
  d_samplings %>%
  mutate(int5 = floor((A - 2) / 5) * 5 + 2) %>%
  select(int5, LOC) %>%
  distinct() %>%
  left_join(d_moths_raw |> 
              # the following are the approximate LOC coordinates
              # the original figure contains the more precise coordinates
              select(LOC = locality, decimalLatitude, decimalLongitude) |> 
              distinct(),
            by = "LOC") |> 
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"),
           crs = 4326) |> 
  st_transform(crs = 2056)

l_plots_maps_CH <- list()
for (int5_i in sort(unique(sf_samplings$int5))){
  l_plots_maps_CH[[as.character(int5_i)]] <- sf_samplings |> 
    filter(int5 == int5_i) |> 
    ggplot() +
    geom_sf(data = sf_ch_border) +
    geom_raster(data = d_relief, inherit.aes = FALSE,
                aes(x = x, y = y, fill = value)) +
    geom_sf(alpha = .5, col = "red3") +
    scale_fill_gradient(low = "grey80", high = "grey10") +
    facet_wrap(~ int5,
               labeller = as_labeller(c("1972" = "1972–1976", "1977" = "1977–1981",
                                        "1982" = "1982–1986", "1987" = "1987–1991",
                                        "1992" = "1992–1996", "1997" = "1997–2001",
                                        "2002" = "2002–2006", "2007" = "2007–2011",
                                        "2012" = "2012–2016", "2017" = "2017–2021"))) +
    theme(axis.text = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          legend.position = "none",
          plot.margin = ggplot2::margin(0, 0, 0, 0, "pt"),
          strip.text = element_text(size = v_textsize["axis.title"]))
  
  if (int5_i == 1972){
    l_plots_maps_CH[[as.character(int5_i)]] <- l_plots_maps_CH[[as.character(int5_i)]] +
      annotation_scale(text_cex = ggplot2:::.pt/v_textsize["additional.text"],
                       location = "tl", pad_y = unit(0.1, "cm"))
  }
  
  
}

p_maps_CH <-plot_grid(plotlist = l_plots_maps_CH, nrow = 2)


# time versus elevation sampling -----------------------------------------------.

p_elev_A <- d_samplings %>%
  select(LOC, A) %>%
  distinct() |> 
  left_join(d_sites %>%
              select(LOC, height) %>%
              distinct(), 
            by = "LOC") %>%
  ggplot(aes(x = A, y = height)) +
  geom_point(colour = "red3", alpha = .5) +
  labs(x = "Year", y = "Elevation") +
  theme(plot.margin = ggplot2::margin(10, 10, 10, 10, "pt"),
        axis.title = element_text(size = v_textsize["axis.title"]),
        axis.text = element_text(size = v_textsize["axis.text"]))

# combined plot ----------------------------------------------------------------.

p_comb <- plot_grid(
  plot_grid(
    p_europe, p_elev_A, 
    nrow = 1, rel_widths = c(1, 3.5),
    labels = c("(a)", "(b)"), label_size = v_textsize["axis.title"]
  ),
  p_maps_CH,
  ncol = 1, rel_heights = c(1, 1.5),
  labels = c("", "(c)"), label_size = v_textsize["axis.title"],
  label_y = 1.072
)

ggsave(p_comb, file = "Output/Figures/Spatiotemporal_coverage2.pdf",
       height = 100, width = 180, units = "mm", dpi = 600)

# ... Figure 2 #################################################################
################################################################################.

plot_grid(
  f_A_height_plot(pred = l_abu_A$l_pred_fe$`A:height`,
                  data_raw = d_mod_z,
                  mean_yday = d_scalings$mean[d_scalings$data == "full" & d_scalings$var == "yday"],
                  response = "abu_tot") +
    scale_y_continuous(breaks = c(0, 10, 100, 1000, 10000), trans = "log1p") +
    labs(x = "Year", y = "Abundance") +
    theme(legend.position = "none",
          axis.text = element_text(size = v_textsize["axis.text"]),
          axis.title = element_text(size = v_textsize["axis.title"])),
  f_A_height_plot(pred = l_ric_A$l_pred_fe$`A:height`,
                  data_raw = d_mod_z,
                  mean_yday = d_scalings$mean[d_scalings$data == "full" & d_scalings$var == "yday"],
                  response = "SCcorr_ric") +
    scale_y_continuous(breaks = c(0, 10, 100, 1000), trans = "log1p") +
    labs(x = "Year", y = "Richness") +
    theme(legend.position = "none",
          axis.text = element_text(size = v_textsize["axis.text"]),
          axis.title = element_text(size = v_textsize["axis.title"])),
  f_A_height_plot(pred = l_mass_A$l_pred_fe$`A:height`,
                  data_raw = d_mod_z,
                  mean_yday = d_scalings$mean[d_scalings$data == "full" & d_scalings$var == "yday"],
                  response = "mass_tot") +
    scale_y_continuous(breaks = c(0, 0.01, 0.1, 1, 10, 100),
                       labels = c(0, 0.01, 0.1, 1, 10, 100),
                       trans = log_plus_trans) +
    labs(x = "Year", y = "Biomass (g)", colour = "Elevation") +
    guides(colour = guide_legend(reverse = T)) +
    theme(axis.text = element_text(size = v_textsize["axis.text"]),
          axis.title = element_text(size = v_textsize["axis.title"]),
          legend.title = element_text(size = v_textsize["legend.title"]),
          legend.text = element_text(size = v_textsize["axis.text"])),
  nrow = 1, rel_widths = c(1.05, 1, 1.6)
) 

ggsave("Output/Figures/Year_trends.pdf", width = 180, height = 65,
       units = "mm", dpi = 600)

# ... Figure 3 / Figures S12 & S13 #############################################
################################################################################.

p_legend <- cowplot::get_legend(f_A_height_plot_comb(list(small = l_abu_A_small$l_pred_fe$`A:height`),
                                                     response = "abu_tot") +
                                  guides(colour = guide_legend(reverse = T)) +
                                  theme(legend.title = element_text(size = v_textsize["legend.title"]),
                                        legend.text = element_text(size = v_textsize["axis.text"])))
p_empty <- ggplot() + theme_nothing()

# abundance --------------------------------------------------------------------.
plot_grid(
  plot_grid(
    plot_grid(
      f_A_height_plot_comb(list(small = l_abu_A_small$l_pred_fe$`A:height`,
                                medium = l_abu_A_medium$l_pred_fe$`A:height`,
                                large = l_abu_A_large$l_pred_fe$`A:height`),
                           response = "abu_tot",
                           name = "Body size") +
        labs(y = "Abundance") +
        theme(axis.title.x = element_blank(),
              legend.position = "none",
              axis.title.y = element_text(size = v_textsize["axis.title"]),
              axis.text = element_text(size = v_textsize["axis.text"]),
              strip.text.x = element_text(size = v_textsize["axis.title"])),
      f_A_height_plot_comb(list(cold = l_abu_A_cold$l_pred_fe$`A:height`,
                                intermediate = l_abu_A_intermediate$l_pred_fe$`A:height`,
                                warm = l_abu_A_warm$l_pred_fe$`A:height`),
                           response = "abu_tot",
                           name = "Temperature niche") +
        labs(y = "Abundance") +
        theme(axis.title.x = element_blank(),
              legend.position = "none",
              axis.title.y = element_text(size = v_textsize["axis.title"]),
              axis.text = element_text(size = v_textsize["axis.text"]),
              strip.text.x = element_text(size = v_textsize["axis.title"])),
      f_A_height_plot_comb(list(monophagous = l_abu_A_spec_m$l_pred_fe$`A:height`,
                                oligophagous = l_abu_A_spec_o$l_pred_fe$`A:height`,
                                polyphagous = l_abu_A_spec_p$l_pred_fe$`A:height`),
                           response = "abu_tot",
                           name = "Food specialisation") +
        labs(y = "Abundance") +
        theme(axis.title.x = element_blank(),
              legend.position = "none",
              axis.title.y = element_text(size = v_textsize["axis.title"]),
              axis.text = element_text(size = v_textsize["axis.text"]),
              strip.text.x = element_text(size = v_textsize["axis.title"])),
      
      ncol = 1, align = "v", 
      labels = c("(a)", "(b)", "(c)"), label_size = v_textsize["axis.title"]),
    plot_grid(p_empty, p_legend, p_empty, rel_heights = c(3, 2, 1), ncol = 1),
    nrow = 1, rel_widths = c(1, .305)),
  f_A_height_plot_comb(list(egg = l_abu_A_egg$l_pred_fe$`A:height`,
                            larva = l_abu_A_larva$l_pred_fe$`A:height`,
                            pupa = l_abu_A_pupa$l_pred_fe$`A:height`,
                            adult = l_abu_A_adult$l_pred_fe$`A:height`),
                       response = "abu_tot",
                       name = "Overwintering stage") +
    labs(y = "Abundance") +
    theme(legend.position = "none",
          axis.title = element_text(size = v_textsize["axis.title"]),
          axis.title.y = element_text(margin = unit(c(0, .4, 0, 0), "cm")),
          axis.text = element_text(size = v_textsize["axis.text"]),
          strip.text.x = element_text(size = v_textsize["axis.title"])),
  ncol = 1, rel_heights = c(3, 1.1),
  labels = c("", "(d)"), label_size = v_textsize["axis.title"])

ggsave("Output/Figures/Trait_Year_trends_abu.pdf", width = 180, height = 200,
       units = "mm", dpi = 600)

# richness ---------------------------------------------------------------------.

plot_grid(
  plot_grid(
    plot_grid(
      f_A_height_plot_comb(list(small = l_ric_A_small$l_pred_fe$`A:height`,
                                medium = l_ric_A_medium$l_pred_fe$`A:height`,
                                large = l_ric_A_large$l_pred_fe$`A:height`),
                           response = "SCcorr_ric",
                           name = "Body size") +
        labs(y = "Richness") +
        theme(axis.title.x = element_blank(),
              legend.position = "none",
              axis.title.y = element_text(size = v_textsize["axis.title"]),
              axis.text = element_text(size = v_textsize["axis.text"]),
              strip.text.x = element_text(size = v_textsize["axis.title"])),
      f_A_height_plot_comb(list(cold = l_ric_A_cold$l_pred_fe$`A:height`,
                                intermediate = l_ric_A_intermediate$l_pred_fe$`A:height`,
                                warm = l_ric_A_warm$l_pred_fe$`A:height`),
                           response = "SCcorr_ric",
                           name = "Temperature niche") +
        labs(y = "Richness") +
        theme(axis.title.x = element_blank(),
              legend.position = "none",
              axis.title.y = element_text(size = v_textsize["axis.title"]),
              axis.text = element_text(size = v_textsize["axis.text"]),
              strip.text.x = element_text(size = v_textsize["axis.title"])),
      f_A_height_plot_comb(list(monophagous = l_ric_A_spec_m$l_pred_fe$`A:height`,
                                oligophagous = l_ric_A_spec_o$l_pred_fe$`A:height`,
                                polyphagous = l_ric_A_spec_p$l_pred_fe$`A:height`),
                           response = "SCcorr_ric",
                           name = "Food specialisation") +
        labs(y = "Richness") +
        theme(axis.title.x = element_blank(),
              legend.position = "none",
              axis.title.y = element_text(size = v_textsize["axis.title"]),
              axis.text = element_text(size = v_textsize["axis.text"]),
              strip.text.x = element_text(size = v_textsize["axis.title"])),
      
      ncol = 1, align = "v", 
      labels = c("(a)", "(b)", "(c)"), label_size = v_textsize["axis.title"]),
    plot_grid(p_empty, p_legend, p_empty, rel_heights = c(3, 2, 1), ncol = 1),
    nrow = 1, rel_widths = c(1, .32)),
  f_A_height_plot_comb(list(egg = l_ric_A_egg$l_pred_fe$`A:height`,
                            larva = l_ric_A_larva$l_pred_fe$`A:height`,
                            pupa = l_ric_A_pupa$l_pred_fe$`A:height`,
                            adult = l_ric_A_adult$l_pred_fe$`A:height`),
                       response = "SCcorr_ric",
                       name = "Overwintering stage") +
    labs(y = "Richness") +
    theme(legend.position = "none",
          axis.title = element_text(size = v_textsize["axis.title"]),
          axis.title.y = element_text(margin = unit(c(0, 0, 0, 0), "cm")),
          axis.text = element_text(size = v_textsize["axis.text"]),
          strip.text.x = element_text(size = v_textsize["axis.title"])),
  ncol = 1, rel_heights = c(3, 1.1),
  labels = c("", "(d)"), label_size = v_textsize["axis.title"])

ggsave("Output/Figures/Trait_Year_trends_ric.pdf", width = 180, height = 200,
       units = "mm", dpi = 600)

# biomass ----------------------------------------------------------------------.

plot_grid(
  plot_grid(
    plot_grid(
      f_A_height_plot_comb(list(small = l_mass_A_small$l_pred_fe$`A:height`,
                                medium = l_mass_A_medium$l_pred_fe$`A:height`,
                                large = l_mass_A_large$l_pred_fe$`A:height`),
                           response = "mass_tot",
                           name = "Body size") +
        labs(y = "Biomass (g)") +
        theme(axis.title.x = element_blank(),
              legend.position = "none",
              axis.title.y = element_text(size = v_textsize["axis.title"]),
              axis.text = element_text(size = v_textsize["axis.text"]),
              strip.text.x = element_text(size = v_textsize["axis.title"])),
      f_A_height_plot_comb(list(cold = l_mass_A_cold$l_pred_fe$`A:height`,
                                intermediate = l_mass_A_intermediate$l_pred_fe$`A:height`,
                                warm = l_mass_A_warm$l_pred_fe$`A:height`),
                           response = "mass_tot",
                           name = "Temperature niche") +
        labs(y = "Biomass (g)") +
        theme(axis.title.x = element_blank(),
              legend.position = "none",
              axis.title.y = element_text(size = v_textsize["axis.title"]),
              axis.text = element_text(size = v_textsize["axis.text"]),
              strip.text.x = element_text(size = v_textsize["axis.title"])),
      f_A_height_plot_comb(list(monophagous = l_mass_A_spec_m$l_pred_fe$`A:height`,
                                oligophagous = l_mass_A_spec_o$l_pred_fe$`A:height`,
                                polyphagous = l_mass_A_spec_p$l_pred_fe$`A:height`),
                           response = "mass_tot",
                           name = "Food specialisation") +
        labs(y = "Biomass (g)") +
        theme(axis.title.x = element_blank(),
              legend.position = "none",
              axis.title.y = element_text(size = v_textsize["axis.title"]),
              axis.text = element_text(size = v_textsize["axis.text"]),
              strip.text.x = element_text(size = v_textsize["axis.title"])),
      
      ncol = 1, align = "v", 
      labels = c("(a)", "(b)", "(c)"), label_size = v_textsize["axis.title"]),
    plot_grid(p_empty, p_legend, p_empty, rel_heights = c(3, 2, 1), ncol = 1),
    nrow = 1, rel_widths = c(1, .305)),
  f_A_height_plot_comb(list(egg = l_mass_A_egg$l_pred_fe$`A:height`,
                            larva = l_mass_A_larva$l_pred_fe$`A:height`,
                            pupa = l_mass_A_pupa$l_pred_fe$`A:height`,
                            adult = l_mass_A_adult$l_pred_fe$`A:height`),
                       response = "mass_tot",
                       name = "Overwintering stage") +
    labs(y = "Biomass (g)") +
    theme(legend.position = "none",
          axis.title = element_text(size = v_textsize["axis.title"]),
          axis.title.y = element_text(margin = unit(c(0, .4, 0, 0), "cm")),
          axis.text = element_text(size = v_textsize["axis.text"]),
          strip.text.x = element_text(size = v_textsize["axis.title"])),
  ncol = 1, rel_heights = c(3, 1.1),
  labels = c("", "(d)"), label_size = v_textsize["axis.title"])

ggsave("Output/Figures/Trait_Year_trends_mass.pdf", width = 180, height = 200,
       units = "mm", dpi = 600)

# ... Figure 4 #################################################################
################################################################################.

diff_h_or <- 1000
diff_h_z <- diff_h_or /
  d_scalings$sd[d_scalings$var == "height" & d_scalings$data == "full"]

diff_A_or <- 10
diff_A_z <- diff_A_or /
  d_scalings$sd[d_scalings$var == "A" & d_scalings$data == "full"]

d_traits_edit <- d_traits |> 
  left_join(d_moths |> 
              group_by(Name_std) |> 
              summarise(grand_total = sum(individualCount),
                        .groups = "drop"),
            by = "Name_std") |> 
  mutate(Spec = tolower(Spec),
         Spec = factor(Spec, levels = c("monophagous", "oligophagous", "polyphagous")),
         overwintering_stage = gsub("adult", "ad.", overwintering_stage),
         overwintering_stage = factor(overwintering_stage, 
                                      levels = c("egg", "larva", "pupa", "ad.")),
         mass_cat = factor(mass_cat, levels = c("small", "medium", "large")),
         Tavg_mean_cat = factor(Tavg_mean_cat, levels = c("cold", "intermediate", "warm")))

d_coefs_mean_edit <- d_coefs_mean |> 
  filter(var == "A:height") |> 
  left_join(d_traits_edit,
            by = "Name_std") |> 
  mutate(fact_mean = exp(diff_h_z * diff_A_z * mean),
         fact_lower = exp(diff_h_z * diff_A_z * CI_low),
         fact_upper = exp(diff_h_z * diff_A_z * CI_high),
         sig = fact_lower > 1 | fact_upper < 1)

d_coefs_sim_edit <- d_coefs_sim |> 
  group_by(Name_std) |> 
  mutate(iter = row_number()) |> 
  ungroup() |> 
  left_join(d_traits_edit,
            by = "Name_std") |> 
  mutate(fact = exp(diff_h_z * diff_A_z * `A:height`))



set.seed(814)
l_intcoef_plots <- list()

# mass_cat-----------------------------------------------------------------.

d_smry_mass_cat <- d_coefs_sim_edit |> 
  group_by(iter, mass_cat) |> 
  summarise(wmean = sum(fact * grand_total) / sum(grand_total),
            .groups = "drop") |> 
  group_by(mass_cat) |> 
  summarise(wmean_mean = mean(wmean),
            wmean_lower = ci(wmean)$CI_low,
            wmean_upper = ci(wmean)$CI_high,
            .groups = "drop")


l_intcoef_plots[["mass_cat"]] <-
  d_coefs_mean_edit |> 
  group_by(mass_cat) |> 
  arrange(mean) |> 
  mutate(index = row_number(),
         median = median(index),
         posneg = which(mean > 0)[1] - .5) |> 
  ungroup() %>%
  ggplot(aes(x = index)) +
  geom_rect(data = function(x) x |>
              group_by(mass_cat) %>%
              group_map(~ data.frame(quantile(sapply(1:10000, f_boot, .), 
                                              probs = c(0.025, 0.975)) %>%
                                       t(), mass_cat = unique(.$mass_cat)),
                        .keep = T) |>
              bind_rows(),
            aes(xmin = X2.5., xmax = X97.5.),
            ymin = -Inf, ymax = Inf, fill = "lightskyblue1", inherit.aes = F) +
  geom_vline(aes(xintercept = median), col = 4, lty = 2) +
  geom_vline(aes(xintercept = posneg), col = 4) +
  geom_hline(data = \(x) x |> 
               select(mass_cat) |> 
               distinct() |> 
               mutate(yintercept = 1),
             aes(yintercept = yintercept), lty = 2) +
  geom_linerange(aes(ymin = fact_lower, ymax = fact_upper), alpha = .5, size = .25) +
  geom_rect(data = d_smry_mass_cat,
            aes(ymin = wmean_lower, ymax = wmean_upper),
            xmin = -Inf, xmax = Inf,
            fill = 2, colour = NA, alpha = .7, inherit.aes = F) +
  geom_hline(data = d_smry_mass_cat,
             aes(yintercept = wmean_mean),
             colour = "red4", size = 1, alpha = .7) +
  geom_point(aes(y = fact_mean), size = .5)  +
  facet_nested(~ "Body size" + mass_cat,
               scales = "free_x", space = "free_x") +
  scale_y_log10(breaks = c(.1, .5, 1, 2, 10)) +
  coord_cartesian(ylim = c(.1, 20)) +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size = v_textsize["axis.title"]),
        axis.text.y = element_text(size = v_textsize["axis.text"]),
        axis.title.y = element_text(size = v_textsize["axis.title"])) +
  labs(y = "Interactive effect")

# Tavg_mean_cat-----------------------------------------------------------------.

d_smry_Tavg_mean_cat <- d_coefs_sim_edit |> 
  group_by(iter, Tavg_mean_cat) |> 
  summarise(wmean = sum(fact * grand_total) / sum(grand_total),
            .groups = "drop") |> 
  group_by(Tavg_mean_cat) |> 
  summarise(wmean_mean = mean(wmean),
            wmean_lower = ci(wmean)$CI_low,
            wmean_upper = ci(wmean)$CI_high,
            .groups = "drop")


l_intcoef_plots[["Tavg_mean_cat"]] <- d_coefs_mean_edit |> 
  group_by(Tavg_mean_cat) |> 
  arrange(mean) |> 
  mutate(index = row_number(),
         median = median(index),
         posneg = which(mean > 0)[1] - .5) |> 
  ungroup() %>%
  ggplot(aes(x = index)) +
  geom_rect(data = function(x) x |> 
              group_by(Tavg_mean_cat) %>%
              group_map(~ data.frame(quantile(sapply(1:10000, f_boot, .), 
                                              probs = c(0.025, 0.975)) %>%
                                       t(), Tavg_mean_cat = unique(.$Tavg_mean_cat)),
                        .keep = T) |> 
              bind_rows(),
            aes(xmin = X2.5., xmax = X97.5.), 
            ymin = -Inf, ymax = Inf, fill = "lightskyblue1", inherit.aes = F) +
  geom_vline(aes(xintercept = median), col = 4, lty = 2) +
  geom_vline(aes(xintercept = posneg), col = 4) +
  geom_hline(data = \(x) x |> 
               select(Tavg_mean_cat) |> 
               distinct() |> 
               mutate(yintercept = 1),
             aes(yintercept = yintercept), lty = 2) +
  geom_linerange(aes(ymin = fact_lower, ymax = fact_upper), alpha = .5, size = .25) +
  geom_rect(data = d_smry_Tavg_mean_cat,
            aes(ymin = wmean_lower, ymax = wmean_upper),
            xmin = -Inf, xmax = Inf,
            fill = 2, colour = NA, alpha = .7, inherit.aes = F) +
  geom_hline(data = d_smry_Tavg_mean_cat,
             aes(yintercept = wmean_mean), 
             colour = "red4", size = 1, alpha = .7) +
  geom_point(aes(y = fact_mean), size = .5)  +
  facet_nested(~ "Temperature niche" + Tavg_mean_cat,
               scales = "free_x", space = "free_x") +
  scale_y_log10(breaks = c(.1, .5, 1, 2, 10)) +
  coord_cartesian(ylim = c(.1, 20)) +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size = v_textsize["axis.title"]),
        axis.text.y = element_text(size = v_textsize["axis.text"]),
        axis.title.y = element_text(size = v_textsize["axis.title"])) +
  labs(y = "Interactive effect")

# Spec-----------------------------------------------------------------.

d_smry_Spec <- d_coefs_sim_edit |>
  filter(!is.na(Spec)) |> 
  group_by(iter, Spec) |> 
  summarise(wmean = sum(fact * grand_total) / sum(grand_total),
            .groups = "drop") |> 
  group_by(Spec) |> 
  summarise(wmean_mean = mean(wmean),
            wmean_lower = ci(wmean)$CI_low,
            wmean_upper = ci(wmean)$CI_high,
            .groups = "drop")


l_intcoef_plots[["Spec"]] <- d_coefs_mean_edit |> 
  filter(!is.na(Spec)) |> 
  group_by(Spec) |> 
  arrange(mean) |> 
  mutate(index = row_number(),
         median = median(index),
         posneg = which(mean > 0)[1] - .5) |> 
  ungroup() %>%
  ggplot(aes(x = index)) +
  geom_rect(data = function(x) x |> 
              group_by(Spec) %>%
              group_map(~ data.frame(quantile(sapply(1:10000, f_boot, .), 
                                              probs = c(0.025, 0.975)) %>%
                                       t(), Spec = unique(.$Spec)),
                        .keep = T) |> 
              bind_rows(),
            aes(xmin = X2.5., xmax = X97.5.), 
            ymin = -Inf, ymax = Inf, fill = "lightskyblue1", inherit.aes = F) +
  geom_vline(aes(xintercept = median), col = 4, lty = 2) +
  geom_vline(aes(xintercept = posneg), col = 4) +
  geom_hline(data = \(x) x |> 
               select(Spec) |> 
               distinct() |> 
               mutate(yintercept = 1),
             aes(yintercept = yintercept), lty = 2) +
  geom_linerange(aes(ymin = fact_lower, ymax = fact_upper), alpha = .5, size = .25) +
  geom_rect(data = d_smry_Spec,
            aes(ymin = wmean_lower, ymax = wmean_upper),
            xmin = -Inf, xmax = Inf,
            fill = 2, colour = NA, alpha = .7, inherit.aes = F) +
  geom_hline(data = d_smry_Spec,
             aes(yintercept = wmean_mean), 
             colour = "red4", size = 1, alpha = .7) +
  geom_point(aes(y = fact_mean), size = .5)  +
  facet_nested(~ "Food specialisation" + Spec,
               scales = "free_x", space = "free_x") +
  scale_y_log10(breaks = c(.1, .5, 1, 2, 10)) +
  coord_cartesian(ylim = c(.1, 20)) +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size = v_textsize["axis.title"]),
        axis.text.y = element_text(size = v_textsize["axis.text"]),
        axis.title.y = element_text(size = v_textsize["axis.title"])) +
  labs(y = "Interactive effect")

# overwintering_stage-----------------------------------------------------------------.

d_smry_overwintering_stage <- d_coefs_sim_edit |> 
  filter(overwintering_stage != "larva/pupa",
         !is.na(overwintering_stage)) |> 
  group_by(iter, overwintering_stage) |> 
  summarise(wmean = sum(fact * grand_total) / sum(grand_total),
            .groups = "drop") |> 
  group_by(overwintering_stage) |> 
  summarise(wmean_mean = mean(wmean),
            wmean_lower = ci(wmean)$CI_low,
            wmean_upper = ci(wmean)$CI_high,
            .groups = "drop")


l_intcoef_plots[["overwintering_stage"]] <-
  d_coefs_mean_edit |> 
  filter(overwintering_stage != "larva/pupa",
         !is.na(overwintering_stage)) |> 
  group_by(overwintering_stage) |> 
  arrange(mean) |> 
  mutate(index = row_number(),
         median = median(index),
         posneg = which(mean > 0)[1] - .5) |> 
  ungroup() %>%
  ggplot(aes(x = index)) +
  geom_rect(data = function(x) x |> 
              group_by(overwintering_stage) %>%
              group_map(~ data.frame(quantile(sapply(1:10000, f_boot, .), 
                                              probs = c(0.025, 0.975)) %>%
                                       t(), overwintering_stage = unique(.$overwintering_stage)),
                        .keep = T) |> 
              bind_rows(),
            aes(xmin = X2.5., xmax = X97.5.), 
            ymin = -Inf, ymax = Inf, fill = "lightskyblue1", inherit.aes = F) +
  geom_vline(aes(xintercept = median), col = 4, lty = 2) +
  geom_vline(aes(xintercept = posneg), col = 4) +
  geom_hline(data = \(x) x |> 
               select(overwintering_stage) |> 
               distinct() |> 
               mutate(yintercept = 1),
             aes(yintercept = yintercept), lty = 2) +
  geom_linerange(aes(ymin = fact_lower, ymax = fact_upper), alpha = .5, size = .25) +
  geom_rect(data = d_smry_overwintering_stage,
            aes(ymin = wmean_lower, ymax = wmean_upper),
            xmin = -Inf, xmax = Inf,
            fill = 2, colour = NA, alpha = .7, inherit.aes = F) +
  geom_hline(data = d_smry_overwintering_stage,
             aes(yintercept = wmean_mean), 
             colour = "red4", size = 1, alpha = .7) +
  geom_point(aes(y = fact_mean), size = .5)  +
  facet_nested(~ "Overwintering stage" + overwintering_stage,
               scales = "free_x", space = "free_x") +
  scale_y_log10(breaks = c(.1, .5, 1, 2, 10)) +
  coord_cartesian(ylim = c(.1, 20)) +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size = v_textsize["axis.title"]),
        axis.text.y = element_text(size = v_textsize["axis.text"]),
        axis.title.y = element_text(size = v_textsize["axis.title"])) +
  labs(y = "Interactive effect")

# combined plot ----------------------------------------------------------------.

plot_grid(plotlist = l_intcoef_plots,
          ncol = 1,
          labels = c("(a)", "(b)", "(c)", "(d)"), 
          label_size = v_textsize["axis.title"])

ggsave(paste0(dat.dir, "Output/Figures/Aheight_interaction.jpg"), 
       width = 180, height = 200,
       units = "mm", dpi = 600)

# ... Figure S7 ################################################################
################################################################################.

l_plots_cov1_A_abu <- c(l_abu_A$l_pred_fe, l_abu_A$l_pred_sm) %>% 
  keep(!names(.) %in% c("A", "height", "A:height")) %>%
  map(f_plot_pred, data = d_mod, response = "abu_tot", line.size = .75, 
      hours_sel = which(d_mod_z$traptype == "p" & d_mod_z$hours_data))

l_plots_cov1_A_abu <- lapply(l_plots_cov1_A_abu,
                             \(x) x +
                               theme(axis.title = element_text(size = v_textsize["axis.title"]),
                                     axis.text = element_text(size = v_textsize["axis.text"])) +
                               labs(y = "Abundance")+ 
                               scale_y_continuous(breaks = c(0, 10, 100, 1000, 10000), trans = "log1p"))

l_plots_cov1_A_abu <- l_plots_cov1_A_abu[c("yday", "P_2day", "T_2day",
                                           "traptype", "bulbtype", "n_trap",
                                           "active_hours", "sample_previous")]

l_plots_cov1_A_abu[c("P_2day", "T_2day",
                     "bulbtype", "n_trap",
                     "sample_previous")] <-
  lapply(l_plots_cov1_A_abu[c("P_2day", "T_2day",
                              "bulbtype", "n_trap",
                              "sample_previous")],
         \(x) x +
           theme(axis.title.y = element_blank(),
                 axis.text.y = element_blank()))

plot_grid(
  plot_grid(plotlist = l_plots_cov1_A_abu[c(1:3)],
            nrow = 1, align = "h", rel_widths = c(1.2, 1, 1)),
  plot_grid(plotlist = l_plots_cov1_A_abu[c(4:6)],
            nrow = 1, align = "h", rel_widths = c(1.2, 1, 1)),
  plot_grid(plotlist = c(l_plots_cov1_A_abu[c(7:8)], list(ggplot() + theme_void())),
            nrow = 1, align = "h", rel_widths = c(1.2, 1, 1)),
  rel_heights = c(1, 1.3, 1), ncol = 1
)

ggsave(file = "Output/Figures/Covariates_abu_A.jpeg",
       width = 180, height = 210, units = "mm", dpi = 600)

# ... Figure S8 ################################################################
################################################################################.

l_plots_cov1_A_ric <- c(l_ric_A$l_pred_fe, l_ric_A$l_pred_sm) %>% 
  keep(!names(.) %in% c("A", "height", "A:height")) %>%
  map(f_plot_pred, data = d_mod, response = "SCcorr_ric", line.size = .75, 
      hours_sel = which(d_mod_z$traptype == "p" & d_mod_z$hours_data))

l_plots_cov1_A_ric <- lapply(l_plots_cov1_A_ric,
                             \(x) x +
                               theme(axis.title = element_text(size = v_textsize["axis.title"]),
                                     axis.text = element_text(size = v_textsize["axis.text"])) +
                               labs(y = "Richness") +
                               scale_y_continuous(breaks = c(0, 10, 100, 1000), 
                                                  trans = "log1p"))

l_plots_cov1_A_ric <- l_plots_cov1_A_ric[c("yday", "P_2day", "T_2day",
                                           "traptype", "bulbtype", "n_trap",
                                           "active_hours", "sample_previous")]

l_plots_cov1_A_ric[c("P_2day", "T_2day",
                     "bulbtype", "n_trap",
                     "sample_previous")] <-
  lapply(l_plots_cov1_A_ric[c("P_2day", "T_2day",
                              "bulbtype", "n_trap",
                              "sample_previous")],
         \(x) x +
           theme(axis.title.y = element_blank(),
                 axis.text.y = element_blank()))

plot_grid(
  plot_grid(plotlist = l_plots_cov1_A_ric[c(1:3)],
            nrow = 1, align = "h", rel_widths = c(1.2, 1, 1)),
  plot_grid(plotlist = l_plots_cov1_A_ric[c(4:6)],
            nrow = 1, align = "h", rel_widths = c(1.2, 1, 1)),
  plot_grid(plotlist = c(l_plots_cov1_A_ric[c(7:8)], list(ggplot() + theme_void())),
            nrow = 1, align = "h", rel_widths = c(1.2, 1, 1)),
  rel_heights = c(1, 1.3, 1), ncol = 1
)

ggsave(file = "Output/Figures/Covariates_ric_A.jpeg",
       width = 180, height = 210, units = "mm", dpi = 600)

# ... Figure S9 ################################################################
################################################################################.

l_plots_cov1_A_mass <- c(l_mass_A$l_pred_fe, l_mass_A$l_pred_sm) %>% 
  keep(!names(.) %in% c("A", "height", "A:height")) %>%
  map(f_plot_pred, data = d_mod, response = "mass_tot", line.size = .75, 
      hours_sel = which(d_mod_z$traptype == "p" & d_mod_z$hours_data))

l_plots_cov1_A_mass <- lapply(l_plots_cov1_A_mass,
                              \(x) x +
                                theme(axis.title = element_text(size = v_textsize["axis.title"]),
                                      axis.text = element_text(size = v_textsize["axis.text"])) +
                                labs(y = "Biomass (g)") + 
                                scale_y_continuous(breaks = c(0, 0.01, 0.1, 1, 10, 100),
                                                   labels = c(0, 0.01, 0.1, 1, 10, 100),
                                                   trans = log_plus_trans))

l_plots_cov1_A_mass <- l_plots_cov1_A_mass[c("yday", "P_2day", "T_2day",
                                             "traptype", "bulbtype", "n_trap",
                                             "active_hours", "sample_previous")]

l_plots_cov1_A_mass[c("P_2day", "T_2day",
                      "bulbtype", "n_trap",
                      "sample_previous")] <-
  lapply(l_plots_cov1_A_mass[c("P_2day", "T_2day",
                               "bulbtype", "n_trap",
                               "sample_previous")],
         \(x) x +
           theme(axis.title.y = element_blank(),
                 axis.text.y = element_blank()))

plot_grid(
  plot_grid(plotlist = l_plots_cov1_A_mass[c(1:3)],
            nrow = 1, align = "h", rel_widths = c(1.2, 1, 1)),
  plot_grid(plotlist = l_plots_cov1_A_mass[c(4:6)],
            nrow = 1, align = "h", rel_widths = c(1.2, 1, 1)),
  plot_grid(plotlist = c(l_plots_cov1_A_mass[c(7:8)], list(ggplot() + theme_void())),
            nrow = 1, align = "h", rel_widths = c(1.2, 1, 1)),
  rel_heights = c(1, 1.3, 1), ncol = 1
)

ggsave(file = "Output/Figures/Covariates_mass_A.jpeg",
       width = 180, height = 210, units = "mm", dpi = 600)

# ... Figure S10 ###############################################################
################################################################################.

plot_grid(
  # Abundance ------------------------------------------------------------------.
  f_A_height_plot(pred = l_abu_A$l_pred_fe$`A:height`,
                  data_raw = d_mod_z,
                  mean_yday = d_scalings$mean[d_scalings$data == "full" & d_scalings$var == "yday"],
                  response = "abu_tot") +
    scale_y_continuous(breaks = c(0, 10, 100, 1000, 10000), trans = "log1p") +
    coord_cartesian(ylim = c(NA, 10000),
                    xlim = as.Date(c("1972-01-01", "2021-12-31"))) +
    labs(x = NULL, y = "Abundance") +
    theme(legend.position = "none",
          axis.text.x = element_blank()) +
    facet_grid(~ "Full model"),
  f_A_height_plot(pred = l_abu_A_lf$l_pred_fe$`A:height`,
                  data_raw = d_mod_lf_z,
                  mean_yday = d_scalings$mean[d_scalings$data == "LF" & d_scalings$var == "yday"],
                  response = "abu_tot") +
    scale_y_continuous(breaks = c(0, 10, 100, 1000, 10000), trans = "log1p") +
    coord_cartesian(ylim = c(NA, 10000),
                    xlim = as.Date(c("1972-01-01", "2021-12-31"))) +
    labs(x = NULL, y = NULL) +
    theme(legend.position = "none",
          axis.text = element_blank()) +
    facet_grid(~ "Fixed only"),
  f_A_height_plot(pred = l_abu_A_p$l_pred_fe$`A:height`,
                  data_raw = d_mod_p_z,
                  mean_yday = d_scalings$mean[d_scalings$data == "p" & d_scalings$var == "yday"],
                  response = "abu_tot") +
    scale_y_continuous(breaks = c(0, 10, 100, 1000, 10000), trans = "log1p") +
    coord_cartesian(ylim = c(NA, 10000),
                    xlim = as.Date(c("1972-01-01", "2021-12-31"))) +
    labs(x = NULL, y = NULL) +
    theme(legend.position = "none",
          axis.text = element_blank()) +
    facet_grid(~ "Manual only"),
  f_A_height_plot(pred = l_abu_A_ne$l_pred_fe$`A:height`,
                  data_raw = d_mod_ne_z,
                  mean_yday = d_scalings$mean[d_scalings$data == "no_estimates" & d_scalings$var == "yday"],
                  response = "abu_tot") +
    scale_y_continuous(breaks = c(0, 10, 100, 1000, 10000), trans = "log1p") +
    coord_cartesian(ylim = c(NA, 10000),
                    xlim = as.Date(c("1972-01-01", "2021-12-31"))) +
    labs(x = NULL, y = NULL) +
    theme(legend.position = "none",
          axis.text = element_blank()) +
    facet_grid(~ "Full data only"),
  
  # Richness -------------------------------------------------------------------.
  f_A_height_plot(pred = l_ric_A$l_pred_fe$`A:height`,
                  data_raw = d_mod_z,
                  mean_yday = d_scalings$mean[d_scalings$data == "full" & d_scalings$var == "yday"],
                  response = "SCcorr_ric") +
    scale_y_continuous(breaks = c(0, 10, 100, 1000), trans = "log1p",
                       labels = c("0", "10", "    100",  "1000")) +
    coord_cartesian(ylim = c(NA, 1200),
                    xlim = as.Date(c("1972-01-01", "2021-12-31"))) +
    labs(x = NULL, y = "Richness") +
    theme(legend.position = "none",
          axis.text.x = element_blank()),
  f_A_height_plot(pred = l_ric_A_lf$l_pred_fe$`A:height`,
                  data_raw = d_mod_lf_z,
                  mean_yday = d_scalings$mean[d_scalings$data == "LF" & d_scalings$var == "yday"],
                  response = "SCcorr_ric") +
    scale_y_continuous(breaks = c(0, 10, 100, 1000), trans = "log1p") +
    coord_cartesian(ylim = c(NA, 1200),
                    xlim = as.Date(c("1972-01-01", "2021-12-31"))) +
    labs(x = NULL, y = NULL) +
    theme(legend.position = "none",
          axis.text = element_blank()),
  f_A_height_plot(pred = l_ric_A_p$l_pred_fe$`A:height`,
                  data_raw = d_mod_p_z,
                  mean_yday = d_scalings$mean[d_scalings$data == "p" & d_scalings$var == "yday"],
                  response = "SCcorr_ric") +
    scale_y_continuous(breaks = c(0, 10, 100, 1000), trans = "log1p") +
    coord_cartesian(ylim = c(NA, 1200),
                    xlim = as.Date(c("1972-01-01", "2021-12-31"))) +
    labs(x = NULL, y = NULL) +
    theme(legend.position = "none",
          axis.text = element_blank()),
  f_A_height_plot(pred = l_ric_A_ne$l_pred_fe$`A:height`,
                  data_raw = d_mod_ne_z,
                  mean_yday = d_scalings$mean[d_scalings$data == "no_estimates" & d_scalings$var == "yday"],
                  response = "SCcorr_ric") +
    scale_y_continuous(breaks = c(0, 10, 100, 1000), trans = "log1p") +
    coord_cartesian(ylim = c(NA, 1200),
                    xlim = as.Date(c("1972-01-01", "2021-12-31"))) +
    labs(x = NULL, y = NULL) +
    theme(legend.position = "none",
          axis.text = element_blank()),
  
  # Biomass --------------------------------------------------------------------.
  f_A_height_plot(pred = l_mass_A$l_pred_fe$`A:height`,
                  data_raw = d_mod_z,
                  mean_yday = d_scalings$mean[d_scalings$data == "full" & d_scalings$var == "yday"],
                  response = "mass_tot") +
    scale_y_continuous(breaks = c(0, 0.01, 0.1, 1, 10, 100),
                       labels = c(0, 0.01, 0.1, 1, 10, "    100"),
                       trans = log_plus_trans) +
    coord_cartesian(ylim = c(NA, 500),
                    xlim = as.Date(c("1972-01-01", "2021-12-31"))) +
    labs(x = NULL, y = "Biomass (g)") +
    theme(legend.position = "none"),
  f_A_height_plot(pred = l_mass_A_lf$l_pred_fe$`A:height`,
                  data_raw = d_mod_lf_z,
                  mean_yday = d_scalings$mean[d_scalings$data == "LF" & d_scalings$var == "yday"],
                  response = "mass_tot") +
    scale_y_continuous(breaks = c(0, 0.01, 0.1, 1, 10, 100),
                       labels = c(0, 0.01, 0.1, 1, 10, "    100"),
                       trans = log_plus_trans) +
    coord_cartesian(ylim = c(NA, 500),
                    xlim = as.Date(c("1972-01-01", "2021-12-31"))) +
    labs(x = NULL, y = NULL) +
    theme(legend.position = "none",
          axis.text.y = element_blank()),
  f_A_height_plot(pred = l_mass_A_p$l_pred_fe$`A:height`,
                  data_raw = d_mod_p_z,
                  mean_yday = d_scalings$mean[d_scalings$data == "p" & d_scalings$var == "yday"],
                  response = "mass_tot") +
    scale_y_continuous(breaks = c(0, 0.01, 0.1, 1, 10, 100),
                       labels = c(0, 0.01, 0.1, 1, 10, "    100"),
                       trans = log_plus_trans) +
    coord_cartesian(ylim = c(NA, 500),
                    xlim = as.Date(c("1972-01-01", "2021-12-31"))) +
    labs(x = NULL, y = NULL) +
    theme(legend.position = "none",
          axis.text.y = element_blank()),
  f_A_height_plot(pred = l_mass_A_ne$l_pred_fe$`A:height`,
                  data_raw = d_mod_ne_z,
                  mean_yday = d_scalings$mean[d_scalings$data == "no_estimates" & d_scalings$var == "yday"],
                  response = "mass_tot") +
    scale_y_continuous(breaks = c(0, 0.01, 0.1, 1, 10, 100),
                       labels = c(0, 0.01, 0.1, 1, 10, "    100"),
                       trans = log_plus_trans) +
    coord_cartesian(ylim = c(NA, 500),
                    xlim = as.Date(c("1972-01-01", "2021-12-31"))) +
    labs(x = NULL, y = NULL) +
    theme(legend.position = "none",
          axis.text.y = element_blank()),
  
  nrow = 3, rel_widths = c(1, .75, .75, .75) )


ggsave(paste0(dat.dir, "Output/Figures/Year_trends_sensibility.jpeg"), width = 180, height = 180,
       units = "mm", dpi = 400)

# ... Figure S11 ###############################################################
################################################################################.

# abundance model coef vs. elevation -------------------------------------------.

d_slope_A_abu_height <- d_slope_A_abu |> 
  group_by(LOC) |> 
  summarise(slope_mean = mean(slope),
            .groups = "drop") |> 
  left_join(d_mod |> 
              group_by(LOC, height) |> 
              summarise(n_years = n_distinct(A),
                        range_years = max(A) - min(A) + 1,
                        .groups = "drop"),
            by = "LOC")

set.seed(41)
mod_slope_abu <- d_slope_A_abu_height |> 
  brm(slope_mean | weights(range_years) ~ height, data = _)

d_pred_slope_abu <- data.frame(height = seq(min(d_slope_A_abu_height$height),
                                            max(d_slope_A_abu_height$height),
                                            length.out = 100))

m_pred_slope_abu <- posterior_epred(mod_slope_abu,
                                    newdata = d_pred_slope_abu)

d_pred_slope_abu$pred <- apply(m_pred_slope_abu, 2, mean)
d_pred_slope_abu$lower <- apply(m_pred_slope_abu, 2, function(x) ci(x)$CI_low)
d_pred_slope_abu$upper <- apply(m_pred_slope_abu, 2, function(x) ci(x)$CI_high)

# richness model coef vs. elevation --------------------------------------------.

d_slope_A_SCcr_height <- d_slope_A_SCcr |> 
  group_by(LOC) |> 
  summarise(slope_mean = mean(slope),
            .groups = "drop") |> 
  left_join(d_mod |> 
              group_by(LOC, height) |> 
              summarise(n_years = n_distinct(A),
                        range_years = max(A) - min(A) + 1,
                        .groups = "drop"),
            by = "LOC")

set.seed(41)
mod_slope_SCcr <- d_slope_A_SCcr_height |> 
  brm(slope_mean | weights(range_years) ~ height, data = _)

d_pred_slope_SCcr <- data.frame(height = seq(min(d_slope_A_SCcr_height$height),
                                             max(d_slope_A_SCcr_height$height),
                                             length.out = 100))

m_pred_slope_SCcr <- posterior_epred(mod_slope_SCcr,
                                     newdata = d_pred_slope_SCcr)

d_pred_slope_SCcr$pred <- apply(m_pred_slope_SCcr, 2, mean)
d_pred_slope_SCcr$lower <- apply(m_pred_slope_SCcr, 2, function(x) ci(x)$CI_low)
d_pred_slope_SCcr$upper <- apply(m_pred_slope_SCcr, 2, function(x) ci(x)$CI_high)

# biomass model coef vs. elevation ---------------------------------------------.

d_slope_A_mass_height <- d_slope_A_mass |> 
  group_by(LOC) |> 
  summarise(slope_mean = mean(slope),
            .groups = "drop") |> 
  left_join(d_mod |> 
              group_by(LOC, height) |> 
              summarise(n_years = n_distinct(A),
                        range_years = max(A) - min(A) + 1,
                        .groups = "drop"),
            by = "LOC")

set.seed(41)
mod_slope_mass <- d_slope_A_mass_height |> 
  brm(slope_mean | weights(range_years) ~ height, data = _)

d_pred_slope_mass <- data.frame(height = seq(min(d_slope_A_mass_height$height),
                                             max(d_slope_A_mass_height$height),
                                             length.out = 100))

m_pred_slope_mass <- posterior_epred(mod_slope_mass,
                                     newdata = d_pred_slope_mass)

d_pred_slope_mass$pred <- apply(m_pred_slope_mass, 2, mean)
d_pred_slope_mass$lower <- apply(m_pred_slope_mass, 2, function(x) ci(x)$CI_low)
d_pred_slope_mass$upper <- apply(m_pred_slope_mass, 2, function(x) ci(x)$CI_high)

# combined plot ----------------------------------------------------------------.

d_slope_A_abu |> 
  group_by(LOC) |> 
  summarise(slope_mean = mean(slope),
            .groups = "drop") |> 
  mutate(par = "Abundance") |> 
  bind_rows(d_slope_A_SCcr |> 
              group_by(LOC) |> 
              summarise(slope_mean = mean(slope),
                        .groups = "drop") |> 
              mutate(par = "Richness")) |> 
  bind_rows(d_slope_A_mass |> 
              group_by(LOC) |> 
              summarise(slope_mean = mean(slope),
                        .groups = "drop") |> 
              mutate(par = "Biomass")) |> 
  mutate(factor_mean = exp(slope_mean * 10 / d_scalings$sd[d_scalings$data == "full" & d_scalings$var == "A"])) |> 
  left_join(d_mod |> 
              group_by(LOC, height) |> 
              summarise(n_years = n_distinct(A),
                        range_years = max(A) - min(A) + 1,
                        .groups = "drop"),
            by = "LOC") |> 
  mutate(par = factor(par, levels = c("Abundance", "Richness", "Biomass"))) |> 
  ggplot(aes(x = height, y = factor_mean)) +
  geom_hline(yintercept = 1, lty = 2) +
  geom_point(aes(size = range_years), alpha = .7) +
  geom_line(data = d_pred_slope_abu |> 
              mutate(par = "Abundance") |> 
              bind_rows(d_pred_slope_SCcr |> 
                          mutate(par = "Richness")) |> 
              bind_rows(d_pred_slope_mass |> 
                          mutate(par = "Biomass")) |> 
              mutate(par = factor(par, levels = c("Abundance", "Richness", "Biomass"))) |> 
              mutate(pred = exp(pred * 10 / d_scalings$sd[d_scalings$data == "full" & d_scalings$var == "A"])), 
            aes(y = pred), colour = "firebrick4") +
  geom_ribbon(data = d_pred_slope_abu |> 
                mutate(par = "Abundance") |> 
                bind_rows(d_pred_slope_SCcr |> 
                            mutate(par = "Richness")) |> 
                bind_rows(d_pred_slope_mass |> 
                            mutate(par = "Biomass")) |> 
                mutate(par = factor(par, levels = c("Abundance", "Richness", "Biomass"))) |> 
                mutate(across(c(lower, upper), ~ exp(. * 10 / d_scalings$sd[d_scalings$data == "full" & d_scalings$var == "A"]))),
              aes(ymin = lower, ymax = upper), y = NA, 
              fill = "firebrick", alpha = .3) +
  labs(x = "Elevation (m asl.)", y = "Year effect", size = "Years range") +
  scale_size_continuous(range = c(1, 4)) +
  scale_y_log10() +
  facet_grid(~ par) +
  theme(strip.text = element_text(size = v_textsize["axis.title"]),
        legend.title = element_text(size = v_textsize["axis.title"]),
        axis.title = element_text(size = v_textsize["axis.title"]),
        axis.text = element_text(size = v_textsize["axis.text"]),
        legend.text = element_text(size = v_textsize["axis.text"]))

ggsave(paste0(dat.dir, "Output/Figures/sitemodels.jpeg"), width = 180, height = 80,
       units = "mm", dpi = 600)


exp(fixef(mod_slope_abu) * 10 / d_scalings$sd[d_scalings$data == "full" & d_scalings$var == "A"] * 1000)
exp(fixef(mod_slope_SCcr) * 10 / d_scalings$sd[d_scalings$data == "full" & d_scalings$var == "A"] * 1000)
exp(fixef(mod_slope_mass) * 10 / d_scalings$sd[d_scalings$data == "full" & d_scalings$var == "A"] * 1000)


# ... Tables S2, S4-S18 ########################################################
################################################################################.

# full models ------------------------------------------------------------------.

wb_out <- createWorkbook()

for (resp_i in c("abu", "ric", "mass")){
  
  mod_i <- eval(parse(text = paste0("l_", resp_i, "_A")))
  
  formula_i <- response ~ s(yday) + P_2day + T_2day +
    C(traptype, "contr.sum") +
    C(bulbtype, "contr.sum") +
    n_trap +
    C(sample_previous, "contr.sum") +
    (1 | spattemp_cluster) +
    (1 | LOC) +
    (1 | night_ID) +
    (1 | trap_ID_A)
  
  
  table_i <- f_summarytable(mod_i$fit, formula_i, v_covlabels_short,
                            data = d_mod_z, model = "year")
  
  addWorksheet(wb_out, resp_i)
  
  writeDataTable(wb_out, resp_i, table_i)
  
  for (c_i in c(1, 2)){
    goon <- T
    start <- 1
    while(goon){
      
      if (n_distinct(table_i[start:nrow(table_i), c_i], na.rm = T) == 1){ # only one value left
        stop <- nrow(table_i)
      } else {
        stop <- which(table_i[start:nrow(table_i), c_i] != table_i[start, c_i])[1] + start - 2
      }
      
      
      if (is.na(stop)) break
      
      if (stop < start) stop <- start
      
      mergeCells(wb_out, resp_i, cols = c_i, rows = seq(start + 1, stop + 1))
      
      if (stop == nrow(table_i)) goon <- F
      
      start <- stop + 1
    }
  }
  
  addStyle(wb_out, sheet = resp_i, rows = seq_len(nrow(table_i) + 1), 
           cols = seq_len(ncol(table_i)), gridExpand = T,
           style = createStyle(fontSize = 9, fontName = "Helvetica", 
                               valign = 'center', border = "TopBottomLeftRight"))
  addStyle(wb_out, sheet = resp_i, rows = 1, cols = 1:100, gridExpand = T,
           style = createStyle(textDecoration = 'bold'), stack = T)
  rows_sig_i <- which(table_i$Parameter == "Fixed effect" &
                        (as.numeric(table_i$`Lower 95%-CI`) > 0 |
                           as.numeric(table_i$`Upper 95%-CI`) < 0))
  addStyle(wb_out, sheet = resp_i, rows = rows_sig_i + 1, cols = 4:6, gridExpand = T,
           style = createStyle(textDecoration = 'bold'), stack = T)
  addStyle(wb_out, sheet = resp_i, rows = 2:100, cols = 4:6, gridExpand = T,
           style = createStyle(halign = 'right'), stack = T)
  addStyle(wb_out, sheet = resp_i, rows = 1, cols = 4:6, gridExpand = T,
           style = createStyle(halign = 'center'), stack = T)
  
  setColWidths(wb_out, resp_i, seq_len(ncol(table_i)), widths = "auto")
}

saveWorkbook(wb_out, file = "Output/Tables/Modsummaries_A_height.xlsx", overwrite = T)

# full observation hours data --------------------------------------------------.

wb_out <- createWorkbook()

for (resp_i in c("abu", "ric", "mass")){
  
  mod_i <- eval(parse(text = paste0("l_", resp_i, "_A_ne")))
  
  formula_i <- response ~ s(yday) + P_2day + T_2day +
    C(traptype, "contr.sum") +
    C(bulbtype, "contr.sum") +
    n_trap +
    C(sample_previous, "contr.sum") +
    (1 | spattemp_cluster) +
    (1 | LOC) +
    (1 | night_ID) +
    (1 | trap_ID_A)
  
  
  table_i <- f_summarytable(mod_i$fit, formula_i, v_covlabels_short,
                            data = d_mod_ne_z, model = "year")
  
  addWorksheet(wb_out, resp_i)
  
  writeDataTable(wb_out, resp_i, table_i)
  
  for (c_i in c(1, 2)){
    goon <- T
    start <- 1
    while(goon){
      
      if (n_distinct(table_i[start:nrow(table_i), c_i], na.rm = T) == 1){ # only one value left
        stop <- nrow(table_i)
      } else {
        stop <- which(table_i[start:nrow(table_i), c_i] != table_i[start, c_i])[1] + start - 2
      }
      
      
      if (is.na(stop)) break
      
      if (stop < start) stop <- start
      
      mergeCells(wb_out, resp_i, cols = c_i, rows = seq(start + 1, stop + 1))
      
      if (stop == nrow(table_i)) goon <- F
      
      start <- stop + 1
    }
  }
  
  
  addStyle(wb_out, sheet = resp_i, rows = seq_len(nrow(table_i) + 1), 
           cols = seq_len(ncol(table_i)), gridExpand = T,
           style = createStyle(fontSize = 8, fontName = "Helvetica", 
                               valign = 'center', border = "TopBottomLeftRight"))
  addStyle(wb_out, sheet = resp_i, rows = 1, cols = 1:100, gridExpand = T,
           style = createStyle(textDecoration = 'bold'), stack = T)
  rows_sig_i <- which(table_i$Parameter == "Fixed effect" &
                        (as.numeric(table_i$`Lower 95%-CI`) > 0 |
                           as.numeric(table_i$`Upper 95%-CI`) < 0))
  addStyle(wb_out, sheet = resp_i, rows = rows_sig_i + 1, cols = 4:6, gridExpand = T,
           style = createStyle(textDecoration = 'bold'), stack = T)
  addStyle(wb_out, sheet = resp_i, rows = 2:100, cols = 4:6, gridExpand = T,
           style = createStyle(halign = 'right'), stack = T)
  addStyle(wb_out, sheet = resp_i, rows = 1, cols = 4:6, gridExpand = T,
           style = createStyle(halign = 'center'), stack = T)
  
  setColWidths(wb_out, resp_i, seq_len(ncol(table_i)), widths = "auto")
}

saveWorkbook(wb_out, file = "Output/Tables/Modsummaries_A_height_ne.xlsx", overwrite = T)

# fixed traps only -------------------------------------------------------------.

wb_out <- createWorkbook()

for (resp_i in c("abu", "ric", "mass")){
  
  mod_i <- eval(parse(text = paste0("l_", resp_i, "_A_lf")))
  
  
  formula_i <- response ~ s(yday) + P_2day + T_2day +
    C(traptype, "contr.sum") +
    C(bulbtype, "contr.sum") +
    C(sample_previous, "contr.sum") +
    (1 | spattemp_cluster) +
    (1 | LOC) +
    (1 | night_ID) +
    (1 | trap_ID_A)
  
  
  table_i <- f_summarytable(mod_i$fit, formula_i, v_covlabels_short,
                            data = d_mod_lf_z, model = "year")
  
  addWorksheet(wb_out, resp_i)
  
  writeDataTable(wb_out, resp_i, table_i)
  
  for (c_i in c(1, 2)){
    goon <- T
    start <- 1
    while(goon){
      
      if (n_distinct(table_i[start:nrow(table_i), c_i], na.rm = T) == 1){ # only one value left
        stop <- nrow(table_i)
      } else {
        stop <- which(table_i[start:nrow(table_i), c_i] != table_i[start, c_i])[1] + start - 2
      }
      
      
      if (is.na(stop)) break
      
      if (stop < start) stop <- start
      
      mergeCells(wb_out, resp_i, cols = c_i, rows = seq(start + 1, stop + 1))
      
      if (stop == nrow(table_i)) goon <- F
      
      start <- stop + 1
    }
  }
  
  
  addStyle(wb_out, sheet = resp_i, rows = seq_len(nrow(table_i) + 1), 
           cols = seq_len(ncol(table_i)), gridExpand = T,
           style = createStyle(fontSize = 8, fontName = "Helvetica", 
                               valign = 'center', border = "TopBottomLeftRight"))
  addStyle(wb_out, sheet = resp_i, rows = 1, cols = 1:100, gridExpand = T,
           style = createStyle(textDecoration = 'bold'), stack = T)
  rows_sig_i <- which(table_i$Parameter == "Fixed effect" &
                        (as.numeric(table_i$`Lower 95%-CI`) > 0 |
                           as.numeric(table_i$`Upper 95%-CI`) < 0))
  addStyle(wb_out, sheet = resp_i, rows = rows_sig_i + 1, cols = 4:6, gridExpand = T,
           style = createStyle(textDecoration = 'bold'), stack = T)
  addStyle(wb_out, sheet = resp_i, rows = 2:100, cols = 4:6, gridExpand = T,
           style = createStyle(halign = 'right'), stack = T)
  addStyle(wb_out, sheet = resp_i, rows = 1, cols = 4:6, gridExpand = T,
           style = createStyle(halign = 'center'), stack = T)
  
  setColWidths(wb_out, resp_i, seq_len(ncol(table_i)), widths = "auto")
}

saveWorkbook(wb_out, file = "Output/Tables/Modsummaries_A_height_LF.xlsx", overwrite = T)


# manual traps only ------------------------------------------------------------.

wb_out <- createWorkbook()

for (resp_i in c("abu", "ric", "mass")){
  
  mod_i <- eval(parse(text = paste0("l_", resp_i, "_A_p")))
  
  formula_i <- response ~ s(yday) + P_2day + T_2day + 
    C(bulbtype, "contr.sum") +
    n_trap +
    C(sample_previous, "contr.sum") + 
    (1 | spattemp_cluster) + 
    (1 | LOC) +
    (1 | night_ID) +
    (1 | trap_ID_A)
  
  table_i <- f_summarytable(mod_i$fit, formula_i, v_covlabels_short,
                            data = d_mod_p_z, model = "year")
  
  addWorksheet(wb_out, resp_i)
  
  writeDataTable(wb_out, resp_i, table_i)
  
  for (c_i in c(1, 2)){
    goon <- T
    start <- 1
    while(goon){
      
      if (n_distinct(table_i[start:nrow(table_i), c_i], na.rm = T) == 1){ # only one value left
        stop <- nrow(table_i)
      } else {
        stop <- which(table_i[start:nrow(table_i), c_i] != table_i[start, c_i])[1] + start - 2
      }
      
      
      if (is.na(stop)) break
      
      if (stop < start) stop <- start
      
      mergeCells(wb_out, resp_i, cols = c_i, rows = seq(start + 1, stop + 1))
      
      if (stop == nrow(table_i)) goon <- F
      
      start <- stop + 1
    }
  }
  
  
  addStyle(wb_out, sheet = resp_i, rows = seq_len(nrow(table_i) + 1), 
           cols = seq_len(ncol(table_i)), gridExpand = T,
           style = createStyle(fontSize = 8, fontName = "Helvetica", 
                               valign = 'center', border = "TopBottomLeftRight"))
  addStyle(wb_out, sheet = resp_i, rows = 1, cols = 1:100, gridExpand = T,
           style = createStyle(textDecoration = 'bold'), stack = T)
  rows_sig_i <- which(table_i$Parameter == "Fixed effect" &
                        (as.numeric(table_i$`Lower 95%-CI`) > 0 |
                           as.numeric(table_i$`Upper 95%-CI`) < 0))
  addStyle(wb_out, sheet = resp_i, rows = rows_sig_i + 1, cols = 4:6, gridExpand = T,
           style = createStyle(textDecoration = 'bold'), stack = T)
  addStyle(wb_out, sheet = resp_i, rows = 2:100, cols = 4:6, gridExpand = T,
           style = createStyle(halign = 'right'), stack = T)
  addStyle(wb_out, sheet = resp_i, rows = 1, cols = 4:6, gridExpand = T,
           style = createStyle(halign = 'center'), stack = T)
  
  setColWidths(wb_out, resp_i, seq_len(ncol(table_i)), widths = "auto")
}

saveWorkbook(wb_out, file = "Output/Tables/Modsummaries_A_height_p.xlsx", overwrite = T)

# ------------------------------------------------------------------------------.
# trait subsets ----------------------------------------------------------------.
# ------------------------------------------------------------------------------.

d_traits <- data.frame(trait = "mass", traitvalue = "small") |> 
  add_row(trait = "mass", traitvalue = "medium") |> 
  add_row(trait = "mass", traitvalue = "large") |> 
  add_row(trait = "Tavg", traitvalue = "cold") |> 
  add_row(trait = "Tavg", traitvalue = "intermediate") |> 
  add_row(trait = "Tavg", traitvalue = "warm") |> 
  add_row(trait = "spec", traitvalue = "m") |> 
  add_row(trait = "spec", traitvalue = "o") |> 
  add_row(trait = "spec", traitvalue = "p") |> 
  add_row(trait = "hib", traitvalue = "egg") |> 
  add_row(trait = "hib", traitvalue = "larva") |> 
  add_row(trait = "hib", traitvalue = "pupa") |> 
  add_row(trait = "hib", traitvalue = "adult")

c_traits <- c(mass = "mass_cat",
              Tavg = "Tavg_mean_cat",
              spec = "Spec",
              hib = "overwintering_stage")
c_traitvalues <- c(m = "Monophagous",
                   o = "Oligophagous",
                   p = "Polyphagous")

for (trait_i in unique(d_traits$trait)){
  wb_out <- createWorkbook()
  for (resp_i in c("abu", "ric", "mass")){
    for (traitvalue_i in d_traits$traitvalue[d_traits$trait == trait_i]){
      
      if (trait_i == "spec"){
        mod_i <- eval(parse(text = paste0("l_", resp_i, "_A_spec_", traitvalue_i)))
      } else {
        mod_i <- eval(parse(text = paste0("l_", resp_i, "_A_", traitvalue_i)))
      }
      
      formula_i <- response ~ s(yday) + P_2day + T_2day +
        C(traptype, "contr.sum") +
        C(bulbtype, "contr.sum") +
        n_trap +
        C(sample_previous, "contr.sum") +
        (1 | spattemp_cluster) +
        (1 | LOC) +
        (1 | night_ID) +
        (1 | trap_ID_A)
      
      data <- eval(parse(text = paste0("d_mod_", trait_i, "_z")))
      
      if (trait_i == "spec"){
        data <- data |> 
          filter(!! sym(c_traits[trait_i]) == c_traitvalues[traitvalue_i])
      } else {
        data <- data |> 
          filter(!! sym(c_traits[trait_i]) == traitvalue_i)
      }
      
      
      table_i <- f_summarytable(mod_i$fit, formula_i, v_covlabels_short,
                                data = data, model = "year")
      
      addWorksheet(wb_out, paste0(resp_i, "_", traitvalue_i))
      
      writeDataTable(wb_out, paste0(resp_i, "_", traitvalue_i), table_i)
      
      for (c_i in c(1, 2)){
        goon <- T
        start <- 1
        while(goon){
          
          if (n_distinct(table_i[start:nrow(table_i), c_i], na.rm = T) == 1){ # only one value left
            stop <- nrow(table_i)
          } else {
            stop <- which(table_i[start:nrow(table_i), c_i] != table_i[start, c_i])[1] + start - 2
          }
          
          
          if (is.na(stop)) break
          
          if (stop < start) stop <- start
          
          mergeCells(wb_out, paste0(resp_i, "_", traitvalue_i), cols = c_i, rows = seq(start + 1, stop + 1))
          
          if (stop == nrow(table_i)) goon <- F
          
          start <- stop + 1
        }
      }
      
      
      addStyle(wb_out, sheet = paste0(resp_i, "_", traitvalue_i), rows = seq_len(nrow(table_i) + 1), 
               cols = seq_len(ncol(table_i)), gridExpand = T,
               style = createStyle(fontSize = 8, fontName = "Helvetica", 
                                   valign = 'center', border = "TopBottomLeftRight"))
      addStyle(wb_out, sheet = paste0(resp_i, "_", traitvalue_i), rows = 1, cols = 1:100, gridExpand = T,
               style = createStyle(textDecoration = 'bold'), stack = T)
      rows_sig_i <- which(table_i$Parameter == "Fixed effect" &
                            (as.numeric(table_i$`Lower 95%-CI`) > 0 |
                               as.numeric(table_i$`Upper 95%-CI`) < 0))
      addStyle(wb_out, sheet = paste0(resp_i, "_", traitvalue_i), rows = rows_sig_i + 1, cols = 4:6, gridExpand = T,
               style = createStyle(textDecoration = 'bold'), stack = T)
      addStyle(wb_out, sheet = paste0(resp_i, "_", traitvalue_i), rows = 2:100, cols = 4:6, gridExpand = T,
               style = createStyle(halign = 'right'), stack = T)
      addStyle(wb_out, sheet = paste0(resp_i, "_", traitvalue_i), rows = 1, cols = 4:6, gridExpand = T,
               style = createStyle(halign = 'center'), stack = T)
      
      setColWidths(wb_out, paste0(resp_i, "_", traitvalue_i), seq_len(ncol(table_i)), widths = "auto")
    }
  }
  
  saveWorkbook(wb_out, file = paste0("Output/Tables/Modsummaries_A_height_", trait_i, ".xlsx"), overwrite = T)
}
