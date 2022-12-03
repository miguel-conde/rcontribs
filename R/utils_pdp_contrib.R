#
# Idea sacada de : https://bradleyboehmke.github.io/HOML/iml.html


# LIBRARIES and SOURCES ---------------------------------------------------

# library(dplyr)
# library(ranger)
# library(ggplot2)
# library(patchwork)
# library(pdp)

get_mode <- function(v, method = c("first", "all")) {

  method <- match.arg(method)

  uniqv <- unique(v)
  tab <- tabulate(match(v, uniqv))

  if (method == "first") out <- uniqv[which.max(tab)]
  else out <- uniqv[tab == max(tab)]

  return(out)
}

ref_var <- function(in_data, pred_vars) {

  # Pone a 0 las pred_vars numéricas y al nivel de referencia las pred_vars
  # factor o character

  out <- in_data %>%
    mutate_at(vars(pred_vars),
              ~ ifelse(is.numeric(.),
                       0,
                       ifelse(is.factor(.),
                              levels(.)[1],
                              ifelse(is.character(.),
                                     levels(factor(.))[1],
                                     .))))

  return(out)
}

avg_mode_var <- function(in_data, pred_var) {

  # Si la columna pred_var es numeric, la cambia por su media
  # Si la columna pred_var es factor o character, la cambia por su moda

  if (is.numeric(in_data[[pred_var]])) {
    avg_x <- in_data %>% pull(!!sym(pred_var)) %>% mean()
    out   <- in_data %>% mutate(!!sym(pred_var) := avg_x)
  } else {
    if (is.factor(in_data[[pred_var]]) | is.character(in_data[[pred_var]])) {
      mode_x <- in_data %>% pull(!!sym(pred_var)) %>% get_mode()
      out    <- in_data %>% mutate(!!sym(pred_var) := mode_x)
    } else
    {
      out <- in_data
    }
  }
  return(out)
}


# AVERAGE PREDICTIONS FOR PDP ---------------------------------------------


pdp_pred_lm <- function(object, newdata)  {
  results <- as.vector(predict(object, newdata))
  return(results)
}

pdp_pred_glm <- function(object, newdata)  {
  results <- as.vector(predict(object, newdata, type = "response"))
  return(results)
}

pdp_pred_rf <- function(object, newdata)  {
  results <- as.vector(predict(object, newdata)$predictions)
  return(results)
}

pdp_pred_xgboost <- function(object, newdata)  {
  results <- as.vector(predict(object, newdata))
  return(results)
}

pdp_pred <- function(object, newdata) {
  results <- as.vector(FUN_PRED(object, newdata))
  return(results)
}

# CONTRIBUTIONS -----------------------------------------------------------

pdp_1_contrib <- function(in_model, in_data, pred_var, pred_fun, grid_resolution = 20) {

  pd_values <- pdp::partial(
    in_model,
    train = in_data,
    pred.var = pred_var,
    pred.fun = pred_fun,
    grid.resolution =  grid_resolution,
    ice = TRUE
  )

  aux_class <- class(pd_values)

  pd_values <- pd_values %>%
    as_tibble() %>%
    group_by(!!sym(pred_var)) %>%
    summarise(yhat = mean(yhat))

  class(pd_values) <- aux_class

  avg_y_hat <- pred_fun(in_model, in_data) %>% mean()
  avg_x <- in_data %>% pull(!!sym(pred_var)) %>% mean()

  X_aux_0 <- in_data %>% mutate_all(~ 0)
  # X_aux_0 <- in_data %>% mutate_if(is.numeric, ~ 0) %>% mutate_if(is.factor, ~ levels(.)[1])
  X_aux_1 <- X_aux_0  %>% mutate(!!sym(pred_var) := avg_x)

  # La distancia es la media de:
  #          fitted - fitted cuando todos valen 0 menos el de interés
  dist <- avg_y_hat - (mean(pred_fun(in_model, X_aux_1)) - mean(pred_fun(in_model, X_aux_0)))

  # aux_class <- class(pd_values)
  # out <- pd_values %>%
  #   as_tibble() %>%
  #   mutate(yhat = yhat - dist) %>%
  #   as.data.frame()
  # class(out) <- aux_class

  out <- pd_values
  out$yhat <- out$yhat - dist

  return(out)
}

pdp_contribs_old <- function(in_model,
                         in_data,
                         pred_vars,
                         pred_fun,
                         grid_resolution = 20) {

  # Lista con las contribuciones de cada variable en su grid
  contrib_grid <- vector(mode = "list", length = length(pred_vars))

  # Lista con las funciones para interpolar la contribucion de cada variable
  # en su grid
  contrib_funs <- vector(mode = "list", length = length(pred_vars))

  # Lista con las contribuciones de cada variable en in_data
  tbl_contribs <- vector(mode = "list", length = length(pred_vars))

  names(contrib_grid) <- names(contrib_funs) <- names(tbl_contribs) <- pred_vars

  for (i in seq_along(pred_vars)) {

    var_i <- pred_vars[i]

    # Contribucion de var_i en su grid
    contrib_grid[[var_i]] <- pdp_1_contrib(in_model, in_data,var_i,
                                           pred_fun, grid_resolution)

    # Función para interpolar la contribución de var_i en su grid
    contrib_funs[[var_i]] <- approxfun(contrib_grid[[var_i]])

    # Tibble de 1 columna con la contribución de var_i en in_data
    contribs_i <- contrib_funs[[var_i]](in_data %>% pull(!!sym(var_i)))
    tbl_contribs[[var_i]] <- tibble(!!sym(var_i) := contribs_i)
  }

  tbl_contribs <- bind_cols(tbl_contribs) %>%
    mutate(baseline = pred_fun(in_model,
                               in_data %>% mutate_all(~ 0)), .before = 1)

  # Regularización para que las contribuciones sumen lo mismo que las
  # predicciones
  # k <- mean(pred_fun(in_model, in_data)) * nrow(in_data) / sum(tbl_contribs)
  k <- pred_fun(in_model, in_data) / apply(tbl_contribs, 1, sum)

  tbl_contribs <- (tbl_contribs * k) %>% as_tibble()

  X_aux_0 <- in_data %>% mutate_all(~ 0)
  tbl_contribs <- tbl_contribs %>%
    rename(y_avg = baseline) %>%
    mutate(baseline = pred_fun(in_model, X_aux_0), .before = 1) %>%
    mutate(a_repartir = y_avg - baseline, .after = 1) %>%
    rowwise() %>%
    mutate(s = sum(c_across(4:ncol(.)))) %>%
    ungroup() %>%
    mutate_at(vars(4:ncol(.)), ~ . + . * a_repartir / s) %>%
    select(-a_repartir, -y_avg, -s)

  out <- list(contribs     = tbl_contribs,
              contrib_grid = contrib_grid,
              contrib_funs = contrib_funs)

  return(out)
}

optim_smooth_contrib <- function(contribs, in_data, seed = NULL) {
  print("Smoothing contribs...")
  the_vars <- setdiff(names(contribs$contribs), "baseline")

  out <- vector(mode = "list", length = length(the_vars))

  for (i in seq_along(out)) {
    print(the_vars[[i]])
    # probe <- in_data %>% select(the_vars[i]) %>%
    #   mutate(yhat = contribs$contribs %>% pull(the_vars[i]))
    probe <- contribs$contrib_grid[[the_vars[[i]]]]
    out[[i]] <- fit_contrib_curve(probe, the_vars[i], seed)
  }

  return(out)

}

#' Title
#'
#' @param in_model
#' @param in_data
#' @param pred_vars
#' @param pred_fun
#' @param grid_resolution
#' @param seed
#'
#' @return
#' @export
#' @import dplyr
#' @import pdp
#'
#' @examples
pdp_contribs <- function(in_model,
                         in_data,
                         pred_vars,
                         pred_fun,
                         grid_resolution = 20,
                         seed = NULL) {

  # Lista con las contribuciones de cada variable en su grid
  contrib_grid <- vector(mode = "list", length = length(pred_vars))

  # Lista con las funciones para interpolar la contribucion de cada variable
  # en su grid
  contrib_funs <- vector(mode = "list", length = length(pred_vars))

  # Lista con las contribuciones de cada variable en in_data
  tbl_contribs <- vector(mode = "list", length = length(pred_vars))

  names(contrib_grid) <- names(contrib_funs) <- names(tbl_contribs) <- pred_vars

  for (i in seq_along(pred_vars)) {

    var_i <- pred_vars[i]

    # Contribucion de var_i en su grid
    contrib_grid[[var_i]] <- pdp_1_contrib(in_model, in_data,var_i,
                                           pred_fun, grid_resolution)

    # Función para interpolar la contribución de var_i en su grid
    contrib_funs[[var_i]] <- approxfun(contrib_grid[[var_i]])

    # Tibble de 1 columna con la contribución de var_i en in_data
    contribs_i <- contrib_funs[[var_i]](in_data %>% pull(!!sym(var_i)))
    tbl_contribs[[var_i]] <- tibble(!!sym(var_i) := contribs_i)
  }

  tbl_contribs <- bind_cols(tbl_contribs) %>%
    mutate(baseline = pred_fun(in_model,
                               in_data %>% mutate_all(~ 0)), .before = 1)
  # tbl_contribs <- bind_cols(tbl_contribs) %>%
  #   mutate(baseline = pred_fun(in_model,
  #                              in_data %>% mutate_if(is.numeric, ~ 0) %>%
  #                                mutate_if(is.factor, ~ levels(.)[1])),
  #          .before = 1)

  # Regularización para que las contribuciones sumen lo mismo que las
  # predicciones
  # k <- mean(pred_fun(in_model, in_data)) * nrow(in_data) / sum(tbl_contribs)
  k <- pred_fun(in_model, in_data) / apply(tbl_contribs, 1, sum)

  # tbl_contribs <- (tbl_contribs * k) %>% as_tibble()
  tbl_contribs <- mutate_all(as_tibble(tbl_contribs), ~. * k)

  X_aux_0 <- in_data %>% mutate_all(~ 0)
  # X_aux_0 <- in_data %>% mutate_if(is.numeric, ~ 0) %>%
  #   mutate_if(is.factor, ~ levels(.)[1])
  tbl_contribs <- tbl_contribs %>%
    rename(y_avg = baseline) %>%
    mutate(baseline = pred_fun(in_model, X_aux_0), .before = 1) %>%
    mutate(a_repartir = y_avg - baseline, .after = 1) %>%
    rowwise() %>%
    mutate(s = sum(c_across(4:ncol(.)))) %>%
    ungroup() %>%
    mutate_at(vars(4:ncol(.)), ~ . + . * a_repartir / s) %>%
    select(-a_repartir, -y_avg, -s)

  # Regularizamos también las grids de pdps acorde con lo anterior
  # Y sus funciones de aproximacion
  for (i in seq_along(pred_vars)) {

    var_i <- pred_vars[i]

    delta <- median((tbl_contribs[[var_i]] - contrib_funs[[var_i]](in_data[[var_i]])))

    contrib_grid[[var_i]]$yhat  <- contrib_grid[[var_i]]$yhat + delta

    contrib_funs[[var_i]] <- approxfun(contrib_grid[[var_i]])
  }


  out <- list(contribs     = tbl_contribs,
              smooth_contribs = NULL,
              contrib_grid = contrib_grid,
              contrib_funs = contrib_funs)

  out$smooth_contribs <- optim_smooth_contrib(out, in_data, seed)
  names(out$smooth_contribs) <- setdiff(names(out$contribs), "baseline")

  return(out)
}


# Plots -------------------------------------------------------------------


ggplot_1_contrib_old <- function(res_all, tgt_var,
                             title = NULL,
                             x_units = "", y_units = "",
                             n_x = 100,
                             the_theme = theme_bw) {

  enquo_tgt_var <- enquo(tgt_var)
  name_tgt_var  <- quo_name(enquo_tgt_var)

  # the_x_data    <- pull(in_data, !!enquo_tgt_var)
  the_x_data    <- pull(res_all$contrib_grid[[name_tgt_var]], !!enquo_tgt_var)

  data_tbl <- tibble(x = seq(from       = min(the_x_data),
                             to         = max(the_x_data),
                             length.out = n_x),
                     y = res_all$contrib_funs[[name_tgt_var]](x))

  the_title <-  ifelse(is.null(title),
                       paste(name_tgt_var, ""),
                       title)
  x_lab <- paste(name_tgt_var,   x_units)
  y_lab <- paste("", y_units)

  p <- ggplot(data = data_tbl, mapping = aes(x, y)) +
    geom_line() +
    geom_hline(yintercept = 0, linetype = 2) +
    labs(title = the_title, x = x_lab, y = y_lab) +
    the_theme()

  return(p)

}

#' Title
#'
#' @param res_all
#' @param tgt_var
#' @param in_data
#' @param pps
#' @param smooth
#' @param title
#' @param x_units
#' @param y_units
#' @param n_x
#' @param the_theme
#'
#' @return
#' @export
#'
#' @examples
ggplot_1_contrib <- function(res_all, tgt_var, in_data,
                             pps = FALSE, smooth = FALSE,
                             title = NULL,
                             x_units = "", y_units = "",
                             n_x = 100,
                             the_theme = theme_bw) {

  enquo_tgt_var <- enquo(tgt_var)
  name_tgt_var  <- quo_name(enquo_tgt_var)

  # the_x_data    <- pull(in_data, !!enquo_tgt_var)
  the_x_data    <- pull(res_all$contrib_grid[[name_tgt_var]], !!enquo_tgt_var)

  data_pdp <- tibble(type = "PDP",
                     x = seq(from       = min(the_x_data),
                             to         = max(the_x_data),
                             length.out = n_x),
                     y = res_all$contrib_funs[[name_tgt_var]](x))

  the_title <-  ifelse(is.null(title),
                       paste(name_tgt_var, ""),
                       title)
  x_lab <- paste(name_tgt_var,   x_units)
  y_lab <- paste("", y_units)

  p <- ggplot(data = data_pdp, mapping = aes(x, y, color = type)) +
    geom_line(size = 2) +
    geom_hline(yintercept = 0, linetype = 2) +
    labs(title = the_title, x = x_lab, y = y_lab,  color = "Contributions") +
    the_theme()

  if (pps == TRUE) {
    data_pps <- tibble(type = "Point Contribs",
                       x = in_data[[name_tgt_var]],
                       y = res_all$contribs[[name_tgt_var]])
                       # y = fix_outliers(res_all$contribs[[name_tgt_var]]))

    p <- p + geom_point(data = data_pps, alpha = .5)
  }

  if (smooth == TRUE) {
    p <- p +
      geom_function(fun = res_all$smooth_contribs[[name_tgt_var]]$fit_fun,
                    args = list(res_all$smooth_contribs[[name_tgt_var]]$pars),
                    size = 1.25, colour = "blue")
  }

  return(p)

}

#' Title
#'
#' @param res_all
#' @param tgt_vars
#' @param y_units
#' @param ...
#'
#' @return
#' @export
#' @import ggplot2
#' @import patchwork
#'
#' @examples
ggplot_contribs <- function(res_all, tgt_vars = NULL, y_units = "", ...) {

  tgt_vars = if(is.null(tgt_vars)) {
    names(res_all$contrib_grid)
  } else {
    tgt_vars
  }

  p_list <- vector(mode = "list", length = length(tgt_vars))
  names(p_list) <- tgt_vars
  for (i in seq_along(p_list)) {
    p_list[[i]] <- ggplot_1_contrib(res_all, !!sym(tgt_vars[[i]]),
                                    y_units = y_units, ...)
  }

  out <- p_list[[1]]
  for (i in 2:length(p_list)) out <- out + p_list[[i]]

  out <- out + guide_area() + plot_layout(guides = "collect")

  return(out)
}


# ROIs --------------------------------------------------------------------

pdp_1_roi_old <- function(in_contribs, in_model, in_data, pred_var, pred_fun, delta = 0.01) {

  k <- 1 + delta

  out <- pdp_1_contrib(in_model, in_data,
                       pred_var = pred_var, pred_fun = pred_fun) %>%
    as_tibble() %>%
    mutate(pred_var_plus = k * !!sym(pred_var),
           y_hat_plus1 = in_contribs$contrib_funs[[pred_var]](pred_var_plus)) %>%
    drop_na() %>%
    mutate(roi = (y_hat_plus1 - yhat) / delta / !!sym(pred_var)) %>%
    mutate(roi = ifelse(is.na(roi), 0, roi)) %>%
    mutate(roi = ifelse(is.infinite(roi), 0, roi)) %>%
    select(-pred_var_plus, -y_hat_plus1, -yhat)

  return(out)
}

#' Title
#'
#' @param x
#' @param in_contribs
#' @param pred_var
#'
#' @return
#' @export
#'
#' @examples
pdp_1_roi <- function(x, in_contribs, pred_var) {

  curve_pars <- as.list(in_contribs$smooth_contribs[[pred_var]]$pars)

  curve_fun <- in_contribs$smooth_contribs[[pred_var]]$roi_fun

  curve_fun(x, curve_pars)
}

pdp_rois_old <- function(in_contribs,
                     in_model,
                     in_data,
                     pred_vars,
                     pred_fun,
                     delta = .01) {

  # Lista con los ROIs de cada variable en su grid
  rois_grid <- vector(mode = "list", length = length(pred_vars))

  # Lista con las funciones para interpolar los ROIs de cada variable
  # en su grid
  rois_funs <- vector(mode = "list", length = length(pred_vars))

  names(rois_grid) <- names(rois_funs) <- pred_vars

  for (i in seq_along(pred_vars)) {

    var_i <- pred_vars[i]

    # Contribucion de var_i en su grid
    rois_grid[[var_i]] <- pdp_1_roi(in_contribs, in_model, in_data, var_i,
                                    pred_fun, delta)

    # Función para interpolar la contribución de var_i en su grid
    rois_funs[[var_i]] <- approxfun(rois_grid[[var_i]])
  }

  out <- list(rois_grid = rois_grid,
              rois_funs = rois_funs)

  return(out)
}

#' Title
#'
#' @param in_contribs
#' @param in_data
#' @param pred_vars
#'
#' @return
#' @export
#'
#' @examples
pdp_rois <- function(in_contribs,
                     in_data,
                     pred_vars = NULL) {

  if (is.null(pred_vars)) {
    pred_vars <- in_contribs$contribs %>% names() %>% setdiff("baseline")
  }

  rois_data <- vector(mode = "list", length = length(pred_vars))

  names(rois_data) <-  pred_vars

  for (i in seq_along(pred_vars)) {

    var_i <- pred_vars[i]

    rois_data[[var_i]] <- pdp_1_roi(in_data[[var_i]], in_contribs, var_i)
  }

  out <- bind_rows(rois_data)

  return(out)
}


# Plots -------------------------------------------------------------------

ggplot_1_roi_old <- function(rois_all, tgt_var,
                         title = NULL,
                         x_units = "", y_units = "",
                         n_x = 100,
                         the_theme = theme_bw) {

  enquo_tgt_var <- enquo(tgt_var)
  name_tgt_var  <- quo_name(enquo_tgt_var)

  # the_x_data    <- pull(in_data, !!enquo_tgt_var)
  the_x_data    <- pull(rois_all$rois_grid[[name_tgt_var]], !!enquo_tgt_var)

  data_tbl <- tibble(x = seq(from       = min(the_x_data),
                             to         = max(the_x_data),
                             length.out = n_x),
                     y = rois_all$rois_funs[[name_tgt_var]](x))

  the_title <-  ifelse(is.null(title),
                       paste(name_tgt_var, ""),
                       title)
  x_lab <- paste(name_tgt_var,   x_units)
  y_lab <- paste("", y_units)

  p <- ggplot(data = data_tbl, mapping = aes(x, y)) +
    geom_line() +
    geom_hline(yintercept = 0, linetype = 2) +
    labs(title = the_title, x = x_lab, y = y_lab) +
    the_theme()

  return(p)

}

#' Title
#'
#' @param in_contribs
#' @param tgt_var
#' @param in_data
#' @param x_min
#' @param x_max
#' @param title
#' @param x_units
#' @param y_units
#' @param n_x
#' @param the_theme
#'
#' @return
#' @export
#'
#' @examples
ggplot_1_roi <- function(in_contribs, tgt_var,
                         in_data = NULL,
                         x_min = ifelse(!is.null(in_data), min(in_data[[tgt_var]]), x_min),
                         x_max = ifelse(!is.null(in_data), max(in_data[[tgt_var]]), x_max),
                         title = NULL,
                         x_units = "", y_units = "",
                         n_x = 100,
                         the_theme = theme_bw) {

  enquo_tgt_var <- enquo(tgt_var)
  name_tgt_var  <- quo_name(enquo_tgt_var)

  if (!is.null(in_data))  the_x_data <- pull(in_data, !!enquo_tgt_var)
  else  the_x_data <- seq(x_min, x_max, length.out = n_x)

  data_tbl <- tibble(x = seq(from       = min(the_x_data),
                             to         = max(the_x_data),
                             length.out = n_x),
                     y = pdp_1_roi(x, in_contribs, name_tgt_var))

  the_title <-  ifelse(is.null(title),
                       paste(name_tgt_var, ""),
                       title)
  x_lab <- paste(name_tgt_var,   x_units)
  y_lab <- paste("", y_units)

  p <- ggplot(data = data_tbl, mapping = aes(x, y)) +
    geom_line() +
    geom_hline(yintercept = 0, linetype = 2) +
    labs(title = the_title, x = x_lab, y = y_lab) +
    the_theme()

  return(p)

}

ggplot_rois_old <- function(rois_all, tgt_vars = NULL, y_units = "", ...) {

  tgt_vars = if(is.null(tgt_vars)) {
    names(rois_all$rois_grid)
  } else {
    tgt_vars
  }

  p_list <- vector(mode = "list", length = length(tgt_vars))
  names(p_list) <- tgt_vars
  for (i in seq_along(p_list)) {
    p_list[[i]] <- ggplot_1_roi(rois_all, !!sym(tgt_vars[[i]]),
                                y_units = y_units, ...)
  }

  out <- p_list[[1]]
  for (i in 2:length(p_list)) out <- out + p_list[[i]]

  return(out)
}

#' Title
#'
#' @param in_contribs
#' @param in_data
#' @param tgt_vars
#' @param y_units
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
ggplot_rois <- function(in_contribs, in_data, tgt_vars = NULL, y_units = "", ...) {

  tgt_vars = if(is.null(tgt_vars)) {
    names(in_contribs$contribs %>% select(-baseline))
  } else {
    tgt_vars
  }

  p_list <- vector(mode = "list", length = length(tgt_vars))
  names(p_list) <- tgt_vars
  for (i in seq_along(p_list)) {
    p_list[[i]] <- ggplot_1_roi(in_contribs, !!sym(tgt_vars[[i]]), in_data,
                                y_units = y_units, ...)
  }

  out <- p_list[[1]]
  for (i in 2:length(p_list)) out <- out + p_list[[i]]

  return(out)
}



# OPTIM CONTRIB CURVES ----------------------------------------------------



# Sigmoid -----------------------------------------------------------------

sigmoid_fun <- function(x, pars) {
  pars[[4]] + pars[[1]] / (1 + exp(-pars[[2]]*(x-pars[[3]])))
}

d_sigmoid_dx_fun <- function(x, pars) {

  the_exp <- exp(pars[[2]]*(x+pars[[3]]))

  pars[[2]] * pars[[1]] * the_exp / (exp(pars[[2]] * pars[[3]]) + exp(pars[[2]]*x))^2
}

estimate_0_sigmoid <- function(contribs, tgt_var) {

  aux <- contribs

  y_inf <- aux[aux[[1]] == max(aux[[1]]), ]$yhat
  y_0   <- mean(aux$yhat[which(aux[[1]] == min(abs(aux[[1]])))])
  y_C   <- (max(aux$yhat) - min(aux$yhat)) / 2
  C_0   <- (max(aux[[1]]) + min(aux[[1]])) / 2
  A_0   <- 2*(y_inf - y_C)
  D_0   <- 2*y_C - y_inf
  k     <- (y_0 - D_0) / A_0 - 1
  B_0   <- ifelse(k > 0, log(k) / C_0, 0)

  out <- list(A_0, B_0, C_0, D_0)

  return(out)
}

# Curve S ---------------------------------------------------------------

curve_s_fun <- function(x, pars) {
  pars[[3]] + exp(pars[[1]] - pars[[2]]/x)
}

# manipulate(curve(curve_s_fun(x, A=A, B=B, C=C), -0, 100, ylim = c(C, C+exp(A))),
#            A = slider(-0, 20, initial = 0),
#            B = slider(-0, 100, initial = 0),
#            C = slider(-0, 20, initial = 0))

curve_ds_dx_fun <- function(x, pars) {
  pars[[2]] / x^2 * exp(pars[[1]] - pars[[2]]/x)
}

# manipulate(curve(curve_ds_dx_fun(x, A=A, B=B, C=C), -0, 100, ylim = c(C, C+exp(A))),
#            A = slider(-0, 20, initial = 0),
#            B = slider(-0, 100, initial = 0),
#            C = slider(-0, 20, initial = 0))

estimate_0_curve_s <- function(contribs, tgt_var) {

  aux <- contribs

  y_inf <- aux[aux[[1]] == max(aux[[1]]), ]$yhat
  y_0   <- mean(aux$yhat[which(aux[[1]] == min(abs(aux[[1]])))])
  B_0   <- max(aux[[1]]) + min(aux[[1]])
  C_0   <- y_0
  k     <- y_inf - C_0
  A_0   <- ifelse( k > 0, log(k), 0)

  out <- c(A_0, B_0, C_0)

  return(out)
}

# Curve tanh ------------------------------------------------------------

curve_tanh_fun <- function(x, pars) {
  pars[[3]] + pars[[1]] * tanh((x - pars[[4]]) / pars[[2]])
}

# manipulate(curve(curve_tanh_fun(x, A=A, B=B, C=C), -10, 10),
#            A = slider(-10,10), B=slider(-10, 10), C= slider(-10,10))

curve_dtanh_dx_fun <- function(x, pars) {
  require(pracma)

  pars[[1]] / pars[[2]] * sech((x - pars[[4]]) / pars[[2]])^2
}

estimate_0_curve_tanh <- function(contribs, tgt_var) {

  aux <- contribs

  y_inf <- aux[aux[[1]] == max(aux[[1]]), ]$yhat
  y_0   <- mean(aux$yhat[which(aux[[1]] == min(abs(aux[[1]])))])
  x_0   <- (max(aux[[1]]) + min(aux[[1]])) / 2

  C_0 <- y_0
  A_0 <- y_inf - C_0
  B_0 <- rnorm(1, 1e-6)

  out <- c(A_0, B_0, C_0, x_0)

  return(out)
}

# Loss Fun ----------------------------------------------------------------

loss_fun <- function(pars, contribs, FUN) {
  out <- contribs %>%
    as_tibble() %>%
    mutate(sigmoid = FUN(contribs[[1]], pars)) %>%
    mutate(d = .[[2]] - sigmoid) %>%
    summarise(sse = sum(d^2)) %>%
    as.numeric()

  return(out)
}

# Optim Fun ---------------------------------------------------------------

set_best_curve <- function(type, res_optim, fit_fun, roi_fun) {
  list(curve_type = type,
       loss       = res_optim$value,
       pars       = res_optim$par,
       fit_fun    = fit_fun,
       res_optim  = res_optim,
       roi_fun    = roi_fun)
}

fit_contrib_curve <- function(contribs, tgt_var, seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  best <- NULL

  res_optim_sigmoid <- tryCatch(
    {
      optim(par = estimate_0_sigmoid(contribs, tgt_var),
            fn = loss_fun,
            contribs = contribs,
            FUN = sigmoid_fun,
            method = "L-BFGS-B")
    },
    error = function(cond) {
      warning(paste(tgt_var, ": Sigmoid failed..."))
      return(NULL)
    }
    )

  if (!is.null(res_optim_sigmoid)) {
  best <- set_best_curve(type      = "sigmoid",
                         res_optim = res_optim_sigmoid,
                         fit_fun   = sigmoid_fun,
                         roi_fun   = d_sigmoid_dx_fun)
  }

  res_optim_s <- tryCatch(
    {
      optim(par = estimate_0_curve_s(contribs, tgt_var),
                       fn = loss_fun,
                       contribs = contribs,
                       FUN = curve_s_fun,
                       method = "L-BFGS-B")
    },
    error = function(cond) {
      warning(paste(tgt_var, ": S Curve failed..."))
      return(NULL)
    }
  )

  if (!is.null(res_optim_s)) {
    if (res_optim_s$value < best$loss) {
      best <- set_best_curve(type      = "s",
                             res_optim = res_optim_s,
                             fit_fun   = curve_s_fun,
                             roi_fun   = curve_ds_dx_fun)
    }
  }

  res_optim_tanh <- tryCatch(
    {
      optim(par = estimate_0_curve_tanh(contribs, tgt_var),
                          fn = loss_fun,
                          contribs = contribs,
                          FUN = curve_tanh_fun,
                          method = "L-BFGS-B")
    },
    error = function(cond) {
      warning(paste(tgt_var, ": TANH Curve failed..."))
      return(NULL)
    }
  )

  if (!is.null(res_optim_tanh)) {
    if (res_optim_tanh$value < best$loss) {
      best <- set_best_curve(type      = "tanh",
                             res_optim = res_optim_tanh,
                             fit_fun   = curve_tanh_fun,
                             roi_fun   = curve_dtanh_dx_fun)
    }
  }

  return(best)
}

fix_outliers <- function(x) {

  bpx <- boxplot(x, plot = FALSE)

  x[x > max(bpx$stats)] <- max(bpx$stats)
  x[x < min(bpx$stats)] <- min(bpx$stats)

  x

}
