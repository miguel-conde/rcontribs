#
# Idea sacada de : https://bradleyboehmke.github.io/HOML/iml.html


# library(dplyr)
# library(ranger)
# library(ggplot2)
# library(patchwork)
# library(pdp)

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
