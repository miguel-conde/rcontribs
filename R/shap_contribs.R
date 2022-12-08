# Idea sacada de https://towardsdatascience.com/explainable-ai-application-of-shapely-values-in-marketing-analytics-57b716fc9d1f

# library(fastshap)


#' Title
#'
#' @param in_model
#' @param X
#' @param nsim
#' @param pred_wrapper
#' @param seed
#'
#' @return
#' @import fastshap
#' @export
#'
#' @examples
get_shap_values <- function(in_model, X, nsim, pred_wrapper, seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  shap_values <- fastshap::explain(in_model,
                               X = X,
                               nsim = nsim,
                               pred_wrapper = pred_wrapper,
                               adjust = TRUE)
  shap_values
}

#' Title
#'
#' @param in_data
#' @param shap_values
#' @param in_model
#' @param pred_wrapper
#' @param pred
#'
#' @return
#' @export
#'
#' @examples
get_shap_contribs <- function(in_data, shap_values, in_model, pred_wrapper, pred = FALSE) {

  delta_vector <- rep(NA, ncol(shap_values)+1)
  names(delta_vector) <- c("baseline", colnames(shap_values))

  delta_vector["baseline"] <-
    pred_wrapper(in_model, in_data %>% mutate_all(~ 0)) %>% .[1]

  for (v in colnames(shap_values)) {

    aux <- in_data %>%
      mutate_at(vars(-all_of(v)), ~ 0) %>%
      mutate_at(all_of(v), ~ mean(.))

    delta_vector[v] <- pred_wrapper(in_model, aux[1,]) - delta_vector["baseline"]
  }

  shap_contribs <- shap_values %>%
    as.data.frame() %>%
    mutate(baseline = 0, .before = 1) %>%
    sweep(2, delta_vector[c("baseline", colnames(shap_values))], "+") %>%
    as_tibble()

  if (pred == TRUE) shap_contribs <- shap_contribs %>%
    mutate(y_hat = rowSums(shap_contribs))

  return(shap_contribs)
}

#' Title
#'
#' @param in_model
#' @param in_data
#' @param pred_vars
#' @param pred_fun
#' @param nsim
#' @param grid_resolution
#' @param seed
#'
#' @return
#' @export
#'
#' @examples
shap_contribs <- function(in_model,
                          in_data,
                          pred_vars,
                          pred_fun,
                          nsim = 10,
                          grid_resolution = 20,
                          seed = NULL) {

  shap_values <- get_shap_values(in_model     = in_model,
                                 X            = in_data,
                                 nsim         = nsim,
                                 pred_wrapper = pred_fun,
                                 seed         = seed)

  shap_contribs <- get_shap_contribs(in_data      = in_data,
                                     shap_values  = shap_values,
                                     in_model     = in_model,
                                     pred_wrapper = pred_fun)

  out <- list(contribs        = shap_contribs,
              smooth_contribs = NULL,
              contrib_grid    = NULL,
              contrib_funs    = NULL)

  return(out)

}
