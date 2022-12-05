# library(tidyverse)
# library(lme4)


#' Title
#'
#' @param x an
#'
#' @return
#' @import lme4
#' @import stringr
#' @import purrr
#' @export
#'
#' @examples
get_lme4_str_formulas <- function(x) {

  if (inherits(x, "lmerMod")) {
    str_f <- paste0(format(summary(x)$call$formula), collapse = "")
  } else {
    if (inherits(x, "formula")) {
      str_f <- paste(format(x), collapse = "")
    } else
    {
      stop("x must be an lmerMod object or a formula")
    }
  }

  str_response <- str_f %>% str_remove("~.*$") %>% str_trim()

  fe_str_f <- str_f %>%
    str_remove("^.*~") %>%
    str_remove("\\(.*$") %>%
    str_remove("[:space:]*\\+[:space:]*$") %>%
    str_squish()
  if (fe_str_f == "") fe_str_f <- NA


  re_str_f <- str_f %>%
    str_extract("\\(.*\\)") %>%
    str_split("\\+[:digit:]*\\(") %>%
    map(~ .x %>% str_remove("\\|.*\\)") %>%
          str_remove("\\(|\\)") %>%
          str_squish()) %>%
    .[[1]]

  groups_f <- str_f %>%
    str_extract_all("\\|[:space:]*[^\\)]+\\)") %>%
    map(~ .x %>%
          stringr::str_remove_all("[\\(|\\||\\)]") %>%
          stringr::str_squish()) %>%
    unlist()

  out <- list(str_response = str_response,
              fe_str_f = fe_str_f,
              re_str_f = re_str_f,
              groups_f = groups_f)

  return(out)
}

#' Title
#'
#' @param x
#' @param type
#'
#' @return
#' @export
#'
#' @examples
get_model_vars <- function(x, type = c("all", "fe", "re")) {

  type <-  match.arg(type, choices = c("all", "fe", "re"))
  f_components <- switch(type,
         all = c("fe_str_f", "re_str_f"),
         fe  = "fe_str_f",
         re  = "re_str_f"
         )

  str_ff <- get_lme4_str_formulas(x)[f_components]

  model_vars <- unique(unlist(map(str_ff, ~ unlist(str_split(.x, "\\+"))))) %>%
    str_trim()

  return(model_vars[!is.na(model_vars)])
}

#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
get_model_response <- function(x) {

  str_response <- get_lme4_str_formulas(x)[["str_response"]]

  return(str_response)
}

#' Title
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
get_model_groups <- function(x) {

  str_response <- get_lme4_str_formulas(x)[["groups_f"]]

  return(str_response)
}


#' Title
#'
#' @param m_h_lm
#' @param base
#'
#' @return
#' @import tidyr
#' @export
#'
#' @examples
calc_contribs_fe <- function(m_h_lm, base) {

  base <- base %>% mutate(contrib_idx = 1:nrow(base), .before = 1)

  model_vars     <- get_model_vars(m_h_lm, type = "fe")
  model_response <- get_model_response(m_h_lm)

  contrib_f <-
    as.formula(paste(model_response, "~ ", paste(model_vars, collapse = "+")))

  vars_in_base <- c(model_response,
                    setdiff(model_vars %>%
                              map(~ str_split(.x, ":")) %>%
                              unlist(),
                            "1"))
  base <- base %>% select(contrib_idx, all_of(vars_in_base)) %>% drop_na()

  mm <- model.matrix(contrib_f, base)
  mm <- mm[, order(colnames(mm))]

  coefs <- fixef(m_h_lm)
  coefs <- coefs[colnames(mm)]

  contribs <- base %>% select(contrib_idx) %>%
    bind_cols(mm %>% sweep(2, coefs, "*")) %>%
    as_tibble() %>%
    arrange(contrib_idx)

  return(contribs)
}


#' Title
#'
#' @param m_h_lm
#' @param base
#' @param groups
#'
#' @return
#' @export
#'
#' @examples
# calc_contribs_re <- function(m_h_lm, base, groups = NULL) {
#
#   base <- base %>% mutate(contrib_idx = 1:nrow(base), .before = 1)
#
#   model_vars     <- get_model_vars(m_h_lm, type = "re")
#   model_response <- get_model_response(m_h_lm)
#   model_groups   <- get_model_groups(m_h_lm)
#   model_groups   <- model_groups[order(model_groups)]
#
#   contrib_f <-
#     as.formula(paste(model_response, "~ ", paste(model_vars, collapse = "+")))
#
#   vars_in_base <- c(model_response,
#                     model_groups,
#                     setdiff(model_vars %>%
#                             map(~ str_split(.x, ":")) %>%
#                             unlist(),
#                           "1"))
#   base <- base %>% select(contrib_idx, all_of(vars_in_base)) %>% drop_na()
#
#   mm <- model.matrix(contrib_f, base)
#   mm <- mm[, order(colnames(mm))]
#
#   if (!is.null(groups)) model_groups <- intersect(model_groups, groups)
#   contribs_g <- vector(mode = "list", length = length(model_groups))
#   names(contribs_g) <- model_groups
#
#   for (g in model_groups) {
#     g_elements <- unique(base[[g]])
#
#     g_elements <- g_elements[order(g_elements)]
#     contribs_lst <- vector(mode = "list", length = length(g_elements))
#     names(contribs_lst) <- g_elements
#
#     for (g_element in g_elements) {
#
#       idx_g_element <- which(base[[g]] == g_element)
#
#       m_matrix_g_element <- mm[idx_g_element, , drop = FALSE]
#       coefs_g_element <- ranef(m_h_lm)[[g]][g_element, colnames(m_matrix_g_element)]
#       coefs_g_element <- coefs_g_element[, order(names(coefs_g_element))] %>%
#         unlist()
#
#       contribs_lst[[g_element]] <-
#         base[idx_g_element,] %>% select(contrib_idx) %>%
#         bind_cols(m_matrix_g_element %>% sweep(2, coefs_g_element, "*")) %>%
#         as_tibble()
#     }
#     contribs_g[[g]] <- bind_rows(contribs_lst)
#   }
#
#   contribs_g <- contribs_g %>% map(~ .x %>% arrange(contrib_idx))
#
#   return(contribs_g)
#
# }

calc_contribs_re <- function(m_h_lm, base, groups = NULL) {
  
  base <- base %>% mutate(contrib_idx = 1:nrow(base), .before = 1)
  
  model_vars     <- get_model_vars(m_h_lm, type = "re")
  model_response <- get_model_response(m_h_lm)
  # model_groups   <- get_model_groups(m_h_lm)
  model_groups   <- names(ranef(m_h_lm))
  model_groups   <- model_groups[order(model_groups)]
  model_groups_vars <- model_groups %>%
    str_split((":")) %>%
    unlist() %>%
    unique()
  
  contrib_f <-
    as.formula(paste(model_response, "~ ", paste(model_vars, collapse = "+")))
  
  vars_in_base <- c(model_response,
                    model_groups_vars,
                    setdiff(model_vars %>%
                              map(~ str_split(.x, ":")) %>%
                              unlist(),
                            "1"))
  base <- base %>% select(contrib_idx, all_of(vars_in_base)) %>% drop_na()
  
  if (!is.null(groups)) model_groups <- intersect(model_groups, groups)
  contribs_g <- vector(mode = "list", length = length(model_groups))
  names(contribs_g) <- model_groups
  
  for (g in model_groups) {
    
    g_vars <- g %>% str_split((":")) %>% unlist() %>%  unique()
    g_elements <- base %>% select(all_of(g_vars)) %>% distinct() %>% arrange()
    
    contribs_lst <- vector(mode = "list", length = nrow(g_elements))
    names(contribs_lst) <- g_elements %>% unite("group", everything()) %>% unlist()
    
    for (idx_g_element in 1:nrow(g_elements)) {
      
      g_element <- g_elements[idx_g_element, , drop = FALSE]
      if (g_element == "susana_sexo") 
      {
        print(g_element)
      }
      aux_base <- base %>% inner_join(g_element, by = names(g_element))
      m_matrix_g_element <- model.matrix(contrib_f, aux_base)
      m_matrix_g_element <- m_matrix_g_element[, order(colnames(m_matrix_g_element)), drop = FALSE]
      
      coefs_g_element <- ranef(m_h_lm)[[g]][unlist(g_element) %>% paste(collapse=":"), colnames(m_matrix_g_element)]
      coefs_g_element <- coefs_g_element[, order(names(coefs_g_element))] %>%
        unlist()
      
      contribs_lst[[idx_g_element]] <-
        aux_base %>% select(contrib_idx) %>%
        bind_cols(m_matrix_g_element %>% sweep(2, coefs_g_element, "*")) %>%
        as_tibble()
    }
    contribs_g[[g]] <- bind_rows(contribs_lst)
  }
  
  contribs_g <- contribs_g %>% map(~ .x %>% arrange(contrib_idx))
  
  return(contribs_g)
  
}

#' Title
#'
#' @param c1
#' @param c2
#'
#' @return
#' @export
#'
#' @examples
sum_contribs <- function(c1, c2) {

  common_vars <- setdiff(intersect(colnames(c1), colnames(c2)), "contrib_idx")
  only_in_c1  <- setdiff(setdiff(colnames(c1), colnames(c2)), "contrib_idx")
  only_in_c2  <- setdiff(setdiff(colnames(c2), colnames(c1)), "contrib_idx")

  out <- NULL

  if (length(common_vars) > 0) {
    out <- bind_cols(c1 %>% select(contrib_idx),
                     c1[, common_vars] + c2[, common_vars])
  }

  if (length(only_in_c1) > 0) {
    if (is.null(out)) out <- c1 %>% select(contrib_idx)
    out <- bind_cols(out, c1[, only_in_c1])
  }

  if (length(only_in_c2) > 0) {
    if (is.null(out)) out <- c2 %>% select(contrib_idx)
    out <- bind_cols(out, c2[, only_in_c1])
  }

  return(out)
}


#' Title
#'
#' @param m_h_lm
#' @param base
#' @param groups
#'
#' @return
#' @export
#'
#' @examples
calc_contribs <- function(m_h_lm, base, groups = NULL) {

  contribs_fe <- calc_contribs_fe(m_h_lm, base)
  contribs_re <- calc_contribs_re(m_h_lm, base, groups = groups)

  contribs <- contribs_fe
  for (i in seq_along(contribs_re)) {
    contribs <- sum_contribs(contribs, contribs_re[[i]])
  }

  return(as_tibble(contribs))

}
