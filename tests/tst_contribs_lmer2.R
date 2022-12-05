
library(tidyverse)
library(lme4)

# source("./Hierarchichal Models/contribs_lmer.R")

data("heights", package = "modelr")
heights

?modelr::heights


h_m1 <- lmer(income ~ 1 + age + sex + age:sex + (1 + age + weight | marital), heights)

summary(h_m1)

coef(h_m1)
fixef(h_m1)
ranef(h_m1)

confint(h_m1, method = "Wald")

X <- getME(h_m1, "X")
colnames(X)
dim(X)

the_coefs <- coef(h_m1)$marital
the_coefs
dim(the_coefs)

get_mlmer_contribs(in_model = h_m1, pred = TRUE)
predict(h_m1) %>% head()


####

get_lme4_str_formulas(h_m1)
get_model_vars(h_m1)

get_model_response(h_m1)
get_model_groups(h_m1)
get_mlmer_contribs(h_m1)


####
f <- income ~ 1 + height + weight + age + education + afqt +
  (1 + height + weight + age + education + afqt | marital) +
  (1 + height + weight + age + education + afqt | sex)

get_lme4_str_formulas(f)
get_model_vars(f)
get_model_vars(f, "all")
get_model_vars(f, "fe")
get_model_vars(f, "re")
get_model_response(f)
get_model_groups(f)

####
h_m2 <- lmer(income ~ 1 + age + weight +
               (1 + age + weight | marital) +
               (1 + age + weight | sex), heights)

calc_contribs_fe(h_m2, heights)

calc_contribs_fe(h_m2, heights) %>%
  dplyr::rowwise() %>%
               dplyr::mutate(pred = sum(dplyr::c_across(-1))) %>%
               dplyr::ungroup()

predict(h_m2, heights, re.form = NA) %>% head(10)

calc_contribs_re(h_m2, heights)
calc_contribs_re(h_m2, heights) %>%
  purrr::map(~ .x %>%
               dplyr::rowwise() %>%
               dplyr::mutate(pred = sum(dplyr::c_across(-1))) %>%
               dplyr::ungroup())

calc_contribs(h_m2, heights) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(pred = sum(dplyr::c_across(-1)))

calc_contribs(h_m2, heights, c("marital", "sex")) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(pred = sum(dplyr::c_across(-1)))

predict(h_m2) %>% head()

predict(h_m2, random.only = TRUE) %>% head()

calc_contribs(h_m2, heights, "marital") %>%
  dplyr::rowwise() %>%
  dplyr::mutate(pred = sum(dplyr::c_across(-1)))

predict(h_m2, re.form = ~ 1 + age + weight +
          (1 + age + weight | marital)) %>% head()

calc_contribs(h_m2, heights, "sex") %>%
  dplyr::rowwise() %>%
  dplyr::mutate(pred = sum(dplyr::c_across(-1)))

predict(h_m2, re.form = ~ 1 + age + weight +
          (1 + age + weight | sex)) %>% head()

####

h_m3 <- lmer(income ~ 1 + age + weight +
               (1 + age + weight | marital:sex)  +
               (1 + age + weight | marital),
             heights)

calc_contribs(h_m3, heights) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(pred = sum(dplyr::c_across(-1)))

predict(h_m3) %>% head()
