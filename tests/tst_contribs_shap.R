library(tidyverse)
library(ranger)

data("Boston", package = "MASS")

x_var <- c("lstat", "rm", "dis", "indus")
y_var <- "medv"

Boston_train <- Boston[-1:-6, c(x_var, y_var)]
Boston_test  <- Boston[1:6, x_var]

f <- as.formula(paste(y_var, "~", paste(x_var, collapse = "+")))
model_lm <- lm(f, data = Boston_train)
model_rf <- ranger(f, data = Boston_train)

## LM
contribs_lm <- model.matrix(f, Boston_train) %>%
  sweep(2, coef(model_lm), "*") %>%
  as_tibble() %>%
  mutate(pred = rowSums(.))
head(contribs_lm)

# Compute approximate Shapley values using 10 Monte Carlo simulations
set.seed(101)  # for reproducibility
shap_lm <- fastshap::explain(model_lm,
                             X = subset(Boston_train, select = -medv),
                             nsim = 10,
                             pred_wrapper = predict,
                             adjust = TRUE)
shap_lm

shap_contribs_lm <- get_shap_contribs(Boston_train, shap_lm, model_lm, predict, pred = TRUE)
shap_contribs_lm
head(contribs_lm)

plot(Boston_train$lstat, shap_contribs_lm$lstat)

## RF

# Compute approximate Shapley values using 10 Monte Carlo simulations
set.seed(101)  # for reproducibility
shap_rf <- fastshap::explain(model_rf,
                             X = subset(Boston_train, select = -medv),
                             nsim = 10,
                             pred_wrapper = function(object, newdata) predict(object, newdata)$predictions,
                             adjust = TRUE)
shap_rf

shap_contribs_rf <- get_shap_contribs(Boston_train, shap_rf,
                                      model_rf,
                                      function(object, newdata) predict(object, newdata)$predictions,
                                      pred = TRUE)
shap_contribs_rf

plot(Boston_train$lstat, shap_contribs_rf$lstat)
plot(Boston_train$rm, shap_contribs_rf$rm)
plot(Boston_train$dis, shap_contribs_rf$dis)
plot(Boston_train$indus, shap_contribs_rf$indus)


shap_contribs(in_model = model_rf,
              in_data = Boston_train,
              pred_vars = NULL,
              pred_fun = function(object, newdata) predict(object, newdata)$predictions,
              nsim = 10,
              grid_resolution = 20,
              seed = NULL)
