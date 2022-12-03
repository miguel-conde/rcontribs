
library(ranger)
library(caret)

library(pdp)

library(ggfortify)
library(tidyverse)

source("interpretability/utils_pdp_contrib.R")

# DATA --------------------------------------------------------------------

data("Boston", package = "MASS")


# MODELS ------------------------------------------------------------------


# Linear model ------------------------------------------------------------

m_lm <- lm(medv ~ ., Boston)

contribs_lm <- model.matrix(medv ~., Boston) %>% sweep(2, coef(m_lm), "*") %>% 
  as_tibble() %>% 
  mutate(y_hat = rowSums(.))

predictions_lm <- fitted(m_lm)
Metrics::rmse(Boston$medv, predictions_lm)
plot(Boston$medv, fitted(m_lm), xlab = "Actual", ylab = "Predicted")
abline(a = 0, b = 1)

# Random Forest -----------------------------------------------------------

m_rf <- ranger(medv ~ ., Boston)

predictions_rf <- predict(m_rf, Boston)$predictions
Metrics::rmse(Boston$medv, predictions_rf)

plot(Boston$medv, predictions_rf, xlab = "Actual", ylab = "Predicted")
abline(a = 0, b = 1)

# Caret - XgBoost ---------------------------------------------------------

the_grid <- data.frame(nrounds = 500,
                       max_depth = 6,
                       eta = 0.3,
                       gamma = 0, # ??
                       colsample_bytree = 1,
                       min_child_weight = 1,
                       subsample = 1)

m_xgboost <- train(medv ~ ., Boston, 
                   method = "xgbTree",
                   tuneGrid = the_grid,
                   trControl = trainControl(method = "none"))

predictions_xgboost <- predict(m_xgboost$finalModel, 
                               Boston %>% select(-medv) %>% as.matrix()) 
Metrics::rmse(Boston$medv, predictions_xgboost)

plot(Boston$medv, predictions_xgboost, xlab = "Actual", ylab = "Predicted")
abline(a = 0, b = 1)


# glm - lognormal ---------------------------------------------------------

m_glm <- glm(medv ~ ., Boston, family = gaussian(link = log))

predictions_glm <- predict(m_glm, Boston, type = "response") 
Metrics::rmse(Boston$medv, predictions_glm)

plot(Boston$medv, predictions_glm, xlab = "Actual", ylab = "Predicted")
abline(a = 0, b = 1)

# PDP ---------------------------------------------------------------------

# https://bradleyboehmke.github.io/HOML/iml.html

# Custom prediction function wrapper
pdp_pred <- function(object, newdata)  {
  results <- mean(as.vector(predict(object, newdata)))
  return(results)
}

# Compute partial dependence values
pd_values <- partial(
  m_lm,
  train = Boston, 
  pred.var = "lstat",
  pred.fun = pdp_pred,
  grid.resolution = 20
)
head(pd_values)  # take a peak

autoplot(pd_values) +
  geom_abline(intercept = 0, slope = as.numeric(coef(m_lm)["lstat"]), color = "red") +
  coord_cartesian(ylim = c(-20, 30))

# Construct c-ICE curves
partial(
  m_lm,
  train = Boston, 
  pred.var = "lstat",
  pred.fun = pdp_pred,
  grid.resolution = 20,
  plot = TRUE,
  center = TRUE,
  plot.engine = "ggplot2"
)  +
  geom_abline(intercept = 0, slope = as.numeric(coef(m_lm)["lstat"]), color = "red") +
  coord_cartesian(ylim = c(-20, 30))

# La distancia es fitted - fitted cuando todos valen 0 menos el de interés
mean(fitted(m_lm) - Boston$lstat * coef(m_lm)["lstat"])
mean(fitted(m_lm)) - mean(Boston$lstat) * coef(m_lm)["lstat"]

# UTILS -------------------------------------------------------------------


# Linear Model ------------------------------------------------------------

res <- pdp_1_contrib(m_lm, Boston, "lstat", pdp_pred_lm)

autoplot(res) +
  geom_abline(intercept = 0, slope = as.numeric(coef(m_lm)["lstat"]), 
              color = "red", linetype = 2) +
  coord_cartesian(ylim = c(-20, 30))

res <- pdp_1_contrib(m_lm, Boston, "black", pdp_pred_lm)

autoplot(res) +
  geom_abline(intercept = 0, slope = as.numeric(coef(m_lm)["black"]), color = "red") +
  coord_cartesian(ylim = c(-20, 30))

res_all_lm <- pdp_contribs(m_lm, Boston, 
                           Boston %>% select(-medv) %>% names(), 
                           pdp_pred_lm)

probe_lm <- res_all_lm$contribs %>% mutate(y_hat = rowSums(.))
probe_lm

plot(predictions_lm, probe_lm$y_hat)
abline(a= 0, b = 1)

probe_lm %>% summarise_all(~ sum(.))
sum(predictions_lm)

contribs_lm %>% summarise_all(~ sum(.))

curve(res_all_lm$contrib_funs$crim(x), from = min(Boston$crim), to = max(Boston$crim))
curve(res_all_lm$contrib_funs$zn(x), from = min(Boston$zn), to = max(Boston$zn))
curve(res_all_lm$contrib_funs$indus(x), from = min(Boston$indus), to = max(Boston$indus))
curve(res_all_lm$contrib_funs$chas(x), from = min(Boston$chas), to = max(Boston$chas))
curve(res_all_lm$contrib_funs$nox(x), from = min(Boston$nox), to = max(Boston$nox))
curve(res_all_lm$contrib_funs$rm(x), from = min(Boston$rm), to = max(Boston$rm))
curve(res_all_lm$contrib_funs$age(x), from = min(Boston$age), to = max(Boston$age))
curve(res_all_lm$contrib_funs$dis(x), from = min(Boston$dis), to = max(Boston$dis))
curve(res_all_lm$contrib_funs$rad(x), from = min(Boston$rad), to = max(Boston$rad))
curve(res_all_lm$contrib_funs$tax(x), from = min(Boston$tax), to = max(Boston$tax))
curve(res_all_lm$contrib_funs$ptratio(x), from = min(Boston$ptratio), to = max(Boston$ptratio))
curve(res_all_lm$contrib_funs$black(x), from = min(Boston$black), to = max(Boston$black))
curve(res_all_lm$contrib_funs$lstat(x), from = min(Boston$lstat), to = max(Boston$lstat))


# Random Forest -----------------------------------------------------------

res_rf <- pdp_1_contrib(m_rf, Boston, "lstat", pdp_pred_rf)

autoplot(res_rf) +
  geom_abline(intercept = 0, slope = as.numeric(coef(m_lm)["lstat"]), color = "red") +
  coord_cartesian(ylim = c(-20, 30))

res_rf <- pdp_1_contrib(m_rf, Boston, "black", pdp_pred_rf)

autoplot(res_rf) +
  geom_abline(intercept = 0, slope = as.numeric(coef(m_lm)["black"]), color = "red") +
  coord_cartesian(ylim = c(-20, 30))

res_rf <- pdp_1_contrib(m_rf, Boston, "dis", pdp_pred_rf)

autoplot(res_rf) +
  geom_abline(intercept = 0, slope = as.numeric(coef(m_lm)["dis"]), color = "red") +
  coord_cartesian(ylim = c(-20, 30))

res_rf <- pdp_1_contrib(m_rf, Boston, "ptratio", pdp_pred_rf)

autoplot(res_rf) +
  geom_abline(intercept = 0, slope = as.numeric(coef(m_lm)["ptratio"]), color = "red") +
  coord_cartesian(ylim = c(-20, 30))


##

res_all_rf <- pdp_contribs(m_rf, Boston, 
                           Boston %>% select(-medv) %>% names(), 
                           pdp_pred_rf)

autoplot(res_all_rf$contrib_grid$crim) +
  geom_abline(intercept = 0, slope = as.numeric(coef(m_lm)["crim"]), color = "red") +
  geom_point(data = tibble(x = seq(0:100), y = res_all_rf$contrib_funs$crim(seq(0:100))),
             aes(x = x, y = y), alpha = .1) +
  coord_cartesian(ylim = c(-20, 30))

autoplot(res_all_rf$contrib_grid$tax) +
  geom_abline(intercept = 0, slope = as.numeric(coef(m_lm)["tax"]), color = "red") +
  coord_cartesian(ylim = c(-20, 30))

autoplot(res_all_rf$contrib_grid$zn) +
  geom_abline(intercept = 0, slope = as.numeric(coef(m_lm)["zn"]), color = "red") +
  coord_cartesian(ylim = c(-20, 30))

autoplot(res_all_rf$contrib_grid$indus) +
  geom_abline(intercept = 0, slope = as.numeric(coef(m_lm)["indus"]), color = "red") +
  coord_cartesian(ylim = c(-20, 30))

autoplot(res_all_rf$contrib_grid$chas) +
  geom_abline(intercept = 0, slope = as.numeric(coef(m_lm)["chas"]), color = "red") +
  coord_cartesian(ylim = c(-20, 30))

probe_rf <- res_all_rf$contribs %>% mutate(y_hat = rowSums(.))
probe_rf

plot(predictions_rf, probe_rf$y_hat)
abline(a= 0, b = 1)

probe_rf %>% summarise_all(~ sum(.))
sum(predictions_rf)

contribs_lm %>% summarise_all(~ sum(.))

curve(res_all_rf$contrib_funs$crim(x), from = min(Boston$crim), to = max(Boston$crim))
curve(res_all_rf$contrib_funs$zn(x), from = min(Boston$zn), to = max(Boston$zn))
curve(res_all_rf$contrib_funs$indus(x), from = min(Boston$indus), to = max(Boston$indus))
curve(res_all_rf$contrib_funs$chas(x), from = min(Boston$chas), to = max(Boston$chas))
curve(res_all_rf$contrib_funs$nox(x), from = min(Boston$nox), to = max(Boston$nox))
curve(res_all_rf$contrib_funs$rm(x), from = min(Boston$rm), to = max(Boston$rm))
curve(res_all_rf$contrib_funs$age(x), from = min(Boston$age), to = max(Boston$age))
curve(res_all_rf$contrib_funs$dis(x), from = min(Boston$dis), to = max(Boston$dis))
curve(res_all_rf$contrib_funs$rad(x), from = min(Boston$rad), to = max(Boston$rad))
curve(res_all_rf$contrib_funs$tax(x), from = min(Boston$tax), to = max(Boston$tax))
curve(res_all_rf$contrib_funs$ptratio(x), from = min(Boston$ptratio), to = max(Boston$ptratio))
curve(res_all_rf$contrib_funs$black(x), from = min(Boston$black), to = max(Boston$black))
curve(res_all_rf$contrib_funs$lstat(x), from = min(Boston$lstat), to = max(Boston$lstat))





# LM vs RF ----------------------------------------------------------------
curve(res_all_rf$contrib_funs$crim(x), from = min(Boston$crim), to = max(Boston$crim),
      ylim = c(-2, 1.5))
curve(res_all_lm$contrib_funs$crim(x), from = min(Boston$crim), to = max(Boston$crim),
      add = TRUE, col = "red")

curve(res_all_rf$contrib_funs$zn(x), from = min(Boston$zn), to = max(Boston$zn),
      ylim = c(-0.2, 4.5))
curve(res_all_lm$contrib_funs$zn(x), from = min(Boston$zn), to = max(Boston$zn),
      add = TRUE, col = "red")

curve(res_all_rf$contrib_funs$indus(x), from = min(Boston$indus), to = max(Boston$indus))
curve(res_all_lm$contrib_funs$indus(x), from = min(Boston$indus), to = max(Boston$indus),
      add = TRUE, col = "red")

curve(res_all_rf$contrib_funs$chas(x), from = min(Boston$chas), to = max(Boston$chas),
      ylim = c(0, 3))
curve(res_all_lm$contrib_funs$chas(x), from = min(Boston$chas), to = max(Boston$chas),
      add = TRUE, col = "red")

curve(res_all_rf$contrib_funs$nox(x), from = min(Boston$nox), to = max(Boston$nox),
      ylim = c(-14, 0.8))
curve(res_all_lm$contrib_funs$nox(x), from = min(Boston$nox), to = max(Boston$nox),
      add = TRUE, col = "red")

curve(res_all_rf$contrib_funs$rm(x), from = min(Boston$rm), to = max(Boston$rm),
      ylim = c(-1, 35))
curve(res_all_lm$contrib_funs$rm(x), from = min(Boston$rm), to = max(Boston$rm),
      add = TRUE, col = "red")

curve(res_all_rf$contrib_funs$age(x), from = min(Boston$age), to = max(Boston$age))
curve(res_all_lm$contrib_funs$age(x), from = min(Boston$age), to = max(Boston$age),
      add = TRUE, col = "red")

curve(res_all_rf$contrib_funs$dis(x), from = min(Boston$dis), to = max(Boston$dis),
      ylim = c(-20, 3))
curve(res_all_lm$contrib_funs$dis(x), from = min(Boston$dis), to = max(Boston$dis),
      add = TRUE, col = "red")

curve(res_all_rf$contrib_funs$rad(x), from = min(Boston$rad), to = max(Boston$rad),
      ylim = c(-0.2, 7))
curve(res_all_lm$contrib_funs$rad(x), from = min(Boston$rad), to = max(Boston$rad),
      add = TRUE, col = "red")

curve(res_all_rf$contrib_funs$tax(x), from = min(Boston$tax), to = max(Boston$tax),
      ylim = c(-9, -1))
curve(res_all_lm$contrib_funs$tax(x), from = min(Boston$tax), to = max(Boston$tax),
      add = TRUE, col = "red")

curve(res_all_rf$contrib_funs$ptratio(x), from = min(Boston$ptratio), to = max(Boston$ptratio),
      ylim = c(-22, 0.5))
curve(res_all_lm$contrib_funs$ptratio(x), from = min(Boston$ptratio), to = max(Boston$ptratio),
      add = TRUE, col = "red")

curve(res_all_rf$contrib_funs$black(x), from = min(Boston$black), to = max(Boston$black),
      ylim = c(0, 4))
curve(res_all_lm$contrib_funs$black(x), from = min(Boston$black), to = max(Boston$black),
      add = TRUE, col = "red")

curve(res_all_rf$contrib_funs$lstat(x), from = min(Boston$lstat), to = max(Boston$lstat),
      ylim = c(-20, 2))
curve(res_all_lm$contrib_funs$lstat(x), from = min(Boston$lstat), to = max(Boston$lstat),
      add = TRUE, col = "red")

# XgBoost -----------------------------------------------------------------

res_all_xgboost <- pdp_contribs(m_xgboost, Boston, 
                                Boston %>% select(-medv) %>% names(), 
                                pdp_pred_xgboost)

res_all_xgboost$contribs

probe_xgboost <- res_all_xgboost$contribs %>% mutate(y_hat = rowSums(.))
probe_xgboost

plot(predictions_xgboost, probe_xgboost$y_hat)
abline(a= 0, b = 1)

probe_xgboost %>% summarise_all(~ sum(.))
sum(predictions_xgboost)

contribs_lm %>% summarise_all(~ sum(.))

curve(res_all_xgboost$contrib_funs$crim(x), from = min(Boston$crim), to = max(Boston$crim))
curve(res_all_xgboost$contrib_funs$zn(x), from = min(Boston$zn), to = max(Boston$zn))
curve(res_all_xgboost$contrib_funs$indus(x), from = min(Boston$indus), to = max(Boston$indus))
curve(res_all_xgboost$contrib_funs$chas(x), from = min(Boston$chas), to = max(Boston$chas))
curve(res_all_xgboost$contrib_funs$nox(x), from = min(Boston$nox), to = max(Boston$nox))
curve(res_all_xgboost$contrib_funs$rm(x), from = min(Boston$rm), to = max(Boston$rm))
curve(res_all_xgboost$contrib_funs$age(x), from = min(Boston$age), to = max(Boston$age))
curve(res_all_xgboost$contrib_funs$dis(x), from = min(Boston$dis), to = max(Boston$dis))
curve(res_all_xgboost$contrib_funs$rad(x), from = min(Boston$rad), to = max(Boston$rad))
curve(res_all_xgboost$contrib_funs$tax(x), from = min(Boston$tax), to = max(Boston$tax))
curve(res_all_xgboost$contrib_funs$ptratio(x), from = min(Boston$ptratio), to = max(Boston$ptratio))
curve(res_all_xgboost$contrib_funs$black(x), from = min(Boston$black), to = max(Boston$black))
curve(res_all_xgboost$contrib_funs$lstat(x), from = min(Boston$lstat), to = max(Boston$lstat))


# Glm ---------------------------------------------------------------------
res_all_glm <- pdp_contribs(m_glm, Boston, 
                                Boston %>% select(-medv) %>% names(), 
                                pdp_pred_glm)

res_all_glm$contribs

probe_glm <- res_all_glm$contribs %>% mutate(y_hat = rowSums(.))
probe_glm

plot(predictions_glm, probe_glm$y_hat)
abline(a= 0, b = 1)

probe_glm %>% summarise_all(~ sum(.))
sum(predictions_glm)

contribs_lm %>% summarise_all(~ sum(.))

curve(res_all_glm$contrib_funs$crim(x), from = min(Boston$crim), to = max(Boston$crim))
curve(res_all_glm$contrib_funs$zn(x), from = min(Boston$zn), to = max(Boston$zn))
curve(res_all_glm$contrib_funs$indus(x), from = min(Boston$indus), to = max(Boston$indus))
curve(res_all_glm$contrib_funs$chas(x), from = min(Boston$chas), to = max(Boston$chas))
curve(res_all_glm$contrib_funs$nox(x), from = min(Boston$nox), to = max(Boston$nox))
curve(res_all_glm$contrib_funs$rm(x), from = min(Boston$rm), to = max(Boston$rm))
curve(res_all_glm$contrib_funs$age(x), from = min(Boston$age), to = max(Boston$age))
curve(res_all_glm$contrib_funs$dis(x), from = min(Boston$dis), to = max(Boston$dis))
curve(res_all_glm$contrib_funs$rad(x), from = min(Boston$rad), to = max(Boston$rad))
curve(res_all_glm$contrib_funs$tax(x), from = min(Boston$tax), to = max(Boston$tax))
curve(res_all_glm$contrib_funs$ptratio(x), from = min(Boston$ptratio), to = max(Boston$ptratio))
curve(res_all_glm$contrib_funs$black(x), from = min(Boston$black), to = max(Boston$black))
curve(res_all_glm$contrib_funs$lstat(x), from = min(Boston$lstat), to = max(Boston$lstat))

# Plots -------------------------------------------------------------------

p1 <- ggplot_1_contrib(res_all_rf, "crim", 
                 x_units = "\n(per capita crime rate by town)",
                 y_units = "Contributions\n(to the median value of owner-occupied homes in $1000s)")
p1

p2 <- ggplot_1_contrib(res_all_rf, "lstat", 
                 x_units = "\n(lower status of the population (percent))",
                 y_units = "Contributions\n(to the median value of owner-occupied homes in $1000s)")
p2

p3 <- ggplot_1_contrib(res_all_rf, "rm", 
                       x_units = "\n(average number of rooms per dwelling)",
                       y_units = "Contributions\n(to the median value of owner-occupied homes in $1000s)")
p3

p <- p1 + p2 + p3 + p2 + p1 + p3 + p1 + p2 + p3 + p2 + p1 + p3 + p1
p

p_rf <- ggplot_contribs(res_all_rf, y_units = "$1000") + 
  plot_annotation(
      title = "Boston - Random Forest Model Contributions",
      subtitle = "(to the median value of owner-occupied homes in $1000s)",
      caption = "Source: mpg dataset in ggplot2"
  )

p_rf

p_xgboost <- ggplot_contribs(res_all_xgboost, y_units = "$1000") + 
  plot_annotation(
    title = "Boston - XgBoost Model Contributions",
    subtitle = "(to the median value of owner-occupied homes in $1000s)",
    caption = "Source: mpg dataset in ggplot2"
  )

p_xgboost

p_lm <- ggplot_contribs(res_all_lm, y_units = "$1000") + 
  plot_annotation(
    title = "Boston - Linear Model Contributions",
    subtitle = "(to the median value of owner-occupied homes in $1000s)",
    caption = "Source: mpg dataset in ggplot2"
  )

p_lm

p_glm <- ggplot_contribs(res_all_glm, y_units = "$1000") + 
  plot_annotation(
    title = "Boston - Generalised Linear Model Contributions",
    subtitle = "(to the median value of owner-occupied homes in $1000s)",
    caption = "Source: mpg dataset in ggplot2"
  )

p_glm

p_rf | p_lm
p_xgboost | p_lm
p_rf | p_xgboost

p_rf + scale_y_continuous(limits = c(-10, 7))

library(RColorBrewer)
colourCount = length(names(Boston))
my_colors = colorRampPalette(brewer.pal(8, "Set2"))(colourCount)
my_colors = colorRampPalette(brewer.pal(9, "Pastel1"))(colourCount)
my_colors = colorRampPalette(brewer.pal(8, "Pastel2"))(colourCount)
my_colors = colorRampPalette(brewer.pal(8, "Accent"))(colourCount)

res_all_rf$contribs %>% 
  mutate(obs = 1:nrow(res_all_rf$contribs), .before = 1) %>% 
  gather(term, contribution, -obs) %>% 
  ggplot(aes(x = obs, y = contribution, fill = term)) +
  # geom_area(colour = "black", size = .2, alpha = .4)  +
  geom_area(position = "fill") +
  scale_fill_manual(values = my_colors) + 
  labs(fill = "Variable")


res_all_lm$contribs %>% 
  mutate(obs = 1:nrow(res_all_rf$contribs), .before = 1) %>% 
  gather(term, contribution, -obs) %>% 
  ggplot(aes(x = obs, y = contribution, fill = term)) +
  # geom_area(colour = "black", size = .2, alpha = .4)  +
  geom_area() +
  scale_fill_manual(values = my_colors) + 
  labs(fill = "Variable")

res_all_xgboost$contribs %>% 
  mutate(obs = 1:nrow(res_all_rf$contribs), .before = 1) %>% 
  gather(term, contribution, -obs) %>% 
  ggplot(aes(x = obs, y = contribution, fill = term)) +
  # geom_area(colour = "black", size = .2, alpha = .4)  +
  geom_area() +
  scale_fill_manual(values = my_colors) + 
  labs(fill = "Variable")


# ROI ---------------------------------------------------------------------

tibble(lstat = Boston$lstat, 
       roi = pdp_1_roi(Boston$lstat, res_all_lm, "lstat")) %>% 
  ggplot(aes(x = lstat, y = roi)) +
  geom_line()

tibble(lstat = Boston$lstat, 
       roi = pdp_1_roi(Boston$lstat, res_all_rf, "lstat")) %>% 
  ggplot(aes(x = lstat, y = roi)) +
  geom_line()

tibble(lstat = Boston$lstat, 
       roi = pdp_1_roi(Boston$lstat, res_all_xgboost, "lstat")) %>% 
  ggplot(aes(x = lstat, y = roi)) +
  geom_line()

tibble(lstat = Boston$lstat, 
       roi = pdp_1_roi(Boston$lstat, res_all_glm, "lstat")) %>% 
  ggplot(aes(x = lstat, y = roi)) +
  geom_line()

## LM
rois_all_lm <- pdp_rois(res_all_lm, Boston)

ggplot_1_roi(res_all_lm, lstat, in_data = Boston,
             x_units = "\nlower status of the population (%)",
             y_units = "ROI ($1000 / 1%)")

ggplot_rois(res_all_lm, in_data = Boston, y_units = "ROI") + 
  plot_annotation(
    title = "Boston - LM Model ROIs",
    # subtitle = "(to the median value of owner-occupied homes in $1000s)",
    caption = "Source: mpg dataset in ggplot2")

## RF
rois_all_rf <- pdp_rois(res_all_rf, Boston)

ggplot_1_roi(res_all_rf, lstat, in_data = Boston, 
             x_units = "\nlower status of the population (%)",
             y_units = "ROI ($1000 / 1%)")

ggplot_rois(res_all_rf, in_data = Boston, y_units = "ROI") + 
  plot_annotation(
    title = "Boston - RF Model ROIs",
    # subtitle = "(to the median value of owner-occupied homes in $1000s)",
    caption = "Source: mpg dataset in ggplot2")

## XGBOOST
rois_all_xgboost <- pdp_rois(res_all_xgboost, Boston)

ggplot_1_roi(res_all_xgboost, lstat, in_data = Boston, 
             x_units = "\nlower status of the population (%)",
             y_units = "ROI ($1000 / 1%)")

ggplot_rois(res_all_xgboost, in_data = Boston, y_units = "ROI") + 
  plot_annotation(
    title = "Boston - XgBoost Model ROIs",
    # subtitle = "(to the median value of owner-occupied homes in $1000s)",
    caption = "Source: mpg dataset in ggplot2")

## GLM
rois_all_glm <- pdp_rois(res_all_glm, Boston)

ggplot_1_roi(res_all_glm, lstat, in_data = Boston, 
             x_units = "\nlower status of the population (%)",
             y_units = "ROI ($1000 / 1%)")

ggplot_rois(res_all_glm, in_data = Boston, y_units = "ROI") + 
  plot_annotation(
    title = "Boston - GLM Model ROIs",
    # subtitle = "(to the median value of owner-occupied homes in $1000s)",
    caption = "Source: mpg dataset in ggplot2")

# Numeric -----------------------------------------------------------------

library(numDeriv)

curve(res_all_glm$contrib_funs$lstat(x), from = min(Boston$lstat), to = max(Boston$lstat))

x <- seq(1.73+.1, 37.97-.1, length.out = 100)
the_grad <- numDeriv::grad(res_all_glm$contrib_funs$lstat, x)

plot(x, the_grad, type = "l")
curve(numDeriv::grad(res_all_glm$contrib_funs$lstat, x), 1.73+.1, 37.97-.10)


# COMPLETE EXAMPLE --------------------------------------------------------


# Setup -------------------------------------------------------------------


# Model
m_glm <- glm(medv ~ ., Boston, family = gaussian(link = log))

# Contributions
contribs_glm <- pdp_contribs(m_glm, Boston, 
                             Boston %>% select(-medv) %>% names(), 
                             pdp_pred_glm)

# ROIs
rois_all_glm <- pdp_rois(contribs_glm, Boston)


# Model Analysis ----------------------------------------------------------

predictions_glm <- predict(m_glm, Boston, type = "response") 
Metrics::rmse(Boston$medv, predictions_glm)

plot(Boston$medv, predictions_glm, xlab = "Actual", ylab = "Predicted")
abline(a = 0, b = 1)

autoplot(m_glm)

# Contributions -----------------------------------------------------------

contribs_glm$contribs

probe_glm <- contribs_glm$contribs %>% mutate(y_hat = rowSums(.))
probe_glm

plot(predictions_glm, probe_glm$y_hat)
abline(a= 0, b = 1)

probe_glm %>% summarise_all(~ sum(.))
sum(predictions_glm)

contribs_lm %>% summarise_all(~ sum(.))


# Plot contributions ------------------------------------------------------

ggplot_1_contrib(contribs_glm, "lstat", 
                 x_units = "\n(lower status of the population (percent))",
                 y_units = "Contributions\n(to the median value of owner-occupied homes in $1000s)")


p_glm <- ggplot_contribs(contribs_glm, y_units = "$1000") + 
  plot_annotation(
    title = "Boston - Generalised Linear Model Contributions",
    subtitle = "(to the median value of owner-occupied homes in $1000s)",
    caption = "Source: mpg dataset in ggplot2"
  )

p_glm


# ROIs --------------------------------------------------------------------

pdp_1_roi(Boston$lstat, contribs_glm, "lstat") 


# Plot ROIs ---------------------------------------------------------------

ggplot_1_roi(contribs_glm, "lstat", Boston,
             x_units = "\nlower status of the population (%)",
             y_units = "ROI ($1000 / 1%)")

ggplot_rois(contribs_glm, Boston, y_units = "ROI") + 
  plot_annotation(
    title = "Boston - GLM Model ROIs",
    # subtitle = "(to the median value of owner-occupied homes in $1000s)",
    caption = "Source: mpg dataset in ggplot2")


# ALL PLOTS ---------------------------------------------------------------

## LM
contribs_lm <- pdp_contribs(m_lm, Boston, 
                            Boston %>% select(-medv) %>% names(), 
                            pdp_pred_lm)


p_lm <- ggplot_contribs(contribs_lm, y_units = "$1000") + 
  plot_annotation(
    title = "Boston - Linear Model Contributions",
    subtitle = "(to the median value of owner-occupied homes in $1000s)",
    caption = "Source: mpg dataset in ggplot2"
  )

rois_lm <- pdp_rois(contribs_lm, Boston)

p_lm_rois <- ggplot_rois(contribs_lm, Boston, y_units = "ROI") + 
  plot_annotation(
    title = "Boston - LM Model ROIs",
    # subtitle = "(to the median value of owner-occupied homes in $1000s)",
    caption = "Source: mpg dataset in ggplot2")

## RF
contribs_rf <- pdp_contribs(m_rf, Boston, 
                            Boston %>% select(-medv) %>% names(), 
                            pdp_pred_rf)

p_rf <- ggplot_contribs(contribs_rf, y_units = "$1000") + 
  plot_annotation(
    title = "Boston - Random Forest Model Contributions",
    subtitle = "(to the median value of owner-occupied homes in $1000s)",
    caption = "Source: mpg dataset in ggplot2"
  )

rois_rf <- pdp_rois(contribs_rf, Boston) 

p_rf_rois <- ggplot_rois(contribs_rf, Boston, y_units = "ROI") + 
  plot_annotation(
    title = "Boston - RF Model ROIs",
    # subtitle = "(to the median value of owner-occupied homes in $1000s)",
    caption = "Source: mpg dataset in ggplot2")

## XgBoost
contribs_xgboost <- pdp_contribs(m_xgboost, Boston, 
                                 Boston %>% select(-medv) %>% names(), 
                                 pdp_pred_xgboost)

p_xgboost <- ggplot_contribs(contribs_xgboost, y_units = "$1000") + 
  plot_annotation(
    title = "Boston - XgBoost Model Contributions",
    subtitle = "(to the median value of owner-occupied homes in $1000s)",
    caption = "Source: mpg dataset in ggplot2"
  )

rois_xgboost <- pdp_rois(contribs_xgboost, Boston)

p_xgboost_rois <- ggplot_rois(contribs_xgboost, Boston, y_units = "ROI") + 
  plot_annotation(
    title = "Boston - XgBoost Model ROIs",
    # subtitle = "(to the median value of owner-occupied homes in $1000s)",
    caption = "Source: mpg dataset in ggplot2")

## GLM
contribs_glm <- pdp_contribs(m_glm, Boston, 
                             Boston %>% select(-medv) %>% names(), 
                             pdp_pred_glm)

p_glm <- ggplot_contribs(contribs_glm, y_units = "$1000") + 
  plot_annotation(
    title = "Boston - Generalised Linear Model Contributions",
    subtitle = "(to the median value of owner-occupied homes in $1000s)",
    caption = "Source: mpg dataset in ggplot2"
  )

rois_glm <- pdp_rois(contribs_glm, Boston)

p_glm_rois <- ggplot_rois(contribs_glm, Boston, y_units = "ROI") + 
  plot_annotation(
    title = "Boston - GLM Model ROIs",
    # subtitle = "(to the median value of owner-occupied homes in $1000s)",
    caption = "Source: mpg dataset in ggplot2")

p_lm | p_glm

p_glm | p_rf

p_glm | p_xgboost

p_rf | p_xgboost


# CURVE FITTING -----------------------------------------------------------

## SIGMOID
estimate_0_sigmoid(contribs_xgboost$contrib_grid$zn, "zn")
res_optim <- optim(par = estimate_0_sigmoid(contribs_xgboost$contrib_grid$zn, "zn"), 
                   fn = loss_fun,
                   contribs = contribs_xgboost$contrib_grid$zn,
                   FUN = sigmoid_fun,
                   method = "L-BFGS-B")

contribs_xgboost$contrib_grid$zn %>% plot(type = "l")
curve(sigmoid_fun(x, as.list(res_optim$par)), 0, 100, col = "blue", add = TRUE)

# 
res_optim <- optim(par = estimate_0_sigmoid(contribs_xgboost$contrib_grid$lstat, "lstat"), 
                   fn = loss_fun,
                   contribs = contribs_xgboost$contrib_grid$lstat,
                   FUN = sigmoid_fun,
                   method = "L-BFGS-B")

contribs_xgboost$contrib_grid$lstat %>% plot(type = "l")
curve(sigmoid_fun(x, as.list(res_optim$par)), 0, 100, col = "blue", add = TRUE)

# 
res_optim <- optim(par = estimate_0_sigmoid(contribs_xgboost$contrib_grid$rm, "rm"), 
                   fn = loss_fun,
                   contribs = contribs_xgboost$contrib_grid$rm,
                   FUN = sigmoid_fun,
                   method = "L-BFGS-B")

contribs_xgboost$contrib_grid$rm %>% plot(type = "l")
curve(sigmoid_fun(x, as.list(res_optim$par)), 0, 100, col = "blue", add = TRUE)

## CURVE S
estimate_0_curve_s(contribs_xgboost$contrib_grid$zn, "zn")
res_optim <- optim(par = estimate_0_curve_s(contribs_xgboost$contrib_grid$zn, "zn"), 
                   fn = loss_fun,
                   contribs = contribs_xgboost$contrib_grid$zn,
                   FUN = curve_s_fun,
                   method = "L-BFGS-B")

contribs_xgboost$contrib_grid$zn %>% plot(type = "l")
curve(curve_s_fun(x, as.list(res_optim$par)), 0, 100, col = "blue", add = TRUE)

# 
res_optim <- optim(par = estimate_0_curve_s(contribs_xgboost$contrib_grid$lstat, "lstat"), 
                   fn = loss_fun,
                   contribs = contribs_xgboost$contrib_grid$lstat,
                   FUN = curve_s_fun,
                   method = "L-BFGS-B")

contribs_xgboost$contrib_grid$lstat %>% plot(type = "l")
curve(curve_s_fun(x, as.list(res_optim$par)), 0, 100, col = "blue", add = TRUE)

# 
res_optim <- optim(par = estimate_0_curve_s(contribs_xgboost$contrib_grid$rm, "rm"), 
                   fn = loss_fun,
                   contribs = contribs_xgboost$contrib_grid$rm,
                   FUN = curve_s_fun,
                   method = "L-BFGS-B")

contribs_xgboost$contrib_grid$rm %>% plot(type = "l")
curve(curve_s_fun(x, as.list(res_optim$par)), 0, 100, col = "blue", add = TRUE)

## TANH
estimate_0_curve_tanh(contribs_xgboost$contrib_grid$zn, "zn")
res_optim <- optim(par = estimate_0_curve_tanh(contribs_xgboost$contrib_grid$zn, "zn"), 
                   fn = loss_fun,
                   contribs = contribs_xgboost$contrib_grid$zn,
                   FUN = curve_tanh_fun,
                   method = "L-BFGS-B")

contribs_xgboost$contrib_grid$zn %>% plot(type = "l")
curve(curve_tanh_fun(x, as.list(res_optim$par)), 0, 100, col = "blue", add = TRUE)

# 
res_optim <- optim(par = estimate_0_curve_tanh(contribs_xgboost$contrib_grid$lstat, "lstat"), 
                   fn = loss_fun,
                   contribs = contribs_xgboost$contrib_grid$lstat,
                   FUN = curve_tanh_fun,
                   method = "L-BFGS-B")

contribs_xgboost$contrib_grid$lstat %>% plot(type = "l")
curve(curve_tanh_fun(x, as.list(res_optim$par)), 0, 100, col = "blue", add = TRUE)

# 
res_optim <- optim(par = estimate_0_curve_tanh(contribs_xgboost$contrib_grid$rm, "rm"), 
                   fn = loss_fun,
                   contribs = contribs_xgboost$contrib_grid$rm,
                   FUN = curve_tanh_fun,
                   method = "L-BFGS-B")

contribs_xgboost$contrib_grid$rm %>% plot(type = "l")
curve(curve_tanh_fun(x, as.list(res_optim$par)), 0, 100, col = "blue", add = TRUE)

## ALL

res_optim <- fit_contrib_curve(contribs_xgboost$contrib_grid$zn, "zn")
res_optim$curve_type

contribs_xgboost$contrib_grid$zn %>% 
  plot(type = "l",  main = paste("zn -", res_optim$curve_type))
curve(res_optim$fit_fun(x, as.list(res_optim$pars)), 0, 100, col = "blue", add = TRUE)

##

res_optim <- fit_contrib_curve(contribs_xgboost$contrib_grid$lstat, "lstat")
res_optim$curve_type

contribs_xgboost$contrib_grid$lstat %>% 
  plot(type = "l",  main = paste("lstat -", res_optim$curve_type))
curve(res_optim$fit_fun(x, as.list(res_optim$pars)), 0, 100, col = "blue", add = TRUE)

##

res_optim <- fit_contrib_curve(contribs_xgboost$contrib_grid$rm, "rm")
res_optim$curve_type

contribs_xgboost$contrib_grid$rm%>% 
  plot(type = "l",  main = paste("rm -", res_optim$curve_type))
curve(res_optim$fit_fun(x, as.list(res_optim$pars)), 0, 100, col = "blue", add = TRUE)

##

res_optim <- fit_contrib_curve(contribs_xgboost$contrib_grid$black, "black")
res_optim$curve_type

contribs_xgboost$contrib_grid$black %>% 
  plot(type = "l",  main = paste("black -", res_optim$curve_type))
curve(res_optim$fit_fun(x, as.list(res_optim$pars)), 0, 400, col = "blue", add = TRUE)

##

res_optim <- fit_contrib_curve(contribs_xgboost$contrib_grid$indus, "indus")
res_optim$curve_type

contribs_xgboost$contrib_grid$indus %>% 
  plot(type = "l",  main = paste("indus -", res_optim$curve_type))
curve(res_optim$fit_fun(x, as.list(res_optim$pars)), 0, 400, col = "blue", add = TRUE)


#####
plot(Boston$zn, contribs_xgboost$contribs$zn)
curve(contribs_xgboost$smooth_contribs$zn$fit_fun(x, as.list(contribs_xgboost$smooth_contribs$zn$pars)), 
      0, 400, col = "blue", add = TRUE)

plot(Boston$lstat, contribs_xgboost$contribs$lstat)
curve(contribs_xgboost$smooth_contribs$lstat$fit_fun(x, as.list(contribs_xgboost$smooth_contribs$lstat$pars)), 
      0, 400, col = "blue", add = TRUE)

plot(Boston$rm, contribs_xgboost$contribs$rm)
curve(contribs_xgboost$smooth_contribs$rm$fit_fun(x, as.list(contribs_xgboost$smooth_contribs$rm$pars)), 
      0, 400, col = "blue", add = TRUE)

plot(Boston$black, contribs_xgboost$contribs$black)
curve(contribs_xgboost$smooth_contribs$black$fit_fun(x, as.list(contribs_xgboost$smooth_contribs$black$pars)), 
      0, 400, col = "blue", add = TRUE)

plot(Boston$indus, contribs_xgboost$contribs$indus)
curve(contribs_xgboost$smooth_contribs$indus$fit_fun(x, as.list(contribs_xgboost$smooth_contribs$indus$pars)), 
      0, 400, col = "blue", add = TRUE)



plot(Boston$lstat, contribs_xgboost$contribs$lstat, type = "p")
lines(contribs_xgboost$contrib_grid$lstat, type = "l", lwd = 2, col = "red")
curve(contribs_xgboost$smooth_contribs$lstat$fit_fun(x, as.list(contribs_xgboost$smooth_contribs$lstat$pars)), 
      0, 400, col = "blue", add = TRUE)

incr <- mean((contribs_xgboost$contribs$lstat - contribs_xgboost$contrib_funs$lstat(Boston$lstat)))
plot(contribs_xgboost$contrib_grid$lstat %>% as_tibble() %>% mutate(yhat = yhat+incr), 
     type = "l", lwd= 2, cex = 2, col = "red", ylim = c(-30, -5))
lines(Boston$lstat, contribs_xgboost$contribs$lstat, type = "p")
curve(contribs_xgboost$smooth_contribs$lstat$fit_fun(x, as.list(contribs_xgboost$smooth_contribs$lstat$pars)), 
        +       0, 400, col = "blue", add = TRUE)

