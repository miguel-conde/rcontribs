source("Hierarchichal Models/contribs_lmer.R", encoding = "utf8")


# MODELS -----------------------------------------------------------------
dt <- read.table("http://bayes.acs.unt.edu:8083/BayesContent/class/Jon/R_SC/Module9/lmm.data.txt",
                 header=TRUE, sep=",", na.strings="NA", dec=".", strip.white=TRUE) %>% 
  as_tibble() 
# data was previously publicly available from
# http://researchsupport.unt.edu/class/Jon/R_SC/Module9/lmm.data.txt
# but the link is now broken

###
dt$classID <- paste(dt$school, dt$class, sep=".")
dt <- dt %>% 
  mutate_if(is.character, ~ factor(.))

# Intercept only ----------------------------------------------------------

m0 <- lmer(extro ~ 1 + (1 | school/classID), data = dt)
summary(m0)

# Varying Intercept only --------------------------------------------------
m2 <- lmer(extro ~ open + agree + social + (1 | school/classID), data = dt)
summary(m2)

m3 <- lmer(extro ~ open + agree + social + (1 | school) + (1 |classID), data = dt)
summary(m3)

# Varying intercept & slope -----------------------------------------------

m_i_s <- lmer(extro ~ open + agree + social + (open + 1 | school/classID), data = dt)
summary(m_i_s)
confint(m_i_s)

# Con interacción
m_int <- lmer(extro ~ open + open:school + agree + social + (open + 1 | school/classID), data = dt)


# TESTS -------------------------------------------------------------------


# Intercept only ----------------------------------------------------------
get_mlmer_contribs2(m0, pred = TRUE)

predict(m0) %>% head()

####

# Con FE, sin RE
get_mlmer_contribs(m0, re.form = NA)

# Con FE, solo RE del modelo a mayor nivel, school
get_mlmer_contribs(m0, re.form = ~ (1 | school))

# Con FE, solo RE del modelo del nivel inferior anidado, classID dentro de school
get_mlmer_contribs(m0, re.form = ~ (1 | classID:school))

# Con FE, todos los RE
get_mlmer_contribs(m0)
get_mlmer_contribs(m0, re.form = NULL)
get_mlmer_contribs(m0, re.form = ~ (1 | school/classID))

####

# Sin FE, sin RE
get_mlmer_contribs(m0, re.form = NA, random.only = TRUE)

# Sin FE, solo RE del modelo a mayor nivel, school
get_mlmer_contribs(m0, re.form = ~ (1 | school), random.only = TRUE)

# Sin FE, solo RE del modelo del nivel inferior anidado, classID dentro de school
get_mlmer_contribs(m0, re.form = ~ (1 | classID:school), random.only = TRUE)

# Sin FE, todos los RE
get_mlmer_contribs(m0, random.only = TRUE)
get_mlmer_contribs(m0, re.form = NULL, random.only = TRUE)
get_mlmer_contribs(m0, re.form = ~ (1 | school/classID), random.only = TRUE)

# Intercept & slope -------------------------------------------------------

get_mlmer_contribs(m_i_s, pred = TRUE)

predict(m_i_s) %>% head()

####

# Con FE, sin RE
get_mlmer_contribs(m_i_s, re.form = NA)

# Con FE, solo RE del modelo a mayor nivel, school
get_mlmer_contribs(m_i_s, re.form = ~ (open + 1 | school))

# Con FE, solo RE del modelo del nivel inferior anidado, classID dentro de school
get_mlmer_contribs(m_i_s, re.form = ~ (open + 1 | classID:school))

# Con FE, todos los RE
get_mlmer_contribs(m_i_s)
get_mlmer_contribs(m_i_s, re.form = NULL)
get_mlmer_contribs(m_i_s, re.form = ~ (open + 1 | school/classID))

####

# Sin FE, sin RE
get_mlmer_contribs(m_i_s, re.form = NA, random.only = TRUE)

# Sin FE, solo RE del modelo a mayor nivel, school
get_mlmer_contribs(m_i_s, re.form = ~ (open + 1 | school), random.only = TRUE)

# Sin FE, solo RE del modelo del nivel inferior anidado, classID dentro de school
get_mlmer_contribs(m_i_s, re.form = ~ (open + 1 | classID:school), random.only = TRUE)

# Sin FE, todos los RE
get_mlmer_contribs(m_i_s, random.only = TRUE)
get_mlmer_contribs(m_i_s, re.form = NULL, random.only = TRUE)
get_mlmer_contribs(m_i_s, re.form = ~ (open + 1 | school/classID), random.only = TRUE)


# Interacciones -----------------------------------------------------------

get_mlmer_contribs(m_int, pred = TRUE)

predict(m_int) %>% head()

####

# Con FE, sin RE
get_mlmer_contribs(m_int, re.form = NA)

# Con FE, solo RE del modelo a mayor nivel, school
get_mlmer_contribs(m_int, re.form = ~ (open + 1 | school))

# Con FE, solo RE del modelo del nivel inferior anidado, classID dentro de school
get_mlmer_contribs(m_int, re.form = ~ (open + 1 | classID:school))

# Con FE, todos los RE
get_mlmer_contribs(m_int)
get_mlmer_contribs(m_int, re.form = NULL)
get_mlmer_contribs(m_int, re.form = ~ (open + 1 | school/classID))

####

# Sin FE, sin RE
get_mlmer_contribs(m_int, re.form = NA, random.only = TRUE)

# Sin FE, solo RE del modelo a mayor nivel, school
get_mlmer_contribs(m_int, re.form = ~ (open + 1 | school), random.only = TRUE)

# Sin FE, solo RE del modelo del nivel inferior anidado, classID dentro de school
get_mlmer_contribs(m_int, re.form = ~ (open + 1 | classID:school), random.only = TRUE)

# Sin FE, todos los RE
get_mlmer_contribs(m_int, random.only = TRUE)
get_mlmer_contribs(m_int, re.form = NULL, random.only = TRUE)
get_mlmer_contribs(m_int, re.form = ~ (open + 1 | school/classID), random.only = TRUE)


# PDPs CONTRIBS -----------------------------------------------------------

source("interpretability/utils_pdp_contrib.R")

pdp_pred_lmer <- function(object, newdata)  {
  results <- as.vector(predict(object, newdata))
  return(results)
}

contribs_lmer <- pdp_contribs(m_int, dt, 
                              c("open", "agree", "social"), 
                              pdp_pred_lmer)

rois_all_lmer <- pdp_rois(contribs_lmer, dt)

# Contributions -----------------------------------------------------------

predictions_lmer <- predict(m_int, dt) 

contribs_lmer$contribs

probe_lmer <- contribs_lmer$contribs %>% mutate(y_hat = rowSums(.))
probe_lmer

plot(predictions_lmer, probe_lmer$y_hat)
abline(a= 0, b = 1)

probe_lmer %>% summarise_all(~ sum(.))
sum(predictions_lmer)

# Plot contributions ------------------------------------------------------

ggplot_1_contrib(contribs_lmer, "open", 
                 x_units = "\n(lower status of the population (percent))",
                 y_units = "Contributions\n(to the median value of owner-occupied homes in $1000s)")


p_lmer <- ggplot_contribs(contribs_lmer, y_units = "$1000") + 
  plot_annotation(
    title = "Boston - Generalised Linear Model Contributions",
    subtitle = "(to the median value of owner-occupied homes in $1000s)",
    caption = "Source: mpg dataset in ggplot2"
  )

p_lmer

# ROIs --------------------------------------------------------------------

pdp_1_roi(dt$open, contribs_lmer, "open") 


# Plot ROIs ---------------------------------------------------------------

ggplot_1_roi(contribs_lmer, "open", dt,
             x_units = "\nlower status of the population (%)",
             y_units = "ROI ($1000 / 1%)")

ggplot_rois(contribs_lmer, dt, y_units = "ROI") + 
  plot_annotation(
    title = "Boston - GLM Model ROIs",
    # subtitle = "(to the median value of owner-occupied homes in $1000s)",
    caption = "Source: mpg dataset in ggplot2")


# Check -------------------------------------------------------------------

get_mlmer_contribs(m_int, pred = TRUE)
contribs_lmer$contribs %>% mutate(pred = rowSums(contribs_lmer$contribs))
