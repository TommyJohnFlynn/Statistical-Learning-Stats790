if(!is.null(dev.list())) dev.off()
rm(list = ls())
cat('\014')

####################
# Data Description # 
####################

# libraries
library(tidymodels)
library(ggcorrplot)
library(vip)
library(pdp)

# seed
set.seed(13)

# load data 
bfat <- read.delim('bodyfat.txt')

########################
# Exploratory Analysis #
########################

# missing values  
any(is.na(bfat))

# summary of data 
summary(bfat)

# box plots 
par(mar = c(7, 5, 5, 3))
boxplot(bfat[6:16], col='lightblue', las=2, ylab='Circumference (cm)',
        main='Box-plots of major body parts', cex.axis=1.5, cex.main=1.5,
        cex.lab=1.5)

# correlation 
ggcorrplot(cor(bfat), type='lower', hc.order=TRUE, lab=TRUE, lab_size=4, 
           colors=c('lightblue', 'white', 'darkorange'), tl.cex=17,
           title='Correlation plot')

# check abdomen vs waist 
bfat$Abdomen / bfat$Waist 

#################
# Preprocessing #
#################

# drop waist and density columns 
bfat <- bfat[, -c(1, 9)]

# drop suspect points 
bfat <- bfat[-c(170, 180, 36, 214), ]

# train/test split
data.split <- initial_split(bfat, prop=0.75)
train <- training(data.split)
test <- testing(data.split)

# 10-fold cross-validation 
cv <- vfold_cv(train, v=10)

################
# Model Choice #
################

# model recipe
model.rec <- recipe(Pct.BF ~ ., data=train) 

# random forest model
model.rf <- rand_forest(
  mtry = tune(),
  trees = tune(),
  min_n = tune()) %>%
  set_engine('ranger', importance = "permutation") %>%
  set_mode('regression')

wf <- workflow() %>%
  add_recipe(model.rec) %>%
  add_model(model.rf) 

################
# Model Tuning #
################

# hyperparameters 
hp.param <- parameters(
  finalize(mtry(), train),
  trees(),
  min_n())

# hyperparameter space 
hp.grid <- grid_max_entropy(
  hp.param, 
  size = 100)

# tune using hp space and cv folds 
hp.tune <- tune_grid(
  object = wf,
  resamples = cv,
  grid = hp.grid,
  metrics = metric_set(rmse),
  control = control_grid(verbose=TRUE))

# display top 5
hp.tune %>% show_best(metric='rmse')

##############
# Best Model #
##############

# select best hyperparameters 
hp.best <- hp.tune %>% select_best("rmse")

# create best model
model.best <- model.rf %>% finalize_model(hp.best)

wf.best <- workflow() %>%
  add_recipe(model.rec) %>%
  add_model(model.best) %>%
  fit(train)

####################
# Model Evaluation #
####################

# RMSE on test set 
predictions <- as.vector(unlist(predict(wf.best, test)))
residuals <- test$Pct.BF - predictions 
RMSE <- sqrt(mean(residuals^2))
RMSE

# variable importance 
wf.best %>% 
  extract_fit_parsnip() %>% 
  vip()

# partial dependence 
wf.best %>% 
  extract_fit_parsnip() %>% 
  partial(pred.var = "Abdomen", plot = TRUE, train = train) 
