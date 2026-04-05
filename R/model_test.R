library(tidymodels)
library(readr)
library(broom.mixed)
library(dotwhisker)

urchins <- read_csv("https://tidymodels.org/start/models/urchins.csv") %>%
  setNames(c("food_regime", "initial_volume", "width")) %>%
  mutate(food_regime = factor(food_regime, levels = c("Initial", "Low", "High")))

ggplot(urchins,
       aes(x = initial_volume, 
           y = width, 
           group = food_regime, 
           col = food_regime)) + 
  geom_point() + 
  geom_smooth(method = lm, se = FALSE) +
  scale_color_viridis_d(option = "plasma", end = .7)



library(tidymodels)  
library(readr)
library(vip)


library(tidymodels)
library(readr)
library(tidyverse)

hotels <- read_csv("https://tidymodels.org/start/case-study/hotels.csv") %>%
  mutate(across(where(is.character), as.factor))

dim(hotels)

hotels %>% 
  count(children) %>% 
  mutate(prop = n/sum(n))


set.seed(123)
splits  <- initial_split(hotels, strata = children)

hotel_other <- training(splits)
hotel_test  <- testing(splits)

hotel_other %>% 
  count(children) %>% 
  mutate(prop = n/sum(n))

hotel_test  %>% 
  count(children) %>% 
  mutate(prop = n/sum(n))

set.seed(234)
val_set <- validation_split(hotel_other, 
                            strata = children, 
                            prop = 0.80)


library(nycflights13)
library(tidyverse)
library(skimr)  
set.seed(123)

flight_data <- nycflights13::flights %>% 
  mutate(
    arr_delay = ifelse(arr_delay >= 30, "late", "on_time"),
    arr_delay = factor(arr_delay),
    date = lubridate::as_date(time_hour)
  ) %>% 
  inner_join(weather, by = c("origin", "time_hour")) %>% 
  dplyr::select(dep_time, flight, origin, dest, air_time, distance, 
         carrier, date, arr_delay, time_hour) %>% 
  na.omit() %>% 
  mutate_if(is.character, as.factor)

glimpse(flight_data)


set.seed(226)
data_split <- initial_split(flight_data, prop = 3/4)
train_data <- training(data_split)
test_data  <- testing(data_split)

## 创建recipe
flights_rec <- recipe(arr_delay ~ ., data = train_data) %>% 
  update_role(flight, time_hour, new_role = "ID") %>% 
  step_date(date, features = c("dow", "month")) %>%               
  step_holiday(
    date, 
    holidays = timeDate::listHolidays("US"), 
    keep_original_cols = FALSE
    ) %>% 
  step_dummy(all_nominal_predictors()) %>% 
  step_zv(all_predictors())

summary(flights_rec)
lr_mod <- logistic_reg() %>% 
  set_engine("glm")

flights_wflow <- 
  workflow() %>% 
  add_model(lr_mod) %>% 
  add_recipe(flights_rec)

flights_fit <- flights_wflow %>% 
  fit(data = train_data)

predict(flights_fit, test_data)

flights_aug <- augment(flights_fit, test_data)

# The data look like: 
flights_aug %>% dplyr::select(arr_delay, time_hour, flight, .pred_class, .pred_on_time)



select_best_vargroup <- list()
library(tidymodels)
library(rpart.plot)
library(vip)
library(MASS)
library(broom)
library(dplyr)
data("Vehicle", package = "mlbench")
Vehicle
Vehicle <- Vehicle %>% 
  mutate(Class = factor(Class))
all_data <- Vehicle
split <- initial_split(Vehicle,strata = Class,prop = 0.7)
train_data <- training(split)
test_data <- testing(split)
rec <- recipe(Class ~ .,data = train_data) %>% 
  step_normalize(all_numeric_predictors())




### lasso 
lasso_spec <- multinom_reg(
  mode = "classification",
  engine = "glmnet",
  penalty = tune(),
  mixture = 1
) %>% 
  set_engine("glmnet")

folds <- vfold_cv(train_data,v = 5,strata = Class)

wf <- workflow() %>% 
  add_recipe(rec) %>% 
  add_model(lasso_spec)

lasso_grid <- grid_regular(
  penalty(),
  levels = 50
)

tune_res <- tune_grid(
  object = wf,
  resamples = folds,
  grid = lasso_grid,
  metrics = metric_set(accuracy)
)

best_lasso <- select_best(
  tune_res,
  metric = "accuracy"
)

final_wf <- finalize_workflow(wf,best_lasso)
final_fit <- fit(final_wf,data = train_data)

lasso_coef <- extract_fit_parsnip(final_fit) %>% 
  pluck("fit") %>% 
  coef(s = best_lasso$penalty)




# 遍历列表，提取所有非零系数
all_features <- purrr::map_dfr(lasso_coef, function(coef_mat) {
  coef_mat %>%
    as.matrix() %>%
    as.data.frame() %>%
    tibble::rownames_to_column("feature") %>%
    dplyr::rename(coef = 2) %>%
    filter(
      coef != 0,
      feature != "(Intercept)"
    )
}, .id = "class")

selected_features <- all_features %>%
  distinct(feature) %>%
  dplyr::pull(feature)

LASSO <- selected_features

select_best_vargroup$LASSO <- LASSO
message("LASSO"," vars selectd!")



### 随机森林筛选
rand_forest_spec <- rand_forest(
  mode = "classification",
  engine = "ranger",
  mtry = tune(),
  trees = 1000,
  min_n = tune()
) %>% 
  set_engine("ranger",importance = "impurity")


folds <- vfold_cv(train_data,v = 5,strata = Class)
wf <- workflow() %>% 
  add_recipe(rec) %>% 
  add_model(rand_forest_spec)

rf_grid <- grid_regular(
  mtry(range = c(1,length(train_data) - 1)),
  min_n(range = c(2,10)),
  levels = 5
)

tune_res <- tune_grid(
  object = wf,
  resamples = folds,
  grid = rf_grid,
  metrics = metric_set(roc_auc,accuracy)
)

best_rf <- select_best(
  tune_res,
  metric = "accuracy"
)

final_wf <- finalize_workflow(wf,best_rf)
final_fit <- fit(final_wf,data = train_data)

rf_vi <- final_fit %>% 
  extract_fit_parsnip() %>% 
  vip::vi()

rf_vi_cum <- rf_vi %>%
  arrange(desc(Importance)) %>%
  mutate(
    Importance_pct = Importance / sum(Importance),
    cum_pct = cumsum(Importance_pct)
  )

selected_features <- rf_vi_cum %>%
  filter(cum_pct <= 0.95) %>%
  pull(Variable)

RF <- selected_features

select_best_vargroup$RF <- RF
message("RF"," vars selectd!")



### xgboost
class_levels <- unique(train_data$Class)
xgb_spec <- boost_tree(
  mode = "classification",
  engine = "xgboost",
  trees = 1000,
  min_n = tune(),
  tree_depth = tune(),
  learn_rate = tune()
) %>% 
  set_engine("xgboost",objective = "multi:softprob")


folds <- vfold_cv(train_data,v = 5,strata = Class)
wf <- workflow() %>% 
  add_recipe(rec) %>% 
  add_model(xgb_spec)

xgb_grid <- grid_regular(
  tree_depth(range = c(2,8)),
  learn_rate(range = c(0.01,0.3)),
  min_n(range = c(2,10)),
  levels = 3
)

tune_res <- tune_grid(
  object = wf,
  resamples = folds,
  grid = xgb_grid,
  metrics = metric_set(roc_auc,accuracy)
)

best_xgboost <- select_best(
  tune_res,
  metric = "accuracy"
)

final_wf <- finalize_workflow(wf,best_xgboost)
final_fit <- fit(final_wf,data = train_data)

xgboost_vi <- final_fit %>% 
  extract_fit_parsnip() %>% 
  vip::vi()

xgboost_vi_cum <- rf_vi %>%
  arrange(desc(Importance)) %>%
  mutate(
    Importance_pct = Importance / sum(Importance),
    cum_pct = cumsum(Importance_pct)
  )

selected_features <- xgboost_vi_cum %>%
  filter(cum_pct <= 0.95) %>%
  pull(Variable)


XGBOOST <- selected_features
select_best_vargroup$XGBOOST <- XGBOOST
message("XGBOOST"," vars selectd!")


### Stepglm both
rec <- rec %>%
  recipes::prep()

train_baked <- recipes::bake(rec, new_data = NULL)

full_model <- nnet::multinom(
  Class ~ .,
  data = train_baked,
  trace = FALSE
)


null_model <- nnet::multinom(
  Class ~ 1,
  data = train_baked,
  trace = FALSE
)


step_forward <- MASS::stepAIC(
  null_model,
  direction = "forward",
  scope = list(lower = ~1, upper = formula(full_model)),
  trace = FALSE
)
message("stepglm[forward]"," vars selectd!")



step_backward <- MASS::stepAIC(
  full_model,
  direction = "backward",
  trace = FALSE
)
message("stepglm[backward]"," vars selectd!")

step_both <- MASS::stepAIC(
  full_model,
  direction = "both",
  trace = FALSE
)
message("stepglm[both]"," vars selectd!")


get_features <- function(model) {
  vars <- all.vars(formula(model))[-1]
  unique(vars)
}

feat_forward  <- get_features(step_forward)
feat_backward <- get_features(step_backward)
feat_both     <- get_features(step_both)

`Stepglm[forward]` <- feat_forward
`Stepglm[backward]` <- feat_backward
`Stepglm[both]` <- feat_both

select_best_vargroup$`Stepglm[forward]` <- feat_forward
select_best_vargroup$`Stepglm[backward]` <- feat_backward
select_best_vargroup$`Stepglm[both]` <- feat_both

return(select_best_vargroup)








select_best_vargroup <- list()
library(tidymodels)
library(rpart.plot)
library(vip)
library(MASS)
library(broom)
library(dplyr)
data("Vehicle", package = "mlbench")
Vehicle
Vehicle <- Vehicle %>% 
  mutate(Class = factor(Class))
all_data <- Vehicle
split <- initial_split(Vehicle,strata = Class,prop = 0.7)
train_data <- training(split)
test_data <- testing(split)
rec <- recipe(Class ~ .,data = train_data) %>% 
  step_normalize(all_numeric_predictors())

folds <- vfold_cv(train_data, v = 5, strata = Class)

lda_rec <- recipe(Class ~ ., data = train_data) %>%
  step_normalize(all_numeric_predictors())

library(discrim)
lda_mod <- discrim_linear() %>%
  set_engine("MASS") %>%
  set_mode("classification")


lda_wf <- workflow() %>%
  add_recipe(lda_rec) %>%
  add_model(lda_mod)

bal_metric <- metric_set(bal_accuracy)
lda_cv <- fit_resamples(lda_wf,resamples = folds,metrics = bal_metric)
lda_final <- last_fit(lda_wf, split, metrics = bal_metric)
final_model <- extract_workflow(lda_final)



### svm
svm_rec <- recipe(Class ~ ., data = train_data) %>%
  step_normalize(all_numeric_predictors())

svm_mod <- svm_rbf(cost = 1, rbf_sigma = 0.1) %>%
  set_engine("kernlab") %>%
  set_mode("classification")

svm_wf <- workflow() %>%
  add_recipe(svm_rec) %>%
  add_model(svm_mod)

# 定义评估指标：仅 Balanced Accuracy
bal_metric <- metric_set(bal_accuracy)

# 执行交叉验证
svm_cv <- fit_resamples(
  object = svm_wf,
  resamples = folds,
  metrics = bal_metric
)

# 在完整训练集上训练，并在测试集上评估
svm_final <- last_fit(svm_wf, split, metrics = bal_metric)
final_svm_model <- extract_workflow(svm_final)





### 随机森林
rf_rec <- recipe(Class ~ ., data = train_data) %>%
  step_normalize(all_numeric_predictors())

rf_mod <- rand_forest(
  trees = 1000,
  mtry = tune(),
  min_n = tune()
) %>%
  set_engine("ranger", importance = "impurity") %>%
  set_mode("classification")


rf_wf <- workflow() %>%
  add_recipe(rf_rec) %>%
  add_model(rf_mod)

rf_grid <- grid_regular(
  mtry(range = c(2, 4)),  # 根据你的特征总数调整范围
  min_n(range = c(5, 20)),
  levels = 5
)

bal_metric <- metric_set(bal_accuracy)

set.seed(2023)
rf_tune <- tune_grid(
  object = rf_wf,
  resamples = folds,
  grid = rf_grid,
  metrics = bal_metric
)
best_rf <- select_best(rf_tune, metric = "bal_accuracy")
final_rf_wf <- finalize_workflow(rf_wf, best_rf)
rf_final <- last_fit(final_rf_wf, split, metrics = bal_metric)
final_rf_model <- extract_workflow(rf_final)



### mlp
mlp_rec <- recipe(Class ~ ., data = train_data) %>%
  step_normalize(all_numeric_predictors())

mlp_mod <- mlp(
  hidden_units = tune(),    # 隐藏层神经元数量，设为待调优
  penalty = tune(),         # L2 正则化系数，设为待调优
  epochs = 100              # 训练轮数
) %>%
  set_engine("nnet") %>%
  set_mode("classification")

mlp_wf <- workflow() %>%
  add_recipe(mlp_rec) %>%
  add_model(mlp_mod)

mlp_grid <- grid_regular(
  hidden_units(range = c(5, 20)),
  penalty(range = c(0.01, 0.1)),
  levels = 5
)
bal_metric <- metric_set(bal_accuracy)

# 执行网格搜索
set.seed(2023)
mlp_tune <- tune_grid(
  object = mlp_wf,
  resamples = folds,
  grid = mlp_grid,
  metrics = bal_metric
)

# 查看调优结果
collect_metrics(mlp_tune)
best_mlp <- select_best(mlp_tune, metric = "bal_accuracy")
final_mlp_wf <- finalize_workflow(mlp_wf, best_mlp)
mlp_final <- last_fit(final_mlp_wf, split, metrics = bal_metric)
final_mlp_model <- extract_workflow(mlp_final)




### gbm
gbm_rec <- recipe(Class ~ ., data = train_data)
gbm_mod <- boost_tree(
  trees = 1000,
  tree_depth = tune(),
  min_n = tune(),
  learn_rate = tune()
) %>%
  set_engine("xgboost") %>%  # 最强最快的 GBM 引擎
  set_mode("classification")

gbm_wf <- workflow() %>%
  add_recipe(gbm_rec) %>%
  add_model(gbm_mod)

bal_metric <- metric_set(bal_accuracy)
gbm_grid <- grid_regular(
  tree_depth(range = c(3, 8)),
  min_n(range = c(5, 20)),
  learn_rate(range = c(0.01, 0.3)),
  levels = 3
)

set.seed(2023)
gbm_tune <- tune_grid(
  gbm_wf,
  resamples = folds,
  grid = gbm_grid,
  metrics = bal_metric
)


best_gbm <- select_best(gbm_tune, metric = "bal_accuracy")
final_gbm_wf <- finalize_workflow(gbm_wf, best_gbm)
gbm_final <- last_fit(final_gbm_wf, split, metrics = bal_metric)
final_gbm_model <- extract_workflow(gbm_final)




## 
nb_rec <- recipe(Class ~ ., data = train_data)
nb_mod <- naive_Bayes() %>%
  set_engine("klaR") %>%
  set_mode("classification")

nb_wf <- workflow() %>%
  add_recipe(nb_rec) %>%
  add_model(nb_mod)

bal_metric <- metric_set(bal_accuracy)
nb_cv <- fit_resamples(
  nb_wf,
  resamples = folds,
  metrics = bal_metric
)

collect_metrics(nb_cv)
nb_final <- last_fit(nb_wf, split, metrics = bal_metric)
final_nb_model <- extract_workflow(nb_final)



### xgboost
xgb_rec <- recipe(Class ~ ., data = train_data)
xgb_mod <- boost_tree(
  trees = 1000,
  tree_depth = tune(),
  min_n = tune(),
  learn_rate = tune(),
  loss_reduction = tune()
) %>%
  set_engine("xgboost") %>%
  set_mode("classification")

xgb_wf <- workflow() %>%
  add_recipe(xgb_rec) %>%
  add_model(xgb_mod)

bal_metric <- metric_set(bal_accuracy)

xgb_grid <- grid_regular(
  tree_depth(range = c(3, 8)),
  min_n(range = c(5, 20)),
  learn_rate(range = c(0.01, 0.3)),
  loss_reduction(range = c(0, 1)),
  levels = 3
)

set.seed(2023)
xgb_tune <- tune_grid(xgb_wf,resamples = folds,grid = xgb_grid,metrics = bal_metric)
best_xgb <- select_best(xgb_tune, metric = "bal_accuracy")
final_xgb_wf <- finalize_workflow(xgb_wf, best_xgb)
xgb_final <- last_fit(final_xgb_wf, split, metrics = bal_metric)
final_xgb_model <- extract_workflow(xgb_final)


### ridge
ridge_rec <- recipe(Class ~ ., data = train_data) %>%
  step_normalize(all_numeric_predictors())


ridge_mod <- multinom_reg(penalty = tune(), mixture = 0) %>%
  set_engine("glmnet") %>%
  set_mode("classification")


ridge_wf <- workflow() %>%
  add_recipe(ridge_rec) %>%
  add_model(ridge_mod)

bal_metric <- metric_set(bal_accuracy)
ridge_grid <- grid_regular(penalty(range = c(-3, 4)),levels = 10)
set.seed(1234)
ridge_tune <- tune_grid(ridge_wf,resamples = folds,grid = ridge_grid,metrics = bal_metric)
best_ridge <- select_best(ridge_tune, metric = "bal_accuracy")
final_ridge_wf <- finalize_workflow(ridge_wf, best_ridge)
ridge_final <- last_fit(final_ridge_wf, split, metrics = bal_metric)
final_ridge_model <- extract_workflow(ridge_final)


### lasso
lasso_rec <- recipe(Class ~ ., data = train_data) %>%
  step_normalize(all_numeric_predictors())

lasso_mod <- multinom_reg(
  penalty = tune(), 
  mixture = 1
  ) %>%
  set_engine("glmnet") %>%
  set_mode("classification")

# 工作流
lasso_wf <- workflow() %>%
  add_recipe(lasso_rec) %>%
  add_model(lasso_mod)

# 调优
lasso_grid <- grid_regular(penalty(range = c(-3, 4)),levels = 10)
set.seed(1234)
lasso_tune <- tune_grid(lasso_wf,resamples = folds,grid = lasso_grid,metrics = bal_metric)

best_lasso <- select_best(lasso_tune, metric = "bal_accuracy")
final_lasso_wf <- finalize_workflow(lasso_wf, best_lasso)
final_lasso_model <- fit(final_lasso_wf,data = train_data)

### Enet
en_rec <- recipe(Class ~ ., data = train_data) %>%
  step_normalize(all_numeric_predictors())

en_mod <- multinom_reg(
  penalty = tune(),
  mixture = tune()
) %>%
  set_engine("glmnet") %>%
  set_mode("classification")

en_wf <- workflow() %>%
  add_recipe(en_rec) %>%
  add_model(en_mod)

en_grid <- grid_regular(penalty(range = c(-3, 4)),mixture(range = c(0.1, 0.9)),levels = c(10, 9))
set.seed(1234)
en_tune <- tune_grid(en_wf,resamples = folds,grid = en_grid,metrics = bal_metric)
best_en <- select_best(en_tune, metric = "bal_accuracy")
final_en_wf <- finalize_workflow(en_wf, best_en)
final_en_model <- fit(final_en_wf,data = train_data)

en_final <- last_fit(final_en_wf, split, metrics = bal_metric)
lasso_final <- last_fit(object = final_lasso_wf, split = split, metrics = bal_metric)
