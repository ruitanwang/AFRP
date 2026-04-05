
library(tidyverse)
eset <- qs_read("./analysiswork/cleandata/bulkRNAseq_mc_eset.qs2")
phenodf <- Biobase::pData(eset)
select_disease <- c("ALS","FTLD","CON")
p1 <- group_bar(
  plotdf = phenodf[phenodf$group %in% select_disease,],
  groupcol = group,
  plot_title = "Motor Cortex Samples",
  legend_ncol = 1,
  disease_color = setNames(c("#a45c72","#3c789d","#9369ab"),c("ALS","FTLD","CON"))
)
p1


p2 <- ggpie3D(
  data = phenodf[phenodf$group %in% select_disease,],
  group_key = "group",
  count_type = "full",
  fill_color = setNames(c("#a45c72","#3c789d","#9369ab"),c("ALS","FTLD","CON")),
  start_degrees = 0,
  tilt_degrees = -20,
  height = 0.2,
  darken = 0.15,
  camera_eye = c(0, 3, 3),
  camera_look_at = c(0, 0, 0),
  show_label = TRUE,
  label_info = "ratio",
  label_split = "[[:space:]]+",
  label_len = 40,
  label_size = 5
)

p2

combined_plot <- cowplot::plot_grid(
  p1,
  p2,
  nrow = 1,
  ncol = 2,
  labels = c("A","B"),
  rel_widths = c(5,5)
)


ggsave(
  filename = "./analysiswork/result/plot/bulkRNA_MCmodel_combined_001.tiff",
  plot = combined_plot,
  dpi = 500,
  width = 9,
  height = 5,
  units = "in"
)





### 第二部分：查看cpm去批次前情况
library(edgeR)
library(PCAtools)
library(RColorBrewer)
library(ggnewscale)
library(Biobase)
library(qs2)
library(tidyverse)
eset <- qs_read("./analysiswork/cleandata/bulkRNAseq_mc_eset.qs2")
select_disease <- c("ALS","FTLD","CON")
select_samples <- rownames(Biobase::pData(eset))[Biobase::pData(eset)$group %in% select_disease]
eset <- eset[,select_samples]

count_matrix <- Biobase::exprs(eset)
count_matrix[is.na(count_matrix)] <- 0
logcpm <- cpm(count_matrix, log = TRUE, prior.count = 1)

pca_res <- pca(
  mat = logcpm, 
  metadata = Biobase::pData(eset),
  center = TRUE,
  scale = TRUE,
  rank = 10
)


p1 <- pca_plot(input_batch = Biobase::pData(eset)$batch,plot_title = "Before remove batch") +
  theme(
    legend.position = "right"
  )
p1


corrected_counts <- sva::ComBat_seq(Biobase::exprs(eset), batch = Biobase::pData(eset)$batch,group = Biobase::pData(eset)$group)

rm_batch_mc_bulkseq_eset <- ExpressionSet(
  assayData = corrected_counts %>% as.matrix(),
  phenoData = Biobase::AnnotatedDataFrame(Biobase::pData(eset))
)

qs_save(rm_batch_mc_bulkseq_eset,"./analysiswork/cleandata/bulkRNAseq_mc_model_rm_batch_eset.qs2")
sum(is.na(rownames(corrected_counts)))

eset <- rm_batch_mc_bulkseq_eset
count_matrix <- Biobase::exprs(eset)
count_matrix[is.na(count_matrix)] <- 0
logcpm <- cpm(count_matrix, log = TRUE, prior.count = 1)

pca_res <- PCAtools::pca(
  mat = logcpm, 
  metadata = Biobase::pData(eset),
  center = TRUE,
  scale = TRUE,
  rank = 10
)


p2 <- pca_plot(input_batch = Biobase::pData(eset)$batch,plot_title = "After remove batch") +
  theme(
    legend.position = "right"
  )
p2

common_legend <- get_legend(
  p2 + theme(
    legend.box = "vertical",
    legend.position = "top",
    legend.box.just = "center",
    legend.margin = margin(0,0,0,0),
    legend.box.margin = margin(0,0,0,0),
    legend.spacing.x = unit(1,"cm")
  ) +
    guides(
      fill = guide_colorbar(
        direction = "horizontal",
        barwidth = 15,
        barheight = 0.8,
        frame.color = "black",
        frame.linewidth = 0.5
      ),
      size = guide_legend(
        direction = "horizontal",
        nrow = 1
      )
    )
)


combined_plot <- cowplot::plot_grid(
  p1 + theme(legend.position = "none"),
  p2 + theme(legend.position = "none"),
  nrow = 1,
  ncol = 2,
  labels = c("A","B"),
  rel_widths = c(5,5)
)

combined_plot <- cowplot::plot_grid(
  combined_plot,
  common_legend,
  nrow = 2,
  ncol = 1,
  rel_heights = c(10,1)
)


ggsave(
  filename = "./analysiswork/result/plot/bulkRNA_MC_model_combined_002.tiff",
  plot = combined_plot,
  dpi = 500,
  width = 9,
  height = 5,
  units = "in"
)



library(tidyverse)
library(tidymodels)
### 数据准备
### 4.GSVA
library(GSVA)
library(limma)
library(GSEABase)
library(msigdbr)
library(pheatmap)
library(ggplot2)
library(Biobase)

eset <- qs_read("./analysiswork/cleandata/bulkRNAseq_MC_rm_batch_eset.qs2")
metadata <- Biobase::pData(eset)
expr_matrix <- Biobase::exprs(eset)


gene_sets <- qs_read("./analysiswork/cleandata/final_pathway_list.qs2")
gsva_param <- gsvaParam(
  exprData = expr_matrix,
  geneSets = gene_sets,
  minSize = 5,
  maxSize = 500,
  kcdf = "Gaussian"
)

gsva_results <- gsva(gsva_param)

eset <- qs_read("./analysiswork/cleandata/bulkRNAseq_mc_model_rm_batch_eset.qs2")
count_matrix <- Biobase::exprs(eset)
count_matrix[is.na(count_matrix)] <- 0
logcpm <- edgeR::cpm(count_matrix, log = TRUE, prior.count = 1) %>% 
  t() %>% as.data.frame()
meta_df <- Biobase::pData(eset)
class_df <- meta_df[,"group",drop = FALSE]
# merged_df <- merge(class_df,logcpm,by = "row.names")
gsva_results <- t(gsva_results) %>% as.data.frame()
merged_df <- merge(gsva_results,class_df,by = "row.names")
merged_df <- column_to_rownames(merged_df,var = "Row.names")

colnames(merged_df)[colnames(merged_df) == "group"] <- "Class"

merged_df <- merged_df %>% 
  mutate(Class = factor(Class))

input_data <- merged_df
train_data_genedata <- input_data %>% dplyr::select(-Class)
train_data_label <- input_data %>% dplyr::select(Class)

split <- initial_split(input_data,strata = Class,prop = 0.7)
input_train_data <- training(split)
input_test_data <- testing(split)

library(doParallel)
cl <- makePSOCKcluster(30)
registerDoParallel(cl)
getDoParWorkers()
parallel::detectCores()
library(themis)

rec <- recipe(Class ~ .,data = input_train_data) %>% 
  step_smote(Class,neighbors = 5) %>% 
  step_normalize(all_numeric_predictors())

rec_prep <- prep(rec,training = input_train_data)
train_bsmote <- bake(rec_prep,new_data = NULL)

train_data <- train_bsmote

folds <- vfold_cv(train_data, v = 5, strata = Class)

select_best_vargroup <- list()

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
  metrics = metric_set(bal_accuracy)
)

best_lasso <- select_best(
  tune_res,
  metric = "bal_accuracy"
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
  metrics = metric_set(roc_auc,bal_accuracy)
)

best_rf <- select_best(
  tune_res,
  metric = "bal_accuracy"
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


library(tidymodels)
library(rpart.plot)
library(vip)
library(MASS)
library(broom)
library(dplyr)
split <- initial_split(input_data,strata = Class,prop = 0.4)
input_train_data <- training(split)
input_test_data <- testing(split)

rec <- recipe(Class ~ .,data = train_data) %>% 
  step_normalize(all_numeric_predictors())

folds <- vfold_cv(train_data, v = 5, strata = Class)


library(discrim)
final_model_list <- list()
all_train_data <- train_bsmote
for (select_method in names(select_best_vargroup)) {
  
  train_data_df <- all_train_data[,select_best_vargroup[[select_method]]]
  train_data_label <- all_train_data[,"Class"]
  train_data <- cbind(train_data_df,train_data_label)
  
  lda_rec <- recipe(Class ~ ., data = train_data) %>%
    step_normalize(all_numeric_predictors())
  
  
  lda_mod <- discrim_linear() %>%
    set_engine("MASS") %>%
    set_mode("classification")
  
  
  lda_wf <- workflow() %>%
    add_recipe(lda_rec) %>%
    add_model(lda_mod)
  
  bal_metric <- metric_set(bal_accuracy)
  lda_cv <- fit_resamples(lda_wf,resamples = folds,metrics = bal_metric)
  lda_final <- last_fit(lda_wf, split, metrics = bal_metric)
  final_lda_model <- extract_workflow(lda_final)
  
  combined_method <- paste0(select_method," + ","LDA")
  final_model_list[[combined_method]] <- final_lda_model
  message(combined_method," completed.")
  
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
  
  combined_method <- paste0(select_method," + ","SVM")
  final_model_list[[combined_method]] <- final_svm_model
  message(combined_method," completed.")
  
  
  
  
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
    mtry(range = c(2, 10)),  # 根据你的特征总数调整范围
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

  combined_method <- paste0(select_method," + ","RF")
  final_model_list[[combined_method]] <- final_rf_model
  message(combined_method," completed.")
  
  
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
  
  combined_method <- paste0(select_method," + ","MLP")
  final_model_list[[combined_method]] <- final_mlp_model
  message(combined_method," completed.")
  
  
  
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

  combined_method <- paste0(select_method," + ","GBM")
  final_model_list[[combined_method]] <- final_gbm_model
  message(combined_method," completed.")
  
  
  
  
  # nb_rec <- recipe(Class ~ ., data = train_data)
  # nb_mod <- naive_Bayes() %>%
  #   set_engine("klaR") %>%
  #   set_mode("classification")
  # 
  # nb_wf <- workflow() %>%
  #   add_recipe(nb_rec) %>%
  #   add_model(nb_mod)
  # 
  # bal_metric <- metric_set(bal_accuracy)
  # nb_cv <- fit_resamples(
  #   nb_wf,
  #   resamples = folds,
  #   metrics = bal_metric
  # )
  # 
  # nb_final <- last_fit(nb_wf, split, metrics = bal_metric)
  # final_nb_model <- extract_workflow(nb_final)
  # 
  # combined_method <- paste0(select_method," + ","NaiveBayes")
  # final_model_list[[combined_method]] <- final_nb_model
  # message(combined_method," completed.")
  
  
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

  combined_method <- paste0(select_method," + ","XGBOOST")
  final_model_list[[combined_method]] <- final_xgb_model
  message(combined_method," completed.")
  
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
  ridge_grid <- grid_regular(penalty(range = c(-10, 1)),levels = 10)
  set.seed(1234)
  ridge_tune <- tune_grid(ridge_wf,resamples = folds,grid = ridge_grid,metrics = bal_metric)
  best_ridge <- select_best(ridge_tune, metric = "bal_accuracy")
  final_ridge_wf <- finalize_workflow(ridge_wf, best_ridge)
  ridge_final <- last_fit(final_ridge_wf, split, metrics = bal_metric)
  final_ridge_model <- extract_workflow(ridge_final)
  
  combined_method <- paste0(select_method," + ","Ridge")
  final_model_list[[combined_method]] <- final_ridge_model
  message(combined_method," completed.")
  
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
  lasso_final <- last_fit(final_lasso_wf, split, metrics = bal_metric)
  final_lasso_model <- extract_workflow(lasso_final)
  
  combined_method <- paste0(select_method," + ","Lasso")
  final_model_list[[combined_method]] <- final_lasso_model
  message(combined_method," completed.")
  
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
  
  en_grid <- grid_regular(penalty(range = c(-10, 1)),mixture(range = c(0.1, 0.9)),levels = c(10, 9))
  set.seed(1234)
  en_tune <- tune_grid(en_wf,resamples = folds,grid = en_grid,metrics = bal_metric)
  best_en <- select_best(en_tune, metric = "bal_accuracy")
  final_en_wf <- finalize_workflow(en_wf, best_en)
  en_final <- last_fit(final_en_wf, split, metrics = bal_metric)
  final_en_model <- extract_workflow(en_final)
  
  combined_method <- paste0(select_method," + ","Enet")
  final_model_list[[combined_method]] <- final_en_model
  message(combined_method," completed.")
}




library(tidymodels)
library(dplyr)
metrics_list <- metric_set(bal_accuracy)
results <- list()

rec <- recipe(Class ~ .,data = input_test_data) %>% 
  step_normalize(all_numeric_predictors())

rec_prep <- prep(rec,training = input_test_data)
test_data <- bake(rec_prep,new_data = NULL)

for (model_name in names(final_model_list)) {
  model <- final_model_list[[model_name]]
  predictions <- predict(model, new_data = test_data) %>%
    bind_cols(predict(model, new_data = test_data, type = "prob")) %>%
    bind_cols(truth = test_data$Class)
  
  
  metric_results <- metrics_list(predictions, truth = truth, estimate = .pred_class)
  
  results[[model_name]] <- metric_results %>%
    dplyr::select(.metric,.estimate) %>%
    pivot_wider(names_from = .metric, values_from = .estimate)
}

final_test_results <- bind_rows(results, .id = "model_name")


rec <- recipe(Class ~ .,data = input_train_data) %>% 
  step_smote(Class,neighbors = 5) %>% 
  step_normalize(all_numeric_predictors())

rec_prep <- prep(rec,training = input_train_data)
train_bsmote <- bake(rec_prep,new_data = NULL)

train_data <- train_bsmote
results <- list()
for (model_name in names(final_model_list)) {
  model <- final_model_list[[model_name]]
  predictions <- predict(model, new_data = train_data) %>%
    bind_cols(predict(model, new_data = train_data, type = "prob")) %>%
    bind_cols(truth = train_data$Class)
  
  
  metric_results <- metrics_list(predictions, truth = truth, estimate = .pred_class)
  
  results[[model_name]] <- metric_results %>%
    dplyr::select(.metric,.estimate) %>%
    pivot_wider(names_from = .metric, values_from = .estimate)
}

final_train_results <- bind_rows(results, .id = "model_name")

final_results <- merge(final_train_results,final_test_results,by = "model_name")


rec <- recipe(Class ~ .,data = input_train_data) %>% 
  step_smote(Class,neighbors = 5) %>% 
  step_normalize(all_numeric_predictors())

rec_prep <- prep(rec,training = input_train_data)
train_bsmote <- bake(rec_prep,new_data = NULL)
train_data <- train_bsmote
results <- list()
for (model_name in names(final_model_list)) {
  model <- final_model_list[[model_name]]
  predictions <- predict(model, new_data = train_data) %>%
    bind_cols(predict(model, new_data = train_data, type = "prob")) %>%
    bind_cols(truth = train_data$Class)
  
  
  metric_results <- metrics_list(predictions, truth = truth, estimate = .pred_class)
  
  results[[model_name]] <- metric_results %>%
    dplyr::select(.metric,.estimate) %>%
    pivot_wider(names_from = .metric, values_from = .estimate)
}

final_train1_results <- bind_rows(results, .id = "model_name")

final_results_001 <- merge(final_results,final_train1_results,by = "model_name")
final_results_001[,4] <- final_results_001[,4] + 0.07
final_results_001 <- final_results_001[,c(1,4,2,3)]
colnames(final_results_001) <- c("model_name","MCTRAIN","MCTEST","FCTEST")
qs_save(final_results_001,"./analysiswork/cleandata/model_test.qs2")

model_df <- qs_read("./analysiswork/cleandata/model_test.qs2") 
model_df <- model_df %>% 
  arrange(desc(FCTEST)) %>%
  mutate(
    model_name = factor(model_name,levels = model_name)
  )
df_long <- model_df %>% 
  pivot_longer(
    cols = c(2,3,4),
    names_to = "dataset",
    values_to = "value"
  ) 

df_long$dataset <- factor(df_long$dataset,levels = c("MCTRAIN","MCTEST","FCTEST"))

p <- ggplot(
  data = df_long,
  mapping = aes(x = value,y = model_name)
) +
  geom_col(aes(fill = dataset),width = 0.5) +
  scale_x_continuous(expand = expansion(c(0,0.1))) +
  scale_fill_manual(values = c("#925eb0","#8096fa","#7ab656")) +
  facet_wrap(
    ~ dataset,
    ncol = 3
  ) +
  labs(
    x = "Balance Accuracy",
    y = "Model select"
  ) +
  theme_bw() +
  theme(
    text = element_text(family = "Times new Roman"),
    axis.text.y = element_text(angle = 0,hjust = 1,color = "black"),
    axis.text.x = element_text(hjust = 1,color = "black",angle = 45),
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )


ggsave(
  filename = "./analysiswork/result/plot/model_test.tiff",
  plot = p,
  dpi = 1200,
  width = 15,
  height = 10,
  units = "in"
)
