
library(tidyverse)
library(caret)
library(pROC)
library(glmnet)
library(DMwR)
library(rmda)
library(ggpubr)
library(ModelGood)
library(Hmisc)
library(htmlTable)
library(rms)
library(mRMRe)
library(DescTools)
library(Publish)
#library(smotefamily)

heatplot <- function(conf_mat)
{
  conf_mat$table %>% 
    heatmap.2(dendrogram = 'none', trace = 'none',
              key = F, cellnote = conf_mat$table,
              Rowv = F, Colv = F, col = 'bluered', 
              xlab = 'Reference', ylab = 'Prediction', 
              notecol = 'yellow')
}


source('assist_file.R')


set.seed(181)
#206 8:2,1817:3
p_thresh <- 0.05
clinics_column <- 9

fpath <- choose.files()
dt_all<- read_csv(fpath)


#dt_all<-dt_all1[seq(2,244,by=2),]



dt_all$Label <- factor(dt_all$Label, ordered = T)

dt <- dt_all[, -c(2:clinics_column)]

#dt <- dt_all

if(!is_empty(nearZeroVar(dt)))
{
  dt <- dt[, -nearZeroVar(dt)]
}

s_pre <- preProcess(dt, method = 'medianImpute')
dt <- predict(s_pre, dt)

# dt <- dt

# clinical data

dt_cli <- dt_all[, c(1:clinics_column)]
s_pre <- preProcess(dt_cli, method = 'medianImpute')
dt_cli <- predict(s_pre, dt_cli)



acc_val_train <- numeric(0)
acc_val_test <- numeric(0)

all_idx <- list()
train_list <- list()
test_list <- list()
cvfit_list <- list()
fit_list <- list()
out_list <- list()

for(numIt  in c(1:100))
{
  idx_train <- createDataPartition(dt$Label, p = 0.7, list = F)
  dt_train <- dt[idx_train, ]
  dt_test <- dt[-idx_train, ]
  
  # browser()
  
  all_idx[[numIt]] <- idx_train
  step_pre <- preProcess(dt_train, method = c('center', 'scale'))
  dt_train_pre <- predict(step_pre, dt_train) %>% as.data.frame
  
  dt_train_pre0 <- dt_train_pre
  dt_train_pre <- SMOTE(Label~., data = dt_train_pre)
  dt_test_pre <- predict(step_pre, dt_test)
  
  dt_mrmr <- mRMR.data(dt_train_pre)
  f_sel <- mRMR.classic(data = dt_mrmr, target_indices = c(1), feature_count = 20)
  
  # browser()
  
  dt_train_pre <- select(dt_train_pre, c('Label', featureNames(f_sel)[unlist(solutions(f_sel))]))
  dt_train_pre0 <- select(dt_train_pre0, c('Label', featureNames(f_sel)[unlist(solutions(f_sel))]))
  dt_test_pre <- select(dt_test_pre, c('Label', featureNames(f_sel)[unlist(solutions(f_sel))]))
  
  x <- as.matrix(dt_train_pre[, -1])
  y <- dt_train_pre$Label
  
  cv.fit <- cv.glmnet(x, y, family = 'binomial')
  fit <- glmnet(x, y, family = 'binomial')
  
  
  train_list[[numIt]] <- dt_train_pre0
  test_list[[numIt]] <- dt_test_pre
  cvfit_list[[numIt]] <- cv.fit
  fit_list[[numIt]] <- fit
  
  # browser()
  pre_res_test <- as.vector(predict(fit, newx = as.matrix(dt_test_pre[, -1]), s = cv.fit$lambda.min))
  roc_res_test <- pROC::roc(dt_test_pre$Label, pre_res_test)
  
  
  
  pre_res_train <- as.vector(predict(fit, newx = x, s = cv.fit$lambda.min))
  roc_res_train <- pROC::roc(dt_train_pre$Label, pre_res_train)
  
  
  dir_sign_test <- roc_res_test$direction
  dir_sign_train <- roc_res_train$direction
  
  if(dir_sign_test == dir_sign_train)
  {
    acc_val_test <- c(acc_val_test, pROC::auc(roc_res_test))
    acc_val_train <- c(acc_val_train, pROC::auc(roc_res_train))
  }
  else
  {
    acc_val_test <- c(acc_val_test, 0)
    acc_val_train <- c(acc_val_train, 0)
  }
}

cv_tbl <- tibble(Training = acc_val_train, Test = acc_val_test)
write_csv(cv_tbl, path = paste0(as.numeric(Sys.time()), '.csv'))

idx_vec <- c(1:length(acc_val_test))
idx_vec <- idx_vec[acc_val_train > acc_val_test]
acc_val <- acc_val_test[acc_val_train > acc_val_test]
init_idx <- which.max(acc_val)
sel_idx <- idx_vec[init_idx]



write.csv(idx_train,"tt.csv")

idx_train <- all_idx[[sel_idx]]
grp_info <- tibble()
grp_info <- tibble(Label = dt$Label, Group = 'Test')
grp_info$Radscore <- 0
grp_info$Group[idx_train] <- 'Training'


dt_train_final <- train_list[[sel_idx]]
dt_test_final <- test_list[[sel_idx]]
cvfit <- cvfit_list[[sel_idx]]
fit <- fit_list[[sel_idx]]

s = cvfit$lambda.min

pre_res_train <- as.vector(predict(fit, newx = as.matrix(dt_train_final[, -1]), s = s))
pre_res_train_prob <- as.vector(predict(fit, newx = as.matrix(dt_train_final[, -1]), s = s, 
                                        type = 'link'))
roc_res_train <- pROC::roc(dt_train_final$Label, pre_res_train, ci = T)

rad_cutoff <- coords(roc_res_train, x = 'best', transpose = F, ret = 'threshold')$threshold[1]




out_res_train <- ifelse(pre_res_train > rad_cutoff, 1, 0)
conf_mat_train <- confusionMatrix(as.factor(out_res_train), as.factor(dt_train_final$Label), 
                                  positive = '1')
rec_train <- c(conf_mat_train$overall[c(1, 3, 4)], conf_mat_train$byClass[c(1:4)])

pre_res_test <- as.vector(predict(fit, newx = as.matrix(dt_test_final[, -1]), s = s))
pre_res_test_prob <- as.vector(predict(fit, newx = as.matrix(dt_test_final[, -1]), s = s, 
                                       type = 'link'))
roc_res_test <- pROC::roc(dt_test_final$Label, pre_res_test, ci = T)

out_res_test <- ifelse(pre_res_test > rad_cutoff, 1, 0)
conf_mat_test <- confusionMatrix(as.factor(out_res_test),as.factor(dt_test_final$Label), 
                                 positive = '1')
rec_test <- c(conf_mat_test$overall[c(1, 3, 4)], conf_mat_test$byClass[c(1:4)])


rec_all_rad <- data.frame(rbind(rec_train, rec_test), row.names = c('Train', 'Test'))

grp_info$Radscore[idx_train] <- pre_res_train
grp_info$Radscore[-idx_train] <- pre_res_test

grp_info <- bind_cols(grp_info, dt_cli[, -1])
for(i_varname in colnames(grp_info)[-c(1:3)])
{
  var_v <- grp_info[[i_varname]]
  if(length(unique(var_v)) < 5)
  {
    grp_info[[i_varname]] <- factor(var_v)
  }
}


formu_sum <- paste('Label', paste(colnames(grp_info)[-c(1, 2, 3)], collapse = '+'), 
                   sep = '~') %>% paste('Q(Radscore)', sep = '+') %>% as.formula()

sum_train <- univariateTable(formu_sum, data = filter(grp_info, Group == 'Training'), 
                             show.totals = F)
sum_test <- univariateTable(formu_sum, data = filter(grp_info, Group == 'Test'), 
                            show.totals = F)

write_csv(grp_info, 'group_info.csv')


formu_sum2<- paste('Group', paste(colnames(grp_info)[-c(1, 2, 3)], collapse = '+'), 
                   sep = '~') %>% paste('Q(Radscore)', sep = '+') %>% as.formula()
sum_all <- univariateTable(formu_sum2, data = grp_info, 
                           show.totals = F)





## rad score
dt_final_test <- tibble(Label = dt_test_final$Label, rad_score = pre_res_test)
dt_final_arr <- arrange(dt_final_test, rad_score)
dt_final_arr$x <- 1:nrow(dt_final_arr)

dt_final_train <- tibble(Label = dt_train_final$Label, rad_score = pre_res_train)

radscore_limit <- c(dt_final_train$rad_score, dt_final_test$rad_score)

limit_y <- c(min(radscore_limit) - 3*sd(radscore_limit), max(radscore_limit) + 3*sd(radscore_limit))



p_train <-  ggboxplot(x = 'Label', y = 'rad_score', data = dt_final_train,
                     add = 'jitter', color = 'Label', palette = 'jco')+
  stat_compare_means(method = 'wilcox.test') + 
 geom_hline(yintercept = rad_cutoff) + theme_bw()
   #geom_hline(yintercept = c(5,10,15)) + theme_bw()
p_test <- ggboxplot(x = 'Label', y = 'rad_score', data = dt_final_test,
                    add = 'jitter', color = 'Label', palette = 'jco')  + 
  stat_compare_means(method = 'wilcox.test') + 
  geom_hline(yintercept = rad_cutoff) + theme_bw()
  #geom_hline(yintercept = c(5,10,15)) + theme_bw()



coefs <- coefficients(fit, s = s)
useful_feature <- unlist(coefs@Dimnames)[coefs@i + 1]
useful_feature <- useful_feature[-1]

dt_coef <- data.frame(Feature = useful_feature, Coef = coefs@x[-1])
dt_coef <- arrange(dt_coef, desc(Coef))
dt_coef$Feature <- factor(dt_coef$Feature, 
                          levels = as.character(dt_coef$Feature))

p_coef <- ggplot(aes(x = Feature, y = Coef), data = dt_coef)
p_coef <- p_coef + geom_col(fill = 'blue', width = 0.7) + coord_flip() + 
  theme_bw() + ylab('Coefficients')

final_data_test <- add_column(dt_test_final, radscore = pre_res_test)
final_data_test <- select(final_data_test, c('Label', useful_feature, 'radscore'))
write_csv(final_data_test, path = 'dataset_test.csv')
final_data_train <- add_column(dt_train_final, radscore = pre_res_train)
final_data_train <- select(final_data_train, c('Label', useful_feature, 'radscore'))
write_csv(final_data_train, path = 'dataset_train.csv')



fit_train <- glm(Label~rad_score, data = dt_final_train, family = 'binomial')

dt_dca_train <- dt_final_train
dt_dca_train$Label <- as.numeric(dt_dca_train$Label) - 1
dca_curve <- decision_curve(Label~rad_score, data = dt_dca_train)

radscore <- paste('Radscore = ', paste(round(coefs@x[-1], 3), 
                                       useful_feature, sep = '*', collapse = '+'), '+', round(coefs@x[1], 3))


# clinics analysis

idx_train <- all_idx[[sel_idx]]
dt_cli_train <- dt_cli[idx_train, ]
dt_cli_test <- dt_cli[-idx_train, ]

p_list <- lapply(dt_cli_train[, -1], inter_test, y = dt_cli_train$Label)


sel_name <- c(names(which(p_list < p_thresh)))
dt_cli_train_1 <- select(dt_cli_train, c('Label', sel_name))

res_ulogit_list <- lapply(colnames(dt_cli_train_1)[-1], ulogit_test, dt = dt_cli_train_1)
res_ulogit <- bind_rows(res_ulogit_list)
res_ulogit <- res_ulogit[-seq(from = 1, to = nrow(res_ulogit), by = 2), ]
res_ulogit$Var <- colnames(dt_cli_train_1)[-1]

res_ulogit_sel <- filter(res_ulogit, p_val < p_thresh)
cli_name <- res_ulogit_sel$Var

dt_cli_train_2 <- select(dt_cli_train_1, c('Label', cli_name))
dt_cli_test_2 <- select(dt_cli_test, c('Label', cli_name))

write.csv(res_ulogit, file = 'ulogit_cli.csv')

pvals_ulogit <- res_ulogit$p_val
names(pvals_ulogit) <- res_ulogit$Var
log_fit <- glm(Label~., data = dt_cli_train_2, family = binomial)
vif_val <- vif(log_fit)
vif_idx <- which(vif_val > 5)

dt_cli_train_i <- dt_cli_train_2
while(!is_empty(vif_idx))
{
  vif_name <- names(vif_idx)
  excl_name <- which.max(pvals_ulogit[vif_name]) %>% names()
  dt_cli_train_i <- select(dt_cli_train_i, -c(excl_name))
  log_fit <- glm(Label~., data = dt_cli_train_i, family = binomial)
  vif_val <- vif(log_fit)
  vif_idx <- which(vif_val > 5)
}
log_fit_final <- step(log_fit)
output_mlogit(log_fit_final, 'mlogit_cli.csv')
dt_cli_train_final <- model.frame(log_fit_final, data = dt_cli_train)
dt_cli_test_final <- model.frame(log_fit_final, data = dt_cli_test)


#output_mlogit(log_fit, 'mlogit_cli.csv')
#dt_cli_train_final <- model.frame(log_fit, data = dt_cli_train)
#dt_cli_test_final <- model.frame(log_fit, data = dt_cli_test)


dt_combined_train <- dt_cli_train_final
dt_combined_train$rad_score <- pre_res_train_prob

dt_combined_test <- dt_cli_test_final
dt_combined_test$rad_score <- pre_res_test_prob


com_form <- paste('Label', paste(colnames(dt_combined_train)[-1], collapse = '+'), sep = '~') %>% as.formula


mod_com_final <- glm(com_form, data = dt_combined_train, family = 'binomial')

output_mlogit(mod_com_final, 'mlogit_com.csv')

dt_combined_train_final <- model.frame(mod_com_final, data = dt_combined_train)
dt_combined_test_final <- model.frame(mod_com_final, data = dt_combined_test)

#com_form <- paste('Label', paste(colnames(dt_combined_train_final)[-1], collapse = '+'), sep = '~') %>% as.formula

#dt_nom <- filter(dt_combined_train_final, rad_score > -8 & rad_score < 5)
dt_nom <- dt_combined_train_final

ddist_train_com <- datadist(dt_nom)
options(datadist = 'ddist_train_com')
mod_train <- lrm(com_form, 
                 data = dt_nom, x = TRUE, y = TRUE)


mod_train_glm <- glm(Label~., 
                     data = dt_nom, family = 'binomial')

coef_nom <- coef(mod_train_glm)
nomo_score <- paste('Nomoscore =', 
                    paste(names(coef_nom), coef_nom, sep = '*') %>% paste(collapse = '+'))
nom_com <- nomogram(mod_train, lp = F, fun = plogis, fun.at = c(0.1, 0.2,0.5,0.8,0.9), 
                    funlabel = 'Risk')
windows()
plot(nom_com)


mod_test <- lrm(com_form, 
                data = dt_combined_test_final, x = TRUE, y = TRUE)



cli_form <- paste('Label', paste(colnames(dt_combined_train_final)[-c(1, ncol(dt_combined_train_final))], collapse = '+'), 
                  sep = '~') %>% as.formula


rad_form <- paste('Label', paste(colnames(dt_combined_train_final)[ncol(dt_combined_train_final)], collapse = '+'), 
                  sep = '~') %>% as.formula



res_train_final <- predict(mod_com_final, newdata = dt_combined_train)
roc_com_train <- pROC::roc(dt_combined_train$Label, res_train_final, 
                           ci = T)

quangege6 <- coords(roc_com_train, x = 'best', transpose = F, ret = 'threshold')$threshold[1]
cutoff<- as.numeric(quangege6)

res_train_final_bin <- as.factor(ifelse(res_train_final > cutoff, 1, 0))

res_conf_com_train <- confusionMatrix(res_train_final_bin, dt_combined_train$Label, positive = '1')
res_test_final <- predict(mod_train, newdata = dt_combined_test)
roc_com_test <- pROC::roc(dt_combined_test$Label, res_test_final,
                          ci = T)
quangege7 <- coords(roc_com_test, x = 'best', transpose = F, ret = 'threshold')$threshold[1]
cutoff<- as.numeric(quangege7)

res_test_final_bin <- as.factor(ifelse(res_test_final > cutoff, 1, 0))
res_conf_com_test <- confusionMatrix(dt_combined_test$Label, res_test_final_bin, 
                                     positive = '1')

rec_train_com <- c(res_conf_com_train$overall[c(1, 3, 4)], res_conf_com_train$byClass[c(1:4)])
rec_test_com <- c(res_conf_com_test$overall[c(1, 3, 4)], res_conf_com_test$byClass[c(1:4)])

rec_all_com <- data.frame(rbind(rec_train_com, rec_test_com), row.names = c('Train', 'Test'))


dt_dca <- dt_combined_train_final
dt_dca$Label <- ifelse(dt_dca$Label == '0', 0, 1)
dca1 <- decision_curve(com_form, 
                       data = dt_dca)
dca2 <- decision_curve(cli_form, 
                       data = dt_dca)

dca3 <- decision_curve(rad_form, 
                       data = dt_dca)




#plot_decision_curve(list(dca1, dca11), confidence.intervals = F, 
                    #col = c('maroon1', 'deepskyblue2', 'springgreen3','black'), 
#curve.names = c('Radiomics nomogram','Model intergrating pathology'),
                    
#legend.position = 'topright', cost.benefit.axis = F, cost.benefits = F,standardize = FALSE)






cli_mod <- glm(cli_form, data = dt_combined_train_final, family = 'binomial')

res_cli_train <- predict(cli_mod, newdata = dt_combined_train_final, type = 'link')
roc_cli_train <- pROC::roc(dt_combined_train_final$Label, res_cli_train,ci=T)
quangege8 <- coords(roc_cli_train, x = 'best', transpose = F, ret = 'threshold')$threshold[1]
cutoff_cli<- as.numeric(quangege8)
res_cli_test <- predict(cli_mod, newdata = dt_combined_test_final, type = 'link')
roc_cli_test <- pROC::roc(dt_combined_test_final$Label, res_cli_test,ci=T)


res_cli_train_bin <- as.factor(ifelse(res_cli_train > cutoff_cli, 1, 0))
res_cli_test_bin <- as.factor(ifelse(res_cli_test > cutoff_cli, 1, 0))

write.csv(cbind(dt_combined_train_final,res_train_final,res_train_final_bin,res_cli_train),"traincl.csv")
write.csv(cbind(dt_combined_test_final,res_test_final,res_test_final_bin,res_cli_test),"testcl.csv")


conf_mat_cli_test <- confusionMatrix(res_cli_test_bin,dt_combined_test_final$Label, positive = '1')

conf_mat_cli_train <- confusionMatrix(res_cli_train_bin, dt_combined_train_final$Label,
                                       positive = '1')


rec_train_cli <- c(conf_mat_cli_train$overall[c(1, 3, 4)], conf_mat_cli_train$byClass[c(1:4)])
rec_test_cli <- c(conf_mat_cli_test$overall[c(1, 3, 4)], conf_mat_cli_test$byClass[c(1:4)])

rec_all_cli <- data.frame(rbind(rec_train_cli, rec_test_cli), row.names = c('Train', 'Test'))

rec_all <- bind_rows(rec_all_rad, rec_all_cli, rec_all_com)
rec_all$Model <- rep(c('Radiomics', 'Clinics', 'Nomogram'), times = rep(2, times = 3))



plot_decision_curve(list(dca1, dca2,dca3), confidence.intervals = F, 
                    col = c('red', 'blue', 'green', 'black'), 
                    curve.names = c('Nomogram','Radiomics', 'Clinics'),
                    cex.axis=1,
                    cex.lab=1,
                    legend.position = 'topright', cost.benefits = FALSE,standardize = FALSE)
#Bar
{
  
  #Data <- test_final
  #Levels = as.factor(response_test)
  Training<-res_train_final 
  Data <- as.data.frame(Training)
  Levels = as.factor(dt_combined_train_final$Label)
  Bar.Chart <- reorder(c(1:length(Levels)),Training)
  adjust<-coords(roc_com_train, x = 'best', transpose = F, ret = 'threshold')$threshold[1]
  Radscore_adjust<-Training-adjust
  p<-ggplot(Data, aes(x = Bar.Chart,
                      xlab ="Bar Chart",
                      y = Training,
                      fill = Levels ,
                      #color = c('red','blue'),
                     width=0.8 ) )+
    geom_bar(stat="identity",color="white", alpha=0.5)
  windows()
  p 
  
  #test
  Test<-res_test_final
  Data <- as.data.frame(Test)
  Levels = as.factor(dt_combined_test_final$Label)
  Bar.Chart <- reorder(c(1:length(Levels)),Test)
  adjust<-coords(roc_com_test, x = 'best', transpose = F, ret = 'threshold')$threshold[1]
  Radscore_adjust<-Test-adjust
  p<-ggplot(Data, aes(x = Bar.Chart,
                      xlab ="Bar Chart",
                      y = Test,
                      fill = Levels ,
                      #color = c('red','blue'),
                      width=0.8 ) )+
    geom_bar(stat="identity",color="white", alpha=0.5)
  windows()
  p 
  
}

cal_train <- calPlot2(mod_train, data = dt_combined_train_final, 
                     legend = F, col = 	'black', lty = 1)
cal_test <- calPlot2(mod_train, data = dt_combined_test_final,
                    legend = F, col = 'black', lty = 1)

HosmerLemeshowTest(cal_train$Frame$lrm, cal_train$Frame$jack)
HosmerLemeshowTest(cal_test$Frame$lrm, cal_test$Frame$jack)

rmarkdown::render('results_final.Rmd', output_file = 'result')

