library(tidyverse)
library(caret)
library(Publish)
library(survival)
library(glmnet)
library(ggpubr)
library(survminer)
library(rolr)
library(survIDINRI)
library(survAUC)
library(rms)
library(DCA)
library(timeROC)




dt.clinics <- read_csv('clinics_rad_DSS.csv') %>% select(-1)
dt.CT <- read_csv('feature_PVP.csv')
dt.PET <- read_csv('feature_AP.csv')

dt.radiomics <- inner_join(dt.CT, dt.PET, by = c('Patient'),
                           suffix = c('.PVP', '.AP')) %>%
  select(-1)


# dt.radiomics <- dt.PET %>% select(-1)
dt.train.list <- list()
dt.test.list <- list()

p.list.train <- list()
p.list.test <- list()

fit.list <- list()
cvfit.list <- list()

set.seed(560)
#65 multi40 1 6 13 17 18 20 #560 dss
for(i.a in 1:10)
{
  idx.all <- 1:nrow(dt.clinics)
  idx.train <- createDataPartition(dt.clinics$Progress, p = 0.7, list = F) %>% c()
  idx.test <- setdiff(idx.all, idx.train)
  
  
  dt.train.clinics <- dt.clinics[idx.train, ]
  dt.train.radiomics <- dt.radiomics[idx.train, ]
  dt.test.clinics <- dt.clinics[idx.test, ]
  dt.test.radiomics <- dt.radiomics[idx.test, ]
  
  
  s_prep <- preProcess(dt.train.radiomics, method = c('medianImpute', 'center', 'scale'))
  
  dt.train.radiomics.1 <- predict(s_prep, dt.train.radiomics)
  dt.test.radiomics.1 <- predict(s_prep, dt.test.radiomics)
  
  
  idx_nz <- nearZeroVar(dt.train.radiomics.1)
  
  if(!is_empty(idx_nz))
  {
    dt.train.radiomics.1 <- dt.train.radiomics.1[, -idx_nz]
    dt.test.radiomics.1 <- dt.test.radiomics.1[, -idx_nz]
  }
  
  dt.temp <- dt.train.radiomics.1 %>% add_column(Progress = dt.train.clinics$Progress, 
                                                 PFS = dt.train.clinics$PFS, 
                                                 .before = 1)
  cox.test <- coxphSeries(Surv(PFS, Progress)~1, 
                          data = dt.temp, 
                          vars = colnames(dt.temp)[-c(1:2)])
  sel.names <- cox.test$Variable[(which(cox.test$Pvalue < 0.1))]
  rm(dt.temp)
  
  dt.train.radiomics.2 <- select(dt.train.radiomics.1, sel.names)
  dt.test.radiomics.2 <- select(dt.test.radiomics.1, sel.names)
  
  idx_exc <- findCorrelation(cor(dt.train.radiomics.2), cutoff = 0.9)
  idx_in <- setdiff(1:ncol(dt.train.radiomics.2), idx_exc)

  dt.train.radiomics.3 <- dt.train.radiomics.2[, idx_in]
  dt.test.radiomics.3 <- dt.test.radiomics.2[, idx_in]
  
  # dt.train.radiomics.3 <- dt.train.radiomics.2
  # dt.test.radiomics.3 <- dt.test.radiomics.2
  
  x <- as.matrix(dt.train.radiomics.3)
  y <- Surv(dt.train.clinics$PFS, dt.train.clinics$Progress)
  
  cv.fit <- cv.glmnet(x, y, family = 'cox', nfolds = 10)
  fit <- glmnet(x, y, family = 'cox')
  
  fit.list[[i.a]] <- fit
  cvfit.list[[i.a]] <- cv.fit
  
  
  
  x.test <- as.matrix(dt.test.radiomics.3)
  radscore.train <- predict(fit, newx = x, s = cv.fit$lambda.min) %>% c()
  radscore.test <- predict(fit, newx = x.test, s = cv.fit$lambda.min) %>% c()
  
  
  dt.train.clinics$Radscore <- radscore.train
  dt.test.clinics$Radscore <- radscore.test
  
  
  dt.os.train <- dt.train.clinics
  dt.os.test <- dt.test.clinics
  
  cat(i.a)
  
  dt.train.list[[i.a]] <- dt.os.train
  dt.test.list[[i.a]] <- dt.os.test
  
  res <- coxph(Surv(dt.test.clinics$PFS, dt.test.clinics$Progress)~radscore.test)
  res <- summary(res)
  p.list.test[[i.a]] <- res$coefficients[length(res$coefficients)]
  
  
  res <- coxph(Surv(dt.train.clinics$PFS, dt.train.clinics$Progress)~radscore.train)
  res <- summary(res)
  p.list.train[[i.a]] <- res$coefficients[length(res$coefficients)]
}

idx.sel <- which.min(p.list.test)
# idx.sel <- 7


fit <- fit.list[[idx.sel]]
cv.fit <- cvfit.list[[idx.sel]]
# LASSO plot --------------------------------------------------------------
pdf('Lasso.pdf', width = 8, height = 9)
oldpar <- par(mfrow = c(2, 1))
plot(cv.fit)
plot(fit, s = cv.fit$lambda.min, xvar = 'lambda')
abline(v = log(cv.fit$lambda.min), lty = 2)
par(oldpar)
dev.off()

coefs <- coefficients(fit, s = cv.fit$lambda.min)

dt.coef <- tibble(Variable = coefs@Dimnames[[1]][coefs@i + 1], 
                  Coefficients = coefs@x) %>% arrange(desc(Coefficients))

Radscore <-paste('Radscore', 
                 paste(dt.coef$Coefficients, dt.coef$Variable, sep = '*') %>% paste(collapse = '+'), 
                 sep = '=')
sink('radscore.txt')
print(Radscore)
sink()
# Coefficients plot -------------------------------------------------------

ggbarplot(data = dt.coef, x = 'Variable', y = 'Coefficients', 
          width = 0.4, fill = 'blue') + 
  coord_flip() + theme_bw() + theme(axis.title = element_text(size = 20), 
                                    axis.text = element_text(size = 20))
ggsave('coefficient.pdf', width = 10, height = 4)

dt.os.train <- dt.train.list[[idx.sel]]
dt.os.test <- dt.test.list[[idx.sel]]

write_csv(dt.os.train, 'data_train.csv')
write_csv(dt.os.test, 'data_test.csv')


