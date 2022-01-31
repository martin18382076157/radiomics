nomo_analysis <- function(dt.os.train, dt.os.test, time, event, TimeInterval, funLabel, output_base)
{
  # ex_event <- setdiff(c('OS', 'Dead', 'PFS', 'Progress'), c(time, event))
  # dt.os.train <- select(dt.os.train, -ex_event)
  # dt.os.test <- select(dt.os.test, -ex_event)
  cox.test.2 <- coxphSeries(Surv(PFS, Progress)~1, 
                            data = dt.os.train, 
                            vars = colnames(dt.os.train)[-c(1:2)])
  
  cox.test.2.test <- coxphSeries(Surv(PFS, Progress)~1, 
                            data = dt.os.test, 
                            vars = colnames(dt.os.test)[-c(1:2)])
  
  write_csv(cox.test.2, sprintf('%s_uni_train.csv', output_base))
  write_csv(cox.test.2.test, sprintf('%s_uni_test.csv', output_base))
  
  sel.names.2 <- cox.test.2$Variable[(which(cox.test.2$Pvalue < 0.05))]
  

  dt.os.train.1 <- select(dt.os.train, c(event, time, sel.names.2))
  dt.os.test.1 <- select(dt.os.test, c(event, time, sel.names.2))
  
  surv.form <- sprintf('Surv(%s, %s)~.', time, event) %>% as.formula()
  cox.os.nom <- coxph(surv.form, data = dt.os.train.1) %>% step()
  
  
  coef.nom <- coefficients(cox.os.nom)
  Nomoscore <- paste('Nomoscore', 
                     paste(coef.nom, names(coef.nom), sep = '*') %>% paste(collapse = '+'), 
                     sep = '=')
  cli_name <- names(coef.nom)
  dt.os.train.2 <<- dt.os.train.1 %>% select(c(event, time, names(coef.nom))) %>% as.data.frame()
  dt.os.test.2 <<- dt.os.test.1 %>% select(c(event, time, names(coef.nom))) %>% as.data.frame()
  
  surv.form <- sprintf('Surv(%s, %s)~%s', time, event, paste(names(coef.nom), 
                                                             collapse = '+')) %>% as.formula()
  
  cox.os.mul.train <- coxph(surv.form, data = dt.os.train.2) %>% publish()
  cox.os.mul.test <- coxph(surv.form, data = dt.os.test.2) %>% publish()

  write_csv(cox.os.mul.train$regressionTable, sprintf('%s_mul_train.csv', output_base))
  write_csv(cox.os.mul.test$regressionTable, sprintf('%s_mul_test.csv', output_base))

  ddist <<- datadist(dt.os.train.2); options(datadist='ddist')
  
  f <- cph(surv.form, data=dt.os.train.2, surv = T)
  med  <- Quantile(f)
  surv <- Survival(f)  # This would also work if f was from cph
  nom <- nomogram(f, fun=list(function(x) surv(TimeInterval[1], x),
                              function(x) surv(TimeInterval[2], x),
                              function(x) surv(TimeInterval[3], x)), funlabel= funLabel, 
                  lp = F, fun.at = seq( 0.1, 0.9, by = 0.1))
  pdf(paste(output_base, 'Nomogram.pdf', sep = '_'), width = 8, height = 5)
  plot(nom)
  dev.off()
  
  nomoscore.train <- predict(cox.os.nom, newdata = dt.os.train.2)
  nomoscore.test <- predict(cox.os.nom, newdata = dt.os.test.2)
  
  return(list(nomoscore.train, nomoscore.test, cli_name))
}
nomo_analysis2 <- function(dt.os.train, dt.os.test, time, event, TimeInterval, funLabel, output_base)
{
  # ex_event <- setdiff(c('OS', 'Dead', 'PFS', 'Progress'), c(time, event))
  # dt.os.train <- select(dt.os.train, -ex_event)
  # dt.os.test <- select(dt.os.test, -ex_event)
  # cox.test.2 <- coxphSeries(Surv(PFS, Progress)~1, 
  #                           data = dt.os.train, 
  #                           vars = colnames(dt.os.train)[-c(1:2)])
  # 
  # cox.test.2.test <- coxphSeries(Surv(PFS, Progress)~1, 
  #                                data = dt.os.test, 
  #                                vars = colnames(dt.os.test)[-c(1:2)])
  # 
  # write_csv(cox.test.2, sprintf('%s_uni_train.csv', output_base))
  # write_csv(cox.test.2.test, sprintf('%s_uni_test.csv', output_base))
  # 
  # sel.names.2 <- cox.test.2$Variable[(which(cox.test.2$Pvalue < 0.05))]
  # 
  # 
  # 
  # 
  # dt.os.train.1 <- select(dt.os.train, c(event, time, sel.names.2))
  # dt.os.test.1 <- select(dt.os.test, c(event, time, sel.names.2))
  
  dt.os.train.1 <- dt.os.train
  dt.os.test.1 <- dt.os.test
  
  surv.form <- sprintf('Surv(%s, %s)~.', time, event) %>% as.formula()
  cox.os.nom <- coxph(surv.form, data = dt.os.train.1)
  
  
  coef.nom <- coefficients(cox.os.nom)
  Nomoscore <- paste('Nomoscore', 
                     paste(coef.nom, names(coef.nom), sep = '*') %>% paste(collapse = '+'), 
                     sep = '=')
  
  dt.os.train.2 <<- dt.os.train.1 %>% select(c(event, time, names(coef.nom))) %>% as.data.frame()
  dt.os.test.2 <<- dt.os.test.1 %>% select(c(event, time, names(coef.nom))) %>% as.data.frame()
  
  surv.form <- sprintf('Surv(%s, %s)~%s', time, event, paste(names(coef.nom), 
                                                             collapse = '+')) %>% as.formula()
  
  cox.os.mul.train <- coxph(surv.form, data = dt.os.train.2) %>% publish()
  cox.os.mul.test <- coxph(surv.form, data = dt.os.test.2) %>% publish()
  
  write_csv(cox.os.mul.train$regressionTable, sprintf('%s_mul_train.csv', output_base))
  write_csv(cox.os.mul.test$regressionTable, sprintf('%s_mul_test.csv', output_base))
  
  ddist <<- datadist(dt.os.train.2); options(datadist='ddist')
  
  f <- cph(surv.form, data=dt.os.train.2, surv = T)
  med  <- Quantile(f)
  surv <- Survival(f)  # This would also work if f was from cph
  nom <- nomogram(f, fun=list(function(x) surv(TimeInterval[1], x),
                              function(x) surv(TimeInterval[2], x), 
                              function(x) surv(TimeInterval[3], x)), funlabel= funLabel, 
                  lp = F, fun.at = seq(0.1, 0.9, by = 0.1))
  pdf(paste(output_base, 'Nomogram.pdf', sep = '_'), width = 8, height = 5)
  plot(nom)
  dev.off()
  
  nomoscore.train <- predict(cox.os.nom, newdata = dt.os.train.2)
  nomoscore.test <- predict(cox.os.nom, newdata = dt.os.test.2)
  
  return(list(nomoscore.train, nomoscore.test))
}

cidx_calc <- function(form1, form2, form3, dt, output)
{
  cox.os.rad.train <- coxph(form1, data = dt)
  cox.os.cli.train <- coxph(form2, data = dt)
  cox.os.com.train <- coxph(form3, data = dt)
  
  sink(output)
  summary(cox.os.rad.train) %>% print()
  summary(cox.os.cli.train) %>% print()
  summary(cox.os.com.train) %>% print()
  sink()
}

km_plot <- function(dt.train, dt.test, time, event, variable, outputBase)
{
  cuts <- surv_cutpoint(dt.train, time = time, event = event, 
                        variables = variable)
  
  cuts <- cuts$cutpoint$cutpoint
  
  dt.train[[variable]] <- ifelse(dt.train[[variable]] < cuts, 0, 1)
  dt.test[[variable]] <- ifelse(dt.test[[variable]] < cuts, 0, 1)
  
  # KM PFS plot Radscore--------------------------------------------------------------
  
  cox.form <- sprintf('Surv(%s, %s)~%s', time, event, variable) %>% as.formula()
  survfit.os.train <- survfit(cox.form, data = dt.train)
  survfit.os.test <- survfit(cox.form, data = dt.test)
  
  survfit.os.train$call$formula <- cox.form
  survfit.os.test$call$formula <- cox.form
  
  
  survplot.os.train <- ggsurvplot(survfit.os.train, data = dt.train, pval = T)
  survplot.os.test <- ggsurvplot(survfit.os.test, data = dt.test, pval = T)
 
  
  ggarrange(survplot.os.train$plot, survplot.os.test$plot, 
            nrow = 1, common.legend = T)
  ggsave(sprintf('%s_KM.pdf', outputBase), width = 12, height = 6)

  cuts
}

calplt <- function(Time, formula, dt, cohort)
{
  cph.model <- cph(formula, data = dt, 
                         surv = T, x = T, y = T, time.inc = Time)
  cal.model <- rms::calibrate(cph.model, u = Time, cmethod = 'KM', m = 10)
  cal.dt <- as_tibble(cal.model[, c(7, 9)]) %>% 
    add_column(Model = sprintf('%s_%02d', cohort, Time))

  cal.dt
}

TimeInterval <- c(12, 24, 36) 
funLabel <- c("1-Years Survival Probability", 
              "2-Years Survival Probability", 
              "3-Years Survival Probability")



dt.os.train <- read_csv('data_train.csv')
dt.os.test <- read_csv('data_test.csv')

dt.os.cli.train <- select(dt.os.train, -c('Radscore'))
dt.os.cli.test <- select(dt.os.test, -c('Radscore'))



output.cli.nomogram <- nomo_analysis(dt.os.cli.train, dt.os.cli.test, 'PFS', 'Progress',TimeInterval, funLabel, 'Clinics')
cli_remain_name <- output.cli.nomogram[[3]]


dt.os.train.com <- select(dt.os.train, c('PFS', 'Progress', cli_remain_name, 'Radscore'))
dt.os.test.com <- select(dt.os.test, c('PFS', 'Progress', cli_remain_name, 'Radscore'))

output.com.nomogram <- nomo_analysis2(dt.os.train.com, dt.os.test.com, 'PFS', 'Progress', 
                                      TimeInterval, funLabel, 'Combined')


nomoscore.com.train <- output.com.nomogram[[1]]
nomoscore.com.test <- output.com.nomogram[[2]]

nomoscore.cli.train <- output.cli.nomogram[[1]]
nomoscore.cli.test <- output.cli.nomogram[[2]]

dt.os.final.train <- dt.os.train %>% select('PFS', 'Progress', 'Radscore') %>% add_column(
  nomo.com = nomoscore.com.train, nomo.cli = nomoscore.cli.train)

dt.os.final.test <- dt.os.test %>% select('PFS', 'Progress', 'Radscore') %>% add_column(
  nomo.com = nomoscore.com.test, nomo.cli = nomoscore.cli.test)

cidx_calc(Surv(PFS,Progress)~Radscore, 
          Surv(PFS,Progress)~nomo.cli,
          Surv(PFS,Progress)~nomo.com, 
          dt = dt.os.final.train, 
          output = 'summary_train.txt')
cidx_calc(Surv(PFS,Progress)~Radscore, 
          Surv(PFS,Progress)~nomo.cli,
          Surv(PFS,Progress)~nomo.com, 
          dt = dt.os.final.test, 
          output = 'summary_test.txt')


km_plot(dt.os.final.train, dt.os.final.test, time = 'PFS', event = 'Progress', 
        variable = 'Radscore', outputBase = 'Radscore')


write_csv( dt.os.final.train, sprintf('%s_value_train.csv', 'Radscore'))
write_csv(dt.os.final.test, sprintf('%s_value_test.csv', 'Radscore'))



km_plot(dt.os.final.train, dt.os.final.test, time = 'PFS', event = 'Progress', 
        variable = 'nomo.cli', outputBase = 'Clinics')

write_csv( dt.os.final.train, sprintf('%s_value_train.csv', 'clinics'))
write_csv(dt.os.final.test, sprintf('%s_value_test.csv', 'clinics'))


km_plot(dt.os.final.train, dt.os.final.test, time = 'PFS', event = 'Progress', 
        variable = 'nomo.com', outputBase = 'Combined')

  
  





# AUC ---------------------------------------------------------------------

auc.os.Combined.train <- AUC.sh(Surv(dt.os.final.train$PFS, dt.os.final.train$Progress), 
                                Surv(dt.os.final.train$PFS, dt.os.final.train$Progress), 
                                dt.os.final.train$nomo.com, dt.os.final.train$nomo.com, 
                                times = seq(1, 8, by = 5))
auc.os.Combined.test <- AUC.sh(Surv(dt.os.final.train$PFS, dt.os.final.train$Progress), 
                               Surv(dt.os.final.test$PFS, dt.os.final.test$Progress), 
                               dt.os.final.train$nomo.com, dt.os.final.test$nomo.com, 
                               times = seq(1, 8, by = 5))

pdf('iAUC_PFS.pdf', width = 8, height = 6)
plot(auc.os.Combined.train, col = 'green', lty = 1)
plot(auc.os.Combined.test, col = 'red', lty = 1, add = T)


legend(x = 1, y = 0.3, col = c('green', 'red'), 
       lty = 1, 
       legend = c(sprintf('iAUC of Combined in training cohort: %.3f', auc.os.Combined.train$iauc), 
                  sprintf('iAUC of Combined in test cohort:        %.3f', auc.os.Combined.test$iauc)), 
       box.col = 'white')
dev.off()


# calibration -------------------------------------------------------------
# 1

cal.train.dt <- lapply(as.list(TimeInterval), calplt, formula = Surv(PFS, Progress)~nomo.com, 
       dt = dt.os.final.train, cohort = 'Train') %>% bind_rows()
cal.test.dt <- lapply(as.list(TimeInterval), calplt, formula = Surv(PFS, Progress)~nomo.com, 
                       dt = dt.os.final.train, cohort = 'Test') %>% bind_rows()

cal.nom.os.dt <- bind_rows(cal.train.dt, cal.test.dt)

cal.p.nom <- ggplot(data = cal.nom.os.dt, aes(x = mean.predicted, y = KM.corrected))
cal.p.nom + geom_line(aes(color = Model)) + geom_point(aes(color = Model)) + 
  xlim(0, 1) + ylim(0, 1) + xlab('Predicted event probability') + ylab('Observed event probability') + 
  geom_abline(slope = 1, intercept = 0, color = 'grey') + theme_bw() + theme(legend.position = c(0.8, 0.2))

ggsave('calibration_PFS_cli.pdf', width = 8, height = 8)

# calibration -------------------------------------------------------------
# 1

cal.train.dt <- lapply(as.list(TimeInterval), calplt, formula = Surv(PFS, Progress)~nomo.cli, 
                       dt = dt.os.final.train, cohort = 'Train') %>% bind_rows()
cal.test.dt <- lapply(as.list(TimeInterval), calplt, formula = Surv(PFS, Progress)~nomo.cli, 
                      dt = dt.os.final.train, cohort = 'Test') %>% bind_rows()

cal.nom.os.dt <- bind_rows(cal.train.dt, cal.test.dt)

cal.p.nom <- ggplot(data = cal.nom.os.dt, aes(x = mean.predicted, y = KM.corrected))
cal.p.nom + geom_line(aes(color = Model)) + geom_point(aes(color = Model)) + 
  xlim(0, 1) + ylim(0, 1) + xlab('Predicted event probability') + ylab('Observed event probability') + 
  geom_abline(slope = 1, intercept = 0, color = 'grey') + theme_bw() + theme(legend.position = c(0.8, 0.2))

ggsave('calibration_PFS_nomogram.pdf', width = 8, height = 8)


# DCA analysis ------------------------------------------------------------

#t.nom.dca <- dt.os.final.train %>% as.data.frame()

#ca.nomo <- stdca(data = dt.nom.dca, outcome = 'Progress', ttoutcome = 'PFS', 
#redictors = c('Radscore', 'nomo.cli', 'nomo.com'), probability = c(F, F, F), timepoint = 4, graph =T)

#t.dca <- tibble(Threshold = dca.nomo$net.benefit$threshold, 
#ll = dca.nomo$net.benefit$all, 
#None = dca.nomo$net.benefit$none,
#adscore = dca.nomo$net.benefit$Radscore, 
#linics = dca.nomo$net.benefit$nomo.cli, 
#ombined = dca.nomo$net.benefit$nomo.com)

#t.dca <- gather(dt.dca, key = 'Model', value = 'Net.Benefit', 2:6)
#t.dca$Model <- factor(dt.dca$Model, levels = c('All', 'None', 'Radscore',
#Clinics', 'Combined'))

#gplot(data = dt.dca, aes(x = Threshold, y = Net.Benefit)) + 
#eom_line(aes(color = Model), size = 1)+
#cale_color_manual(values = c('black', 'grey', 'green', 'red', 'blue')) + 
#heme_bw() + theme(legend.position = c(0.8, 0.8), 
#xis.text = element_text(size = 18), 
#xis.title = element_text(size = 18)) + ylim(-0.1, 0.5)

#gsave('decision_curve.pdf', width = 8, height = 6)

# IDI NRI -----------------------------------------------------------------
# 
# dt.idi.nri.train <- dt.os.train %>% select('Time', 'Event', 'PFS', 'Progress', 'UISS', 'Combined')
# dt.idi.nri.test <- dt.os.test %>% select('Time', 'Event', 'PFS', 'Progress', 'UISS', 'Combined')
# # sink('NRI.txt')
# # nricens(time = dt.idi.nri.train$Time, event = dt.idi.nri.train$Event, 
# #         p.std = dt.idi.nri.train$UISS, p.new = dt.idi.nri.train$Combined, cut = c(0.5, 1.5), 
# #         t0 = 4)
# # nricens(time = dt.idi.nri.train$PFS, event = dt.idi.nri.train$Progress, 
# #         p.std = dt.idi.nri.train$UISS, p.new = dt.idi.nri.train$Combined, cut = c(0.5, 1.5), 
# #         t0 = 4)
# # sink()
# 
# sink('nri_idi.txt')
# cat('Time train')
# dt.idi.nri.train$Event <- ifelse(dt.idi.nri.train$Time > 4, 0, dt.idi.nri.train$Event)
# 
# reclassification(data = as.data.frame(dt.idi.nri.train), cOutcome = 2, predrisk1 = dt.idi.nri.train$UISS/2,
#                  predrisk2 = dt.idi.nri.train$Combined/2, 
#                  cutoff = c(0, 0.25, 0.75, 1))
# lp.os.risk.uiss.test <- get.risk.coxph(cox.os.uiss.test, t0 = 4)
# lp.os.risk.Combined.test <- get.risk.coxph(cox.os.Combined.test, t0 = 4)
# 
# cat('Time test')
# dt.idi.nri.test$Event <- ifelse(dt.idi.nri.test$Time > 4, 0, dt.idi.nri.test$Event)
# 
# reclassification(data = as.data.frame(dt.idi.nri.test), cOutcome = 2, predrisk1 = dt.idi.nri.test$UISS/2,
#                  predrisk2 = dt.idi.nri.test$Combined/2, 
#                  cutoff = c(0, 0.25, 0.75, 1))
# 
# 
# dt.idi.nri.train$PFS <- ifelse(dt.idi.nri.train$PFS > 4, 0, dt.idi.nri.train$Progress)
# 
# cat('PFS train')
# reclassification(data = as.data.frame(dt.idi.nri.train), cOutcome = 4, predrisk1 = dt.idi.nri.train$UISS/2,
#                  predrisk2 = dt.idi.nri.train$Combined/2, 
#                  cutoff = c(0, 0.25, 0.75, 1))
# lp.pfs.risk.uiss.test <- get.risk.coxph(cox.pfs.uiss.test, t0 = 4)
# lp.pfs.risk.Combined.test <- get.risk.coxph(cox.pfs.Combined.test, t0 = 4)
# 
# cat('PFS test')
# 
# dt.idi.nri.test$Progress <- ifelse(dt.idi.nri.test$PFS > 4, 0, dt.idi.nri.test$Progress)
# 
# reclassification(data = as.data.frame(dt.idi.nri.test), cOutcome = 4, predrisk1 = dt.idi.nri.test$UISS/2,
#                  predrisk2 = dt.idi.nri.test$Combined/2, 
#                  cutoff = c(0, 0.25, 0.75, 1))
# 
# sink()

# AUC for cli---------------------------------------------------------------------

auc.os.Cli.train <- AUC.sh(Surv(dt.os.final.train$PFS, dt.os.final.train$Progress), 
                           Surv(dt.os.final.train$PFS, dt.os.final.train$Progress), 
                           dt.os.final.train$nomo.cli, dt.os.final.train$nomo.cli, 
                           times = seq(1, 8, by = 5))
auc.os.Cli.test <- AUC.sh(Surv(dt.os.final.train$PFS, dt.os.final.train$Progress), 
                          Surv(dt.os.final.test$PFS, dt.os.final.test$Progress), 
                          dt.os.final.train$nomo.cli, dt.os.final.test$nomo.cli, 
                          times = seq(1, 8, by = 5))

pdf('iAUC_PFS_cli.pdf', width = 8, height = 6)
plot(auc.os.Cli.train, col = 'green', lty = 1)
plot(auc.os.Cli.test, col = 'red', lty = 1, add = T)


legend(x = 1, y = 0.3, col = c('green', 'red'), 
       lty = 1, 
       legend = c(sprintf('iAUC of Rad_score in training cohort: %.3f', auc.os.Cli.train$iauc), 
                  sprintf('iAUC of Rad_score in test cohort:        %.3f', auc.os.Cli.test$iauc)), 
       box.col = 'white')
dev.off()
cuts_cli<-km_plot(dt.os.final.train, dt.os.final.test, time = 'PFS', event = 'Progress', 
                  variable = 'nomo.cli', outputBase = 'cli')


