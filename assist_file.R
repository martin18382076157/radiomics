inter_test <- function(x, y)
{
  # browser()
  if(is.factor(x) | is.character(x))
  {
    x <- as.factor(x)
    res <- fisher.test(x = x, y = y)
  }
  else
  {
    x1 <- x[which(y == 0)]
    x2 <- x[which(y == 1)]
    
    res <- wilcox.test(x1, x2)
  }
  
  res$p.value
}

ulogit_test <- function(var_name, dt)
{
  fml <- as.formula(paste('Label', var_name, sep = '~'))
  fit <- glm(fml, data = dt, family = binomial)
  fit_sum <- summary(fit)
  tbl_tmp <- as_tibble(confint(fit))
  tbl_tmp$OR <- coefficients(fit)
  # browser()
  tbl_tmp <- exp(tbl_tmp)
  tbl_tmp$p_val <- fit_sum$coefficients[, 4]
  tbl_tmp
}


output_mlogit <- function(mod_com_final, filename)
{
  res_dt <- publish(mod_com_final)
  write_csv(res_dt$regressionTable, path = filename)
}
