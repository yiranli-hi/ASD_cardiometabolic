###Part 1. Basic description
###Part 2.Cox regression main analysis with ASD as timevarying predictor

library(dplyr)
library(data.table)
library(lubridate)
library(gtsummary)
library(flextable)
library(splines)
library(survival)
library(IRanges)

###Part 1: Basic description
vars_median <- c("age_onset")         # non-normal continuous vars
vars_factor <- c("ASD", "sex", "SEC", "birth_year_cat", "depression", "anxiety", "personality", "ADHD", "SUD", "schizophrenia", "bipolar", "eating", "intellectual") # categorical variables

# Standardize variable types
participant <- participant %>%
  mutate(
    across(all_of(vars_factor), ~ as.factor(.x))
  )
tbl1 <- participant %>%
  tbl_summary(
    by = ASD,
    type = list(
      all_categorical() ~ "categorical",
      all_continuous()  ~ "continuous",
      any_of(vars_median) ~ "continuous2"
    ),
    statistic = list(
      all_continuous()    ~ "{mean} ({sd})",
      any_of(vars_median) ~ "{median} ({p25}, {p75})",
      all_categorical()   ~ "{n} / {N} ({p}%)"
    ),
    digits = list(
      all_continuous()    ~ 2,
      any_of(vars_median) ~ 2,
      all_categorical()   ~ 0   
    ),
    missing = "ifany"
  ) %>%
  add_overall()

save_as_docx(
  table = tbl1 %>% as_flex_table(),
  path  = "Results/Table1_Descriptives.docx"
)
###Part 2: Main analyses_Cox regression(timevarying ASD)
# load data
data_timevy<-c("CM_timevy","diabetes_timevy","hypertension_timevy","dyslipidemia_timevy","stroke_timevy","AP_timevy","MI_timevy","HF_timevy")
for (i in 1:8) {
  assign(data_timevy[i],readRDS(paste0("Processed_data/",data_timevy[i],".rds")))
}
# define results format
myformat <- function(x, digit = 3, sep = TRUE) {
  format(round(as.numeric(x), digit), nsmall = digit, trim = TRUE, big.mark = ifelse(sep, ",", ""))}

# functions for Models 1, 2, and 3
# Model 1: adjusted sex, birth_year
adjusted_model1 <- function(df){
  results <- summary(coxph(Surv(tstart, tstop, event) ~ ASD_timevarying + sex + birth_year, data = df))
  p.value<- myformat(results$coefficients[,"Pr(>|z|)"])
  HR <- myformat(results$conf.int[,"exp(coef)"])
  HR.confint.lower <- myformat(results$conf.int[,"lower .95"])
  HR.confint.upper <- myformat(results$conf.int[,"upper .95"])
  HR <- paste0(HR, " (",HR.confint.lower, "-", HR.confint.upper, ")")
  res<-as.data.frame(cbind(HR, p.value))
  names(res)<-c("HR (95% CI)", "p value")
  # Select the ASD row by variable name
  res <- res[rownames(res) == "ASD_timevarying", , drop = FALSE]
  return(res)}

# Model 2: adjusted sex, birth_year, SEC, income
adjusted_model2 <- function(df){
  results <- summary(coxph(Surv(tstart, tstop, event) ~ ASD_timevarying + sex + birth_year + SEC +income, data = df))
  p.value<- myformat(results$coefficients[,"Pr(>|z|)"])
  HR <- myformat(results$conf.int[,"exp(coef)"])
  HR.confint.lower <- myformat(results$conf.int[,"lower .95"])
  HR.confint.upper <- myformat(results$conf.int[,"upper .95"])
  HR <- paste0(HR, " (",HR.confint.lower, "-", HR.confint.upper, ")")
  res<-as.data.frame(cbind(HR, p.value))
  names(res)<-c("HR (95% CI for HR)", "p.value")
  # Select the ASD row by variable name
  res <- res[rownames(res) == "ASD_timevarying", , drop = FALSE]
  return(res)}

##Model 3 (eTable 3 in paper), adjusting each psychiatric comorbidity (here, two of the psychiatric comorbidities were shown as examples)
# Function for model 3a: adjusted sex, birth_year,  SEC, income, depressive disorder
adjusted_model3 <- function(df){
  results <- summary(coxph(Surv(tstart, tstop, event) ~ ASD_timevarying + sex + birth_year + SEC + income + depression, data = df))
  p.value<- myformat(results$coefficients[,"Pr(>|z|)"])
  HR <- myformat(results$conf.int[,"exp(coef)"])
  HR.confint.lower <- myformat(results$conf.int[,"lower .95"])
  HR.confint.upper <- myformat(results$conf.int[,"upper .95"])
  HR <- paste0(HR, " (",HR.confint.lower, "-", HR.confint.upper, ")")
  res<-as.data.frame(cbind(HR, p.value))
  names(res)<-c("HR (95% CI for HR)", "p.value")
  # Select the ASD row by variable name
  res <- res[rownames(res) == "ASD_timevarying", , drop = FALSE]
  return(res)}

# Function for model 3b: adjusted sex, birth_year,  SEC, income, anxiety
adjusted_model4 <- function(df){
  results <- summary(coxph(Surv(tstart, tstop, event) ~ ASD_timevarying + sex + birth_year + SEC + income + anxiety, data = df))
  p.value<- myformat(results$coefficients[,"Pr(>|z|)"])
  HR <- myformat(results$conf.int[,"exp(coef)"])
  HR.confint.lower <- myformat(results$conf.int[,"lower .95"])
  HR.confint.upper <- myformat(results$conf.int[,"upper .95"])
  HR <- paste0(HR, " (",HR.confint.lower, "-", HR.confint.upper, ")")
  res<-as.data.frame(cbind(HR, p.value))
  names(res)<-c("HR (95% CI for HR)", "p.value")
  # Select the ASD row by variable name
  res <- res[rownames(res) == "ASD_timevarying", , drop = FALSE]
  return(res)}

# Get results_main analysis
outcome_timevy<-list(CM_timevy=CM_timevy,diabetes_timevy=diabetes_timevy,hypertension_timevy=hypertension_timevy,dyslipidemia_timevy=dyslipidemia_timevy,stroke_timevy=stroke_timevy,AP_timevy=AP_timevy,MI_timevy=MI_timevy,HF_timevy=HF_timevy)
#Apply the cox model functions to different outcomes
model1_hr<-sapply(outcome_timevy,adjusted_model1)
model2_hr<-sapply(outcome_timevy,adjusted_model2)
cox_results <- rbind(model1_hr,model2_hr)
colnames(cox_results)<-c("any cardiometabolic disease","diabetes","hypertension","dyslipidemia","stroke","angina pectoris","myocardial infarction","heart failure")
t_cox_results<-t(cox_results)
colnames(t_cox_results)<-c("model1_HR (95% CI)","model1_p","model2_HR (95% CI)","model2_p")
write.csv(t_cox_results,"Main_results.csv",row.names=T,quote=F)

# Get results_adjusting each psychiatry comorbidity
#Apply the cox model functions to different outcomes
model3_hr<-sapply(outcome_timevy,adjusted_model3)
model4_hr<-sapply(outcome_timevy,adjusted_model4)

cox_results <- rbind(model3_hr,model4_hr)
colnames(cox_results)<-c("any cardiometabolic disease","diabetes","hypertension","dyslipidemia","stroke","angina pectoris","myocardial infarction","heart failure")
t_cox_results<-t(cox_results)
colnames(t_cox_results)<-c("model3_HR (95% CI)","model3_p","model4_HR (95% CI)","model4_p")
write.csv(t_cox_results,"Model3_results.csv",row.names=T,quote=F)
########################END##############################
