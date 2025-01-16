#Part 1. Basic description
#Part 2.Cox regression main analysis with ASD as timevarying predictor
library(dplyr)
library(data.table)
library(lubridate)
library(gtsummary)
library(flextable)
library(splines)
library(survival)
library(IRanges)

###Part 1: Basic description
#1.1 for continuous variables with normal distribution
#prepare dataset: participant: consists of the continuous variables for calculating mean (SD)
#subgroup of ASD and non-ASD
table1<-participant%>%
  tbl_summary(by = ASD,digits=list(all_categorical() ~c(0,2),all_continuous()~2),statistic=list(all_continuous()~"{mean} ({sd})"))
#overall
table1_all<-participant%>%
  tbl_summary(digits=list(all_categorical() ~c(0,2),all_continuous()~2),statistic=list(all_continuous()~"{mean} ({sd})"))
table<-tbl_merge(tbls=list(table1_all,table1),tab_spanner = c("Total","Subgroup"))
table %>% as_flex_table() 

#1.2 for continuous variables with non-normal distribution
#prepare dataset: participant: consists of the continuous variables for calculating median (IQR): e.g., age onset
#subgroup of ASD and non-ASD
table1<-participant%>%
  tbl_summary(by = ASD)

#overall
table1_all<-participant%>%
  tbl_summary()
table<-tbl_merge(tbls=list(table1_all,table1),tab_spanner = c("Total","Subgroup"))
table %>% as_flex_table() %>%
  save_as_docx(table,path="Table1_Descriptives_median.docx")

###Part 2: Main analyses_Cox regression(timevarying ASD)
#2.1 prepare data for time-varying ASD 

#function to prepare data for Cox analysis
Cox_prep<-function(outcome_start){
  #Define event and survival time 
  outcome_data<-ASD_cardio%>% mutate(last_fup_date=pmin(.data[[outcome_start]],end_follow,na.rm=T))%>%
    mutate(outcome=(ifelse(is.na(.data[[outcome_start]]),0,1))) %>%
    mutate(os_days=as.numeric(difftime(last_fup_date,as.Date("2014-01-01"),units="days")))
#1: prepare dataset of time-independent covariates,event, and survival time
  indepdt_var<-outcome_data%>%
    select(RINPERSOON:birth_year,SEC,income,ADHD:any_psychiatry,outcome,os_days)
#2:Prepare dataset including ASD (time-dependent covariate),event, and survival time
  #if ASD was diagnosis before baseline, the ASD_time is 0; or else ASD_time is the days after baseline;
  #if no ASD dianosis, ASD_time is survival time.
  depdt_var<-outcome_data%>% mutate(ASD=ifelse(is.na(ASD_diag_start),0,1)) %>%
    mutate(ASD_time=(ifelse(is.na(ASD_diag_start),os_days,ifelse(ASD_diag_start<=as.Date("2014-01-01"),0,
                                                                 as.numeric(difftime(ASD_diag_start,as.Date("2014-01-01"),units="days")))))) %>%
    select(RINPERSOON,outcome,os_days,ASD,ASD_time)
  
#3: merge the datasets
  cox<-tmerge(data1=indepdt_var,data2=depdt_var,id=RINPERSOON,
                       event=event(os_days,outcome),exposure=tdc(ASD_time,ASD))
  cox$ASD_timevarying<-0
  cox$ASD_timevarying[!is.na(cox$exposure)]<-1
  cox<-cox %>% select(-exposure,-outcome) %>% mutate (tstart=tstart+(2014-birth_year)*365,tstop=tstop+(2014-birth_year)*365)
  return(cox)
}

#use function Cox_prep to get dataset for Cox regression with time-varying ASD for each outcome

CM_timevy<-Cox_prep("CM_diag_start")
diabetes_timevy<-Cox_prep("diabetes_start")
hypertension_timevy<-Cox_prep("hypertension_start")
dyslipidemia_timevy<-Cox_prep("dyslipidemia_start")
stroke_timevy<-Cox_prep("stroke_diag_start")
AP_timevy<-Cox_prep("AP_diag_start")
MI_timevy<-Cox_prep("MI_diag_start")
HF_timevy<-Cox_prep("HF_diag_start")

#save datasets
outcome_timevy<-list(CM_timevy=CM_timevy,diabetes_timevy=diabetes_timevy,hypertension_timevy=hypertension_timevy,dyslipidemia_timevy=dyslipidemia_timevy,stroke_timevy=stroke_timevy,AP_timevy=AP_timevy,MI_timevy=MI_timevy,HF_timevy=HF_timevy)
for (i in names(outcome_timevy)) {
  saveRDS(outcome_timevy[[i]],paste0("Processed_data/",i,".rds"))
}

#2.2 Cox regression analysis
#load data
data_timevy<-c("CM_timevy","diabetes_timevy","hypertension_timevy","dyslipidemia_timevy","stroke_timevy","AP_timevy","MI_timevy","HF_timevy")
for (i in 1:8) {
  assign(data_timevy[i],readRDS(paste0("Processed_data/",data_timevy[i],".rds")))
}
# functions for Model 1 and 2
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
  res <- res[c("1"),]
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
  res <- res[c("1"),]
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
  res <- res[c("1"),]
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
  res <- res[c("1"),]
  return(res)}

# Get results_main analysis
outcome_timevy<-list(CM_timevy=CM_timevy,diabetes_timevy=diabetes_timevy,hypertension_timevy=hypertension_timevy,dyslipidemia_timevy=dyslipidemia_timevy,stroke_timevy=stroke_timevy,AP_timevy=AP_timevy,MI_timevy=MI_timevy,HFAC_timevy=HFAC_timevy)
#Apply the cox model functions to different outcomes
model1_hr<-sapply(outcome_timevy,adjusted_model1)
model2_hr<-sapply(outcome_timevy,adjusted_model2)
cox_results <- rbind(model1_hr,model2_hr)
colnames(cox_results)<-c("any cardiometabolic disease","diabetes","hypertension","dyslipidemia","stroke","angina pectoris","myocardial infarction","heart failure")
t_cox_results<-t(cox_results)
colnames(t_cox_results)<-c("model1_HR (95% CI)","model1_p","model2_HR (95% CI)","model2_p")
write.csv(t_cox_results,"/Main_results.csv",row.names=T,quote=F)

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

