################################################################################
######Body-weight gain after a diagnosis of ANCA vasculitis (2024)
#
##Descriptive statistics, linear regression, logistic regression and survival ax 
#
#Thomas French
################################################################################
#packages
library(readxl) #to read excel workbooks into R
library(dplyr) #for %>%  
library(plyr) #for ddply() function
library(pROC) #ROC curves 
library(CIplot) #plotting CI intervals 
library(survival) #for survival analysis 
library(ggsurvfit) #for kaplan meir plots 
library(ggplot2) #plotting data 
options(scipen=999) #turning off scientific notation
################################################################################
#import data & set seed 

Final_spreadsheet <- read_excel("Final spreadsheet.xlsx")
VSC <- Final_spreadsheet #rename for ease of access 

set.seed(123) #set seed for reproducibility
################################################################################
#preparing data for analysis 

#adjusting steroids at 0 months data as appropriate 
VSC$Steroids.at.0.months[VSC$Steroids.at.0.months == '?60'] <- 60
VSC$Steroids.at.0.months[VSC$Steroids.at.0.months == '?'] <- NA
VSC$Steroids.at.0.months[VSC$Steroids.at.0.months == '30(?)'] <- 30
VSC$Steroids.at.0.months[VSC$Steroids.at.0.months == '30 (?)'] <- 30
VSC$Steroids.at.0.months[VSC$Steroids.at.0.months == 'None (frailty)'] <- 0
VSC$Steroids.at.0.months[VSC$Steroids.at.0.months == '60(?)'] <- 60
VSC$Steroids.at.0.months[VSC$Steroids.at.0.months == 'NA (Died)'] <- 0
VSC$Steroids.at.0.months[VSC$Steroids.at.0.months == 'NA (relapses prior to edinburgh service)'] <- NA
VSC$Steroids.at.0.months[is.na(VSC$Steroids.at.0.months)] <- NA

#adjusting steroids at 3 months data as appropriate 
VSC$Steroid.at.3.months[VSC$Steroid.at.3.months == '?'] <- NA
VSC$Steroid.at.3.months[VSC$Steroid.at.3.months == '10(?)'] <- 10
VSC$Steroid.at.3.months[VSC$Steroid.at.3.months == 'Deceased'] <- 0
VSC$Steroid.at.3.months[VSC$Steroid.at.3.months == 'Died? (No follow-up)'] <- NA
VSC$Steroid.at.3.months[VSC$Steroid.at.3.months == 'NA (Died)'] <- NA
VSC$Steroid.at.3.months[VSC$Steroid.at.3.months == 'NA (lost to follow-up)'] <- NA
VSC$Steroid.at.3.months[VSC$Steroid.at.3.months == 'None (frailty)'] <- 0
VSC$Steroid.at.3.months[is.na(VSC$Steroid.at.3.months)] <- NA
table(VSC$Steroid.at.3.months)

#adjusting induction therapy regime as appropriate 
VSC$Induction.agent[VSC$Induction.agent == 'Both (MMF)'] <- 'Both'
VSC$Induction.agent[VSC$Induction.agent == 'Cyclophosphamide (MMF)'] <- "Cyclophosphamide"
VSC$Induction.agent[VSC$Induction.agent == 'Rituximab (MMF)'] <- 'Rituximab'
VSC$Induction.agent[VSC$Induction.agent == 'NA'] <- NA

################################################################################
#examining missing data
md.pattern(VSC) #examine missing data patterns

p_missing <- unlist(lapply(VSC, function(x) sum(is.na(x))))/nrow(VSC) #look at proportion of data missing in each category
sort(p_missing[p_missing > 0], decreasing = TRUE)




################################################################################
#store variables as factors or as continuous 

VSC$Steroids.at.0.months <- as.numeric(VSC$Steroids.at.0.months)
VSC$Steroid.at.3.months <- as.numeric(VSC$Steroid.at.3.months)
VSC$age_at_presentation <- as.numeric(VSC$age_at_presentation)
VSC$eGFR_CKDEPI <- as.numeric((VSC$eGFR_CKDEPI))
VSC$uPCR_presentation <- as.numeric(VSC$uPCR_presentation)

VSC$sex <- as.factor(VSC$sex)
VSC$ANCA_group <- as.factor(VSC$ANCA_group)
VSC$baseline_hyperlipidaemia <- as.factor(VSC$baseline_hyperlipidaemia)
VSC$baseline_DM <- as.factor(VSC$baseline_DM)
VSC$Induction.agent <- as.factor(VSC$Induction.agent)

################################################################################
#dummy coding 
VSC$sex <- ifelse(VSC$sex=='Male', 1, 0)
VSC$baseline_hyperlipidaemia <- ifelse(VSC$baseline_hyperlipidaemia=='yes', 1, 0)
VSC$baseline_DM <- ifelse(VSC$baseline_DM=='yes', 1, 0)
VSC$Rituximab <- ifelse(VSC$Induction.agent == 'Rituximab', 1, 0)
VSC$Cyclophosphamide <- ifelse(VSC$Induction.agent == 'Cyclophosphamide', 1, 0)
VSC$MMF_alone <- ifelse(VSC$Induction.agent == 'MMF', 1, 0)
VSC$Both <- ifelse(VSC$Induction.agent == 'Both', 1, 0)
VSC$MPO <- ifelse(VSC$ANCA_group == 'MPO', 1, 0)
VSC$PR3 <- ifelse(VSC$ANCA_group == 'PR3', 1, 0)
VSC$dual_pos <- ifelse(VSC$ANCA_group == 'dual_pos', 1, 0)
VSC$ANCA_neg <- ifelse(VSC$ANCA_group == 'ANCA_neg', 1, 0)
VSC$age_deciles <- VSC$age_at_presentation/10
VSC$eGFR_10mls <- VSC$eGFR_CKDEPI/10
VSC$censor_date <- as.Date("2024-01-01")
VSC$date_index <- as.Date(VSC$date_index)
VSC$BMI_change <- VSC$BMI_at_6m - VSC$BMI_baseline 
VSC$BMI_pct_change <- ((VSC$BMI_change/VSC$BMI_baseline)*100)
VSC$pct_GT5 <- ifelse(VSC$BMI_pct_change >= 5, 1, 0)
table(VSC$pct_GT5)
VSC$years_prior_to_censor <- (2024-VSC$year_presentation)
################################################################################
#basic descriptive statistics
library(arsenal)
mycontrols  <- tableby.control(test=TRUE, total=FALSE,
                               numeric.test="kwt", cat.test="chisq",
                               numeric.stats=c("N", "median", "q1q3"),
                               cat.stats=c("N","countpct"),
                               stats.labels=list(N='Count', median='Median', q1q3='Q1,Q3'))

VSC$sex <- as.factor(VSC$sex)
VSC$baseline_hyperlipidaemia <- as.factor(VSC$baseline_hyperlipidaemia)
VSC$baseline_DM <- as.factor(VSC$baseline_DM)

tab <- tableby(pct_GT5 ~age_deciles + sex +
                 ANCA_group + eGFR_10mls 
               + BMI_baseline +   
                 Induction.agent + years_prior_to_censor + 
                 baseline_DM + baseline_hyperlipidaemia + 
                 Steroid.at.3.months + Steroids.at.0.months + 
                 uPCR_presentation, data=VSC, control=mycontrols)

summary(tab)
################################################################################
#checking distributions, correlations
ggplot(VSC, aes(x=age_at_presentation)) + 
  geom_density() #left-tail 

ggplot(VSC, aes( y=age_at_presentation)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=2, notch=FALSE)
sort(VSC$age_at_presentation)


ggplot(VSC, aes(x=BMI_at_6m)) + 
  geom_density() #right-tail 
ggplot(VSC, aes( y=BMI_at_6m)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=2, notch=FALSE)


ggplot(VSC, aes(x=Steroids.at.0.months)) + 
  geom_density() 

ggplot(VSC, aes( y=Steroids.at.0.months)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=2, notch=FALSE)


ggplot(VSC, aes(x=Steroid.at.3.months)) + 
  geom_density() 
ggplot(VSC, aes( y=Steroid.at.3.months)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=2, notch=FALSE)


ggplot(VSC, aes(x=BMI_pct_change)) + 
  geom_density()

ggplot(VSC, aes( y=BMI_pct_change)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=2, notch=FALSE)

ggplot(VSC, aes(x=eGFR_CKDEPI)) + 
  geom_density()  
ggplot(VSC, aes( y=eGFR_CKDEPI)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=2, notch=FALSE)

ggplot(VSC, aes(x=years_prior_to_censor, y=BMI_pct_change, color=pct_GT5)) + 
  geom_point(size=6) + scale_x_reverse()
abline(lm(pct_GT5~years_prior_to_censor, data=VSC), col='red')

#relationships between numeric variables 
pairs <- VSC %>% select(age_at_presentation, eGFR_CKDEPI, Steroids.at.0.months,
                        Steroid.at.3.months, BMI_pct_change, years_prior_to_censor,
                        BMI_baseline)

pairs(pairs, panel=function(x,y){
  points(x,y)
  abline(lm(y~x), col='red')})

################################################################################
#describing relationship between variables
VSC$age_2 <- VSC$age_at_presentation *VSC$age_at_presentation
agemodel <- lm(BMI_pct_change~age_at_presentation, data=VSC)
age_quadratic_model <- lm(BMI_pct_change~age_at_presentation+age_2, data=VSC)
summary(age_quadratic_model)
agevalues <- seq(0,90, 0.1)
agepredict <- predict(age_quadratic_model,list(age_at_presentation=agevalues, age_2=agevalues^2))
plot(VSC2$age_at_presentation, VSC2$BMI_pct_change,pch=16)
lines(agevalues, agepredict, col='blue')

plot(VSC$years_prior_to_censor, VSC$BMI_pct_change)
abline(lm(BMI_pct_change~years_prior_to_censor, data=VSC))
years_model <- lm(BMI_pct_change~years_prior_to_censor, data=VSC)
summary(years_model)

plot(VSC$BMI_baseline, VSC$BMI_pct_change)
abline(lm(BMI_pct_change~BMI_baseline, data=VSC))
################################################################################
#removing outliers from the outcome variable in order to approximate a normal distribution
testing <- VSC[complete.cases(VSC[ , c('BMI_pct_change')]), ]
Q1 <- quantile(testing$BMI_pct_change, .25)
Q3 <- quantile(testing$BMI_pct_change, .75)
IQR <- IQR(testing$BMI_pct_change)
outliers <- subset(testing, testing$BMI_pct_change<(Q1 - 1.5*IQR) | testing$BMI_pct_change>(Q3 + 1.5*IQR))
outliers <- sort(outliers$BMI_pct_change) #identified seven 
VSC2 = subset(VSC, !(BMI_pct_change %in% outliers))
################################################################################
#univariate linear regression 
df <- data.frame("Characteristic", "coefficient", "2.5%", "97.5%", "p_value")


list  <- c('MPO', 'PR3', 'dual_pos', 'eGFR_10mls', 'baseline_hyperlipidaemia', 'Steroid.at.3.months', 'MMF_alone', 'Cyclophosphamide', 
           'baseline_DM', 'BMI_baseline', 'Steroids.at.0.months', 'Both','sex','uPCR_presentation', 'years_prior_to_censor')

for (var in list) {
  
  formula <- as.formula(paste("BMI_pct_change ~", var))
  model1 <- lm(formula, data=VSC2)
  output = c(var, model1$coefficients[2], confint(model1, level=0.95)[2], confint(model1, level=0.95)[4], coef(summary(model1))[,4][2])
  df = rbind(df, output)
  
  
  
}
df
View(df)

################################################################################
#backwards stepwise linear regression 

all <- lm(BMI_pct_change~ age_at_presentation + age_2 + eGFR_10mls + 
            years_prior_to_censor + Steroids.at.0.months + Steroid.at.3.months + 
            sex + MPO + PR3 + dual_pos + BMI_baseline + Cyclophosphamide + 
            MMF_alone + Both, data=VSC2)

backward <- step(all, direction='backward', scope=formula(all), trace=0)
backward$anova
backward$coefficients 

#adjusted model 
adjusted_model <- lm (BMI_pct_change ~ age_at_presentation + age_2 + eGFR_10mls +
                        years_prior_to_censor + Steroids.at.0.months + BMI_baseline + 
                        Cyclophosphamide + MMF_alone + Both, data=VSC2)
summary(adjusted_model)

string <- summary(adjusted_model)
coefs <- string$coefficients[2:10]
p_vals <- string$coefficients[32:40]
lower_CI <- confint(adjusted_model, level=0.95)[2:10]
upper_CI <- confint(adjusted_model, level=0.95)[12:20]
results <- data.frame(coefs, lower_CI, upper_CI, p_vals)
results 


df <- data.frame("Characteristic", "coefficient", "2.5%", "97.5%", "p_value")


list  <- c('MPO', 'PR3', 'dual_pos', 'baseline_hyperlipidaemia', 'Steroid.at.3.months', 
           'baseline_DM','sex','uPCR_presentation')

for (var in list) {
  
  formula <- as.formula(paste("BMI_pct_change ~ age_at_presentation + age_2 + eGFR_10mls +
                        years_prior_to_censor + Steroids.at.0.months + BMI_baseline + 
                        Cyclophosphamide + MMF_alone + Both +", var))
  model1 <- lm(formula, data=VSC2)
  str <- summary(model1)
  output = c(var, str$coefficients[11], confint(model1, level=0.95)[11], confint(model1, level=0.95)[22], str$coefficients[44])
  df = rbind(df, output)
  
  
  
}
df
View(df)

################################################################################
plot(adjusted_model) #zred/zpred plot, Q-Q plot of residuals etc 
################################################################################
#logistic regression 
################################################################################
#univariate logistic regression 

df <- data.frame("Characteristic", "Odds_ratio", "2.5%", "97.5%", "p_value")


list  <- c('MPO', 'PR3', 'dual_pos', 'eGFR_10mls', 'baseline_hyperlipidaemia', 'Steroid.at.3.months', 'MMF_alone', 'Cyclophosphamide', 
           'baseline_DM', 'BMI_baseline', 'age_deciles', 'Steroids.at.0.months', 'Both','sex','uPCR_presentation', 'years_prior_to_censor')

for (var in list) {
  
  formula <- as.formula(paste("obese_at6months ~", var))
  model1 <- glm(formula, data=VSC, family='binomial')
  output = c(var, exp(model1$coefficients[2]), exp(confint.default(model1))[2], exp(confint.default(model1))[4], coef(summary(model1))[,4][2])
  df = rbind(df, output)
  
  
  
}
df
View(df)


################################################################################
#multivariable logistic regression model 

fit <- glm(obese_at6months ~ age_deciles + sex +
             MPO + PR3 + dual_pos + eGFR_10mls 
           + BMI_baseline +   
             Cyclophosphamide + MMF_alone + Both + years_prior_to_censor + 
             Steroids.at.0.months + Steroid.at.3.months, family=binomial(link="logit"), data=VSC)
summary(fit)
coefs <- c(exp(fit$coefficients)[2:14])

lower_CI <- c(exp(confint.default(fit))[2:14])
upper_CI <- c(exp(confint.default(fit))[16:28])
p_values <- c(coef(summary(fit))[,4][2:14])
results <- data.frame(coefs, lower_CI, upper_CI, p_values)
View(results)
results
################################################################################
#Assessing logistic regression model 
#calculating ROC curve 
library(pROC)
prob=predict(fit, type="response")

VSC2 <- VSC[complete.cases(VSC[ , c('obese_at6months', 'sex', 'MPO', 'PR3', 'dual_pos', 'eGFR_10mls', 'Steroid.at.3.months', 'MMF_alone', 'Cyclophosphamide', 
                                    'BMI_baseline', 'Steroids.at.0.months', 'Both', 'years_prior_to_censor', 'age_deciles')]), ]

VSC2$prob=prob
g <- roc(obese_at6months ~ prob, data = VSC2)

ROC <- plot(g)

ROC
auc(g)

#calculating variance inflation factor for independent variables 
vif_values <- vif(fit)
barplot(vif_values, main = "VIF Values", horiz = TRUE, col = "#8B1A1A")
################################################################################
#calculating ORs for variables not included in adjusted logistic regression model 

df <- data.frame("Characteristic", "Odds_ratio", "2.5%", "97.5%", "p_value")
list  <- c('uPCR_presentation', 'uACR_presentation', 'baseline_hyperlipidaemia', 
           'baseline_DM')

fit <- glm(pct_GT5 ~ age_deciles + sex +
             MPO + PR3 + dual_pos + eGFR_10mls 
           + BMI_baseline +   
             Cyclophosphamide + MMF_alone + Both + years_prior_to_censor + 
             Steroids.at.0.months + Steroid.at.3.months + uPCR_presentation, family=binomial(link="logit"), data=VSC)
summary(fit)

output <- c("uPCR_presentation", exp(fit$coefficients)[15], exp(confint.default(fit))[15], exp(confint.default(fit))[30], coef(summary(fit))[,4][15])
df = rbind(df, output)
View(df)

################################################################################
#Survival analysis 
################################################################################
#Computing time to event variable and coding 'event' outcome 
VSC$TTE_unrelapsed <- VSC$censor_date - VSC$date_index
VSC$TTE_unrelapsed[VSC$date_first_relapse!='None'] <- 0
VSC$TTE_dead <- VSC$days_to_death
VSC$TTE_dead[VSC$TTE_dead=='N/A'] <- 0
VSC$TTE_unrelapsed[VSC$TTE_dead != 0] <- 0 #if they died make that the censor date 
VSC$TTE_relapse <- VSC$days_to_first_relapse
VSC$TTE_relapse[VSC$TTE_relapse=='N/A'] <- 0
VSC$TTE_dead[VSC$TTE_relapse != 0] <- 0
VSC$TTE_total <- as.numeric(VSC$TTE_dead) + as.numeric(VSC$TTE_relapse) + as.numeric(VSC$TTE_unrelapsed)
VSC$TTE_total_months <- VSC$TTE_total/30
VSC$TTE_total_years <- VSC$TTE_total_months/12
VSC$event <- ifelse(VSC$date_first_relapse != 'None', 1, 0)
VSC$BMI_GET_30 <- ifelse(VSC$BMI_baseline >= 30, "Obese", "Not obese")

################################################################################
#cox ph regression 
#First, visually checking proportional hazards assumption for each variable 
ggsurvfit(survfit(Surv(VSC$TTE_total, VSC$event == 1)~ VSC$Induction.agent)) + coord_cartesian(xlim=c(0, 3000))
ggsurvfit(survfit(Surv(VSC$TTE_total, VSC$event == 1)~ VSC$ANCA_group)) + coord_cartesian(xlim=c(0, 3000))
ggsurvfit(survfit(Surv(VSC$TTE_total, VSC$event == 1)~ VSC$age_categories)) + coord_cartesian(xlim=c(0, 3000))
ggsurvfit(survfit(Surv(VSC$TTE_total, VSC$event == 1)~ VSC$GFR_categories)) + coord_cartesian(xlim=c(0, 3000))
ggsurvfit(survfit(Surv(VSC$TTE_total, VSC$event == 1)~ VSC$Steroids.at.0.months)) + coord_cartesian(xlim=c(0, 3000))


#univariate analysis first 
covariates <- c('MPO', 'PR3', 'dual_pos', 'eGFR_10mls', 'baseline_hyperlipidaemia', 'Steroid.at.3.months', 'MMF_alone', 'Cyclophosphamide', 
                'baseline_DM', 'BMI_GET_30', 'age_deciles', 'Steroids.at.0.months', 'Both', 'years_prior_to_censor')

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(TTE_total, event)~', x)))
univ_models <- lapply( univ_formulas, function(x){coxph(x, data = VSC)})
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         wald.test<-signif(x$wald["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, wald.test, p.value)
                         names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                       "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res <- t(as.data.frame(univ_results, check.names = FALSE))
as.data.frame(res)

################################################################################
#multivariate analysis 
res.cox <- coxph(Surv(TTE_total, event) ~ age_deciles + sex + eGFR_10mls + 
                   + BMI_GET_30 + MPO + PR3 + dual_pos + MMF_alone + Cyclophosphamide + Both + 
                   years_prior_to_censor      
                 , data =  VSC)
summary(res.cox)
################################################################################
#kaplan meier plots 
Y = Surv(VSC$TTE_total_years, VSC$event == 1)
kmfit = survfit(Y ~ VSC$BMI_GET_30)
p1 <- survdiff(formula = Surv(TTE_total, event) ~ BMI_GET_30, data=VSC)
p1
km <- ggsurvfit(kmfit, lwd=1.5) + add_confidence_interval(type='lines', alpha=0.2, lwd=1) + 
  scale_ggsurvfit() +
  add_legend_title(title="Obesity status at baseline: ") + 
  coord_cartesian(xlim = c(0, 5), ylim =c(0.5, 1)) +
  labs(x="Time from index date (years)",
       y="Freedom from relapse",
       caption="Log-rank p-value: 0.900") + 
  scale_color_manual(values=c("red", "black"), labels=c("Not obese", "Obese")) + 
  scale_fill_manual(values=c("red", "black"), labels=c("Not obese", "Obese")) 

km + theme(
  plot.title = element_text(size=20, hjust=0.5),
  axis.text = element_text(size=12),
  axis.title =element_text(size=15),
  legend.position="bottom",
  legend.key.size = unit(1, "cm"),
  plot.margin=unit(c(1,1,0.5,1), 'cm'), 
  axis.line = element_line(color='black'),
  panel.border = element_blank(),
  plot.caption = element_text(hjust = 0.2, size=15, vjust=90, face="italic")
)  + add_censor_mark(shape="|", size=0.01, stroke=10)  



km2 <- ggsurvfit(kmfit, lwd=1.2, type="cumhaz") + add_confidence_interval(type='lines', alpha=0.2, lwd=1) + 
  scale_ggsurvfit() +
  add_legend_title(title="Obesity status at baseline: ") + 
  coord_cartesian(xlim = c(0, 5), ylim =c(0, 1)) +
  labs(,
       x="Time from index date (years)",
       y="Cumulative risk of relapse",
       caption="Log-rank p-value: 0.900") + 
  scale_color_manual(values=c("red", "black"), labels=c("Not obese", "Obese")) + 
  scale_fill_manual(values=c("red", "black"), labels=c("Not obese", "Obese")) 

km2 + theme(
  plot.title = element_text(size=20, hjust=0.5),
  axis.text = element_text(size=12),
  axis.title =element_text(size=15),
  legend.position="bottom",
  legend.key.size = unit(1, "cm"),
  plot.margin=unit(c(1,1,0.5,1), 'cm'), 
  axis.line = element_line(color='black'),
  panel.border = element_blank(),
  plot.caption = element_text(hjust = 0.2, size=15, vjust=90, face="italic")
)  + add_censor_mark(shape="|", size=0.01, stroke=10)  





