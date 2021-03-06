# -------------------------------------- #
#           Spat21/Mozzie Study          #
#     Run final computational model      #
#               CSP Target               # 
#                 Aim 2                  #
#            Mozzie Phase 1              #
#               K. Sumner                #
#           February 18, 2020            #
# -------------------------------------- #

# good resource for trouble shooting convergence problems
# https://rstudio-pubs-static.s3.amazonaws.com/33653_57fc7b8e5d484c909b615d8633c01d51.html


#### -------- load packages ------------ ####

# load in the packages of interest
library(dplyr)
library(readr)
library(tidyr)
library(betareg)
library(lme4)
library(ggplot2)
library(sjstats)
library(lmerTest)
library(glmmTMB)


#### ----- read in the data sets ----- ####

# read in the combined ama and csp data set for mosquito abdomens
model_data = read_rds("Desktop/Dissertation Materials/SpatialR21 Grant/Final Dissertation Materials/Aim 2/clean_ids_haplotype_results/AMA_and_CSP/final/model data/final_model_data/spat21_aim2_merged_data_with_weights_5MAR2020.rds")



#### ------ check covariate coding ------- ####

# rescale csp and ama moi
model_data$csp_moi_rescaled = scale(model_data$csp_moi)
model_data$ama_moi_rescaled = scale(model_data$ama_moi)

# subset the data set to samples that passed pfcsp sequencing only
model_data = model_data %>%
  filter(!(is.na(csp_haps_shared)))

# check the covariates
str(model_data$sample_id_human)
str(model_data$HH_ID_human)
str(model_data$unq_memID)
str(model_data$p_te_all_csp)
str(model_data$age_cat_baseline)
model_data$age_cat_baseline = as.factor(model_data$age_cat_baseline)
str(model_data$village_name)
model_data$village_name = as.factor(model_data$village_name)
model_data$village_name = relevel(model_data$village_name,ref = "Maruti")
str(model_data$mosquito_week_count_cat)
model_data$mosquito_week_count_cat = as.factor(model_data$mosquito_week_count_cat)
str(model_data$pfr364Q_std_combined_rescaled)
str(model_data$aim2_exposure)
model_data$aim2_exposure = as.factor(model_data$aim2_exposure)
model_data$aim2_exposure = relevel(model_data$aim2_exposure,ref = "symptomatic infection")
str(model_data$mean_moi_category)



#### ------ make some plots of covariates ------- ####

# make a plot of p_te_all over the exposure
ggplot(model_data, aes(x = p_te_all_csp)) + geom_density() + facet_wrap(~aim2_exposure)
ggplot(model_data, aes(x = pfr364Q_std_combined_rescaled)) + geom_density() + facet_wrap(~aim2_exposure)
ggplot(model_data, aes(x = age_cat_baseline)) + geom_density() + facet_wrap(~aim2_exposure)
ggplot(model_data, aes(x = csp_moi_rescaled)) + geom_density() + facet_wrap(~aim2_exposure)
ggplot(model_data, aes(x = mosquito_week_count_cat)) + geom_density() + facet_wrap(~aim2_exposure)
ggplot(model_data, aes(x = village_name)) + geom_density() + facet_wrap(~aim2_exposure)
ggplot(model_data, aes(x = HH_ID_human)) + geom_density() + facet_wrap(~aim2_exposure)
table(model_data$HH_ID_human,model_data$aim2_exposure)
table(model_data$unq_memID,model_data$aim2_exposure)
table(model_data$sample_id_human,model_data$aim2_exposure)
ggplot(model_data, aes(x=human_date,y = mosquito_week_count)) + geom_line()



#### ------ run the final models and do model selection ------- ####

# run the original multi-level model with all covariates and interaction term
model1 <- glmmTMB(p_te_all_csp~aim2_exposure+age_cat_baseline+mosquito_week_count_cat +pfr364Q_std_combined_rescaled+village_name+csp_moi_rescaled+aim2_exposure*age_cat_baseline+(1|HH_ID_human/unq_memID),family=binomial(link = "logit"), data = model_data)
summary(model1)
exp(-0.01347)
exp(0.20624)

# run the model with all covariates but interaction removed
model2 <- glmmTMB(p_te_all_csp~aim2_exposure+age_cat_baseline+mosquito_week_count_cat +pfr364Q_std_combined_rescaled+village_name+csp_moi_rescaled+(1|HH_ID_human/unq_memID),family=binomial(link = "logit"), data = model_data)
summary(model2)
exp(0.37386)
anova(model1,model2)

# now run the model removing village
model3 <- glmmTMB(p_te_all_csp~aim2_exposure+age_cat_baseline+mosquito_week_count_cat +pfr364Q_std_combined_rescaled+csp_moi_rescaled+(1|HH_ID_human/unq_memID),family=binomial(link = "logit"), data = model_data)
summary(model3)
exp(0.36189)
anova(model2,model3) 
anova(model3,model1)
(exp(0.37386)-exp(0.36189))/exp(0.37386)

# now run the model removing parasite density but adding back in village
model4 <- glmmTMB(p_te_all_csp~aim2_exposure+age_cat_baseline+mosquito_week_count_cat +village_name+csp_moi_rescaled+(1|HH_ID_human/unq_memID),family=binomial(link = "logit"), data = model_data)
summary(model4)
anova(model2,model4) 
anova(model3,model4)
anova(model1,model4)
exp(0.36558)
(exp(0.37386)-exp(0.36558))/exp(0.37386)

# now run the model removing age
model5 <- glmmTMB(p_te_all_csp~aim2_exposure+mosquito_week_count_cat +village_name+csp_moi_rescaled+(1|HH_ID_human/unq_memID),family=binomial(link = "logit"), data = model_data)
summary(model5)
exp(0.35012)
anova(model2,model5) 
(exp(0.37386)-exp(0.35012))/exp(0.37386)

# now run the model removing mosquito_week_count
model6 <- glmmTMB(p_te_all_csp~aim2_exposure+age_cat_baseline +village_name+csp_moi_rescaled+(1|HH_ID_human/unq_memID),family=binomial(link = "logit"), data = model_data)
summary(model6)
exp(0.32707)
anova(model2,model6)
(exp(0.37386)-exp(0.32707))/exp(0.37386)

# now run the model removing mean moi
model_data_subset = model_data %>%
  filter(!(is.na(csp_moi_rescaled)))
model7 <- glmmTMB(p_te_all_csp~aim2_exposure+age_cat_baseline+mosquito_week_count_cat +village_name+pfr364Q_std_combined_rescaled+(1|HH_ID_human/unq_memID),family=binomial(link = "logit"), data = model_data_subset)
summary(model7)
exp(0.50159)
anova(model2,model7) # model 2 is better
(exp(0.37386)-exp(0.50159))/exp(0.37386)

# now run the crude model with no covariates
model8 <- glmmTMB(p_te_all_csp~aim2_exposure+(1|HH_ID_human/unq_memID),family=binomial(link = "logit"), data = model_data_subset)
summary(model8)
exp(0.4708)
anova(model2,model8) # model 2 is better
# model 2 is better but this one didn't have convergence issues
 
# a few more model comparisons
anova(model4,model5) # model 4 better than 5
anova(model3,model5)
anova(model5,model6) # model 5 better
anova(model5,model7) # model 5 better
anova(model5,model8) # model 5 better


# summary:
# deciding to go with model 2 but will remove csp moi


#### ------ create a forest plot of the final model output ------ ####

# run model 2 output without csp moi
model2 <- glmmTMB(p_te_all_csp~aim2_exposure+pfr364Q_std_combined_rescaled+age_cat_baseline+mosquito_week_count_cat+village_name+(1|HH_ID_human/unq_memID),family=binomial(link = "logit"), data = model_data)
summary(model2)
performance::icc(model2)
exp(confint(model2,method="Wald"))

# run model without being multilevel
model2_all_nmlm <- glm(p_te_all_csp~aim2_exposure+age_cat_baseline+mosquito_week_count_cat+pfr364Q_std_combined_rescaled+village_name,family=binomial(link = "logit"), data = model_data)
summary(model2_all_nmlm)
exp(0.432712)
exp(confint(model2_all_nmlm,method="Wald"))

table1 = exp(confint(model2,method="Wald"))
estimates = c(table1[2,3],NA,table1[3,3],NA,NA,NA,table1[5,3],table1[4,3],NA,table1[6,3],NA,NA,NA,table1[7,3],table1[8,3])
lower_ci = c(table1[2,1],NA,table1[3,1],NA,NA,NA,table1[5,1],table1[4,1],NA,table1[6,1],NA,NA,NA,table1[7,1],table1[8,1])
upper_ci = c(table1[2,2],NA,table1[3,2],NA,NA,NA,table1[5,2],table1[4,2],NA,table1[6,2],NA,NA,NA,table1[7,2],table1[8,2])
names = c("Asymptomatic infection","","Parasite density (parasite/uL whole blood)"," ","Participant age           ","<5 years (REF)","5-15 years",">15 years","  ","High mosquito abundance","   ","Village                ","Maruti (REF)","Kinesamo","Sitabicha")
forest_plot_df = data.frame(names,estimates,lower_ci,upper_ci)
forest_plot_df$names = factor(forest_plot_df$names, levels = c("Asymptomatic infection","","Parasite density (parasite/uL whole blood)"," ","Participant age           ","<5 years (REF)","5-15 years",">15 years","  ","High mosquito abundance","   ","Village                ","Maruti (REF)","Kinesamo","Sitabicha"))
forest_plot_df$names = ordered(forest_plot_df$names, levels = c("Asymptomatic infection","","Parasite density (parasite/uL whole blood)"," ","Participant age           ","<5 years (REF)","5-15 years",">15 years","  ","High mosquito abundance","   ","Village                ","Maruti (REF)","Kinesamo","Sitabicha"))

# create a forest plot
library(forcats)
fp <- ggplot(data=forest_plot_df, aes(x=fct_rev(names), y=estimates, ymin=lower_ci, ymax=upper_ci)) +
  geom_pointrange(size=c(3,1,1,1,1,1,1,1,1,1,1,1,1,1,1),colour=c("#E1AF00","#969696","#969696","#969696","#969696","#969696","#969696","#969696","#969696","#969696","#969696","#969696","#969696","#969696","#969696")) + 
  geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("") + ylab("Odds ratio (95% CI)") +
  scale_y_continuous(trans="log10") +
  theme_bw() +
  theme(text = element_text(size=25)) 
fp

# export the plot
ggsave(fp, filename="/Users/kelseysumner/Desktop/forest_plot_aim2_model_continuous_outcome.png", device="png",
       height=10, width=12.5, units="in", dpi=400)



#### ------- make a plot of p_te_all stratified -------- ####

# make a density plot of p_te_all 
model_data$aim2_exposure = relevel(model_data$aim2_exposure, ref="asymptomatic infection")
p_te_all_plot = ggplot(data=model_data,aes(x=p_te_all_csp,fill=aim2_exposure)) +
  geom_density(alpha=0.6) + 
  scale_fill_manual(values=c("#E1AF00","#3B9AB2")) + 
  labs(fill="Symptomatic status") +
  theme_bw() + 
  xlab("Probability of transmission across all variables") +
  ylab("Density") +
  theme(plot.title = element_text(size = 26, face = "bold", hjust = 0.5), text = element_text(size=25), legend.position = c(0.7, 0.77),legend.box.background = element_rect(colour = "black"), legend.text = element_text(size=40))
p_te_all_plot
ggsave(p_te_all_plot, filename="/Users/kelseysumner/Desktop/p_te_all_plot_density.png", device="png",
       height=8, width=14, units="in", dpi=500)





#### ------ make a plot of the odds ratios of p_te_all coded binary  -------- ####

# look at a summary of the outcome variable
summary(model_data$p_te_all_csp)
length(which(is.na(model_data$p_te_all_csp)))
hist(model_data$p_te_all_csp)

# make a binary variable for 0 or >0
model_data$outcome_binary_lessthan0 = ifelse(model_data$p_te_all_csp > 0,"greater than 0.00","equal to 0.00")
table(model_data$outcome_binary_lessthan0,model_data$p_te_all_csp,useNA = "always")
table(model_data$outcome_binary_lessthan0, useNA = "always")
model_data$outcome_binary_lessthan0 = factor(model_data$outcome_binary_lessthan0)
levels(model_data$outcome_binary_lessthan0)
model_data$outcome_binary_lessthan0 = relevel(model_data$outcome_binary_lessthan0,ref = "equal to 0.00")

# make a binary variable for <0.05 or >= 0.05
model_data$outcome_binary_lessthan0.05 = ifelse(model_data$p_te_all_csp < 0.05,"less than 0.05","greater than 0.05")
table(model_data$outcome_binary_lessthan0.05,model_data$p_te_all_csp,useNA = "always")
table(model_data$outcome_binary_lessthan0.05, useNA = "always")
model_data$outcome_binary_lessthan0.05 = factor(model_data$outcome_binary_lessthan0.05)
levels(model_data$outcome_binary_lessthan0.05)
model_data$outcome_binary_lessthan0.05 = relevel(model_data$outcome_binary_lessthan0.05,ref = "less than 0.05")

# make a binary variable for <0.1 or >= 0.1
model_data$outcome_binary_lessthan0.1 = ifelse(model_data$p_te_all_csp < 0.1,"less than 0.1","greater than 0.1")
table(model_data$outcome_binary_lessthan0.1,model_data$p_te_all_csp,useNA = "always")
table(model_data$outcome_binary_lessthan0.1, useNA = "always")
model_data$outcome_binary_lessthan0.1 = factor(model_data$outcome_binary_lessthan0.1)
levels(model_data$outcome_binary_lessthan0.1)
model_data$outcome_binary_lessthan0.1 = relevel(model_data$outcome_binary_lessthan0.1,ref = "less than 0.1")

# make a binary variable for <0.15 or >= 0.15
model_data$outcome_binary_lessthan0.15 = ifelse(model_data$p_te_all_csp < 0.15,"less than 0.15","greater than 0.15")
table(model_data$outcome_binary_lessthan0.15,model_data$p_te_all_csp,useNA = "always")
table(model_data$outcome_binary_lessthan0.15, useNA = "always")
model_data$outcome_binary_lessthan0.15 = factor(model_data$outcome_binary_lessthan0.15)
levels(model_data$outcome_binary_lessthan0.15)
model_data$outcome_binary_lessthan0.15 = relevel(model_data$outcome_binary_lessthan0.15,ref = "less than 0.15")

# make a binary variable for <0.2 or >= 0.2
model_data$outcome_binary_lessthan0.2 = ifelse(model_data$p_te_all_csp < 0.2,"less than 0.2","greater than 0.2")
table(model_data$outcome_binary_lessthan0.2,model_data$p_te_all_csp,useNA = "always")
table(model_data$outcome_binary_lessthan0.2, useNA = "always")
model_data$outcome_binary_lessthan0.2 = factor(model_data$outcome_binary_lessthan0.2)
levels(model_data$outcome_binary_lessthan0.2)
model_data$outcome_binary_lessthan0.2 = relevel(model_data$outcome_binary_lessthan0.2,ref = "less than 0.2")

# make a binary variable for <0.25 or >= 0.25
model_data$outcome_binary_lessthan0.25 = ifelse(model_data$p_te_all_csp < 0.25,"less than 0.25","greater than 0.25")
table(model_data$outcome_binary_lessthan0.25,model_data$p_te_all_csp,useNA = "always")
table(model_data$outcome_binary_lessthan0.25, useNA = "always")
model_data$outcome_binary_lessthan0.25 = factor(model_data$outcome_binary_lessthan0.25)
levels(model_data$outcome_binary_lessthan0.25)
model_data$outcome_binary_lessthan0.25 = relevel(model_data$outcome_binary_lessthan0.25,ref = "less than 0.25")

# make a binary variable for <0.3 or >= 0.3
model_data$outcome_binary_lessthan0.3 = ifelse(model_data$p_te_all_csp < 0.3,"less than 0.3","greater than 0.3")
table(model_data$outcome_binary_lessthan0.3,model_data$p_te_all_csp,useNA = "always")
table(model_data$outcome_binary_lessthan0.3, useNA = "always")
model_data$outcome_binary_lessthan0.3 = factor(model_data$outcome_binary_lessthan0.3)
levels(model_data$outcome_binary_lessthan0.3)
model_data$outcome_binary_lessthan0.3 = relevel(model_data$outcome_binary_lessthan0.3,ref = "less than 0.3")

# make a binary variable for <0.35 or >= 0.35
model_data$outcome_binary_lessthan0.35 = ifelse(model_data$p_te_all_csp < 0.35,"less than 0.35","greater than 0.35")
table(model_data$outcome_binary_lessthan0.35,model_data$p_te_all_csp,useNA = "always")
table(model_data$outcome_binary_lessthan0.35, useNA = "always")
model_data$outcome_binary_lessthan0.35 = factor(model_data$outcome_binary_lessthan0.35)
levels(model_data$outcome_binary_lessthan0.35)
model_data$outcome_binary_lessthan0.35 = relevel(model_data$outcome_binary_lessthan0.35,ref = "less than 0.35")

# make a binary variable for <0.4 or >= 0.4
model_data$outcome_binary_lessthan0.4 = ifelse(model_data$p_te_all_csp < 0.4,"less than 0.4","greater than 0.4")
table(model_data$outcome_binary_lessthan0.4,model_data$p_te_all_csp,useNA = "always")
table(model_data$outcome_binary_lessthan0.4, useNA = "always")
model_data$outcome_binary_lessthan0.4 = factor(model_data$outcome_binary_lessthan0.4)
levels(model_data$outcome_binary_lessthan0.4)
model_data$outcome_binary_lessthan0.4 = relevel(model_data$outcome_binary_lessthan0.4,ref = "less than 0.4")

# make a binary variable for <0.45 or >= 0.45
model_data$outcome_binary_lessthan0.45 = ifelse(model_data$p_te_all_csp < 0.45,"less than 0.45","greater than 0.45")
table(model_data$outcome_binary_lessthan0.45,model_data$p_te_all_csp,useNA = "always")
table(model_data$outcome_binary_lessthan0.45, useNA = "always")
model_data$outcome_binary_lessthan0.45 = factor(model_data$outcome_binary_lessthan0.45)
levels(model_data$outcome_binary_lessthan0.45)
model_data$outcome_binary_lessthan0.45 = relevel(model_data$outcome_binary_lessthan0.45,ref = "less than 0.45")

# make a binary variable for <0.5 or >= 0.5
model_data$outcome_binary_lessthan0.5 = ifelse(model_data$p_te_all_csp < 0.5,"less than 0.5","greater than 0.5")
table(model_data$outcome_binary_lessthan0.5,model_data$p_te_all_csp,useNA = "always")
table(model_data$outcome_binary_lessthan0.5, useNA = "always")
model_data$outcome_binary_lessthan0.5 = factor(model_data$outcome_binary_lessthan0.5)
levels(model_data$outcome_binary_lessthan0.5)
model_data$outcome_binary_lessthan0.5 = relevel(model_data$outcome_binary_lessthan0.5,ref = "less than 0.5")

# make a binary variable for <0.55 or >= 0.55
model_data$outcome_binary_lessthan0.55 = ifelse(model_data$p_te_all_csp < 0.55,"less than 0.55","greater than 0.55")
table(model_data$outcome_binary_lessthan0.55,model_data$p_te_all_csp,useNA = "always")
table(model_data$outcome_binary_lessthan0.55, useNA = "always")
model_data$outcome_binary_lessthan0.55 = factor(model_data$outcome_binary_lessthan0.55)
levels(model_data$outcome_binary_lessthan0.55)
model_data$outcome_binary_lessthan0.55 = relevel(model_data$outcome_binary_lessthan0.55,ref = "less than 0.55")

# make a binary variable for <0.6 or >= 0.6
model_data$outcome_binary_lessthan0.6 = ifelse(model_data$p_te_all_csp < 0.6,"less than 0.6","greater than 0.6")
table(model_data$outcome_binary_lessthan0.6,model_data$p_te_all_csp,useNA = "always")
table(model_data$outcome_binary_lessthan0.6, useNA = "always")
model_data$outcome_binary_lessthan0.6 = factor(model_data$outcome_binary_lessthan0.6)
levels(model_data$outcome_binary_lessthan0.6)
model_data$outcome_binary_lessthan0.6 = relevel(model_data$outcome_binary_lessthan0.6,ref = "less than 0.6")

# make a binary variable for <0.65 or >= 0.65
model_data$outcome_binary_lessthan0.65 = ifelse(model_data$p_te_all_csp < 0.65,"less than 0.65","greater than 0.65")
table(model_data$outcome_binary_lessthan0.65,model_data$p_te_all_csp,useNA = "always")
table(model_data$outcome_binary_lessthan0.65, useNA = "always")
model_data$outcome_binary_lessthan0.65 = factor(model_data$outcome_binary_lessthan0.65)
levels(model_data$outcome_binary_lessthan0.65)
model_data$outcome_binary_lessthan0.65 = relevel(model_data$outcome_binary_lessthan0.65,ref = "less than 0.65")

# make a binary variable for <0.7 or >= 0.7
model_data$outcome_binary_lessthan0.7 = ifelse(model_data$p_te_all_csp < 0.7,"less than 0.7","greater than 0.7")
table(model_data$outcome_binary_lessthan0.7,model_data$p_te_all_csp,useNA = "always")
table(model_data$outcome_binary_lessthan0.7, useNA = "always")
model_data$outcome_binary_lessthan0.7 = factor(model_data$outcome_binary_lessthan0.7)
levels(model_data$outcome_binary_lessthan0.7)
model_data$outcome_binary_lessthan0.7 = relevel(model_data$outcome_binary_lessthan0.7,ref = "less than 0.7")

# make a binary variable for <0.75 or >= 0.75
model_data$outcome_binary_lessthan0.75 = ifelse(model_data$p_te_all_csp < 0.75,"less than 0.75","greater than 0.75")
table(model_data$outcome_binary_lessthan0.75,model_data$p_te_all_csp,useNA = "always")
table(model_data$outcome_binary_lessthan0.75, useNA = "always")
model_data$outcome_binary_lessthan0.75 = factor(model_data$outcome_binary_lessthan0.75)
levels(model_data$outcome_binary_lessthan0.75)
model_data$outcome_binary_lessthan0.75 = relevel(model_data$outcome_binary_lessthan0.75,ref = "less than 0.75")

# make a binary variable for <0.8 or >= 0.8
model_data$outcome_binary_lessthan0.8 = ifelse(model_data$p_te_all_csp < 0.8,"less than 0.8","greater than 0.8")
table(model_data$outcome_binary_lessthan0.8,model_data$p_te_all_csp,useNA = "always")
table(model_data$outcome_binary_lessthan0.8, useNA = "always")
model_data$outcome_binary_lessthan0.8 = factor(model_data$outcome_binary_lessthan0.8)
levels(model_data$outcome_binary_lessthan0.8)
model_data$outcome_binary_lessthan0.8 = relevel(model_data$outcome_binary_lessthan0.8,ref = "less than 0.8")

# make a binary variable for <0.85 or >= 0.85
model_data$outcome_binary_lessthan0.85 = ifelse(model_data$p_te_all_csp < 0.85,"less than 0.85","greater than 0.85")
table(model_data$outcome_binary_lessthan0.85,model_data$p_te_all_csp,useNA = "always")
table(model_data$outcome_binary_lessthan0.85, useNA = "always")
model_data$outcome_binary_lessthan0.85 = factor(model_data$outcome_binary_lessthan0.85)
levels(model_data$outcome_binary_lessthan0.85)
model_data$outcome_binary_lessthan0.85 = relevel(model_data$outcome_binary_lessthan0.85,ref = "less than 0.85")

# make a binary variable for <0.9 or >= 0.9
model_data$outcome_binary_lessthan0.9 = ifelse(model_data$p_te_all_csp < 0.9,"less than 0.9","greater than 0.9")
table(model_data$outcome_binary_lessthan0.9,model_data$p_te_all_csp,useNA = "always")
table(model_data$outcome_binary_lessthan0.9, useNA = "always")
model_data$outcome_binary_lessthan0.9 = factor(model_data$outcome_binary_lessthan0.9)
levels(model_data$outcome_binary_lessthan0.9)
model_data$outcome_binary_lessthan0.9 = relevel(model_data$outcome_binary_lessthan0.9,ref = "less than 0.9")

# make a binary variable for <0.95 or >= 0.95
model_data$outcome_binary_lessthan0.95 = ifelse(model_data$p_te_all_csp < 0.95,"less than 0.95","greater than 0.95")
table(model_data$outcome_binary_lessthan0.95,model_data$p_te_all_csp,useNA = "always")
table(model_data$outcome_binary_lessthan0.95, useNA = "always")
model_data$outcome_binary_lessthan0.95 = factor(model_data$outcome_binary_lessthan0.95)
levels(model_data$outcome_binary_lessthan0.95)
model_data$outcome_binary_lessthan0.95 = relevel(model_data$outcome_binary_lessthan0.95,ref = "less than 0.95")

# binary outcome 0 with a logistic model - this is basically a hurdle model
model0 <- glmmTMB(outcome_binary_lessthan0~aim2_exposure+pfr364Q_std_combined_rescaled+age_cat_baseline+mosquito_week_count_cat +village_name+(1|HH_ID_human/unq_memID),family=binomial(link = "logit"), data = model_data)
summary(model0)
exp(confint(model0,method="Wald"))

# binary outcome <0.05 with a logistic model
model.05 <- glmmTMB(outcome_binary_lessthan0.05~aim2_exposure+pfr364Q_std_combined_rescaled+age_cat_baseline+mosquito_week_count_cat +village_name+(1|HH_ID_human/unq_memID),family=binomial(link = "logit"), data = model_data)
summary(model.05)
exp(confint(model.05,method="Wald"))
# converged

# binary outcome <0.1 with a logistic model
model.1 <- glmmTMB(outcome_binary_lessthan0.1~aim2_exposure+pfr364Q_std_combined_rescaled+age_cat_baseline+mosquito_week_count_cat +village_name+(1|HH_ID_human/unq_memID),family=binomial(link = "logit"), data = model_data)
summary(model.1)
exp(confint(model.1,method="Wald"))
# converged

# binary outcome <0.15 with a logistic model
model.15 <- glmmTMB(outcome_binary_lessthan0.15~aim2_exposure+pfr364Q_std_combined_rescaled+age_cat_baseline+mosquito_week_count_cat +village_name+(1|HH_ID_human/unq_memID),family=binomial(link = "logit"), data = model_data)
summary(model.15)
exp(confint(model.15,method="Wald"))
# converged

# binary outcome <0.2 with a logistic model
model.2 <- glmmTMB(outcome_binary_lessthan0.2~aim2_exposure+pfr364Q_std_combined_rescaled+age_cat_baseline+mosquito_week_count_cat+village_name+(1|HH_ID_human/unq_memID),family=binomial(link = "logit"), data = model_data)
summary(model.2)
exp(confint(model.2, method="Wald"))
# converged

# binary outcome <0.25 with a logistic model
model.25 <- glmmTMB(outcome_binary_lessthan0.25~aim2_exposure+pfr364Q_std_combined_rescaled+age_cat_baseline+mosquito_week_count_cat+village_name+(1|HH_ID_human/unq_memID),family=binomial(link = "logit"), data = model_data)
summary(model.25)
exp(confint(model.25, method="Wald"))
# converged

# binary outcome <0.3 with a logistic model
model.3 <-  glmmTMB(outcome_binary_lessthan0.3~aim2_exposure+pfr364Q_std_combined_rescaled+age_cat_baseline+mosquito_week_count_cat+village_name+(1|HH_ID_human/unq_memID),family=binomial(link = "logit"), data = model_data)
summary(model.3)
exp(confint(model.3, method="Wald"))
# converged

# binary outcome <0.35 with a logistic model
model.35 <- glmmTMB(outcome_binary_lessthan0.35~aim2_exposure+pfr364Q_std_combined_rescaled+age_cat_baseline+mosquito_week_count_cat+village_name+(1|HH_ID_human/unq_memID),family=binomial(link = "logit"), data = model_data)
summary(model.35)
exp(confint(model.35, method="Wald"))
# converged

# binary outcome <0.4 with a logistic model
model.4 <- glmmTMB(outcome_binary_lessthan0.4~aim2_exposure+pfr364Q_std_combined_rescaled+age_cat_baseline+mosquito_week_count_cat+village_name+(1|HH_ID_human/unq_memID),family=binomial(link = "logit"), data = model_data)
summary(model.4)
exp(confint(model.4, method="Wald")) 
# converged

# binary outcome <0.45 with a logistic model
model.45 <- glmmTMB(outcome_binary_lessthan0.45~aim2_exposure+pfr364Q_std_combined_rescaled+age_cat_baseline+mosquito_week_count_cat+village_name+(1|HH_ID_human/unq_memID),family=binomial(link = "logit"), data = model_data)
summary(model.45)
exp(confint(model.45, method="Wald")) 
# converged

# binary outcome <0.5 with a logistic model
model.5 <- glmmTMB(outcome_binary_lessthan0.5~aim2_exposure+pfr364Q_std_combined_rescaled+age_cat_baseline+mosquito_week_count_cat +village_name+(1|HH_ID_human/unq_memID),family=binomial(link = "logit"), data = model_data)
summary(model.5)
exp(confint(model.5, method="Wald")) 
# converged

# binary outcome <0.55 with a logistic model
model.55 <- glmmTMB(outcome_binary_lessthan0.55~aim2_exposure+pfr364Q_std_combined_rescaled+age_cat_baseline+mosquito_week_count_cat +village_name+(1|HH_ID_human/unq_memID),family=binomial(link = "logit"), data = model_data)
summary(model.55)
exp(confint(model.55, method="Wald")) 
# converged

# binary outcome <0.6 with a logistic model
model.6 <- glmmTMB(outcome_binary_lessthan0.6~aim2_exposure+pfr364Q_std_combined_rescaled+age_cat_baseline+mosquito_week_count_cat +village_name+(1|HH_ID_human/unq_memID),family=binomial(link = "logit"), data = model_data)
summary(model.6)
exp(confint(model.6, method="Wald")) # can't compute confidence interval
# converged but producing weird results

# binary outcome <0.65 with a logistic model
model.65 <- glmmTMB(outcome_binary_lessthan0.65~aim2_exposure+pfr364Q_std_combined_rescaled+age_cat_baseline+mosquito_week_count_cat +village_name+(1|HH_ID_human/unq_memID),family=binomial(link = "logit"), data = model_data)
summary(model.65)
exp(confint(model.65, method="Wald")) # can't compute confidence interval
# converged but producing weird results

# binary outcome <0.7 with a logistic model
model.7 <- glmmTMB(outcome_binary_lessthan0.7~aim2_exposure+pfr364Q_std_combined_rescaled+age_cat_baseline+mosquito_week_count_cat +village_name+(1|HH_ID_human/unq_memID),family=binomial(link = "logit"), data = model_data)
summary(model.7)
exp(confint(model.7, method="Wald")) # can't compute confidence interval
# converged but weird results

# read in the model results
model_results = read_csv("Desktop/Dissertation Materials/SpatialR21 Grant/Final Dissertation Materials/Aim 2/computational_model_materials/aim2_binary_outcome_final.csv")
model_results = model_results[-c(13),]

# make the plot
model_plot = ggplot(data=model_results,aes(x=binary_outcome,y=estimate,group=1),cex=1.5,col="#E1AF00") +
  geom_line(data=model_results,aes(x=binary_outcome,y=estimate,group=1),cex=1.5,col="#E1AF00") +
  geom_ribbon(data=model_results,aes(x=1:length(binary_outcome),ymin = lower_ci, ymax = upper_ci),alpha=0.2,fill="#E1AF00") +
  theme_bw() +
  xlab("Cutoff for what is a transmission event") + 
  ylab("Odds ratio (95% CI)") + 
  scale_y_continuous(breaks=c(0,1,2,3),trans="log10") +
  geom_hline(yintercept=1,linetype="dashed") + 
  coord_flip() +
  theme(text = element_text(size=16)) 
model_plot
ggsave(model_plot, filename="/Users/kelseysumner/Desktop/binary_coding_model_plot.png", device="png",
       height=7, width=5, units="in", dpi=500)



