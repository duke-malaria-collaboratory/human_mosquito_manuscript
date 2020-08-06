# -------------------------------------- #
#           Spat21/Mozzie Study          #
#  Use Amy's asymp/symp matching method  #
#                 Aim 2                  #
#            Mozzie Phase 1              #
#               K. Sumner                #
#            April 7, 2020               #
# -------------------------------------- #


#### -------- load packages ------------ ####

# load in the packages of interest
library(tidyverse)
library(BSDA)
library(glmmTMB)



#### ----- read in the data sets ----- ####

# read in the combined ama and csp data set for mosquito abdomens
model_data = read_rds("Desktop/Dissertation Materials/SpatialR21 Grant/Final Dissertation Materials/Aim 2/clean_ids_haplotype_results/AMA_and_CSP/final/model data/final_model_data/spat21_aim2_merged_data_with_weights_5MAR2020.rds")

# subset the data set to samples that passed pfcsp sequencing only
csp_data = model_data %>%
  filter(!(is.na(csp_haps_shared)))

# subset the data set to samples that passed pfcsp sequencing only
ama_data = model_data %>%
  filter(!(is.na(ama_haps_shared)))

# read in the full human demographic data set
final_data = read_rds("Desktop/Dissertation Materials/SpatialR21 Grant/Final Dissertation Materials/Final Data Sets/Final Cohort data June 2017 to July 2018/Human data/spat21_clean_human_files/merged_files/final merged data/final_recoded_data_set/spat21_human_final_censored_data_for_dissertation_with_exposure_outcome_1MAR2020.rds")
final_data = final_data %>%
  select(c(unq_memID,age_cat_baseline))



#### ------ look at pairs of infections -------- ####

# first subset to only the individuals who have at least 1 symptomatic and 1 asymptomatic infection
# for csp - end up with 65 unique participants and 1565 obs
csp_infxns_to_include_p1 = csp_data %>%
  group_by(unq_memID,aim2_exposure) %>%
  summarize(n=n())
csp_infxns_to_include_p2 = csp_infxns_to_include_p1 %>%
  group_by(unq_memID) %>%
  summarize(n=n())
csp_infxns_to_include_p3 = csp_infxns_to_include_p2 %>%
  filter(n>1)
csp_subset_data = csp_data %>%
  filter(csp_data$unq_memID %in% csp_infxns_to_include_p3$unq_memID)
# for ama - end up with 56 unique participants and 1197 obs
ama_infxns_to_include_p1 = ama_data %>%
  group_by(unq_memID,aim2_exposure) %>%
  summarize(n=n())
ama_infxns_to_include_p2 = ama_infxns_to_include_p1 %>%
  group_by(unq_memID) %>%
  summarize(n=n())
ama_infxns_to_include_p3 = ama_infxns_to_include_p2 %>%
  filter(n>1)
ama_subset_data = ama_data %>%
  filter(ama_data$unq_memID %in% ama_infxns_to_include_p3$unq_memID)

# now calculate the average p_te_all within asymptomatic and symptomatic infections
# for csp
csp_average_data = csp_subset_data %>%
  group_by(village_name,unq_memID,aim2_exposure) %>%
  summarize(average_p_te_all_csp=mean(p_te_all_csp,na.rm=T), median_p_te_all_csp = median(p_te_all_csp,na.rm = T))
# for ama
ama_average_data = ama_subset_data %>%
  group_by(village_name,unq_memID,aim2_exposure) %>%
  summarize(average_p_te_all_ama=mean(p_te_all_ama,na.rm=T), median_p_te_all_ama = median(p_te_all_ama,na.rm = T))  

# now split up the data sets so have separate columns for asymptomatic and symptomatic infections
# for csp
csp_asymp = csp_average_data %>%
  filter(aim2_exposure=="asymptomatic infection") %>%
  select(-c(aim2_exposure)) %>%
  rename("asymp_average_p_te_all_csp" = average_p_te_all_csp,"asymp_median_p_te_all_csp" = median_p_te_all_csp)
csp_symp = csp_average_data %>%
  filter(aim2_exposure == "symptomatic infection") %>%
  select(-c(village_name,aim2_exposure)) %>%
  rename("symp_average_p_te_all_csp" = average_p_te_all_csp,"symp_median_p_te_all_csp" = median_p_te_all_csp)
csp_all = left_join(csp_asymp,csp_symp,by=c("unq_memID"))
csp_all$village_name.y <- NULL
csp_all = rename(csp_all,village_name = village_name.x)
# ama
ama_asymp = ama_average_data %>%
  filter(aim2_exposure=="asymptomatic infection") %>%
  select(-c(aim2_exposure)) %>%
  rename("asymp_average_p_te_all_ama" = average_p_te_all_ama,"asymp_median_p_te_all_ama" = median_p_te_all_ama)
ama_symp = ama_average_data %>%
  filter(aim2_exposure == "symptomatic infection") %>%
  select(-c(village_name,aim2_exposure)) %>%
  rename("symp_average_p_te_all_ama" = average_p_te_all_ama,"asymp_median_p_te_all_ama" = median_p_te_all_ama)
ama_all = left_join(ama_asymp,ama_symp,by=c("unq_memID"))
ama_all$village_name.y <- NULL
ama_all = rename(ama_all,village_name = village_name.x)

# check the mean of means by symptomatic status
# for csp
mean(csp_all$asymp_average_p_te_all_csp)
mean(csp_all$symp_average_p_te_all_csp)
# for ama
mean(ama_all$asymp_average_p_te_all_ama)
mean(ama_all$symp_average_p_te_all_ama)

# check the median of means by symptomatic status
# for csp
median(csp_all$asymp_average_p_te_all_csp)
median(csp_all$symp_average_p_te_all_csp)
# for ama
median(ama_all$asymp_average_p_te_all_ama)
median(ama_all$symp_average_p_te_all_ama)

# check normality (but N>30 so not a large issue)
# for csp
d <- csp_all$asymp_average_p_te_all_csp-csp_all$symp_average_p_te_all_csp
shapiro.test(d) # normality not an issue
hist(d)
# for ama
d <- ama_all$asymp_average_p_te_all_ama-ama_all$symp_average_p_te_all_ama
shapiro.test(d) # normality could be an issue but N>30 so central limit theorem applies
hist(d)

# paired t-test
t.test(csp_all$asymp_average_p_te_all_csp, csp_all$symp_average_p_te_all_csp, paired = TRUE, alternative = "greater")
t.test(ama_all$asymp_average_p_te_all_ama, ama_all$symp_average_p_te_all_ama, paired = TRUE, alternative = "greater")

# sign test just to be careful
SIGN.test(csp_all$asymp_average_p_te_all_csp, csp_all$symp_average_p_te_all_csp,md=0,alternative = "greater")
SIGN.test(ama_all$asymp_average_p_te_all_ama, ama_all$symp_average_p_te_all_ama,md=0,alternative = "greater")


#### -------- for csp: compute proportion nonzero pairings with a mosquito ------ ####

# what this does: for each person's infection, compute the proportion of nonzero pairings with mosquito
# then take the median of that proportion across symptomatic status

# only for csp right now

# create a variable that says whether or not p_te_all_csp is non-zero or not
csp_subset_data$nonzero = ifelse(csp_subset_data$p_te_all_csp > 0,"nonzero","zero")
table(csp_subset_data$nonzero,useNA = "always")
length(which(csp_subset_data$p_te_all_csp == 0))

# for each person's infection, count the number of nonzero pairings
csp_part1_data = csp_subset_data %>%
  filter(nonzero == "nonzero") %>%
  group_by(village_name,unq_memID,aim2_exposure,sample_id_human) %>%
  summarize(numerator=n())
csp_nonzero_data = csp_subset_data %>%
  group_by(village_name,unq_memID,aim2_exposure,sample_id_human) %>%
  summarize(denominator = n())

# join the data frames
csp_all_nonzero_data = left_join(csp_nonzero_data,csp_part1_data,by=c("sample_id_human"))

# clean up the data frames
csp_all_nonzero_data$village_name.y <- NULL
csp_all_nonzero_data$aim2_exposure.y <- NULL
csp_all_nonzero_data$unq_memID.y <- NULL
csp_all_nonzero_data = csp_all_nonzero_data %>% rename(village_name = village_name.x,unq_memID=unq_memID.x,aim2_exposure=aim2_exposure.x)

# change the numerators that are NA to 0
csp_all_nonzero_data$numerator[which(is.na(csp_all_nonzero_data$numerator))] = 0

# calculate the proportion of infections that are nonzero
csp_all_nonzero_data$prop_nonzero = csp_all_nonzero_data$numerator/csp_all_nonzero_data$denominator

# create a variable for household ID
csp_all_nonzero_data$HH_ID_human = rep(NA,nrow(csp_all_nonzero_data))
for (i in 1:nrow(csp_all_nonzero_data)) {
  csp_all_nonzero_data$HH_ID_human[i] = str_split(csp_all_nonzero_data$unq_memID[i],"_")[[1]][1]
}
table(csp_all_nonzero_data$HH_ID_human, useNA = "always")

# add covarites to the data set for the regression
colnames(csp_data)
csp_to_merge_data = csp_data %>% select(sample_id_human,pfr364Q_std_combined_rescaled,age_cat_baseline,mosquito_week_count_cat)
csp_all_nonzero_data = left_join(csp_all_nonzero_data,csp_to_merge_data,by="sample_id_human")
csp_all_nonzero_data = unique(csp_all_nonzero_data)

# check the covariates
str(csp_all_nonzero_data$sample_id_human)
str(csp_all_nonzero_data$HH_ID_human)
str(csp_all_nonzero_data$unq_memID)
str(csp_all_nonzero_data$age_cat_baseline)
str(csp_all_nonzero_data$village_name)
csp_all_nonzero_data$village_name = relevel(csp_all_nonzero_data$village_name,ref = "Maruti")
str(csp_all_nonzero_data$mosquito_week_count_cat)
str(csp_all_nonzero_data$pfr364Q_std_combined_rescaled)
str(csp_all_nonzero_data$aim2_exposure)

# do a multi-level logistic regression across symptomatic status
csp_model = glmmTMB(prop_nonzero ~ aim2_exposure + pfr364Q_std_combined_rescaled+mosquito_week_count_cat+(1|HH_ID_human/unq_memID),family=binomial(link = "logit"), data = csp_all_nonzero_data)
summary(csp_model)
exp(confint(csp_model,method="Wald"))
# OR: 2.56 (95% CI: 1.36 to 4.81)

table1 = exp(confint(csp_model,method="Wald"))
estimates = c(table1[2,3],table1[3,3],table1[4,3])
lower_ci = c(table1[2,1],table1[3,1],table1[4,1])
upper_ci = c(table1[2,2],table1[3,2],table1[4,2])
names = c("Asymptomatic infection","Parasite density (parasite/uL whole blood)","High mosquito abundance")
forest_plot_df = data.frame(names,estimates,lower_ci,upper_ci)
forest_plot_df$names = factor(forest_plot_df$names, levels = c("Asymptomatic infection","Parasite density (parasite/uL whole blood)","High mosquito abundance"))
forest_plot_df$names = ordered(forest_plot_df$names, levels = c("Asymptomatic infection","Parasite density (parasite/uL whole blood)","High mosquito abundance"))

# create a forest plot
library(forcats)
fp <- ggplot(data=forest_plot_df, aes(x=fct_rev(names), y=estimates, ymin=lower_ci, ymax=upper_ci)) +
  geom_pointrange(size=c(3,1,1),colour=c("#006d2c","#969696","#969696")) + 
  geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("") + ylab("Odds ratio (95% CI)") +
  scale_y_continuous(trans="log10") +
  theme_bw() +
  theme(text = element_text(size=14))
fp

# export the plot
ggsave(fp, filename="/Users/kelseysumner/Desktop/forest_plot_aim2_model_continuous_outcome_figure3coding.png", device="png",
       height=3, width=7, units="in", dpi=400)


# now take the median of that proportion by symptomatic status
csp_final_data = csp_all_nonzero_data %>%
  group_by(village_name,unq_memID,aim2_exposure) %>%
  summarize(median_prop_nonzero = median(prop_nonzero,na.rm = T),num_infections = n())

# now split up the data sets so have separate columns for asymptomatic and symptomatic infections
# for csp
csp_asymp = csp_final_data %>%
  filter(aim2_exposure=="asymptomatic infection") %>%
  select(-c(aim2_exposure)) %>%
  rename("asymp_median_prop_nonzero" = median_prop_nonzero)
csp_symp = csp_final_data %>%
  filter(aim2_exposure=="symptomatic infection") %>%
  select(-c(aim2_exposure)) %>%
  rename("symp_median_prop_nonzero" = median_prop_nonzero)
csp_all = left_join(csp_asymp,csp_symp,by=c("unq_memID"))
csp_all$village_name.y <- NULL
csp_all = rename(csp_all,village_name = village_name.x)
csp_all = data.frame(csp_all)
csp_all$total_infections = csp_all$num_infections.x + csp_all$num_infections.y
csp_all$num_infections.x <- NULL
csp_all$num_infections.y <- NULL

# now add a variable for the age categories
final_data = final_data %>% group_by(unq_memID,age_cat_baseline) 
final_data = data.frame(final_data)
csp_all = left_join(csp_all,final_data,by="unq_memID")
csp_all = unique(csp_all)

# now make a plot of the average probability of p_te_all by symptomatic status
csp_plot = ggplot(csp_all, aes(x=symp_median_prop_nonzero, y=asymp_median_prop_nonzero)) +
  geom_point(aes(size=total_infections),pch=21,alpha=0.75,color="black",fill="#8AAF9D") + 
  theme_bw() +
  geom_abline(intercept = 0, slope = 1,linetype="dashed") +
  xlab("Median proportion infected mosquitoes with parasites matching a symptomatic infection") +
  ylab("Median proportion infected mosquitoes with \n parasites matching an asymptomatic infection") +
  labs(size="Number of malaria infections") +
  scale_x_continuous(limits=c(0,1)) +
  scale_y_continuous(limits=c(0,1)) +
  theme(text = element_text(size=14), legend.position = c(0.85,0.25), legend.box.background = element_rect(colour = "black")) 
ggsave(csp_plot, filename="/Users/kelseysumner/Desktop/csp_matched_prob_plot_alt3.png", device="png",
       height=5, width=11, units="in", dpi=500)



#### -------- now do for ama: compute proportion nonzero pairings with a mosquito ------ ####

# what this does: for each person's infection, compute the proportion of nonzero pairings with mosquito
# then take the median of that proportion across symptomatic status

# only for ama right now

# create a variable that says whether or not p_te_all_csp is non-zero or not
ama_subset_data$nonzero = ifelse(ama_subset_data$p_te_all_ama > 0,"nonzero","zero")
table(ama_subset_data$nonzero,useNA = "always")
length(which(ama_subset_data$p_te_all_ama == 0))

# for each person's infection, count the number of nonzero pairings
ama_part1_data = ama_subset_data %>%
  filter(nonzero == "nonzero") %>%
  group_by(village_name,unq_memID,aim2_exposure,sample_id_human) %>%
  summarize(numerator=n())
ama_nonzero_data = ama_subset_data %>%
  group_by(village_name,unq_memID,aim2_exposure,sample_id_human) %>%
  summarize(denominator = n())

# join the data frames
ama_all_nonzero_data = left_join(ama_nonzero_data,ama_part1_data,by=c("sample_id_human"))

# clean up the data frames
ama_all_nonzero_data$village_name.y <- NULL
ama_all_nonzero_data$aim2_exposure.y <- NULL
ama_all_nonzero_data$unq_memID.y <- NULL
ama_all_nonzero_data = ama_all_nonzero_data %>% rename(village_name = village_name.x,unq_memID=unq_memID.x,aim2_exposure=aim2_exposure.x)

# change the numerators that are NA to 0
ama_all_nonzero_data$numerator[which(is.na(ama_all_nonzero_data$numerator))] = 0

# calculate the proportion of infections that are nonzero
ama_all_nonzero_data$prop_nonzero= ama_all_nonzero_data$numerator/ama_all_nonzero_data$denominator

# now take the median of that proportion by symptomatic status
ama_final_data = ama_all_nonzero_data %>%
  group_by(village_name,unq_memID,aim2_exposure) %>%
  summarize(median_prop_nonzero = median(prop_nonzero,na.rm = T),num_infections = n())

# now split up the data sets so have separate columns for asymptomatic and symptomatic infections
# for ama
ama_asymp = ama_final_data %>%
  filter(aim2_exposure=="asymptomatic infection") %>%
  select(-c(aim2_exposure)) %>%
  rename("asymp_median_prop_nonzero" = median_prop_nonzero)
ama_symp = ama_final_data %>%
  filter(aim2_exposure=="symptomatic infection") %>%
  select(-c(aim2_exposure)) %>%
  rename("symp_median_prop_nonzero" = median_prop_nonzero)
ama_all = left_join(ama_asymp,ama_symp,by=c("unq_memID"))
ama_all$village_name.y <- NULL
ama_all = rename(ama_all,village_name = village_name.x)
ama_all = data.frame(ama_all)
ama_all$total_infections = ama_all$num_infections.x + ama_all$num_infections.y
ama_all$num_infections.x <- NULL
ama_all$num_infections.y <- NULL

# now add a variable for the age categories
ama_all = left_join(ama_all,final_data,by="unq_memID")
ama_all = unique(ama_all)

# now make a plot of the average probability of p_te_all by symptomatic status
ama_plot = ggplot(ama_all, aes(x=symp_median_prop_nonzero, y=asymp_median_prop_nonzero)) +
  geom_point(aes(size=total_infections),pch=21,alpha=0.75,color="black",fill="#8AAF9D") + 
  theme_bw() +
  geom_abline(intercept = 0, slope = 1,linetype="dashed") +
  xlab("Median proportion infected mosquitoes for symptomatic infections") +
  ylab("Median proportion infected mosquitoes for asymptomatic infections") +
  labs(size="Number of malaria infections") +
  scale_x_continuous(limits=c(0,1)) +
  scale_y_continuous(limits=c(0,1)) +
  theme(text = element_text(size=12)) 
ggsave(ama_plot, filename="/Users/kelseysumner/Desktop/ama_matched_prob_plot_alt3.png", device="png",
       height=5.5, width=8, units="in", dpi=500)

# check the median of medians by symptomatic status
# for ama
median(ama_all$asymp_median_prop_nonzero)
median(ama_all$symp_median_prop_nonzero)
IQR(ama_all$asymp_median_prop_nonzero)
IQR(ama_all$symp_median_prop_nonzero)

# sign test
SIGN.test(ama_all$asymp_median_prop_nonzero, ama_all$symp_median_prop_nonzero,md=0,alternative = "greater")

# add covariates to the data set for the regression
colnames(ama_data)
ama_to_merge_data = ama_data %>% select(sample_id_human,pfr364Q_std_combined_rescaled,age_cat_baseline,mosquito_week_count_cat)
ama_all_nonzero_data = left_join(ama_all_nonzero_data,ama_to_merge_data,by="sample_id_human")
ama_all_nonzero_data = unique(ama_all_nonzero_data)

# create a variable for household ID
ama_all_nonzero_data$HH_ID_human = rep(NA,nrow(ama_all_nonzero_data))
for (i in 1:nrow(ama_all_nonzero_data)) {
  ama_all_nonzero_data$HH_ID_human[i] = str_split(ama_all_nonzero_data$unq_memID[i],"_")[[1]][1]
}
table(ama_all_nonzero_data$HH_ID_human, useNA = "always")

# check the covariates
str(ama_all_nonzero_data$sample_id_human)
str(ama_all_nonzero_data$HH_ID_human)
str(ama_all_nonzero_data$unq_memID)
str(ama_all_nonzero_data$age_cat_baseline)
str(ama_all_nonzero_data$village_name)
ama_all_nonzero_data$village_name = relevel(ama_all_nonzero_data$village_name,ref = "Maruti")
str(ama_all_nonzero_data$mosquito_week_count_cat)
str(ama_all_nonzero_data$pfr364Q_std_combined_rescaled)
str(ama_all_nonzero_data$aim2_exposure)

# do a multi-level logistic regression across symptomatic status
ama_model = glmmTMB(prop_nonzero ~ aim2_exposure + pfr364Q_std_combined_rescaled+mosquito_week_count_cat+ (1|HH_ID_human/unq_memID),family=binomial(link = "logit"), data = ama_all_nonzero_data)
summary(ama_model)
exp(confint(ama_model,method="Wald"))
# OR: 1.30 (95% CI: 0.63 to 2.69)

table1 = exp(confint(ama_model,method="Wald"))
estimates = c(table1[2,3],NA,table1[3,3],NA,NA,NA,table1[4,3])
lower_ci = c(table1[2,1],NA,table1[3,1],NA,NA,NA,table1[4,1])
upper_ci = c(table1[2,2],NA,table1[3,2],NA,NA,NA,table1[4,2])
names = c("Asymptomatic infection","","Parasite density (parasite/uL whole blood)"," ","Mosquito abundance          ","Low (REF)","High")
forest_plot_df = data.frame(names,estimates,lower_ci,upper_ci)
forest_plot_df$names = factor(forest_plot_df$names, levels = c("Asymptomatic infection","","Parasite density (parasite/uL whole blood)"," ","Mosquito abundance          ","Low (REF)","High"))
forest_plot_df$names = ordered(forest_plot_df$names, levels = c("Asymptomatic infection","","Parasite density (parasite/uL whole blood)"," ","Mosquito abundance          ","Low (REF)","High"))

# create a forest plot
library(forcats)
fp <- ggplot(data=forest_plot_df, aes(x=fct_rev(names), y=estimates, ymin=lower_ci, ymax=upper_ci)) +
  geom_pointrange(size=c(3,1,1,1,1,1,1),colour=c("#006d2c","#969696","#969696","#969696","#969696","#969696","#969696")) + 
  geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("") + ylab("Odds ratio (95% CI)") +
  scale_y_continuous(trans="log10") +
  theme_bw() +
  theme(text = element_text(size=14)) 
fp

# export the plot
ggsave(fp, filename="/Users/kelseysumner/Desktop/forest_plot_aim2_model_continuous_outcome_figure3coding_ama_supplement.png", device="png",
       height=5, width=8, units="in", dpi=400)
