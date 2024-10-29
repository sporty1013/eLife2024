
library(sjstats)
library(lme4)
library(lmerTest)
library(MuMIn)
library(effsize)
library(Rmisc)
library(smplot2)

#Import data
eLife2024<-read.csv("eLife2024.csv", header = TRUE, na.strings = "") 

#Preparing dataset for Table 2 and Figure 4.
#Data Frame for the difference
id<-eLife2024$id[eLife2024$Task =="STW"]
DLPFC_diff<-ifelse(is.na(eLife2024$DLPFC[eLife2024$Task == "STUP"]) | is.na(eLife2024$DLPFC[eLife2024$Task == "STW"]), NA, eLife2024$DLPFC[eLife2024$Task == "STUP"]-eLife2024$DLPFC[eLife2024$Task == "STW"])
VLPFC_diff<-ifelse(is.na(eLife2024$VLPFC[eLife2024$Task == "STUP"]) | is.na(eLife2024$VLPFC[eLife2024$Task == "STW"]), NA, eLife2024$VLPFC[eLife2024$Task == "STUP"]-eLife2024$VLPFC[eLife2024$Task == "STW"])
DMPFC_diff<-ifelse(is.na(eLife2024$DMPFC[eLife2024$Task == "STUP"]) | is.na(eLife2024$DMPFC[eLife2024$Task == "STW"]), NA, eLife2024$DMPFC[eLife2024$Task == "STUP"]-eLife2024$DMPFC[eLife2024$Task == "STW"])
VMPFC_diff<-ifelse(is.na(eLife2024$VMPFC[eLife2024$Task == "STUP"]) | is.na(eLife2024$VMPFC[eLife2024$Task == "STW"]), NA, eLife2024$VMPFC[eLife2024$Task == "STUP"]-eLife2024$VMPFC[eLife2024$Task == "STW"])

age<-eLife2024$age[eLife2024$Task =="STW"]
sex<-eLife2024$sex[eLife2024$Task =="STW"]
bmi<-eLife2024$bmi[eLife2024$Task =="STW"]
step_number_ave<-eLife2024$step_number_ave[eLife2024$Task =="STUP"]

stw_mep_post<-eLife2024$stw_mep_post[eLife2024$Task =="STW"]
stup_mep_post<-eLife2024$stup_mep_post[eLife2024$Task =="STUP"]

phq8_total<-eLife2024$phq8_total[eLife2024$Task =="STW"]
pcs_total<-eLife2024$pcs_total[eLife2024$Task =="STW"]
paindetect_total<-eLife2024$paindetect_total[eLife2024$Task =="STW"]
koos_pain<-eLife2024$koos_pain[eLife2024$Task =="STW"]
peg_total<-eLife2024$peg_total[eLife2024$Task =="STW"]
fabq_total<-eLife2024$fabq_total[eLife2024$Task =="STW"]

diff<-data.frame(id,VLPFC_diff, VMPFC_diff, DLPFC_diff,DMPFC_diff, age, sex, bmi, stw_mep_post, stup_mep_post, 
                 phq8_total, pcs_total, paindetect_total, koos_pain, peg_total,fabq_total, step_number_ave)


###############################################################################################
#Table 2
a<-lm(DLPFC_diff~koos_pain + age+ bmi+ sex, data = diff)
summary(a)
sd <- sd(diff$koos_pain, na.rm = TRUE)
coef <- coef(a)["koos_pain"]
estimate <- coef * sd
confint <- confint(a)["koos_pain",]
confint_sd<-confint*sd
cat(round(estimate, 2),"(",round(confint_sd,2),")\n")


a<-lm(DLPFC_diff~stup_mep_post + age+ bmi+ sex, data = diff)
summary(a)
sd <- sd(diff$stup_mep_post, na.rm = TRUE)
coef <- coef(a)["stup_mep_post"]
estimate <- coef * sd
confint <- confint(a)["stup_mep_post",]
confint_sd<-confint*sd
cat(round(estimate, 2),"(",round(confint_sd,2),")\n")


a<-lm(DLPFC_diff~peg_total + age+ bmi+ sex, data = diff)
summary(a)
sd <- sd(diff$peg_total, na.rm = TRUE)
coef <- coef(a)["peg_total"]
estimate <- coef * sd
confint <- confint(a)["peg_total",]
confint_sd<-confint*sd
cat(round(estimate, 2),"(",round(confint_sd,2),")\n")


a<-lm(DLPFC_diff~phq8_total + age+ bmi+ sex, data = diff)
summary(a)
sd <- sd(diff$phq8_total, na.rm = TRUE)
coef <- coef(a)["phq8_total"]
estimate <- coef * sd
confint <- confint(a)["phq8_total",]
confint_sd<-confint*sd
cat(round(estimate, 2),"(",round(confint_sd,2),")\n")


a<-lm(DLPFC_diff~pcs_total + age+ bmi+ sex, data = diff)
summary(a)
sd <- sd(diff$pcs_total, na.rm = TRUE)
coef <- coef(a)["pcs_total"]
estimate <- coef * sd
confint <- confint(a)["pcs_total",]
confint_sd<-confint*sd
cat(round(estimate, 2),"(",round(confint_sd,2),")\n")

a<-lm(DLPFC_diff~fabq_total + age+ bmi+ sex, data = diff)
summary(a)
sd <- sd(diff$fabq_total, na.rm = TRUE)
coef <- coef(a)["fabq_total"]
estimate <- coef * sd
confint <- confint(a)["fabq_total",]
confint_sd<-confint*sd
cat(round(estimate, 2),"(",round(confint_sd,2),")\n")

################################################################################################3
#Figure 4

DLPFC_koos<-ggplot(data = diff, mapping = aes(x = koos_pain, y = DLPFC_diff)) +
  geom_point(shape = 17, color = "#0f993d", size = 0.5) +
  sm_statCorr(color = "darkgrey", corr_method = "spearman", linetype = "solid", text_size = 2, size = 0.5)+
  theme_gray(base_size = 8)+
  theme(panel.grid.minor = element_blank())

DLPFC_mep<-ggplot(data = diff, mapping = aes(x = stup_mep_post, y = DLPFC_diff)) +
  geom_point(shape = 17, color = "#0f993d", size = 0.5) +
  sm_statCorr(color = "darkgrey", corr_method = "spearman",linetype = "solid", text_size = 2, size = 0.5)+
  theme_gray(base_size = 8)+
  theme(panel.grid.minor = element_blank())

DLPFC_peg<-ggplot(data = diff, mapping = aes(x = peg_total, y = DLPFC_diff)) +
  geom_point(shape = 17, color = "#0f993d", size = 0.5) +
  sm_statCorr(color = "darkgrey", corr_method = "spearman",linetype = "solid", text_size = 2, size = 0.5)+
  theme_gray(base_size = 8)+
  theme(panel.grid.minor = element_blank())

DLPFC_phq<-ggplot(data = diff, mapping = aes(x = phq8_total, y = DLPFC_diff)) +
  geom_point(shape = 17,color = "#0f993d",size = 0.5) +
  sm_statCorr(color = "darkgrey", corr_method = "spearman",linetype = "solid", text_size = 2, size = 0.5)+
  theme_gray(base_size = 8)+
  theme(panel.grid.minor = element_blank())

DLPFC_pcs<-ggplot(data = diff, mapping = aes(x = pcs_total, y = DLPFC_diff)) +
  geom_point(shape = 17,color = "#0f993d",size =0.5) +
  sm_statCorr(color = "darkgrey", corr_method = "spearman",linetype = "solid", text_size = 2, size = 0.5)+
  theme_gray(base_size = 8)+
  theme(panel.grid.minor = element_blank())

DLPFC_fabq<-ggplot(data = diff, mapping = aes(x = fabq_total, y = DLPFC_diff)) +
  geom_point(shape = 17, color = "#0f993d", size = 0.5) +
  sm_statCorr(color = "darkgrey", corr_method = "spearman",linetype = "solid", text_size = 2, size = 0.5)+
  theme_gray(base_size = 8)+
  theme(panel.grid.minor = element_blank())

all<-DLPFC_koos + DLPFC_mep + DLPFC_peg+DLPFC_phq + DLPFC_pcs+ DLPFC_fabq+ plot_layout(ncol = 3)
all