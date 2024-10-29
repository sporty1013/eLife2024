library(sjstats)
library(lme4)
library(lmerTest)
library(MuMIn)
library(effsize)
library(Rmisc)
library(smplot2)

eLife2024<-read.csv("eLife2024.csv", header = TRUE, na.strings = "") 

# Figure 3
# Repeated ANCOVA
# Mean difference (95% CI) & Cohen's d
model<-lmer(VLPFC ~ Task + age +sex + bmi + (1|id), data = eLife2024)
estimate<-round(fixef(model)[2],2)
confint<-round(confint(model)[4,],2)
cat(estimate,"(",confint,")\n")
cohen<-cohens_d(VLPFC ~ Task, data = eLife2024)
cohend<-round(cohen$Cohens_d,2)
cohen_low<-round(cohen$CI_low,2)
cohen_high<-round(cohen$CI_high,2)
cat(cohend,"(",cohen_low,",",cohen_high,")\n")

model<-lmer(VMPFC ~ Task + age +sex + bmi + (1|id), data = eLife2024)
estimate<-round(fixef(model)[2],2)
confint<-round(confint(model)[4,],2)
cat(estimate,"(",confint,")\n")
cohen<-cohens_d(VMPFC ~ Task, data = eLife2024)
cohend<-round(cohen$Cohens_d,2)
cohen_low<-round(cohen$CI_low,2)
cohen_high<-round(cohen$CI_high,2)
cat(cohend,"(",cohen_low,",",cohen_high,")\n")


model<-lmer(DLPFC ~ Task + age +sex + bmi + (1|id), data = eLife2024)
estimate<-round(fixef(model)[2],2)
confint<-round(confint(model)[4,],2)
cat(estimate,"(",confint,")\n")
cohen<-cohens_d(DLPFC ~ Task, data = eLife2024)
cohend<-round(cohen$Cohens_d,2)
cohen_low<-round(cohen$CI_low,2)
cohen_high<-round(cohen$CI_high,2)
cat(cohend,"(",cohen_low,",",cohen_high,")\n")

model<-lmer(DMPFC ~ Task + age +sex + bmi + (1|id), data = eLife2024)
estimate<-round(fixef(model)[2],2)
confint<-round(confint(model)[4,],2)
cat(estimate,"(",confint,")\n")
cohen<-cohens_d(DMPFC ~ Task, data = eLife2024)
cohend<-round(cohen$Cohens_d,2)
cohen_low<-round(cohen$CI_low,2)
cohen_high<-round(cohen$CI_high,2)
cat(cohend,"(",cohen_low,",",cohen_high,")\n")

# boxplots
df<-subset(eLife2024, select = c('id','Task', "VLPFC", "VMPFC", "DLPFC", "DMPFC"))

df_stacked<-df %>% pivot_longer(cols=c("VLPFC", "VMPFC", "DLPFC", "DMPFC"),
                                names_to='ROI',
                                values_to='HbO')
df_stacked

HbO<-ggplot(df_stacked, aes(x = ROI, y = HbO)) +
  geom_boxplot(aes(fill = Task), show.legend=T) +
  labs(x = "",
       y = "HbO (µmol/L mm)") +
  scale_x_discrete(limits=c("VLPFC", "VMPFC", "DLPFC", "DMPFC")) +
  theme_hc(base_size = 14)+
  scale_fill_manual(values=c("#f1eef6", "#69b3a2"))

HbO

