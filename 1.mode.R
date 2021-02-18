library(ggplot2)
library(reshape2)
library(car)

# import data set
dat=read.csv("./merged_sorted_cut.csv")

# create metadata (met) ----
met = data.frame(sample = names(dat[,2:ncol(dat)]))
met$organ = met$sample
met$organ = sub(".*_", "", met$organ)
met$year = met$sample
met$year = sub("TK.*_Y", "", met$year)
met$year = sub("_L", "", met$year)
met$year = sub("_M", "", met$year)
met$TK= met$sample
met$TK = sub("_.*", "", met$TK)

# import another data set
spec_exam=read.csv("./specimens_examined.csv")

# merging met and spec_exam
met=merge(met, spec_exam, by="sample")

# prepare data frame "max" ----####
max = c()
for(i in 2:ncol(dat)){
  s = max(grep(max(dat[,i]), dat[,i]))
  m = dat$Size..bp.[s]
  n = names(dat[i])
  temp = data.frame(sample = n, maxi = m)
  max = rbind(max, temp)
  print(paste(i))
}

max = merge(met, max, by = "sample")
max$organ=sub("L", "Liver", max$organ)
max$organ=sub("M", "Muscle", max$organ)
names(max)[2]<-"Tissue_Type"

boxplot(maxi~year, max)
which(max$maxi>150000)
max = max[-78,]
boxplot(maxi~year, max)

# Check normality for each Tissue_Type*Year combination
# Liver
tempo=subset(max, Tissue_Type=="Liver" & year=="1986", select = "maxi")
shapiro.test(tempo$maxi)

tempo=subset(max, Tissue_Type=="Liver" & year=="1996", select = "maxi")
shapiro.test(tempo$maxi)

tempo=subset(max, Tissue_Type=="Liver" & year=="2006", select = "maxi")
shapiro.test(tempo$maxi)

tempo=subset(max, Tissue_Type=="Liver" & year=="2016", select = "maxi")
shapiro.test(tempo$maxi)

# Muscle
tempo=subset(max, Tissue_Type=="Muscle" & year=="1986", select = "maxi")
shapiro.test(tempo$maxi)

tempo=subset(max, Tissue_Type=="Muscle" & year=="1996", select = "maxi")
shapiro.test(tempo$maxi)

tempo=subset(max, Tissue_Type=="Muscle" & year=="2006", select = "maxi")
shapiro.test(tempo$maxi)

tempo=subset(max, Tissue_Type=="Muscle" & year=="2016", select = "maxi")
shapiro.test(tempo$maxi)

# Check homogenity
max$ty = paste0(max$Tissue_Type, max$year)
max$ty = sub("Liver", "L", max$ty)
max$ty = sub("Muscle", "M", max$ty)

library(car)
leveneTest(maxi~ty, data=max) # Not Homogen
bartlett.test(maxi~ty, data=max) # Homogen
fligner.test(maxi~ty, data=max) # Homogen

# anova for mode
anova1 = aov(maxi ~ Genus + Tissue_Type/year, data = max)
summary(anova1)
summary(lm(anova1))

r2 = summary(anova1)[[1]][,2]
r2 = r2/sum(r2)
r2

# checking homogenity variance on model anova1
plot(anova1, 1)
# 3 different tests were run above, 2 homogen, 1 is not homogen
# check normality
plot(anova1, 2) 
# Extract the residuals
anova1_res = residuals(object = anova1)
# Run Shapiro-Wilk test
shapiro.test(x = anova1_res) # the data is normal 

TukeyHSD(anova1)
anova1_Tukey = TukeyHSD(anova1)
temp = data.frame(anova1_Tukey$`Tissue_Type:year`)
temp$Comparison = rownames(temp)
temp$Comparison = sub(":","", temp$Comparison)
temp$Comparison = sub(":","", temp$Comparison)
library(rcompanion)
cldList(p.adj ~ Comparison, data = temp, threshold = 0.05)

# plotting mode
temp5=max
temp5$year2 = temp5$year
temp5$year2 = sub("2016", "4", temp5$year2)
temp5$year2 = sub("2006", "3", temp5$year2)
temp5$year2 = sub("1996", "2", temp5$year2)
temp5$year2 = sub("1986", "1", temp5$year2)
temp5$year2 = as.numeric(temp5$year2)

temp5L = temp5[temp5[,"Tissue_Type"]=="Liver",]
summary(lm(maxi~year2, data=temp5L))

temp5M = temp5[temp5[,"Tissue_Type"]=="Muscle",]
summary(lm(maxi~ year2 + Genus, data=temp5M))

library(ggplot2)
p = ggplot(temp5, aes(year, maxi)) +
  theme_bw(base_size = 15) +
  geom_boxplot(aes(year, maxi, fill=Tissue_Type), width = 0.8, outlier.shape = NA) +
  geom_point(aes(year, maxi, fill=Tissue_Type),position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8, seed = 25), size=1, shape=19) +
  geom_abline(intercept = lm(maxi~year2, data=temp5L)$coefficients[1], slope= lm(maxi~year2, data=temp5L)$coefficients[2], linetype="longdash") +
  geom_abline(intercept = lm(maxi~year2, data=temp5M)$coefficients[1], slope= lm(maxi~year2, data=temp5M)$coefficients[2]) +
  scale_fill_manual(values=c("gray100", "gray81")) +
  ylim(0, 65000) +
  labs(x="Year", y="Mode (bp)", fill="Tissue Type") +
  theme (axis.text  = element_text(size=12, colour="black"),
         text = element_text(colour="black"),
         axis.title = element_text(size = 12),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.background = element_rect(fill=NA),
         axis.line=element_line(colour="black", size=0.6, lineend = "square")) +
  annotate("text", x = 1.8, y = 65000, label = "Muscle = 13787 + 7322 year") +
  annotate("text", x = 1.8, y = 60000, parse = TRUE, label = as.character(expression(R^{2}*" = "*"0.214"*","*" "*"p<0.001"))) + 
  theme(legend.position="none")
p
