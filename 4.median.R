# ---- Median ----
# Check normality for each Tissue_Type*Year combination
# Liver
tempo=subset(all_data, Tissue_Type=="Liver" & year=="1986", select = "median")
shapiro.test(tempo$median)

tempo=subset(all_data, Tissue_Type=="Liver" & year=="1996", select = "median")
shapiro.test(tempo$median)

tempo=subset(all_data, Tissue_Type=="Liver" & year=="2006", select = "median")
shapiro.test(tempo$median)

tempo=subset(all_data, Tissue_Type=="Liver" & year=="2016", select = "median")
shapiro.test(tempo$median)

# Muscle
tempo=subset(all_data, Tissue_Type=="Muscle" & year=="1986", select = "median")
shapiro.test(tempo$median)

tempo=subset(all_data, Tissue_Type=="Muscle" & year=="1996", select = "median")
shapiro.test(tempo$median)

tempo=subset(all_data, Tissue_Type=="Muscle" & year=="2006", select = "median")
shapiro.test(tempo$median)

tempo=subset(all_data, Tissue_Type=="Muscle" & year=="2016", select = "median")
shapiro.test(tempo$median)

# Check homogenity
leveneTest(median~Tissue_Type:year, data=all_data) 
bartlett.test(median~ty, data=all_data)
fligner.test(median~ty, data=all_data)

# anova for median
anova4 = aov(median ~ Genus + Tissue_Type/year, data=all_data)
summary(anova4)
summary(lm(anova4))

r2 = summary(anova4)[[1]][,2]
r2 = r2/sum(r2)
r2

# check the homogenity
plot(anova4, 1)

# check normality
plot(anova4, 2) 
# Extract the residuals
anova4_res = residuals(object = anova4)
# Run Shapiro-Wilk test
shapiro.test(x = anova4_res)

anova4_Tukey = TukeyHSD(anova4)
temp = data.frame(anova4_Tukey$`Tissue_Type:year`)
temp$Comparison = rownames(temp)
temp$Comparison = sub(":","", temp$Comparison)
temp$Comparison = sub(":","", temp$Comparison)
rcompanion::cldList(p.adj ~ Comparison, data = temp, threshold = 0.05)

# plotting median fragment size ----
temp5=all_data
temp5$year2 = temp5$year
temp5$year2 = sub("2016", "4", temp5$year2)
temp5$year2 = sub("2006", "3", temp5$year2)
temp5$year2 = sub("1996", "2", temp5$year2)
temp5$year2 = sub("1986", "1", temp5$year2)
temp5$year2 = as.numeric(temp5$year2)

temp5L = temp5[temp5[,"Tissue_Type"]=="Liver",]
summary(lm(median~year2, data=temp5L))

temp5M = temp5[temp5[,"Tissue_Type"]=="Muscle",]
summary(lm(median~year2, data=temp5M))

p4 = ggplot(temp5, aes(year, median, fill=Tissue_Type)) +
  theme_bw(base_size = 15) +
  geom_boxplot(width = 0.8, outlier.shape = NA) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8, seed = 25), size=1, shape=19) +
  geom_abline(intercept = lm(median~year2, data=temp5L)$coefficients[1], slope= lm(median~year2, data=temp5L)$coefficients[2], linetype="longdash") +
  geom_abline(intercept = lm(median~year2, data=temp5M)$coefficients[1], slope= lm(median~year2, data=temp5M)$coefficients[2]) +
  scale_fill_manual(values=c("gray100", "gray81")) +
  ylim(0,56000) +
  labs(x="Year", y="Median (bp)", fill="Tissue Type") +
  theme(axis.text  = element_text(size=12, colour="black"),
        text = element_text(colour="black"),
        axis.title = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill=NA),
        axis.line=element_line(colour="black", size=0.6, lineend = "square")) +
  annotate("text", x = 3.2, y = 4500, label = "Muscle = 23350 + 3752 year") +
  annotate("text", x = 3.2, y = 1000, parse = TRUE, label = as.character(expression(R^{2}*" = "*"0.14"*","*" "*"p<0.001"))) +
  theme(legend.position="none")
p4
