# ---- weighted mean ####

rata =c()
for(i in 2:ncol(dat)){
  n = names(dat[i])
  w = weighted.mean(dat[,1], dat[,i])
  temp = data.frame(sample=n, w.mean=w)
  rata = rbind(rata, temp)
  print(paste(i))
}

head(rata)
str(rata)
head(slope_df)
max_weighted = merge(slope_df, rata, by="sample")

library(moments)
md = c()
av = c()
sk = c()
kr = c()
vr = c()
for(j in 2:ncol(dat)){
  l = vector("list", nrow(dat))
  for(i in 1:length(l)){
    l[[i]] = dat[i,1]
    l[[i]] = rep(l[[i]], dat[i,j])
  }
  l = unlist(l)
  
  m = data.frame(sample = names(dat)[j], median = median(l))
  a = data.frame(sample = names(dat)[j], meand = mean(l))
  s = data.frame(sample = names(dat)[j], skew = skewness(l))
  k = data.frame(sample = names(dat)[j], kurtosis = kurtosis(l))
  v = data.frame(sample = names(dat)[j], variance = var(l))
  
  md = rbind(md, m)
  av = rbind(av, a)
  sk = rbind(sk, s)
  kr = rbind(kr, k)
  vr = rbind(vr, v)

    print(paste(j, " of ", ncol(dat)))
}


max_weighted = merge(max_weighted, md, by = "sample")
max_weighted = merge(max_weighted, sk, by = "sample")
max_weighted = merge(max_weighted, kr, by = "sample")
max_weighted = merge(max_weighted, vr, by = "sample")

all_data = max_weighted

# plotting weighted.mean
boxplot(w.mean ~ year, max_weighted) # one outlier on 2006 > 80000, 3 outliers on 1996 below 20000
qqnorm(max_weighted$w.mean)
qqline(max_weighted$w.mean)
which(max_weighted$w.mean > 80000)
max_weighted2 = max_weighted[-18,]

# Check normality for each Tissue_Type*Year combination
# Liver
tempo=subset(max_weighted2, Tissue_Type=="Liver" & year=="1986", select = "w.mean") 
shapiro.test(tempo$w.mean)

tempo=subset(max_weighted2, Tissue_Type=="Liver" & year=="1996", select = "w.mean")
shapiro.test(tempo$w.mean)

tempo=subset(max_weighted2, Tissue_Type=="Liver" & year=="2006", select = "w.mean")
shapiro.test(tempo$w.mean)

tempo=subset(max_weighted2, Tissue_Type=="Liver" & year=="2016", select = "w.mean")
shapiro.test(tempo$w.mean)

# Muscle
tempo=subset(max_weighted2, Tissue_Type=="Muscle" & year=="1986", select = "w.mean")
shapiro.test(tempo$w.mean)

tempo=subset(max_weighted2, Tissue_Type=="Muscle" & year=="1996", select = "w.mean")
shapiro.test(tempo$w.mean)

tempo=subset(max_weighted2, Tissue_Type=="Muscle" & year=="2006", select = "w.mean")
shapiro.test(tempo$w.mean)

tempo=subset(max_weighted2, Tissue_Type=="Muscle" & year=="2016", select = "w.mean")
shapiro.test(tempo$w.mean)

car::leveneTest(w.mean~ty, data=max_weighted2)
bartlett.test(w.mean~ty, data=max_weighted2)
fligner.test(w.mean~ty, data=max_weighted2)


# anova for mean
anova3 = aov(w.mean ~ Genus + Tissue_Type/year, data=max_weighted2)
summary(anova3)
summary(lm(anova3))

r2 = summary(anova3)[[1]][,2]
r2 = r2/sum(r2)
r2

# check the homogenity
plot(anova3, 1)

# check normality
plot(anova3, 2) 
# Extract the residuals
anova3_res = residuals(object = anova3)
# Run Shapiro-Wilk test
shapiro.test(x = anova3_res) # the data is NOT normal

anova3_Tukey = TukeyHSD(anova3)
temp = data.frame(anova3_Tukey$`Tissue_Type:year`)
temp$Comparison = rownames(temp)
temp$Comparison = sub(":","", temp$Comparison)
temp$Comparison = sub(":","", temp$Comparison)
rcompanion::cldList(p.adj ~ Comparison, data = temp, threshold = 0.05)

# plotting average fragment size ----
temp5=max_weighted2
temp5$year2 = temp5$year
temp5$year2 = sub("2016", "4", temp5$year2)
temp5$year2 = sub("2006", "3", temp5$year2)
temp5$year2 = sub("1996", "2", temp5$year2)
temp5$year2 = sub("1986", "1", temp5$year2)
temp5$year2 = as.numeric(temp5$year2)

temp5L = temp5[temp5[,"Tissue_Type"]=="Liver",]
summary(lm(w.mean~year2, data=temp5L))

temp5M = temp5[temp5[,"Tissue_Type"]=="Muscle",]
summary(lm(w.mean~year2, data=temp5M))

p3 = ggplot(temp5, aes(year, w.mean, fill=Tissue_Type)) +
  theme_bw(base_size = 15) +
  geom_boxplot(width = 0.8, outlier.shape = NA) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8, seed = 25), size=1, shape=19) +
  geom_abline(intercept = lm(w.mean~year2, data=temp5L)$coefficients[1], slope= lm(w.mean~year2, data=temp5L)$coefficients[2], linetype="longdash") +
  geom_abline(intercept = lm(w.mean~year2, data=temp5M)$coefficients[1], slope= lm(w.mean~year2, data=temp5M)$coefficients[2]) +
  scale_fill_manual(values=c("gray100", "gray81")) +
  ylim(0,65000) +
  labs(x="Year", y="Mean (bp)", fill="Tissue Type") +
  theme(axis.text  = element_text(size=12, colour="black"),
        text = element_text(colour="black"),
        axis.title = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill=NA),
        axis.line=element_line(colour="black", size=0.6, lineend = "square")) +
  annotate("text", x = 3.2, y = 4500, label = "Muscle = 29008.4 + 3225.7 year") +
  annotate("text", x = 3.2, y = 1000, parse = TRUE, label = as.character(expression(R^{2}*" = "*"0.125"*","*" "*"p<0.01"))) +
  theme(legend.position=c(.2, .8))
p3
