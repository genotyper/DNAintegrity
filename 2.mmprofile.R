# proportion
# ---- using functions to create prop & slope_df objects ##### 
over= c() # created empty object for over data frame
for(i in 2: ncol(dat)){
  temp = c()
  for (j in c(seq(0, 50000, 5000))){
    n = names(dat[i])
    kol=colSums(dat[i])
    a = sum(dat[which(dat$Size..bp.>j), i])
    b = a/kol
    temp1 = data.frame(sample=n, b=b)
    temp = cbind(temp, temp1$b)
  }
  temp2= data.frame(sample=n, temp)
  over= rbind(over, temp2)
  print(paste(i))
}

nama = c()
for (j in c(seq(0, 50000, 5000))){
  nama_temp = paste("over", j/1000, "k", sep="")
  nama = c(nama, nama_temp)
}

colnames(over)
colnames(over) = c("sample", paste(nama))

prop = merge(max, over, by= "sample")
prop = prop[,c(-2:-7)]

library(reshape2)
prop = melt(prop, id.vars="sample")
prop = merge(max, prop, by="sample")
prop$variable2 = prop$variable

replaces <- data.frame(from = c(nama), to = c(1:length(nama)))

prop = DataCombine::FindReplace(prop, Var = "variable2", replaceData = replaces, from = "from", to = "to", exact = FALSE)
prop$variable2=as.numeric(prop$variable2)

rownames(over)
row.names(over) = over$sample
over = over[,-1] 
over_t = as.data.frame(t(over))
over_t = cbind(thres=replaces$to, over_t) # adding thres column to over_t from replaces$to

y_col = c(paste(colnames(over_t[,-1])))
x_col = "thres"

library(dplyr)
slope_df = as.data.frame(expand.grid(y=y_col, x=x_col, stringsAsFactors = F) %>%
                           mutate(formula = paste("I", "(", y, "-1", ")", "~", x, "-1")) %>%
                           group_by(formula) %>%
                           mutate(slope = summary(lm(formula, data=over_t))$coefficients[1,1]) %>%
                           ungroup())

slope_df = slope_df[,c(-2,-3)]
names(slope_df)[1] = "sample"

slope_df = merge(max, slope_df, by="sample")
slope_df$yen = 1+slope_df$slope*11 # add yen column for multiple geom_segment
temp = slope_df[,c(1,8,9)]

prop = merge(prop, temp, by="sample")

q = ggplot(prop, aes(NULL)) + 
  geom_point(aes(x=variable, y=value), position = position_jitter(width=.35, height=0, seed=25), shape=19) +
  theme_bw(base_size = 15) +
  labs(x="Base Pair Threshold", y="Proportion of DNA Fragments") +
  geom_segment(x = 0, y = 1, xend=11, yend=prop$yen, size=0.4) +
  facet_grid(Tissue_Type~year) +
  theme(legend.position="none", axis.text.x = element_text(angle=35, hjust=0.8, vjust=0.8))
q
ggsave("./figures/Percentage.png", plot = q, device = "png", width = 14, height = 8, units = "in", dpi = 200) # Figure 2

# ---- Molecular mass profile ----
# Check normality for each Tissue_Type*Year combination
# Liver
tempo=subset(slope_df, Tissue_Type=="Liver" & year=="1986", select = "slope")
shapiro.test(tempo$slope)

tempo=subset(slope_df, Tissue_Type=="Liver" & year=="1996", select = "slope") #
shapiro.test(tempo$slope)

tempo=subset(slope_df, Tissue_Type=="Liver" & year=="2006", select = "slope")
shapiro.test(tempo$slope)

tempo=subset(slope_df, Tissue_Type=="Liver" & year=="2016", select = "slope") #
shapiro.test(tempo$slope)

# Muscle
tempo=subset(slope_df, Tissue_Type=="Muscle" & year=="1986", select = "slope")
shapiro.test(tempo$slope)

tempo=subset(slope_df, Tissue_Type=="Muscle" & year=="1996", select = "slope")
shapiro.test(tempo$slope)

tempo=subset(slope_df, Tissue_Type=="Muscle" & year=="2006", select = "slope")
shapiro.test(tempo$slope)

tempo=subset(slope_df, Tissue_Type=="Muscle" & year=="2016", select = "slope")
shapiro.test(tempo$slope)

car::leveneTest(slope~ty, data=slope_df)
bartlett.test(slope~ty, data=slope_df)
fligner.test(slope~ty, data=slope_df)

# anova for slope
anova2 = aov(slope ~ Genus + Tissue_Type/year, data=slope_df) 
summary(anova2)
summary(lm(anova2))

r2 = summary(anova2)[[1]][,2]
r2 = r2/sum(r2)
r2

# check the homogenity
plot(anova2, 1)

# check normality
plot(anova2, 2) 
# Extract the residuals
anova2_res = residuals(object = anova2)
# Run Shapiro-Wilk test
shapiro.test(x = anova2_res) # the data is NOT normal

TukeyHSD(anova2)
anova2_Tukey = TukeyHSD(anova2)
temp = data.frame(anova2_Tukey$`Tissue_Type:year`)
temp$Comparison = rownames(temp)
temp$Comparison = sub(":","", temp$Comparison)
temp$Comparison = sub(":","", temp$Comparison)
rcompanion::cldList(p.adj ~ Comparison, data = temp, threshold = 0.05)

# plotting regression coefficient/Molecular Mass Profile ----
temp5=slope_df
temp5$year2 = temp5$year
temp5$year2 = sub("2016", "4", temp5$year2)
temp5$year2 = sub("2006", "3", temp5$year2)
temp5$year2 = sub("1996", "2", temp5$year2)
temp5$year2 = sub("1986", "1", temp5$year2)
temp5$year2 = as.numeric(temp5$year2)

temp5L = temp5[temp5[,"Tissue_Type"]=="Liver",]
summary(lm(slope~year2, data=temp5L))

temp5M = temp5[temp5[,"Tissue_Type"]=="Muscle",]
summary(lm(slope~year2, data=temp5M))

p2 = ggplot(temp5, aes(x=year, y=slope, fill=Tissue_Type)) +
  theme_bw(base_size = 15) +
  geom_boxplot(width = 0.8, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8, seed = 25), size=1, shape=19) +
  geom_abline(intercept = lm(slope~year2, data=temp5L)$coefficients[1], slope= lm(slope~year2, data=temp5L)$coefficients[2], linetype="longdash") +
  geom_abline(intercept = lm(slope~year2, data=temp5M)$coefficients[1], slope= lm(slope~year2, data=temp5M)$coefficients[2]) +
  scale_fill_manual(values=c("gray100", "gray81")) +
  labs(x="Year", y="Molecular Mass Profile", fill="Tissue Type") +
  theme (axis.text  = element_text(size=12, colour="black"),
         text = element_text(colour="black"),
         axis.title = element_text(size = 12),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.background = element_rect(fill=NA),
         axis.line=element_line(colour="black", size=0.6, lineend = "square")) +
  annotate("text", x = 3.1, y = -0.12, label = "Muscle = -0.08327 + 0.00699 year") +
  annotate("text", x = 3.1, y = -0.125, parse = TRUE, label = as.character(expression(R^{2}*" = "*"0.14"*","*" "*"p<0.001"))) +
  theme(legend.position="none")
p2
