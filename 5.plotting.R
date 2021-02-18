# plotting peak ----
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
ggsave("./Peak_trendline3.png", plot = p, device = "png", width = 6, height = 4, units = "in", dpi = 300) # Figure 1
# end of plotting peak ====

# plotting regression coefficient ----
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
ggsave("./Regression_Coefficient3.png", plot = p2, device = "png", width = 6, height = 4, units = "in", dpi = 300) # Figure 2.1
# end of plotting regression coefficient ====

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
ggsave("./Average_Fragment4.png", plot = p3, device = "png", width = 6, height = 4, units = "in", dpi = 300) # Figure 4
# end of plotting average fragment size ====

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
# ggsave("./median_Fragment4.png", plot = p4, device = "png", width = 6, height = 4, units = "in", dpi = 300)

library(patchwork)
png(file = "figures/combined.png", height = 8.5, width = 10, res = 600, units = "in")
print(p3 + p4 + p + p2 + plot_annotation(tag_levels = 'A'))
dev.off()