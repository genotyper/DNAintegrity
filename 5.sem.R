all_data = max_weighted

all_data$Genus = as.character(all_data$Genus)

#assess correlations
cor(all_data[,c(6,8,10:14)])
plot(all_data$maxi, all_data$median)

all_data$skew = -all_data$skew
all_data$kurtosis = -all_data$kurtosis


plot(all_data[,6])
plot(all_data[,8])
plot(all_data[,10]) #
plot(all_data[,11])
plot(all_data[,12]) #
plot(all_data[,13]) #
plot(all_data[,14]) #
#treat outliers as NA and fiml estimate

all_data$w.mean[which(all_data$w.mean == max(all_data$w.mean))] = NA
all_data$skew[which(all_data$skew == min(all_data$skew))] = NA
all_data$kurtosis[which(all_data$kurtosis == min(all_data$kurtosis))] = NA
all_data$variance[which(all_data$variance == max(all_data$variance))] = NA

#scale setting
pom = function(var, add.one = FALSE){
  var = var - min(var, na.rm = TRUE)
  var = var/max(var, na.rm = TRUE)
  if(add.one == TRUE){
    var = var + 1
  }
  return(var)
}

all_data$maxi = pom(all_data$maxi)
all_data$slope = pom(all_data$slope)
all_data$w.mean = pom(all_data$w.mean)
all_data$median = pom(all_data$median)
all_data$variance = pom(all_data$variance)
all_data$skew = pom(all_data$skew)
all_data$kurtosis = pom(all_data$kurtosis)

dim(all_data)
table(all_data$year, all_data$Tissue_Type)
cor(all_data[,c(6,8,10:14)], use = "pairwise.complete.obs")

#####
#options
#####
# #create variable that is orthogonal to tissue type and year - this is the interaction of tissue type and year
# 
# #make tissue type and year integers
# all_data$ttn = as.numeric(as.factor(all_data$Tissue_Type))
# all_data$yn = as.numeric(as.factor(all_data$year))
# #multiple them - highly correlated with both
# all_data$ttyn = all_data$ttn * all_data$yn
# #regress this vector onto tissue type and year
# m = lm(ttyn ~ ttn + yn, all_data)
# #take residual
# all_data$ortho_interaction = resid(m)
# 
# #bin genera with small sample sizes
# all_data$gn = NA
# all_data$gn[all_data$Genus == "Peromyscus"] = 0
# all_data$gn[all_data$Genus == "Reithrodontomys"] = 1
# all_data$gn[all_data$Genus == "Oryzomys"] = 2
# all_data$gn[all_data$Genus == "Sigmodon"] = 3
# all_data$gn[all_data$Genus == "Onychomys"] = 3
# all_data$gn[all_data$Genus == "Neotoma"] = 3
#####
names(all_data)[6] = "mode"

all_data$year = as.numeric(all_data$year)

#####
#models
#####

t(combn(c("mode", "slope", "w.mean", "median", "skew", "kurtosis", "variance"), 3))

m1 = '
dna =~ mode + slope + w.mean
dna ~ year + Genus
'

m2 = '
dna =~ mode + slope + median
dna ~ year + Genus
'

m3 = '
dna =~ mode + slope + skew
dna ~ year + Genus
'

m4 = '
dna =~ mode + slope + kurtosis
dna ~ year + Genus
'

m5 = '
dna =~ mode + slope + variance
dna ~ year + Genus
'

m6 = '
dna =~ mode + w.mean + median
dna ~ year + Genus
'

m7 = '
dna =~ mode + w.mean + skew
dna ~ year + Genus
'

m8 = '
dna =~ mode + w.mean + kurtosis
dna ~ year  + Genus
'

m9 = '
dna =~ mode + w.mean + variance
dna ~ year + Genus
'

m10 = '
dna =~ mode + median + skew
dna ~ year + Genus
'

m11 = '
dna =~ mode + median + kurtosis
dna ~ year + Genus
'

m12 = '
dna =~ mode + median + variance
dna ~ year + Genus
'

m13 = '
dna =~ mode + skew + kurtosis
dna ~ year + Genus
'

m14 = '
dna =~ mode + skew + variance
dna ~ year + Genus
'

m15 = '
dna =~ mode + kurtosis + variance
dna ~ year + Genus
'

m16 = '
dna =~ slope + w.mean + median
dna ~ year + Genus
'

m17 = '
dna =~ slope + w.mean + skew
dna ~ year + Genus
'

m18 = '
dna =~ slope + w.mean + kurtosis
dna ~ year + Genus
'

m19 = '
dna =~ slope + w.mean + variance
dna ~ year + Genus
'

m20 = '
dna =~ slope + median + skew
dna ~ year + Genus
'

m21 = '
dna =~ slope + median + kurtosis
dna ~ year + Genus
'

m22 = '
dna =~ slope + median + variance
dna ~ year + Genus
'

m23 = '
dna =~ slope + skew + kurtosis
dna ~ year + Genus
'

m24 = '
dna =~ slope + skew + variance
dna ~ year + Genus
'

m25 = '
dna =~ slope + kurtosis + variance
dna ~ year + Genus
'

m26 = '
dna =~ w.mean + median + skew
dna ~ year + Genus
'

m27 = '
dna =~ w.mean + median + kurtosis
dna ~ year + Genus
'

m28 = '
dna =~ w.mean + median + variance
dna ~ year + Genus
'

m29 = '
dna =~ w.mean + skew + kurtosis
dna ~ year + Genus
'

m30 = '
dna =~ w.mean + skew + variance
dna ~ year + Genus
'

m31 = '
dna =~ w.mean + kurtosis + variance
dna ~ year + Genus
'

m32 = '
dna =~ median + skew + kurtosis
dna ~ year + Genus
'

m33 = '
dna =~ median + skew + variance
dna ~ year + Genus
'

m34 = '
dna =~ median + kurtosis + variance
dna ~ year + Genus
'

m35 = '
dna =~ skew + kurtosis + variance
dna ~ year + Genus
'
m36 = '
dna =~ mode + slope
dna ~ year + Genus
'

m37 = '
dna =~ mode + w.mean
dna ~ year + Genus
'

m38 = '
dna =~ mode + median
dna ~ year + Genus
'

m39 = '
dna =~ mode + skew
dna ~ year + Genus
'

m40 = '
dna =~ mode + kurtosis
dna ~ year + Genus
'

m41 = '
dna =~ mode + variance
dna ~ year + Genus
'

m42 = '
dna =~ slope + w.mean
dna ~ year + Genus
'

m43 = '
dna =~ slope + median
dna ~ year + Genus
'

m44 = '
dna =~ slope + skew
dna ~ year + Genus
'

m45 = '
dna =~ slope + kurtosis
dna ~ year + Genus
'

m46 = '
dna =~ slope + variance
dna ~ year + Genus
'

m47 = '
dna =~ w.mean + median
dna ~ year + Genus
'

m48 = '
dna =~ w.mean + skew
dna ~ year + Genus
'

m49 = '
dna =~ w.mean + kurtosis
dna ~ year + Genus
'

m50 = '
dna =~ w.mean + variance
dna ~ year + Genus
'

m51 = '
dna =~ median + skew
dna ~ year + Genus
'

m52 = '
dna =~ median + kurtosis
dna ~ year + Genus
'

m53 = '
dna =~ median + variance
dna ~ year + Genus
'

m54 = '
dna =~ skew + kurtosis
dna ~ year + Genus
'

m55 = '
dna =~ skew + variance
dna ~ year + Genus
'

m56 = '
dna =~ kurtosis + variance
dna ~ year + Genus
'

models = vector("list", 56)
models[[1]] = m1
models[[2]] = m2
models[[3]] = m3
models[[4]] = m4
models[[5]] = m5
models[[6]] = m6
models[[7]] = m7
models[[8]] = m8
models[[9]] = m9
models[[10]] = m10
models[[11]] = m11
models[[12]] = m12
models[[13]] = m13
models[[14]] = m14
models[[15]] = m15
models[[16]] = m16
models[[17]] = m17
models[[18]] = m18
models[[19]] = m19
models[[20]] = m20
models[[21]] = m21
models[[22]] = m22
models[[23]] = m23
models[[24]] = m24
models[[25]] = m25
models[[26]] = m26
models[[27]] = m27
models[[28]] = m28
models[[29]] = m29
models[[30]] = m30
models[[31]] = m31
models[[32]] = m32
models[[33]] = m33
models[[34]] = m34
models[[35]] = m35
models[[36]] = m36
models[[37]] = m37
models[[38]] = m38
models[[39]] = m39
models[[40]] = m40
models[[41]] = m41
models[[42]] = m42
models[[43]] = m43
models[[44]] = m44
models[[45]] = m45
models[[46]] = m46
models[[47]] = m47
models[[48]] = m48
models[[49]] = m49
models[[50]] = m50
models[[51]] = m51
models[[52]] = m52
models[[53]] = m53
models[[54]] = m54
models[[55]] = m55
models[[56]] = m56
#####
library(lavaan)

r2l = vector()
r2m = vector()
for(i in 1:56){
  x = sem(models[[i]], all_data, std.lv = TRUE, missing = "fiml", group = "Tissue_Type")
  grab = grep("dna", names(inspect(x, "r2")$Liver))
  l = inspect(x, "r2")$Liver[grab][[1]]
  m = inspect(x, "r2")$Muscle[grab][[1]]
  r2l = c(r2l, l)
  r2m = c(r2m, m)
  print(i)
}

plot(r2l, r2m)

#model 40 and 45 are the top two with the largest R2 and significant regressions
x = sem(models[[45]], all_data, std.lv = TRUE, missing = "fiml", group = "Tissue_Type")
summary(x)
inspect(x, "r2")

#impose weak factorial invariance across groups- so the construct is the same thing (comparable) across groups
#here mode and kurtosis is best performer for maximizing R2
m = '
dna =~ c(v3,v3)*slope + c(v4,v4)*kurtosis
dna ~ year + Genus
'
x =sem(m, all_data, std.lv = TRUE, missing = "fiml", group = "Tissue_Type")
summary(x)
inspect(x, "r2")


m = '
dna =~ c(v3,v3)*mode + c(v4,v4)*kurtosis
dna ~ year + Genus
'
x =sem(m, all_data, std.lv = TRUE, missing = "fiml", group = "Tissue_Type")
summary(x)
inspect(x, "r2")

library(semPlot)
lbls = c("Mode", "Kurtosis", "Year", "Genus", "DNA")

# layout(t(1:2))
f = semPaths(x,"std", residuals = T,edge.label.cex=3, what = "path",
         sizeMan=15, sizeMan2 = 8, nodeLabels=lbls,
         intercepts = F, edge.color = "grey30", title = F, DoNotPlot = TRUE)
png("./both_sem.png", width = 8, height = 4, units = "in", res = 600)
par(mfrow=c(1,2))
plot(f[[1]])
plot(f[[2]])
dev.off()

ggplot(max_weighted, aes(maxi, kurtosis)) + geom_point() + ylim(0,15)

temp5=max_weighted

temp5L = temp5[temp5[,"Tissue_Type"]=="Liver",]
summary(lm(kurtosis~maxi, data=temp5L))

temp5M = temp5[temp5[,"Tissue_Type"]=="Muscle",]
summary(lm(kurtosis~maxi, data=temp5M))

p = ggplot(temp5, aes(maxi, kurtosis)) +
  theme_bw(base_size = 20) + 
  geom_point(size = 3, alpha = .4) +
  geom_abline(intercept = lm(kurtosis~maxi, data=temp5L)$coefficients[1], slope= lm(kurtosis~maxi, data=temp5L)$coefficients[2], linetype="longdash") + 
  geom_abline(intercept = lm(kurtosis~maxi, data=temp5M)$coefficients[1], slope= lm(kurtosis~maxi, data=temp5M)$coefficients[2]) +
  ylim(0,15) + labs(x = "Mode (bp)", y = "Kurtosis")
ggsave("./kurtosis_mode.png", plot = p, device = "png", width = 6, height = 4, units = "in", dpi = 600) # Figure 4

#effect of 1986
all_data = all_data[all_data$year != "1986",]
x = sem(m, all_data, std.lv = TRUE, missing = "fiml", group = "Tissue_Type")
summary(x)
inspect(x, "r2")




