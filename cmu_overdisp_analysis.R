#script for looking at overdispersion in genetic, ct, and wifi data
# also includes isolation period wifi contact counts;
setwd("~/Desktop/projects/covid data analysis/cmu")
library(data.table)
library(fitdistrplus)
library(showtext) 
library(svglite)
library(stringr)
library(readxl)
font_add_google("Montserrat", "montserrat") #load montserrat
showtext_auto()
par(family = 'montserrat')

##### using cluster size data #####
cmu_clusters_data = read.delim("cluster_attribute_metrics.tsv")
cmu_clusters_new = cmu_clusters_data[,c(1,4)]
cmu_clusters_new = data.table(cmu_clusters_new)
colnames(cmu_clusters_new) = c("cluster_num","N")
cluster_sizes_new = cmu_clusters_new
cluster_sizes_new$N = cluster_sizes_new$N-1 #reducing by one to include 0 sized clusters (i.e. no offspring)

svglite("cluster_sizes.svg", width = 8, height = 6)
par(family = 'montserrat')
hist(cluster_sizes_new$N, xlim=c(0,40), breaks=seq(0,40,1), freq=FALSE, ylim = c(0,0.8),
     main = "Distribution of Cluster Offspring", ylab = "Density", xlab = "Offspring", yaxt= "n", 
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5)  
axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8),
     labels = c("0%","20%","40%","60%","80%"),cex.axis=1.5)
cluster_size_fit = fitdist(cluster_sizes_new$N, "nbinom") 
cluster_size_fitpoints = dnbinom(0:38, size=cluster_size_fit$estimate[1], mu=cluster_size_fit$estimate[2])
lines(seq(0.5, 38.5,1), cluster_size_fitpoints, lwd="1", type='b',pch=20,cex=1)
text(20,0.4, paste("k=", round(cluster_size_fit$estimate[1],2), sep=""),cex=1.5)
dev.off()

cluster_sizes_new = cluster_sizes_new[order(-N)]
cluster_pct_80 = sum(cluster_sizes_new$N)*0.8
running_tot = 0
i = 0
while(running_tot < cluster_pct_80){
  i = i+1
  running_tot = running_tot + cluster_sizes_new$N[i]
}
pct_resp_cluster_80 = 100*i/dim(cluster_sizes_new)[1]

##### using contact tracing data #####
ct_data <- read_excel("Broad_data_positives_with_contacts_DEIDENTIFIED.xlsx", sheet = "Sheet2_reformatted")
cmu_ct = data.table(ct_data)

svglite("reported_contacts.svg", width = 8, height = 6)
par(family = 'montserrat')
hist(cmu_ct$total_contacts,xlim=c(0,40), breaks=seq(0,30,1), freq=FALSE, ylim = c(0,0.8),
     main = "Distribution of Reported Contacts", ylab = "Density", xlab = "Contacts", yaxt= "n", 
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8),
     labels = c("0%","20%","40%","60%","80%"),cex.axis=1.5)
ct_fit = fitdist(cmu_ct$total_contacts, "nbinom") 
ct_fitpoints = dnbinom(0:38, size=ct_fit$estimate[1], mu=ct_fit$estimate[2])
lines(seq(0.5, 38.5,1), ct_fitpoints, lwd="1", type='b',pch=20,cex=1)
text(20,0.4, paste("k=", round(ct_fit$estimate[1],2), sep=""), cex=1.5)
dev.off()

cmu_ct = cmu_ct[order(-total_contacts)]
ct_pct_80 = sum(cmu_ct$total_contacts)*0.8
running_tot = 0
i = 1
while(running_tot < ct_pct_80){
  running_tot = running_tot + cmu_ct$total_contacts[i]
  i = i+1
}
pct_resp_ct_80 = 100*i/dim(cmu_ct)[1]

##### using wifi data #####
cmu_meta = read.csv("metadata_all.csv")
cmu_meta_table = data.table(cmu_meta)
cmu_meta_table[Symptom.Onset == 'unknown']$Symptom.Onset = NA
cmu_meta_table[Symptom.Onset == 'asymptomatic']$Symptom.Onset = NA
barcode_dates = data.table(cmu_meta_table[,.(barcode,as.Date(Test.Date), as.Date(Symptom.Onset))])
barcode_dates$min = NA
for(i in 1:dim(barcode_dates)[1]){
  barcode_dates$min[i] = min(barcode_dates$V2[i], barcode_dates$V3[i], na.rm=TRUE)
}
barcode_dates = data.table(barcode_dates$barcode, as.Date(barcode_dates$min, "1970/1/1"))
colnames(barcode_dates) = c("indA", "trace_date")

fall_wifi_data = read.csv("final_fall.csv")
spring_wifi_data = read.csv("final_spring.csv")

fall_wifi_table = data.table(fall_wifi_data)
fall_wifi_table$day = as.Date(fall_wifi_table$day)
fall_wifi_table$pair = str_remove_all(fall_wifi_table$pair, "[\\['\\]]")
fall_wifi_table = fall_wifi_table[, c("indA","indB") := tstrsplit(pair, ", ")]
spring_wifi_table = data.table(spring_wifi_data)
spring_wifi_table$day = as.Date(spring_wifi_table$day)
spring_wifi_table$pair = str_remove_all(spring_wifi_table$pair, "[\\['\\]]")
spring_wifi_table = spring_wifi_table[, c("indA","indB") := tstrsplit(pair, ", ")]
year_wifi_table = rbind(fall_wifi_table, spring_wifi_table)

year_wifi_table_dated = barcode_dates[year_wifi_table, on=c(indA = "indA")]
year_wifi_table_dated = barcode_dates[year_wifi_table_dated, on=c(indA = "indB")]
colnames(year_wifi_table_dated) = c("indA", "testA", "indB", "testB", "day", "pair", "sum","median", "count")
traced_interactions = year_wifi_table_dated[(day <= testB & day >= testB-2) | (day <= testA & day >= testA-2),]

traced_interactions_a = traced_interactions[,c(1,2,5)]
traced_interactions_b = traced_interactions[,c(3,4,5)]

all_traced = rbind(traced_interactions_a, traced_interactions_b, use.names = FALSE)
all_traced = all_traced[day <= testA & day >= testA-2, ]
interaction_table = all_traced[,.N,indA]

svglite("wifi_contacts.svg", width = 8, height = 6)
par(family = 'montserrat')
hist(interaction_table$N,xlim=c(0,1200), breaks=seq(0,1200,15), freq=FALSE, ylim = c(0,0.01), 
     main = "Distribution of WiFi Contacts", ylab = "Density", xlab = "Contacts", yaxt= "n", 
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
axis(2, at = c(0, 0.0025, 0.005, 0.0075, 0.01),
     labels = c("0%","0.25%","0.5%","0.75%","1%"),cex.axis=1.5)
traced_fit = fitdist(interaction_table$N, "nbinom") 
traced_fitpoints = dnbinom(seq(0, 1200,20), size=traced_fit$estimate[1], mu=traced_fit$estimate[2])
lines(seq(0.5, 1200.5,20), traced_fitpoints, lwd="1", type='b',pch=20,cex=1)
text(600,0.005, paste("k=", round(traced_fit$estimate[1],2), sep=""), cex=1.5)
dev.off()

interaction_table = interaction_table[order(-N)]
wifi_pct_80 = sum(interaction_table$N)*0.8
running_tot = 0
i = 1
while(running_tot < wifi_pct_80){
  running_tot = running_tot + interaction_table$N[i]
  i = i+1
}
pct_resp_wifi_80 = 100*i/dim(interaction_table)[1]


##### sup: isolation period contacts from wifi data #####
pre_interval = 10
post_interval = 20
daily_table = year_wifi_table_dated[(day >= testB-pre_interval & day <= testB+post_interval) 
                                         | (day >= testA-pre_interval & day <= testA+post_interval),]
daily_table$offsetA = as.numeric(daily_table$day-daily_table$testA)
daily_table$offsetB = as.numeric(daily_table$day-daily_table$testB)
daily_table$pospos = as.numeric(daily_table$testA-daily_table$testB)
daily_table_a = daily_table[,c(1,2,5,10,12)]
daily_table_b = daily_table[,c(3,4,5,11,12)]
daily_ints = rbind(daily_table_a,daily_table_b, use.names = FALSE)
daily_ints = daily_ints[offsetA <= post_interval & offsetA >= -pre_interval]
daily_ints[10 >= pospos & pospos >= -10,]$pospos = 1
daily_ints[pospos != 1 | is.na(pospos),]$pospos = 0
daily_ints_user = daily_ints[,list(.N, (sum(pospos)/length(pospos))),.(indA,offsetA)]
daily_meds = daily_ints_user[,.(median(N), sd(N)),offsetA][order(offsetA)]

pre_mean = mean(daily_meds[1:11,]$V1)
iso_mean = mean(daily_meds[11:21,]$V1)
post_mean = mean(daily_meds[21:31,]$V1)

svglite("isolation_period_wifi.svg", width = 11, height = 4)
par(family = 'montserrat')
plot((daily_meds$V1-pre_mean)/pre_mean, type='b',xaxt="n",xlab="Days Post-Positive",ylab="Median Interactions",
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
axis(1, at=seq(1,diff(range(daily_meds$offsetA))+1,2),
     labels=daily_meds$offsetA[seq(1,diff(range(daily_meds$offsetA))+1,2)],cex.axis=1.5)
lines(seq(1,11,1),rep(0,11),lwd=2,col='darkgrey')
dif = -(pre_mean-iso_mean)/pre_mean
lines(seq(11,21,1), rep(dif,11),lwd=2,col='darkgrey')
abline(v=11,lty=3,col='red')
post_dif = -(pre_mean-post_mean)/pre_mean
lines(seq(21,31,1), rep(post_dif,11),lwd=2,col='darkgrey')
abline(v=21,lty=3,col='red')
dev.off()