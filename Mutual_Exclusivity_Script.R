#PROP_Test
setwd("")
library(reshape2)


DF1 = read.delim("RANDOM_Nets_TCGA_mutational_data.txt", header = T)
DF2 = read.delim("CM_Nets_TCGA_mutational_data.txt", header = T)
AR_Net_nam = t(DF2[,1])
Random = DF1[,1]
row.names(DF1) = DF1[,1]
DF1 = DF1[,-1]
row.names(DF2) = DF2[,1]
DF2 = DF2[,-1]




######### PROPORTION TEST for enrichment  analysis ##########################################
#                                                                                           #
# This code calculatess the statistical significance of the enrichment of mutated genes in  #
# the CM-GRNs versus randomly generated ones by using the proportion test.                  #
#                                                                                           #
#############################################################################################




test <- list()
k=1
for (i in 1:dim(DF2)[1]){
  for(j in 1:dim(DF1)[1]){
    test[[k]] <- prop.test(rbind(as.numeric(DF1[j,]), as.numeric(DF2[i,])), correct=TRUE, alternative="two.sided", conf.level=.95)
    k=k+1
  }
}




test <- data.frame(matrix(unlist(test), nrow=48000, byrow=T))
  test = as.data.frame(test[,3])
   colnames(test) = c("p_Value")
    index = as.data.frame(rep(1:48, each = 1000, len = 48000))
      colnames(index) = c("Index")
        cbnd = as.data.frame(cbind(index, test))
          spl = as.data.frame(split(cbnd, cbnd$Index))
            Ind = paste("Index", 1:48, sep = "")
          P_Val = paste("P_Val", 1:48, sep = "")
        cb = as.data.frame(as.vector(t(cbind(Ind, P_Val))))
      colnames(cb) = c("colnames")
    cb = t(cb)
  colnames(spl) = cb
spl <- spl[ , -grep("\\Index.", colnames(spl))]
colnames(spl) = AR_Net_nam

chi_sq_p_val_0.95_conf_0_ALL = cbind(Random, spl)
  row.names(chi_sq_p_val_0.95_conf_0_ALL) = chi_sq_p_val_0.95_conf_0_ALL[,1]
  write.table(chi_sq_p_val_0.95_conf_0_ALL, "chi_sq_p_val_0.95_conf_0_ALL.txt", col.names = T, row.names = F, sep = "\t")

    chi_sq_p_val_0.95_conf_0_ALL = chi_sq_p_val_0.95_conf_0_ALL[,-1]
      chi_sq_p_val_0.95_conf_0_ALL <- as.data.frame(lapply(chi_sq_p_val_0.95_conf_0_ALL, function(x) {as.numeric(as.character(x))}))
        BH_correction_0.95_conf_0_ALL <- as.data.frame(lapply(chi_sq_p_val_0.95_conf_0_ALL, p.adjust, method="BH", n = 1000))
          BH_correction_0.95_conf_0_ALL = cbind(Random, BH_correction_0.95_conf_0_ALL)
          write.table(BH_correction_0.95_conf_0_ALL, "BH_correction_0.95_conf_0_ALL.txt", col.names = T, row.names = F, sep = "\t")


tmp1 = read.delim("BH_correction_0.95_conf_0_ALL.txt", header = T) 
  tmp1 = as.data.frame(t(tmp1))
   cln = as.data.frame(rownames(tmp1))
    cb1 = cbind(cln, tmp1)
     cln =sub("Random", "Net", cb1[,1]) 
      cb1 = cb1[,-1]
       cb1 = cbind(cln, tmp1)
        Ar_Net_names = as.matrix(cb1[,1])
        Ar_Net_names = as.data.frame(Ar_Net_names[2:49])
         colnames(Ar_Net_names) = c("Net")

         cb1 = cb1[,-1]
          Rnd_Net_names = as.character(as.matrix(cb1[1,]))
           cb1 = cb1[-1,]
            cb1 <- as.data.frame(lapply(cb1, function(x) {as.numeric(as.character(x))}))

below= c("0.1")
  below = as.matrix(as.data.frame(rep(below, 48)))
      below = as.data.frame(as.numeric(as.character(below)))
        colnames(below) = c("below")


colnames(cb1) = Rnd_Net_names
cb = cbind(Ar_Net_names, below, cb1)


x = within(cb, {
below = rowSums(cb[-c(1:2)] < cb$below)})

write.table(x, "BH_correction_0.95_conf_0_ALL_0.1_below.txt", col.names = T, row.names = F, sep = "\t")

# FREQUENCIES DIFFERENCES COMPUTING
# First, you have to compute the diferences of frequencies of 0, 1, 2 mutations and so on.
DF1 = read.delim("RANDOM_Net_for_chi_sq_0_ALL.txt", header = T)
DF2 = read.delim("ARACNe_Net_for_chi_sq_0_ALL.txt", header = T)
AR_Net_nam = t(DF2[,1])
Random = DF1[,1]
  row.names(DF1) = DF1[,1]
   DF1 = as.data.frame(DF1[, -c(1,3)])
    colnames(DF1) = c("#1#_freq")
     row.names(DF2) = DF2[,1]
      DF2 = as.data.frame(DF2[,-c(1,3)])
       colnames(DF2) = c("#1#_freq")

test <- list()
k=1
for (i in 1:dim(DF2)[1]){
  for(j in 1:dim(DF1)[1]){
    test[[k]] <- as.numeric(DF2[i,]) - as.numeric(DF1[j,])
    k=k+1
  }
}
test = as.character(test)

index = as.data.frame(rep(1:48, each = 1000, len = 48000))
  colnames(index) = c("Index")
    cbnd = as.data.frame(cbind(index, test))
      spl = as.data.frame(split(cbnd, cbnd$Index))
        Ind = paste("Index", 1:48, sep = "")
          Diff = paste("Diff", 1:48, sep = "")
            cb = as.data.frame(as.vector(t(cbind(Ind, Diff))))
              colnames(cb) = c("colnames")
               cb = t(cb)
                 colnames(spl) = cb
             spl <- spl[ , -grep("\\Index.", colnames(spl))]
          colnames(spl) = AR_Net_nam
     Diff_0_ARACNe_RANDOM = cbind(Random, spl)

write.table(Diff_0_ARACNe_RANDOM, "Diff_0_ARACNe_RANDOM.txt", col.names = T, row.names = F, sep = "\t")



DF1 = read.delim("BH_correction_0.95_conf_0_ALL.txt", header = T)
  DF2 = read.delim("Diff_0_ARACNe_RANDOM.txt", header = T)
    DF2 = DF2[,-1]
      colnames(DF2)= as.character(paste("Diff", colnames(DF2), sep=""))
        DFF = cbind(DF1, DF2)


dat.c <- melt(DFF, 
              id.var='Random', 
              measure.var=grep('ARACNeNet', names(DFF), value=TRUE),
              variable.name='aRACNeNet',
              value.name='aRACNeNet.val')
dat.c$idx <- gsub('ARACNeNet', '', dat.c$aRACNeNet)
dat.s <- melt(DFF, 
              id.var='Random', 
              measure.var=grep('DiffARACNeNet', names(DFF), value=TRUE),
              variable.name='diffARACNeNet',
              value.name='diffARACNeNet.val')
dat.s$idx <- gsub('DiffARACNeNet', '', dat.s$diffARACNeNet)
dat <- merge(dat.c, dat.s)

out <- dat[dat$aRACNeNet.val < 0.1, ]
  out.c <- dcast(out, Random ~ aRACNeNet, value.var='aRACNeNet.val')
    out.s <- dcast(out, Random ~ diffARACNeNet, value.var='diffARACNeNet.val')
      cbind_ARACNe_RANDOM_0_ALL_Freq= merge(out.c, out.s)

write.table(cbind_ARACNe_RANDOM_0_ALL_Freq, "cbind_ARACNe_RANDOM_0_ALL_Freq_0.1.txt", col.names = T, row.names = F, sep = "\t")



tmp1 = read.delim("cbind_ARACNe_RANDOM_0_ALL_Freq_0.1.txt", header = T) 
  tmp1 = as.data.frame(t(tmp1))
   cln = as.data.frame(rownames(tmp1))
    cb1 = cbind(cln, tmp1)
     cln =sub("Random", "Net", cb1[,1]) 
      cb1 = cb1[,-1]
       cb1 = cbind(cln, tmp1)
        Ar_Net_names = as.matrix(cb1[,1])
        Ar_Net_names = as.data.frame(Ar_Net_names[2:49])
         colnames(Ar_Net_names) = c("Net")

         cb1 = cb1[,-1]
          Rnd_Net_names = as.character(as.matrix(cb1[1,]))
           cb1 = cb1[-1,]
            cb1 <- as.data.frame(lapply(cb1, function(x) {as.numeric(as.character(x))}))

below= c("0")
  below = as.matrix(as.data.frame(rep(below, 96)))
      below = as.data.frame(as.numeric(as.character(below)))
        colnames(below) = c("below")

up= c("0")
  up = as.matrix(as.data.frame(rep(up, 96)))
     up = as.data.frame(as.numeric(as.character(up)))
        colnames(up) = c("up")

colnames(cb1) = Rnd_Net_names
cb = cbind(Ar_Net_names, below, up, cb1)


x = within(cb, {
below = rowSums(cb[-c(1:3)] < cb$below, na.rm = T)})

y = within(cb, {
up = rowSums(cb[-c(1:3)] > cb$up, na.rm = T)})


write.table(x, "BH_correction_0.95_conf_0_ALL_0.1_below_counting.txt", col.names = T, row.names = F, sep = "\t")
write.table(y, "BH_correction_0.95_conf_0_ALL_0.1_up_counting.txt", col.names = T, row.names = F, sep = "\t")














#final_table
final_tab_below = x[, c(1:2)]
final_tab1_below = as.data.frame(final_tab_below[1:48,])
final_tab2_below = as.data.frame(final_tab_below[49:96,])
final_tab_below = cbind(final_tab1_below, final_tab2_below)
colnames(final_tab_below) = c("ARACNe_Net", "Not_sig", "ARACNe_Net",  "Neg_freq")



final_tab2_below = final_tab2_below[, -1]
cb1 = cbind(final_tab1_below, final_tab2_below)
colnames(cb1) = c("ARACNe_Net", "Not_sig", "Neg_freq")


final_tab_up = y[, c(1:2)]
final_tab1_up = as.data.frame(final_tab_up[1:48,])
final_tab2_up = as.data.frame(final_tab_up[49:96,])
final_tab2_up = final_tab2_up[, -1]
cb2 = cbind(final_tab1_up, final_tab2_up)
colnames(cb2) = c("ARACNe_Net", "Sig", "Pos_freq")

cb_final = cbind(cb1, cb2)
row.names(cb_final) = cb_final[,1]
cb_final = cb_final[, -4]


#plot
tmp1 = read.delim("BH_correction_0.95.txt", header = T)
row.names(tmp1) = tmp1[,1]
tmp1 = tmp1[,-1]
den<-apply(tmp1, 2, density)
pdf("p_val_BH_0.95_corrected_plots.pdf")
for(i in 1:length(den)){
  plot(den[[i]], 
       main=paste('density ', i), 
       ylim=c(0, 800), xlim =c(0, 1),  
       xlab = "BH corrected p val")
}

dev.off()












