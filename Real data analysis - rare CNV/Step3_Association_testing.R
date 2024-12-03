# Chi-square test was used to test the association between treatment and CNV status(whether a CNV exist)
two.prop.test<-function(res,
                        CNVR.start,
                        CNVR.end){
  
  samples<-names(table(res$Sample_ID))
  # length(samples)
  n.treated<-length(which(substr(samples,1,1) %in% c("t","T")))
  n.untreated<-length(which(substr(samples,1,1) %in% c("u","U")))
  
  
  overlap.seq<-which(!(res$Start > CNVR.end | res$End < CNVR.start))
  
  tmp.res<-res[overlap.seq,]
  # dulication proportion
  dup.prop<-sum(tmp.res$Event=="Duplication" | tmp.res$Event=="Gain")/nrow(tmp.res)
  
  # duplications in each CNVR
  n.dup.treated.sample <- tmp.res$Sample_ID[which((tmp.res$Event=="Duplication" | tmp.res$Event=="Gain") 
                                                  & substr(tmp.res$Sample_ID,1,1) %in% c("t","T"))]
  n.dup.treated <- length(names(table(n.dup.treated.sample)))
  
  n.dup.untreated.sample <- tmp.res$Sample_ID[which((tmp.res$Event=="Duplication" | tmp.res$Event=="Gain") 
                                                    & substr(tmp.res$Sample_ID,1,1) %in% c("u","U"))]
  n.dup.untreated <- length(names(table(n.dup.untreated.sample)))
  
  # deletions in each CNVR
  n.del.treated.sample <- tmp.res$Sample_ID[which((tmp.res$Event=="Loss") 
                                                  & substr(tmp.res$Sample_ID,1,1) %in% c("t","T"))]
  n.del.treated <- length(names(table(n.del.treated.sample)))
  
  n.del.untreated.sample <- tmp.res$Sample_ID[which((tmp.res$Event=="Loss") 
                                                    & substr(tmp.res$Sample_ID,1,1) %in% c("u","U"))]
  n.del.untreated <- length(names(table(n.del.untreated.sample)))
  
  
  
  tmp.sample<-names(table(tmp.res$Sample_ID))
  A<-length(which(substr(tmp.sample,1,1) %in% c("t","T")))
  B<-n.treated-A
  
  C<-length(which(substr(tmp.sample,1,1) %in% c("u","U")))
  D<-n.untreated-C
  
  tab <- matrix(c(A, B, C, D), ncol = 2, byrow = T)
  test.res<-prop.test(tab)
  
  
  
  tmp.data<-data.frame(start = CNVR.start,
                       end   = CNVR.end,
                       width = CNVR.end-CNVR.start,
                       n.dup.treated   = n.dup.treated,
                       n.dup.untreated = n.dup.untreated,
                       n.del.treated   = n.del.treated,
                       n.del.untreated = n.del.untreated,
                       n.CNV.treated   = A,
                       n.CNV.untreated = C,
                       N.total         = A+C,
                       treated.prop = test.res$estimate[1],
                       untreated.prop = test.res$estimate[2],
                       prop.diff      = paste0(round((test.res$estimate[1]-test.res$estimate[2]),2), " (",
                                               round(test.res$conf.int[1],2),", ",
                                               round(test.res$conf.int[2],2),")"),
                       pval=test.res$p.value)
  
  return(test.res=tmp.data)
}


# ======================================================================================== #
res<-readRDS("CBS.res.t1.RDS")

# step 1: construct CNVR
# use RO method by CNVRuler

# step 2: for each CNVR: (1) state the consist of deletion or duplication, 
# (2) perform a Chi-square test to find the CNVRs that are differently distributed in two groups (treated/ untreated cells).
#                        
#                  yes            no
# treated           A              B
# untreated         C              D
CNVRs<-read.delim("CNVR.txt",row.names = 1)
CNVRs$Type=rownames(CNVRs)
colnames(CNVRs)<-colnames(CNVRs)[c(2,3,4,5,1)]

chisq.data<-data.frame()
for (i in 1:nrow(CNVRs)){
  print(i)
  tmp.data<-two.prop.test(res=res,CNVR.start=CNVRs$Start[i], CNVR.end=CNVRs$End[i])
  chisq.data<-rbind(chisq.data,tmp.data)
}


# add position info
norl.res<-get(load(file="norl.res.DC.map.Rdata"))
pos.data<-norl.res$RC_norm[,c(1,2,3)]

pos.data2<-data.frame()
for (i in 1:nrow(chisq.data)){
  print(i)
  tmp.pos.data<-data.frame(chr       = pos.data$seqnames[CNVRs$Start[i]],
                           raw.start = pos.data$start[CNVRs$Start[i]],
                           raw.end   = pos.data$end[CNVRs$End[i]])
  pos.data2<-rbind(pos.data2,tmp.pos.data)
}


chisq.data<-cbind(pos.data2,
                  data.frame(Type=CNVRs$Type,CNVR.ID=CNVRs$CNVR.ID),
                  chisq.data)
rownames(chisq.data)= chisq.data$CNVR.ID

chisq.data$treated.prop<-format(chisq.data$treated.prop,digits=3)
chisq.data$untreated.prop<-format(chisq.data$untreated.prop,digits=3)
chisq.data$fort.pval<-format(chisq.data$pval, scientific = TRUE,digits = 3)

chisq.data<-chisq.data[,c(1:4,13:17,19,18,5:12)]
chisq.data<-chisq.data[order(chisq.data$pval),]
write.csv(chisq.data,"chisq.data.res.csv",row.names = FALSE)
