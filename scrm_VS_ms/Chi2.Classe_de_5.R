#Comparaison sorties de ms et scrm


Tvalue = c(0.1, "1.0", "2.0")
alphavalue = c(0.1, 2, 10, 50)
Mvalue = c(1, 0.1, 50)
nvalue = c(10, 2, 100)
long_locusvalue = c(100, 200)
simul = c(0:99)
ms=list()
ms1=list()
scrm=list()
scrm1=list()
COMP = list()
COMP1 =list()
Results <- matrix(ncol=length(Tvalue)*length(alphavalue)*length(long_locusvalue)+length(nvalue)*length(Mvalue)*length(long_locusvalue),nrow=length(simul),NA)
#n-island models
x=1
for (g in 1:length(long_locusvalue)){
        for (h in 1:length(Mvalue)){
                for (k in 1:length(nvalue)){
                        for (l in simul){
                                ms[[x]] <- as.data.frame(read.table(paste0("~/Documents/Mes Documents/Cours et Erasmus/MASTER/M2RBEE/Stage/Stage.data/Test.comparaisonMS.SCRM/sim_indep_10000frag_size",long_locusvalue[g],"_n=",nvalue[k],"_M=",Mvalue[h],"/mssim_indep_histo_10000frag_size",long_locusvalue[g],"_n=",nvalue[k],"_M=",Mvalue[h],"simul_",l,".txt"), quote="\"", stringsAsFactors=FALSE))[,2]
                                scrm[[x]] <- as.data.frame(read.table(paste0("~/Documents/Mes Documents/Cours et Erasmus/MASTER/M2RBEE/Stage/Stage.data/Test.comparaisonMS.SCRM/sim_indep_10000frag_size",long_locusvalue[g],"_n=",nvalue[k],"_M=",Mvalue[h],"/scrmsim_indep_histo_10000frag_size",long_locusvalue[g],"_n=",nvalue[k],"_M=",Mvalue[h],"simul_",l,".txt"), quote="\"", stringsAsFactors=FALSE))[,2]
                                j = 1
                                if (length(ms[[x]])<length(scrm[[x]])){
                                        b=length(scrm[[x]])-length(ms[[x]])  #equalize length of our vector
                                        ms[[x]]<-c(ms[[x]],rep(0,b))
                                }
                                else if (length(ms[[x]])>length(scrm[[x]])){
                                        b=length(ms[[x]])-length(scrm[[x]])
                                        scrm[[x]]<-c(scrm[[x]],rep(0,b))
                                }

                                for (j in 1:length(ms[[x]])){                
                                        if ((ms[[x]][j]>=10) & (scrm[[x]][j]>=10) & (j+1<=length(ms[[x]])) & (j+1<=length(scrm[[x]]))) {
                                                ms[[x]][j] <- ms[[x]][j] 
                                                scrm[[x]][j] <- scrm[[x]][j] 
                                                }
                                        else {
                                                while ((ms[[x]][j]<10 | scrm[[x]][j]<10) & j+1<=length(ms[[x]]) & j+1<=length(scrm[[x]]))   {
                                                        i=1
                                                        ms[[x]][j] <- ms[[x]][j]+ms[[x]][i+j] 
                                                        ms[[x]] <- ms[[x]][-(i+j)]
                                                        scrm[[x]][j] <- scrm[[x]][j]+scrm[[x]][i+j] 
                                                        scrm[[x]] <- scrm[[x]][-(i+j)]
                                                        i=i+1
                                                        }
                                                }
                                        j=j+1
                                }
                        COMP<-chisq.test(scrm[[x]],ms[[x]])   #Chisq test
                        Results[l+1,x] <- COMP[[3]]
                        l=l+1
                        }
                x=x+1
                }                                        
        }
}



x=1
for (g in 1:length(long_locusvalue)){
        for (h in 1:length(alphavalue)){
                for (k in 1:length(Tvalue)){
                        for (l in simul){
                                ms1[[x]] <- as.data.frame(read.table(paste0("~/Documents/Mes Documents/Cours et Erasmus/MASTER/M2RBEE/Stage/Stage.data/Test.comparaisonMS.SCRM/sim_indep_10000frag_size",long_locusvalue[g],"_T=",Tvalue[k],"_alpha=",alphavalue[h],"/mssim_indep_histo_10000frag_size",long_locusvalue[g],"_T=",Tvalue[k],"_alpha=",alphavalue[h],"simul_",l,".txt"), quote="\"", stringsAsFactors=FALSE))[,2]
                                scrm1[[x]] <- as.data.frame(read.table(paste0("~/Documents/Mes Documents/Cours et Erasmus/MASTER/M2RBEE/Stage/Stage.data/Test.comparaisonMS.SCRM/sim_indep_10000frag_size",long_locusvalue[g],"_T=",Tvalue[k],"_alpha=",alphavalue[h],"/scrmsim_indep_histo_10000frag_size",long_locusvalue[g],"_T=",Tvalue[k],"_alpha=",alphavalue[h],"simul_",l,".txt"), quote="\"", stringsAsFactors=FALSE))[,2]                                  
                                j = 1
                                if (length(ms1[[x]])<length(scrm1[[x]])){
                                        b=length(scrm1[[x]])-length(ms1[[x]])  #equalize length of our vector
                                        ms1[[x]]<-c(ms1[[x]],rep(0,b))
                                }
                                else if (length(ms1[[x]])>length(scrm1[[x]])){
                                        b=length(ms1[[x]])-length(scrm1[[x]])
                                        scrm1[[x]]<-c(scrm1[[x]],rep(0,b))
                                }
                                
                                for (j in 1:length(ms1[[x]])){                
                                        if ((ms1[[x]][j]>=10) & (scrm1[[x]][j]>=10) & (j+1<=length(ms1[[x]])) & (j+1<=length(scrm1[[x]]))) {
                                                ms1[[x]][j] <- ms1[[x]][j] 
                                                scrm1[[x]][j] <- scrm1[[x]][j] 
                                        }
                                        else {
                                                while ((ms1[[x]][j]<10 | scrm1[[x]][j]<10) & j+1<=length(ms1[[x]]) & j+1<=length(scrm1[[x]]))   {
                                                        i=1
                                                        ms1[[x]][j] <- ms1[[x]][j]+ms1[[x]][i+j] 
                                                        ms1[[x]] <- ms1[[x]][-(i+j)]
                                                        scrm1[[x]][j] <- scrm1[[x]][j]+scrm1[[x]][i+j] 
                                                        scrm1[[x]] <- scrm1[[x]][-(i+j)]
                                                        i=i+1
                                                }
                                        }
                                        j=j+1
                                }
                                COMP<-chisq.test(scrm1[[x]],ms1[[x]])   #Chisq test
                                Results[l+1,x+18] <- COMP[[3]]
                                l=l+1
                        }
                        x=x+1
                }                                        
        }
}

Results  
sampAll <- expand.grid(Tvalue,alphavalue,long_locusvalue,KEEP.OUT.ATTRS = FALSE)
sampAll2<-expand.grid(nvalue,Mvalue,long_locusvalue,KEEP.OUT.ATTRS = FALSE)
names = c()
for (i in 1:c(length(nvalue)*length(Mvalue)*length(long_locusvalue))){
        names[i]<-paste0("Locus length=",sampAll2[i,3],", M=",sampAll2[i,2],", n=",sampAll2[i,1])      
}
for (i in 1:c(length(Tvalue)*length(alphavalue)*length(long_locusvalue))){
        names[i+length(nvalue)*length(Mvalue)*length(long_locusvalue)]<-paste0("Locus length=",sampAll[i,3],", T=",sampAll[i,1],", alpha=",sampAll[i,2])
}

colnames(Results) <- names
Counts <- Results<0.05
Counts=apply(Counts,2,sum)
Mean_stats=apply(Results,2,mean)
Final=cbind(as.matrix(Mean_stats),as.matrix(Counts))
colnames(Final)<-c("Chisq mean P-value","Reject H0 counts")
View(Final)
#hyp nulle homo = variable aleatoire suivent la même loi
#if < 0.05 => rejet de hyp nulle et donc variable suivent loi différentes.
write.table(Final, "~/Documents/Mes Documents/Cours et Erasmus/MASTER/M2RBEE/Stage/Stage.data/Test.comparaisonMS.SCRM/Chisq_10classes.csv", quote = FALSE, sep=";",row.names = TRUE,col.names=NA)                                
     

####Test pour scénario extrême


Mvalue = c(1, 0.1, 50)
nvalue = c(10, 2, 100)
long_locusvalue = c(100, 200)
ms=list()
scrm=list()
ms[[x]] <- as.data.frame(read.table(paste0("~/Documents/Mes Documents/Cours et Erasmus/MASTER/M2RBEE/Stage/Stage.data/Test.comparaisonMS.SCRM/sim_indep_10000frag_size",long_locusvalue[1],"_n=",nvalue[3],"_M=",Mvalue[1],"/mssim_indep_histo_10000frag_size",long_locusvalue[1],"_n=",nvalue[3],"_M=",Mvalue[1],"simul_",67,".txt"), quote="\"", stringsAsFactors=FALSE))[,2]
scrm[[x]] <- as.data.frame(read.table(paste0("~/Documents/Mes Documents/Cours et Erasmus/MASTER/M2RBEE/Stage/Stage.data/Test.comparaisonMS.SCRM/sim_indep_10000frag_size",long_locusvalue[1],"_n=",nvalue[3],"_M=",Mvalue[1],"/scrmsim_indep_histo_10000frag_size",long_locusvalue[1],"_n=",nvalue[3],"_M=",Mvalue[1],"simul_",67,".txt"), quote="\"", stringsAsFactors=FALSE))[,2]
x = 1
j = 1
if (length(ms[[x]])<length(scrm[[x]])){
        b=length(scrm[[x]])-length(ms[[x]])  #equalize length of our vector
        ms[[x]]<-c(ms[[x]],rep(0,b))
}
if (length(ms[[x]])>length(scrm[[x]])){
        b=length(ms[[x]])-length(scrm[[x]])
        scrm[[x]]<-c(scrm[[x]],rep(0,b))
}

for (j in 1:length(ms[[x]])){                
        if ((ms[[x]][j]>=23) & (scrm[[x]][j]>=23) & (j+1<=length(ms[[x]])) & (j+1<=length(scrm[[x]]))) {
                ms[[x]][j] <- ms[[x]][j] 
                scrm[[x]][j] <- scrm[[x]][j] 
        }
        else {
                while ((ms[[x]][j]<23 | scrm[[x]][j]<23) & j+1<=length(ms[[x]]) & j+1<=length(scrm[[x]]))   {
                        i=1
                        ms[[x]][j] <- ms[[x]][j]+ms[[x]][i+j] 
                        ms[[x]] <- ms[[x]][-(i+j)]
                        scrm[[x]][j] <- scrm[[x]][j]+scrm[[x]][i+j] 
                        scrm[[x]] <- scrm[[x]][-(i+j)]
                        i=i+1
                }
        }
        j=j+1
}

ms
scrm
chisq.test(ms[[1]],scrm[[1]],rescale.p=TRUE)
chisq.test(ms[[1]],simulate.p.value = TRUE,B=100)
?chisq.test

barplot(ms[[1]],xlab="Mutations count",ylab="Freq",main=paste0("ms ",names[x]),las=1)                          
barplot(scrm[[1]],xlab="Mutations count",ylab="Freq",main=paste0("ms ",names[x]),las=1)                          



