#Script pour comparer les outputs de scrm et ms

#Entrer ces vecteurs de valeurs pour les scénarios
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
#n-islands model
x=1
for (i in 1:length(long_locusvalue)){
        for (j in 1:length(Mvalue)){
                for (k in 1:length(nvalue)){
                        for (l in simul){
                                ms[[x]] <- as.data.frame(read.table(paste0("~/Documents/Mes Documents/Cours et Erasmus/MASTER/M2RBEE/Stage/Stage.data/Test/sim_indep_10000frag_size",long_locusvalue[i],"_n=",nvalue[k],"_M=",Mvalue[j],"/mssim_indep_histo_10000frag_size",long_locusvalue[i],"_n=",nvalue[k],"_M=",Mvalue[j],"simul_",l,".txt"), quote="\"", stringsAsFactors=FALSE))[,2]
                                scrm[[x]] <- as.data.frame(read.table(paste0("~/Documents/Mes Documents/Cours et Erasmus/MASTER/M2RBEE/Stage/Stage.data/Test/sim_indep_10000frag_size",long_locusvalue[i],"_n=",nvalue[k],"_M=",Mvalue[j],"/scrmsim_indep_histo_10000frag_size",long_locusvalue[i],"_n=",nvalue[k],"_M=",Mvalue[j],"simul_",l,".txt"), quote="\"", stringsAsFactors=FALSE))[,2]                                  
                                if (length(ms[[x]])<length(scrm[[x]])){
                                        b=length(scrm[[x]])-length(ms[[x]])  #equalize length of our vector
                                        ms[[x]]<-c(ms[[x]],rep(0,b))
                                }
                                else if (length(ms[[x]])>length(scrm[[x]])){
                                        b=length(ms[[x]])-length(scrm[[x]])
                                        scrm[[x]]<-c(scrm[[x]],rep(0,b))
                                }
                                if (length(ms[[x]]>10)) {                    #Second : group when mutation's count is above 10
                                        scrm[[x]] <- c(scrm[[x]][1:10],sum(scrm[[x]][11:length(scrm[[x]])]))
                                        ms[[x]] <- c(ms[[x]][1:10],sum(ms[[x]][11:length(ms[[x]])]))      
                                }
                                
                                COMP<-chisq.test(scrm[[x]],ms[[x]])   #Chisq test
                                Results[l+1,x] <- COMP[[3]]
                                l=l+1
                        }  
                        x=x+1
                }                                        
        }
}
#SSPSC
x=1
for (i in 1:length(long_locusvalue)){
        for (j in 1:length(alphavalue)){
                for (k in 1:length(Tvalue)){
                        for (l in simul){
                                ms1[[x]] <- as.data.frame(read.table(paste0("~/Documents/Mes Documents/Cours et Erasmus/MASTER/M2RBEE/Stage/Stage.data/Test/sim_indep_10000frag_size",long_locusvalue[i],"_T=",Tvalue[k],"_alpha=",alphavalue[j],"/mssim_indep_histo_10000frag_size",long_locusvalue[i],"_T=",Tvalue[k],"_alpha=",alphavalue[j],"simul_",l,".txt"), quote="\"", stringsAsFactors=FALSE))[,2]
                                scrm1[[x]] <- as.data.frame(read.table(paste0("~/Documents/Mes Documents/Cours et Erasmus/MASTER/M2RBEE/Stage/Stage.data/Test/sim_indep_10000frag_size",long_locusvalue[i],"_T=",Tvalue[k],"_alpha=",alphavalue[j],"/scrmsim_indep_histo_10000frag_size",long_locusvalue[i],"_T=",Tvalue[k],"_alpha=",alphavalue[j],"simul_",l,".txt"), quote="\"", stringsAsFactors=FALSE))[,2]                                  
                                if (length(ms1[[x]])<length(scrm1[[x]])){
                                        b=length(scrm1[[x]])-length(ms1[[x]])
                                        ms1[[x]]<-c(ms1[[x]],rep(0,b))
                                }
                                else if (length(ms1[[x]])>length(scrm1[[x]])){
                                        b=length(ms1[[x]])-length(scrm1[[x]])
                                        scrm1[[x]]<-c(scrm1[[x]],rep(0,b))
                                }
                                if (length(ms1[[x]]>10)) {
                                        scrm1[[x]] <- c(scrm1[[x]][1:10],sum(scrm1[[x]][11:length(scrm1[[x]])]))
                                        ms1[[x]] <- c(ms1[[x]][1:10],sum(ms1[[x]][11:length(ms1[[x]])]))      
                                }
                                
                                COMP<-chisq.test(scrm1[[x]],ms1[[x]])
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
write.table(Final, "~/Documents/Mes Documents/Cours et Erasmus/MASTER/M2RBEE/Stage/Stage.data/Test/ChisqGroupe.csv", quote = FALSE, sep=";",row.names = TRUE,col.names=NA)

#Ploter certains histogrammes.

Testouille2 <- as.data.frame(read.table(paste0("~/Documents/Mes Documents/Cours et Erasmus/MASTER/M2RBEE/Stage/Stage.data/Test/sim_indep_10000frag_size",long_locusvalue[1],"_n=",nvalue[3],"_M=",Mvalue[2],"/scrmsim_indep_histo_10000frag_size",long_locusvalue[1],"_n=",nvalue[3],"_M=",Mvalue[2],"simul_",3,".txt"), quote="\"", stringsAsFactors=FALSE))[,2]                                  
Testouille2
Testouille3 <- c(Testouille2[1:10],sum(Testouille2[11:length(Testouille2)]))
Testouille3
barplot(Testouille3,ylim=c(0,10000),names.arg=c(1:10,"11+"),xlab="Mutations count",ylab="Freq",main="Locus length=100, M=0.1, n=100",las=1)
