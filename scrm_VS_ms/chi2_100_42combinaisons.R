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
Results <- matrix(ncol=42,nrow=100,NA)
Results
colnames=c(Tvalue)
x=1
for (i in 1:length(long_locusvalue)){
        for (j in 1:length(Mvalue)){
                for (k in 1:length(nvalue)){
                        for (l in simul){
                                
                                ms[[x]] <- as.data.frame(read.table(paste0("~/Documents/Mes Documents/Cours et Erasmus/MASTER/M2RBEE/Stage/Stage.data/Test/sim_indep_10000frag_size",long_locusvalue[i],"_n=",nvalue[k],"_M=",Mvalue[j],"/mssim_indep_histo_10000frag_size",long_locusvalue[i],"_n=",nvalue[k],"_M=",Mvalue[j],"simul_",l,".txt"), quote="\"", stringsAsFactors=FALSE))[,2]
                                scrm[[x]] <- as.data.frame(read.table(paste0("~/Documents/Mes Documents/Cours et Erasmus/MASTER/M2RBEE/Stage/Stage.data/Test/sim_indep_10000frag_size",long_locusvalue[i],"_n=",nvalue[k],"_M=",Mvalue[j],"/scrmsim_indep_histo_10000frag_size",long_locusvalue[i],"_n=",nvalue[k],"_M=",Mvalue[j],"simul_",l,".txt"), quote="\"", stringsAsFactors=FALSE))[,2]                                  
                                if (length(ms[[x]])<length(scrm[[x]])){
                                        b=length(scrm[[x]])-length(ms[[x]])
                                        ms[[x]]<-c(ms[[x]],rep(0,b))
                                        COMP<-chisq.test(ms[[x]],scrm[[x]])
                                }
                                else if (length(ms[[x]])>length(scrm[[x]])){
                                        b=length(ms[[x]])-length(scrm[[x]])
                                        scrm[[x]]<-c(scrm[[x]],rep(0,b))
                                        COMP<-chisq.test(ms[[x]],scrm[[x]])
                                }
                                else {
                                        COMP<-chisq.test(scrm[[x]],ms[[x]])
                                }
                                Results[l+1,x] <- COMP[[3]]
                                l=l+1
                        }  
                        x=x+1
                }                                        
        }
}
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
                                        COMP<-chisq.test(ms1[[x]],scrm1[[x]])
                                }
                                else if (length(ms1[[x]])>length(scrm1[[x]])){
                                        b=length(ms1[[x]])-length(scrm1[[x]])
                                        scrm1[[x]]<-c(scrm1[[x]],rep(0,b))
                                        COMP<-chisq.test(ms1[[x]],scrm1[[x]])
                                }
                                else {
                                        COMP<-chisq.test(scrm1[[x]],ms1[[x]])
                                }
                                Results[l+1,x+18] <- COMP[[3]]
                                l=l+1
                        }  
                        x=x+1
                }                                        
        }
}

sampAll <- expand.grid(Tvalue,alphavalue,long_locusvalue,KEEP.OUT.ATTRS = FALSE)
sampAll2<-expand.grid(Mvalue,nvalue,long_locusvalue,KEEP.OUT.ATTRS = FALSE)

names = c()
for (i in 1:24){
        names[i]<-paste0("T=",sampAll[i,1],"_alpha=",sampAll[i,2],"_fraglength=",sampAll[i,3])
}
for (i in 1:18){
        names[i+24]<-paste0("M=",sampAll2[i,1],"_n=",sampAll2[i,2],"_fraglength=",sampAll2[i,3])
}
colnames(Results) <- names
Results


Results
V <- Results[,2]<0.05
sum(V)

?read.data





