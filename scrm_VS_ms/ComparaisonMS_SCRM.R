Tvalue = c(0.1, "1.0", "2.0")
alphavalue = c(0.1, 2, 10, 50)
Mvalue = c(1, 0.1, 50)
nvalue = c(10, 2, 100)
long_locusvalue = c(100, 200)

ms=list()
ms1=list()
scrm=list()
scrm1=list()
COMP = list()
COMP1 =list()
x=1
for (i in 1:length(long_locusvalue)){
        for (j in 1:length(Mvalue)){
                for (k in 1:length(nvalue)){
                                        ms[[x]] <- as.data.frame(read.table(paste0("~/Documents/Mes Documents/Cours et Erasmus/MASTER/M2RBEE/Stage/Stage.data/Test/sim_indep_10000frag_size",long_locusvalue[i],"_n=",nvalue[k],"_M=",Mvalue[j],"/mssim_indep_histo_10000frag_size",long_locusvalue[i],"_n=",nvalue[k],"_M=",Mvalue[j],".txt"), quote="\"", stringsAsFactors=FALSE))
                                        scrm[[x]] <- as.data.frame(read.table(paste0("~/Documents/Mes Documents/Cours et Erasmus/MASTER/M2RBEE/Stage/Stage.data/Test/sim_indep_10000frag_size",long_locusvalue[i],"_n=",nvalue[k],"_M=",Mvalue[j],"/scrmsim_indep_histo_10000frag_size",long_locusvalue[i],"_n=",nvalue[k],"_M=",Mvalue[j],".txt"), quote="\"", stringsAsFactors=FALSE))                                  
                                        COMP[[x]]<-cbind(c(ms[[x]][2],scrm[[x]][2]))
                                        x=x+1
                }                                        
        }
}
x=1
for (i in 1:length(long_locusvalue)){
        for (j in 1:length(alphavalue)){
                for (k in 1:length(Tvalue)){
                        ms1[[x]] <- as.data.frame(read.table(paste0("~/Documents/Mes Documents/Cours et Erasmus/MASTER/M2RBEE/Stage/Stage.data/Test/sim_indep_10000frag_size",long_locusvalue[i],"_T=",Tvalue[k],"_alpha=",alphavalue[j],"/mssim_indep_histo_10000frag_size",long_locusvalue[i],"_T=",Tvalue[k],"_alpha=",alphavalue[j],".txt"), quote="\"", stringsAsFactors=FALSE))
                        scrm1[[x]] <- as.data.frame(read.table(paste0("~/Documents/Mes Documents/Cours et Erasmus/MASTER/M2RBEE/Stage/Stage.data/Test/sim_indep_10000frag_size",long_locusvalue[i],"_T=",Tvalue[k],"_alpha=",alphavalue[j],"/scrmsim_indep_histo_10000frag_size",long_locusvalue[i],"_T=",Tvalue[k],"_alpha=",alphavalue[j],".txt"), quote="\"", stringsAsFactors=FALSE))                                  
                        COMP1[[x]]<-cbind(c(ms1[[x]][2],scrm1[[x]][2]))
                        x=x+1
                }                                        
        }
}

COMP
COMP1
ms1
C=list()
D=list()
Results = matrix(ncol=2,nrow=42)

for (i in 1:length(ms)){
        A<-as.vector(ms[[i]][2])
        B<-as.vector(scrm[[i]][2])
        if (length(A[[1]])<length(B[[1]])){
                C<-chisq.test(A[[1]],B[[1]][1:length(A[[1]])])
        }
        else if (length(B[[1]])<length(A[[1]])){
                C<-chisq.test(B[[1]],A[[1]][1:length(B[[1]])])
        }
        else {
                C<-chisq.test(B[[1]],A[[1]])
        }
        C<-unlist(C)
        Results[i,1] <- C[[1]]
        Results[i,2] <- C[[3]]
        i = i+1
}

Results

for (j in 1:length(ms1)){
        E<-as.vector(ms1[[j]][2])
        G<-as.vector(scrm1[[j]][2])
        if (length(E[[1]])<length(G[[1]])){
                D<-chisq.test(E[[1]],G[[1]][1:length(E[[1]])])
        }
        else if (length(G[[1]])<length(E[[1]])){
                D<-chisq.test(G[[1]],E[[1]][1:length(G[[1]])])
        }
        else {
                D<-chisq.test(G[[1]],E[[1]])
        }
        D<-unlist(D)
        Results[j+18,1] <- D[[1]]
        Results[j+18,2] <- D[[3]]
        j = j+1
}

Results
V <- Results[,2]<0.05
sum(V)


test = as.data.frame(read.table(paste0("~/Documents/Mes Documents/Cours et Erasmus/MASTER/M2RBEE/Stage/Stage.data/Test/sim_indep_10000frag_size",long_locusvalue[2],"_n=",nvalue[3],"_M=",Mvalue[2],"/mssim_indep_histo_10000frag_size",long_locusvalue[2],"_n=",nvalue[3],"_M=",Mvalue[2],".txt"), quote="\"", stringsAsFactors=FALSE))
test2 = test[,1]*test[,2]
test3 = sum(test2)
test3
test4 = test3/200
test4

?read.data
