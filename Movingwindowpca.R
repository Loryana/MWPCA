#MWPCA -----

rm(list=ls())

window_pca <- function(H,datagen_table=data.frame(datagen)){
  comp=1
  tcarre_val=c()
  sspe_res=c()
  #faire un test de validation!!!!!!!
  datagen_table_val=datagen_table[sample(rownames(datagen_table),20),]
  
  for(i in 1:(26-H)){
    choice=c()
    for(j in comp:(comp+H)){
      choice=c(choice,which(date_gen==j))
    }
    datagen_table_train=datagen_table[choice,]
    datagen_table_test=datagen_table[max(choice)+1,]
    res_rnirs=pca(datagen_table_train,datagen_table_test[1,],ncomp=10)
    pred_testn=res_rnirs$Tu
    
    tcarre=c()
    
    for(w in 1:nrow(datagen_table_test)){
      tcarre=c(tcarre,t(t(as.vector(datagen_table_test[w,])))%*%
                 res_rnirs$P%*%
                 pinv(t(res_rnirs$Tr)%*%res_rnirs$Tr)$Xplus%*%
                 t(res_rnirs$P)%*%
                 t(as.vector(datagen_table_test[w,])))
      
      
      
      rstat=t(t(as.vector(datagen_table_test[w,])))-t(t(as.vector(datagen_table_test[w,])))%*%res_rnirs$P%*%t(res_rnirs$P)
      
      qstat=t(rstat)%*%rstat
    }
    
    
    sspe=c()
    
    #utilisation du test de validation ici
    for(k in 1:nrow(datagen_table_val)){
      #sspe=t(t(as.vector(datagen_table_val[1,])))-t(t(as.vector(datagen_table_val[1,])))%*%res_rnirs$P%*%t(res_rnirs$P)
      
      sspe=c(sspe,abs(t(t(as.vector(datagen_table_val[k,])))-t(t(as.vector(datagen_table_val[k,])))%*%res_rnirs$P%*%t(res_rnirs$P))^2)
      
    }
    
    sspe_res=c(sspe_res,sum(sspe))
    
    tcarre_val=c(tcarre_val,mean(tcarre))
    
    comp=comp+1
  }
  res=list(tcarre_val,qstat,sspe_res)
  return(res)
}