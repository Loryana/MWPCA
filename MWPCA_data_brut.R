#MWPCA -----
rm(list=ls())
setwd("~/Thèse/Script")

source("Movingwindowpca_modif_save.R")
source("asd_read.R")

library(rnirs)
library(ggplot2)

#fonction pour récupérer le nom spécifiques (et date) des fichiers / dossiers
infofile <- function(x,sta,sto){substr(x,sta,sto)}

setwd("~/Thèse/Données blé")

filename=dir()

#fonction pour récupérer le nom spécifiques (et date) des fichiers / dossiers
infofile <- function(x,sta,sto){substr(x,sta,sto)}

#récupération des nom de dossier
filedate=sapply(filename,infofile,sta=1,sto=4)
#récupération des noms de dossier qui nous intéressent uniquement
filename=filename[which(filedate=="2011")]


load("Données.RData")

setwd("~/Thèse/Script")


for(i in 1:length(data)){
  print(which(sapply(rownames(data[[i]]),infofile,sta=3,sto=5)=="90_"))
}

valeurN=read.csv("valeursN_2011.csv",h=T,sep=';')
valeurN2=split(valeurN,valeurN$CodePlante)

# pour un génotype ----

# genotyp=c('C_136','C_141','C_165','C_167','C_258','C_264','C_315','C_346','C_39_','C_408','C_444','C_492','C_51_','C_577','C_589','C_668','C_720','C_77_','C_82_','C_90_',
#           'D_136','D_141','D_165','D_167','D_258','D_264','D_315','D_346','D_39_','D_408','D_444','D_492','D_51_','D_577','D_589','D_668','D_720','D_77_','D_82_','D_90_') #07/04
# gen2=c('C_136','C_141','C_165','C_167','C_258','C_264','C_315','C_346','C_39','C_408','C_444','C_492','C_51','C_577','C_589','C_668','C_720','C_77','C_82_','C_90',
#           'D_136','D_141','D_165','D_167','D_258','D_264','D_315','D_346','D_39','D_408','D_444','D_492','D_51','D_577','D_589','D_668','D_720','D_77','D_82_','D_90') #07/04
# setwd("~/Thèse/Script/Groupes CD") #pour l'enregistrement des figures

# 
# genotyp=c('A_264','A_316','A_321','A_437','A_526','A_539','A_674','A_90_',
#           'B_264','B_316','B_321','B_437','B_526','B_539','B_674','B_90_') #07/04
# gen2=c('A_264','A_316','A_321','A_437','A_526','A_539','A_674','A_90',
#        'B_264','B_316','B_321','B_437','B_526','B_539','B_674','B_90') #07/04


genotyp=c('A_674','B_90_')
gen2=c('A_674','B_90')



setwd("~/Thèse/Script/Groupes AB brut") #pour l'enregistrement des figures

# genotyp=c("A_90_")
# gen2=c('A_90')
# setwd("~/Thèse/Script/A90") #pour l'enregistrement des figures



res=list()
res2=vector("list",length(genotyp))
names(res2)=genotyp


for(g in 1:length(genotyp)){
  svMisc::progress(g/length(genotyp)*100)
  res[[g]]=matrix(NA,ncol=2,nrow=8)
  gen=genotyp[g]
  datagen=c() #matrice contenant les spectres pour le génotype 131
  for(i in 1:length(data)){
    row_num=c()
    for( k in 1:length(gen)){
      row_num=c(row_num,which(sapply(rownames(data[[i]]),infofile,sta=1,sto=5)==gen[k]))
      #print(rownames(data[[i]][row_num,]))
    }
    if(length(row_num)>1){
      datagen=rbind(datagen,data[[i]][row_num,])
    }
  }
  
  mm=valeurN2[[which(names(valeurN2)==gen2[g])]]
  date_plot=mm$degresjour
  
  # plot des spectres
  for(k in 1:nrow(datagen)){
    if(k==1){
      plot(colnames(datagen),datagen[k,],type='l',ylim=c(min(datagen),max(datagen)),xlab="Longueur d'ondes",ylab="")
    }else{
      lines(colnames(datagen),datagen[k,],col=k)
    }
  }
  # dev.print(device = png, file = paste(genotyp[g],"_spectresbruts.png",sep=''), width = 600)
  
  
  res2[[g]]=vector('list',8)
  names(res2[[g]])=c('w2','w3','w4','w5','w6','w7','w8','w9')
  
  for(h in 2:8){
    date_gen=sort(c(rep((1:length(date_plot)),2)))
    
    setwd("~/Thèse/Script/Groupes AB brut/MWPCA")
    
    oh=window_pca(H=h,genotyp=gen,datagen_table=data.frame(datagen),datagen_test=data.frame(datagen),nco=5,date_gen=date_gen,date_gen_test=date_gen,datemax=max(date_gen),rnirmeth=TRUE)
    
    setwd("~/Thèse/Script/Groupes AB brut")
    
    bidule=data.frame(date=mm$degresjour[-c(1:(h+1))],Tcarre=oh[[1]])
    machin=data.frame(Qstat=oh[[2]])
    
    hist(oh[[1]],breaks=15,main=paste(genotyp[g],": Taille de la fenêtre =",h + 1))
    dev.print(device = png, file = paste(genotyp[g],"_hist_T-w_",h+1,".png",sep=''), width = 600)
    hist(oh[[2]],breaks=15,main=paste(genotyp[g],": Taille de la fenêtre =",h + 1))
    dev.print(device = png, file = paste(genotyp[g],"_hist_Q_w_",h+1,".png",sep=''), width = 600)
    
    
    df <- cbind(bidule,machin)
    
    p <- ggplot()
    p <- p + geom_point(data=df,aes(x=date, y=Tcarre,color="#F8766D"),size=4)
    p <- p + geom_line(data=df,aes(x=date, y=Tcarre,color="#F8766D"))
    p <- p + geom_point(data=df,aes(x=date, y=Qstat,color="#00BFC4"),size=4)
    p <- p + geom_line(data=df,aes(x=date, y=Qstat,color="#00BFC4"))
    p <- p + ggtitle(paste(genotyp[g],": Taille de la fenêtre =",h + 1))
    p <- p + xlab("Degrés jour")
    p <- p + ylab("Statistiques")
    p <- p + scale_colour_manual(name = "Statistiques", values = c("#F8766D", "#00BFC4"),labels=c("T²","Q"))
    # p <- p + theme(legend.justification = c(0, 1),legend.position = c(0.03, 0.97),text = element_text(size = 20))
    p <- p + theme(plot.title = element_text(hjust = 0.5))
    p
    ggsave(paste(genotyp[g],"_w_",h+1,".png",sep=''))
    
    # plot(mm$degresjour[-c(1:(h+1))],oh[[1]],type='b',ylab="T²",xlab="Pas de temps",ylim=c(0,max(c(oh[[1]],oh[[2]]))),main=paste(genotyp[g],": Taille de la fenêtre =",h + 1))
    # lines(mm$degresjour[-c(1:(h+1))],oh[[2]],col='red',type='b')
    # dev.print(device = png, file = paste(genotyp[g],"_w_",h,".png",sep=''), width = 600)
    # boxplot(oh[[1]],title=paste(genotyp[g],": Taille de la fenêtre =",h + 1))
    # dev.print(device = png, file = paste(genotyp[g],"_box_T_",h,".png",sep=''), width = 600)
    # boxplot(oh[[2]],title=paste(genotyp[g],": Taille de la fenêtre =",h + 1))
    # dev.print(device = png, file = paste(genotyp[g],"_box_Q_",h,".png",sep=''), width = 600)
    res[[g]][h,1]=date_plot[which.max(oh[[1]])+h]
    res[[g]][h,2]=date_plot[which.max(oh[[2]])+h]
    supra=c()
    for(k in 1:length(boxplot(oh[[1]])$out)){
      supra=c(supra,date_plot[which(oh[[1]]==boxplot(oh[[1]])$out[k])+h])
    }
    mana=c()
    for(l in 1:length(boxplot(oh[[2]])$out)){
      mana=c(mana,date_plot[which(oh[[2]]==boxplot(oh[[2]])$out[l])+h])
    }
    res2[[g]][[h]]=vector('list',2)
    names(res2[[g]][[h]])=c('T','Q')
    res2[[g]][[h]][[1]]=supra
    res2[[g]][[h]][[2]]=mana
  }
}

save(res,file='GroupeAB_res_2gen.RData')
save(res2,file='GroupeAB_res2_2gen.RData')

# save(res,file='GroupeCD_res.RData')

# save(res,file='A90.RData')



setwd("~/Thèse/Script")

arg=10 # 7 : azote, 10: chloro, 13: NDVI, 22: eau

load(paste("Tommy",arg-1,".RData",sep='_')) #mm2011 azote

modele.ingrid <- function(a,b,c,d,x){  # calcul la valeur du point y pour la valeur x
  res <- d+ a/(1+b*exp(c*x))
  return(res)
}

# genotyp=c('A_264','A_316','A_321','A_437','A_526','A_539','A_674','A_90',
#           'B_264','B_316','B_321','B_437','B_526','B_539','B_674','B_90') #07/04 ####"problème des 90_ pour valeurN2


genotyp=c('A_674','B_90')

resupoint=vector("list",length(genotyp))
resupoint2=vector("list",length(genotyp))
nom_gen=c()



for(g in 1:length(genotyp)){
  troup=c()
  soldier=c()
  for(i in 1:8){
    troup=c(troup,res2[[g]][[i]][[1]])
    soldier=c(soldier,res2[[g]][[i]][[2]])
  }
  
  
  if(any(which(summary(as.factor(troup))>2))){
    saveme=summary(as.factor(troup))[which(summary(as.factor(troup))>2)]
  }
  if(any(which(summary(as.factor(soldier))>2))){
    savemeleretour=summary(as.factor(soldier))[which(summary(as.factor(soldier))>2)]
  }
  
  gen=genotyp[g]
  date_num=valeurN2[[which(names(valeurN2)==gen)]]$date_bis
  
  data_gen=valeurN2[[which(names(valeurN2)==gen)]]
  # candy=c()
  # for(i in 1:length(date_num)){
  #   crush=which(date_num[i]==data_gen$date_bis)
  #   candy=c(candy,crush)
  # }
  
  #plot(strptime(date_num,"%d/%m/%Y"),data_gen$azmoy2011,type='b',xlab='',ylab='',main=paste(gen,' azmoy'))
  p <- ggplot()
  if(arg==7){
    p <- p + geom_point(aes(x=valeurN$degresjour[which(valeurN$CodePlante==gen)], y=data_gen$azmoy2011),color="#F8766D",size=4)
    p <- p + geom_line(aes(x=valeurN$degresjour[which(valeurN$CodePlante==gen)], y=data_gen$azmoy2011),color="#F8766D")
    p <- p + ylab("Teneur en Azote")
  }else if(arg==10){
    p <- p + geom_point(aes(x=valeurN$degresjour[which(valeurN$CodePlante==gen)], y=data_gen$chlomoy2011),color="#F8766D",size=4)
    p <- p + geom_line(aes(x=valeurN$degresjour[which(valeurN$CodePlante==gen)], y=data_gen$chlomoy2011),color="#F8766D")
    p <- p + ylab("Teneur en Clorophylle")
  }else if(arg==13){
    p <- p + geom_point(aes(x=valeurN$degresjour[which(valeurN$CodePlante==gen)], y=data_gen$NDVImoy2011),color="#F8766D",size=4)
    p <- p + geom_line(aes(x=valeurN$degresjour[which(valeurN$CodePlante==gen)], y=data_gen$NDVImoy2011),color="#F8766D")
    p <- p + ylab("NDVI")
  }else if(arg==22){
    p <- p + geom_point(aes(x=valeurN$degresjour[which(valeurN$CodePlante==gen)], y=data_gen$EAUmoy2011),color="#F8766D",size=4)
    p <- p + geom_line(aes(x=valeurN$degresjour[which(valeurN$CodePlante==gen)], y=data_gen$EAUmoy2011),color="#F8766D")
    p <- p + ylab("NDVI")
  }
  p <- p + xlab("Degrés jour")
  
  
  # plot(valeurN$degresjour[which(valeurN$CodePlante==gen)],data_gen$azmoy2011,type='b',xlab='',ylab='',main=paste(gen,' azmoy'))
  # plot(valeurN$degresjour[which(valeurN$CodePlante==gen)],data_gen$chlomoy2011,type='b',xlab='',ylab='',main=paste(gen,' azmoy'))
  # plot(valeurN$degresjour[which(valeurN$CodePlante==gen)],data_gen$NDVImoy2011,type='b',xlab='',ylab='',main=paste(gen,' azmoy'))
  # plot(valeurN$degresjour[which(valeurN$CodePlante==gen)],data_gen$EAUmoy2011,type='b',xlab='',ylab='',main=paste(gen,' azmoy'))
  
  
  pointfin=c()
  pointfin2=c()
  
  if(any(is.element(names(saveme),names(savemeleretour)))){
    nom_gen=c(nom_gen,genotyp[[g]])
    for(k in 1:length(names(saveme[is.element(names(saveme),names(savemeleretour))]))){
      obj_date=names(saveme[is.element(names(saveme),names(savemeleretour))])[k]
      Y_deg=valeurN$degresjour[which(valeurN$CodePlante==gen)][which(obj_date==valeurN$degresjour[which(valeurN$CodePlante==gen)])]
      pointfin=c(pointfin,Y_deg)
    }
    for(t in 1:length(names(savemeleretour[is.element(names(savemeleretour),names(saveme))]))){
      obj_date=names(savemeleretour[is.element(names(savemeleretour),names(saveme))])[t]
      Y_deg=valeurN$degresjour[which(valeurN$CodePlante==gen)][which(obj_date==valeurN$degresjour[which(valeurN$CodePlante==gen)])]
      pointfin2=c(pointfin2,Y_deg)
    }
  }
  
  resupoint[[g]]=pointfin
  resupoint2[[g]]=pointfin2
  
}


resupoint
resupoint2

names(resupoint)=nom_gen
names(resupoint2)=nom_gen

setwd("~/Thèse/Script/Groupes AB brut/Résultat")


for(i in 1:length(genotyp)){
  plot(valeurN$degresjour[which(valeurN$CodePlante==genotyp[i])],valeurN2[[which(names(valeurN2)==genotyp[i])]][,arg],type='b',xlab='',ylab='',main=paste(genotyp[i],colnames(valeurN2[[i]])[arg]))
  for(k in 1:length(resupoint[[i]])){
    points(resupoint[[i]][k],
           valeurN2[[which(names(valeurN2)==genotyp[i])]][which(valeurN2[[which(names(valeurN2)==genotyp[i])]]$degresjour==resupoint[[i]][k]),arg],
           col='purple',cex=3)
  }
  # points(resupoint[[i]][1],
  #        valeurN2[[which(names(valeurN2)==genotyp[i])]][which(valeurN2[[which(names(valeurN2)==genotyp[i])]]$degresjour==resupoint[[i]][1]),arg],
  #        col='purple')
  # points(resupoint[[i]][2],
  #        valeurN2[[which(names(valeurN2)==genotyp[i])]][which(valeurN2[[which(names(valeurN2)==genotyp[i])]]$degresjour==resupoint[[i]][2]),arg],
  #        col='purple')
  points(mm2011$t0[which(mm2011$CodePlante==genotyp[i])], modele.ingrid(a = mm2011$a[which(mm2011$CodePlante==genotyp[i])] ,b = mm2011$b[which(mm2011$CodePlante==genotyp[i])],c = mm2011$c[which(mm2011$CodePlante==genotyp[i])] ,
                                                                        d = mm2011$d[which(mm2011$CodePlante==genotyp[i])],x=mm2011$t0[which(mm2011$CodePlante==genotyp[i])]),col='red',cex=3,pch=2)
  points(mm2011$t1[which(mm2011$CodePlante==genotyp[i])], modele.ingrid(a = mm2011$a[which(mm2011$CodePlante==genotyp[i])] ,b = mm2011$b[which(mm2011$CodePlante==genotyp[i])],c = mm2011$c[which(mm2011$CodePlante==genotyp[i])] ,
                                                                        d = mm2011$d[which(mm2011$CodePlante==genotyp[i])],x=mm2011$t1[which(mm2011$CodePlante==genotyp[i])]),col='blue',cex=3,pch=2)
  points(mm2011$tinf[which(mm2011$CodePlante==genotyp[i])], modele.ingrid(a = mm2011$a[which(mm2011$CodePlante==genotyp[i])] ,b = mm2011$b[which(mm2011$CodePlante==genotyp[i])],c = mm2011$c[which(mm2011$CodePlante==genotyp[i])] ,
                                                                          d = mm2011$d[which(mm2011$CodePlante==genotyp[i])],x=mm2011$tinf[which(mm2011$CodePlante==genotyp[i])]),col='green',cex=3,pch=2)
  dev.print(device = png, file = paste(genotyp[i],"_resupoint_nco=5_",colnames(valeurN2[[i]])[arg],".png",sep=''), width = 600)
}




