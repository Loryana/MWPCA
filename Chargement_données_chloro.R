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
#chargement des sepctres bruts
load("Données.RData")


################################################################################pour la chlorophylle

data_chloro=data

for(i in 1:length(data_chloro)){
  data_chloro[[i]]=data_chloro[[i]][,which(as.numeric(colnames(data_chloro[[1]]))<=1350 & as.numeric(colnames(data_chloro[[1]]))>=450)]
}

for(i in 1:length(data_chloro)){
  data_chloro[[i]]=snv(data_chloro[[i]])
}

for(i in 1:length(data_chloro)){
  data_chloro[[i]]=savgol(data_chloro[[i]],n=9,p=2,m=1)
}
plot(as.numeric(colnames(data_chloro[[1]])),data_chloro[[1]][1,],type='l')


save(data_chloro,file='Données_chloro.RData')
