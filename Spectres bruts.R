#Spectres bruts

rm(list=ls())

setwd("~/Thèse/Script/")

asd <- source("asd_read.R")

#fonction pour récupérer le nom spécifiques (et date) des fichiers / dossiers
infofile <- function(x,sta,sto){substr(x,sta,sto)}

setwd("~/Thèse/Données Blé/")
path_cwd <- "C:/Users/heloise.villesseche/Documents/Thèse/Données Blé"
path_cwd

filename=dir(path = path_cwd)
#récupération des nom de dossier
filedate=sapply(filename,infofile,sta=1,sto=4)
#récupération des noms de dossier qui nous intéressent uniquement
filename=filename[which(filedate=="2011")]

# #retrait des dossier contenant des fichirs vides (0 Ko)
# # filename=filename[-c(which(filename=="20110506"),which(filename=="20110509")
# #                      ,which(filename=="20110513"))]
# 
# setwd("~/Thèse/Script/")
#ouverture des dossier contenant les spectres
data=list()
for(i in 1:length(filename)){
  data[[i]]=SIGNE_load(filename[i])
}
plot(as.numeric(colnames(data[[1]])),data[[1]][1,],type='l')


for(i in 1:length(data)){
  data[[i]]=data[[i]][,which(as.numeric(colnames((data[[i]])))>=450)]
}
plot(as.numeric(colnames(data[[1]])),data[[1]][1,],type='l')


for(i in 1:length(data)){
   data[[i]]=adj_asd(data[[i]],c(602, 1402))
}
plot(data[[1]][1,],type='l')
plot(data[[5]][1,],type='l')


save(data,file="Données.RData")
