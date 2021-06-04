library(tidyverse)
library(ggsci)
library(RColorBrewer)

Fonction_cluster = function(My_data,dist_meth,clust_meth,n_clust){
  
  
  expMatrix = read.table(My_data, header = T, row.names = 1)
  
  if(dist_meth =="eucli"){
    matDist = dist(expMatrix) 
  }else if(dist_meth =="corr"){
    matDist = as.dist(1 - cor(t(expMatrix)))
  }else if(dist_meth =="Max"){
    matDist = dist(expMatrix, method = "maximum")
  }
  
  if(clust_meth=="km"){
    res = kmeans(matDist, n_clust)
    vecCluster = res$cluster
  }else if(clust_meth=="hcl"){
    res = hclust(matDist)
    vecCluster = cutree(res, n_clust)
  }
  


  #fusionner les cluster avec les données de départ
  expMatrix = cbind(as.data.frame(vecCluster),data.frame(expMatrix))
  #mettre les noms de lignes dans la première colonne (pour transformation plus tard)
  cluster <- cbind(rownames(expMatrix), data.frame(expMatrix, row.names=NULL))
  #changer le nom des colonnes pour avoir juste le numéro des points de temps
  colnames(cluster) <- gsub("^[^-]+T", "", colnames(cluster))
  #mise en forme des données pour être utilisble par ggplot2
  cluster_mod = cluster %>% pivot_longer(cols = -c("rownames(expMatrix)","vecCluster"))
  #changement des noms de colonnes
  colnames(cluster_mod) = c("gene","cluster","time_point","value")
  #définir les points de temps comme valeur numérique et les cluster comme caractère (pour ggplot)
  cluster_mod <- transform(cluster_mod, time_point = as.numeric(time_point))
  cluster_mod2 <- transform(cluster_mod, cluster = as.character(cluster))
  

  
 
  
  #plot gene cluster
 graph <- ggplot(cluster_mod2, aes(x=time_point,y=value,color=cluster,shape=gene))+
    geom_line(alpha=0.7, size=0.5)+
    facet_wrap(~cluster)+
    scale_color_futurama()+
    theme_minimal()+
    theme(legend.position = "none")+xlab("Points de temps")+ylab("Expression")
 
 ggsave("Graph.png", graph)
 
 pdf("heatmap.pdf")
  for (i in c(1:n_clust)){
      geneCluster = names(which(vecCluster == i))
    cluster = expMatrix[geneCluster,]
    name <- paste("Cluster", i,"_",dist_meth,"_",clust_meth, sep=" ")
    cluster_mat = as.matrix(cluster)
    
    coul <- colorRampPalette(brewer.pal(8, "Blues"))(25)
    heatmap(cluster_mat,col = coul,
            main= name, xlab = "time points", ylab ="genes",
            cexRow = 0.5, cexCol = 0.5)
  }
  dev.off()
  
  #extraire liste de gènes
  gene_list <- NULL
  for (i in c(1:n_clust)){
    geneCluster <- as.data.frame(names(which(vecCluster == i)))
    geneCluster <- cbind(geneCluster,i)
    gene_list <- rbind(gene_list,geneCluster)
    
  }
  colnames(gene_list) <- c("Genes","Cluster")
  write.table(gene_list, "gene_list.txt")
}

#My_data = nom du fichier
#dist_meth = "eucli" , "corr" ou "Max" - methode calcul des distances
#clust_meth = "km" ou "hcl" - methode de clustering
#n_clust = nombre de cluster


Fonction_cluster("Mito_Genes.txt","Max","km",4)








gene_list = read.table("gene_list.txt", header = T, row.names = 1)

#séparer nom du gène et numéro
gene_list <- gene_list %>% separate(c("Genes"), sep="([|])", into= c("number","name"))

#Cluster 1
write.table(gene_list %>% filter(Cluster==1) %>% select(number), "cluster1.txt", quote=F, row.names = F, col.names = F)
#Cluster 2
write.table(gene_list %>% filter(Cluster==2) %>% select(number), "cluster2.txt", quote=F, row.names = F, col.names = F)
#Cluster 3
write.table(gene_list %>% filter(Cluster==3) %>% select(number), "cluster3.txt", quote=F, row.names = F, col.names = F)
#Cluster 4
write.table(gene_list %>% filter(Cluster == 4) %>% select(number), "cluster4.txt", quote=F, row.names = F, col.names = F)


