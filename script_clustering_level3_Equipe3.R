# This function is useful to draw gene expression profiles
plotGenes <- function(expData, title = "", yMax = NULL, meanProfile = TRUE){
  
  # Check function parameters
  if(is.null(yMax)){
    
    print("You must specify a maximal value for Y axis")
    
  }else{
    
    # Representation of the first expression profile
    plot(1:ncol(expData), expData[1,], col = "grey", type = "l",
         ylim = c(0, ceiling(yMax)),
         xlab = "Time point", ylab = "Gene expression level",
         main = title)
    
    # Add expression profile for other genes
    for(i in 2:nrow(expData)){
      
      lines(1:ncol(expData), expData[i,], col = "grey")
      
      # end of for()  
    }
    
    # Average expression profile
    if(meanProfile == TRUE){
      expMean = apply(expData, 2, mean)
      lines(1:ncol(expData), expMean, col = "red", 
            lwd = 1.5, lty = "dashed")
    }
    
    # end of else()   
  }
  
  # end of function plotGenes()  
}


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
  

  
  pdf("Plot_cluster.pdf")
  for (i in c(1:n_clust)){
    
    geneCluster = names(which(vecCluster == i))
    cluster = expMatrix[geneCluster,]
    name <- paste("Cluster", i,"_",dist_meth,"_",clust_meth, sep=" ")
    cluster_mat = as.matrix(cluster)
    
    plotGenes(cluster, title = name, yMax = max(expMatrix))
    heatmap(cluster_mat, main= name)
  }
  dev.off()
  
}


#My_data = nom du fichier
#dist_meth = "eucli" , "corr" ou "Max" - methode calcul des distances
#clust_meth = "km" ou "hcl" - methode de clustering
#n_clust = nombre de cluster
Fonction_cluster("Mito_Genes.txt","Max","km",6)
