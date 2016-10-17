#---------------------------------------------
# Program to study longiudinal structures
# among partons in AMPT d+Au events
#
# JDOK - 10/14/16
#---------------------------------------------

visualize_partons <- function(collsyst)
{
  library(scatterplot3d)
  library(car)
  
  # Read the initial nucleon information from file
  if(collsyst == "dau")
  {
    initialNucleonInfo <- read.table("npart-xy.dat", skip = 1)
  }
  else if(collsyst == "auau")
  {
    initialNucleonInfo <- read.table("npart-xy_auau.dat", skip = 1)
  }  
  colnames(initialNucleonInfo) <- c("x","y","seq","status","znuc","junk1","junk2")

  # Determine the number of inelastically wounded nucleons (status = 3)
  npart <- sum(initialNucleonInfo$status == "3")

  # Read in initial parton information from file
  if(collsyst == "dau")
  {
    initialPartonInfo <- read.table("parton-initial-afterPropagation.dat", skip = 1)
  }
  else if(collsyst == "auau")
  {
    initialPartonInfo <- read.table("parton-initial-afterPropagation_AuAu.dat", skip = 1)
  }
  colnames(initialPartonInfo) <- c("pid","px","py","pz","e","x","y","z","t")
  
  # Impose a cut on the formation time t < 0.5 fm/c
  initialPartonInfo <- initialPartonInfo[which(initialPartonInfo$t < 0.5),]

  # Get the polar coordinates of partons and add to existing dataframe
  initialPartonInfo <- transform(initialPartonInfo, r=sqrt(x*x + y*y))
  initialPartonInfo <- transform(initialPartonInfo, phi=atan2(y,x))
  initialPartonInfo <- transform(initialPartonInfo, eta=atanh(r/sqrt(x*x + y*y + z*z)))
  
  # We now want to use the K-means algorithm to cluster the partons intro strings
  # or flux tubes based on the difference in transverse position of parton pairs.
  # The algorithm is seeded with the number of wounded nucleons in the event
  kMeansObject <- kmeans(initialPartonInfo[,c("r","phi")], npart, 500, npart, "Hartigan-Wong")
  print(kMeansObject)
  
  # Assign a cluster number to each parton
  initialPartonInfo$clus <- kMeansObject$cluster
  initialPartonInfo$clus <- factor(initialPartonInfo$clus)

  # Draw 3D visualization of initial partons and the clusters they belong to
  colors  <- rainbow(npart)
  
  if(collsyst == "dau")
  {
    markers <- as.numeric(initialPartonInfo$clus)
    s3d <- scatterplot3d(initialPartonInfo$x, initialPartonInfo$y, initialPartonInfo$z, 
                         color = colors[as.numeric(initialPartonInfo$clus)], 
                         pch = markers, 
                         xlab="x [fm]", ylab="y [fm]", zlab="z [fm]")
    
    plot(initialPartonInfo$r, initialPartonInfo$phi, 
         col = colors[as.numeric(initialPartonInfo$clus)], 
         xlab="r [fm]", ylab="phi [rad]", pch = markers)
  }
  else if(collsyst == "auau")
  {
    s3d <- scatterplot3d(initialPartonInfo$x, initialPartonInfo$y, initialPartonInfo$z, 
                         color = colors[as.numeric(initialPartonInfo$clus)], 
                         xlab="x [fm]", ylab="y [fm]", zlab="z [fm]")
    
    plot(initialPartonInfo$r, initialPartonInfo$phi, 
         col = colors[as.numeric(initialPartonInfo$clus)], xlab="r [fm]", ylab="phi [rad]")
  }
  
  #legend(s3d$xyz.convert(0.4, 3.7, -1.85), yjust=0, pch = markers,
  #       legend = levels(initialPartonInfo$clus), col = colors)
  
  #Determine, for each cluster the min and max z
  clusters <- levels(initialPartonInfo$clus)
  minz  <- aggregate(initialPartonInfo$z, by=list(initialPartonInfo$clus), FUN=min)[2]
  maxz  <- aggregate(initialPartonInfo$z, by=list(initialPartonInfo$clus), FUN=max)[2]
  meanr <- aggregate(initialPartonInfo$r, by=list(initialPartonInfo$clus), FUN=mean)[2]
  
  clusterInfo <- data.frame(clusters, minz, maxz, meanr)
  colnames(clusterInfo) <- c("clus","minz","maxz","meanr")
  
  #Determine the extent in pseudorapidity of each cluster
  clusterInfo <- transform(clusterInfo, etamin=atanh(meanr/sqrt(meanr*meanr + minz*minz)))
  clusterInfo <- transform(clusterInfo, etamax=atanh(meanr/sqrt(meanr*meanr + maxz*maxz)))
  
  if(collsyst == "dau")
  {
    plot(clusterInfo$etamin, clusterInfo$meanr, col = colors, xlab = "Eta", ylab = "r", xlim = c(-4,6), pch = markers)
    points(clusterInfo$etamax, clusterInfo$meanr, col = colors, pch = markers)
  }
  else if(collsyst == "auau")
  {
    plot(clusterInfo$etamin, clusterInfo$meanr, col = colors, xlab = "Eta", ylab = "r", xlim = c(-4,6))
    points(clusterInfo$etamax, clusterInfo$meanr, col = colors)
  }
  segments(clusterInfo$etamin, clusterInfo$meanr, clusterInfo$etamax, clusterInfo$meanr, col = colors, lwd=2)
}
