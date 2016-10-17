#---------------------------------------------
# Program to study longiudinal structures
# among partons in AMPT d+Au events
#
# JDOK - 10/14/16
#---------------------------------------------

visualize_partons <- function(collsyst,b)
{
  library(scatterplot3d)
  library(car)
  
  # Read the initial nucleon information from file
  if(collsyst == "dau")
  {
    initialNucleonInfo <- read.table(paste("npart-xy_dAu_", b, "fm.dat", sep=""), skip = 1)
  }
  else if(collsyst == "auau")
  {
    initialNucleonInfo <- read.table(paste("npart-xy_AuAu_", b, "fm.dat", sep=""), skip = 1)
  }  
  colnames(initialNucleonInfo) <- c("x","y","seq","status","znuc","junk1","junk2")

  # Determine the number of inelastically wounded nucleons (status = 3)
  npart <- sum(initialNucleonInfo$status == "3")

  # Read in initial parton information from file
  if(collsyst == "dau")
  {
    initialPartonInfo <- read.table(paste("parton-initial-afterPropagation_dAu_", b, "fm.dat", sep=""), skip = 1)
  }
  else if(collsyst == "auau")
  {
    initialPartonInfo <- read.table(paste("parton-initial-afterPropagation_AuAu_", b, "fm.dat", sep=""), skip = 1)
  }
  colnames(initialPartonInfo) <- c("pid","px","py","pz","m","x","y","z","t")
  
  # Plot the distribution of formation times
  #timeCounts <- table(initialPartonInfo$t)
  #barplot(timeCounts, main="Formation Time", xlab="Formation Time")
  
  # Impose a cut on the formation time t < 0.5 fm/c
  initialPartonInfo <- initialPartonInfo[which(initialPartonInfo$t < 0.5),]

  # Get the polar coordinates of partons and add to existing dataframe
  initialPartonInfo <- transform(initialPartonInfo, r=sqrt(x*x + y*y))
  initialPartonInfo <- transform(initialPartonInfo, phi=atan2(y,x))
  initialPartonInfo <- transform(initialPartonInfo, gamma=1/sqrt(1-(px*px+py*py+pz*pz)/(px*px+py*py+pz*pz+m*m)))
  initialPartonInfo <- transform(initialPartonInfo, tau=t/gamma)
  initialPartonInfo <- transform(initialPartonInfo, rap=asinh(z/tau))
  
  # We now want to use the K-means algorithm to cluster the partons intro strings
  # or flux tubes based on the difference in transverse position of parton pairs.
  # The algorithm is seeded with the number of wounded nucleons in the event
  kMeansObject <- kmeans(initialPartonInfo[,c("r","phi")], npart, 200, npart, "Hartigan-Wong")

  # Assign a cluster number to each parton
  initialPartonInfo$clus <- kMeansObject$cluster
  initialPartonInfo$clus <- factor(initialPartonInfo$clus)

  # Draw 3D visualization of initial partons and the clusters they belong to
  colors  <- rainbow(npart)
  markers <- as.numeric(initialPartonInfo$clus)
  
  if(collsyst == "dau")
  {
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
  minrap  <- aggregate(initialPartonInfo$rap, by=list(initialPartonInfo$clus), FUN=min)[2]
  maxrap  <- aggregate(initialPartonInfo$rap, by=list(initialPartonInfo$clus), FUN=max)[2]
  meanr <- aggregate(initialPartonInfo$r, by=list(initialPartonInfo$clus), FUN=mean)[2]
  meanx <- aggregate(initialPartonInfo$x, by=list(initialPartonInfo$clus), FUN=mean)[2]
  meany <- aggregate(initialPartonInfo$y, by=list(initialPartonInfo$clus), FUN=mean)[2]

  clusterInfo <- data.frame(clusters, minrap, maxrap, meanr, meanx, meany)
  colnames(clusterInfo) <- c("clus","minrap","maxrap","meanr", "meanx", "meany")

  if(collsyst == "dau")
  {
    plot(clusterInfo$minrap, clusterInfo$meanx, col = colors, xlab = "Eta", ylab = "r", xlim = c(-8,8), pch = markers)
    points(clusterInfo$maxrap, clusterInfo$meanx, col = colors, pch = markers)
  }
  else if(collsyst == "auau")
  {
    plot(clusterInfo$minrap, clusterInfo$meanx, col = colors, xlab = "Eta", ylab = "r", xlim = c(-8,8), ylim = c(-8,8))
    points(clusterInfo$maxrap, clusterInfo$meanx, col = colors)
  }
  
  segments(clusterInfo$minrap, clusterInfo$meanx, clusterInfo$maxrap, clusterInfo$meanx, col = colors, lwd=2)
}
