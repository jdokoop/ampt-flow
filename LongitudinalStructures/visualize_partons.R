#---------------------------------------------
# Program to study longiudinal structures
# among partons in AMPT d+Au events
#
# JDOK - 10/14/16
#---------------------------------------------

visualize_partons <- function()
{
  library(scatterplot3d)
  library(car)
  
  # Read the initial nucleon information from file
  initialNucleonInfo <- read.table("npart-xy.dat", skip = 1)
  colnames(initialNucleonInfo) <- c("x","y","seq","status","znuc","junk1","junk2")

  # Determine the number of inelastically wounded nucleons (status = 3)
  npart <- sum(initialNucleonInfo$status == "3")

  # Read in initial parton information from file
  initialPartonInfo <- read.table("parton-initial-afterPropagation.dat", skip = 1)
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

  # Assign a cluster number to each parton
  initialPartonInfo$clus <- kMeansObject$cluster
  initialPartonInfo$clus <- factor(initialPartonInfo$clus)
  levels(initialPartonInfo$clus)
  
  # Draw 3D visualization of initial partons and the clusters they belong to
  colors  <- rainbow(npart)
  markers <- as.numeric(initialPartonInfo$clus)
  s3d <- scatterplot3d(initialPartonInfo$x, initialPartonInfo$y, initialPartonInfo$z, 
                       color = colors[as.numeric(initialPartonInfo$clus)], 
                       pch = markers, 
                       xlab="x [fm]", ylab="y [fm]", zlab="z [fm]")
  
  legend(s3d$xyz.convert(0.4, 3.7, -1.85), yjust=0, pch = markers,
         legend = levels(initialPartonInfo$clus), col = colors)
  
  # Draw a 2D scatterplot of polar parton coordinates
  plot(initialPartonInfo$r, initialPartonInfo$phi, col = colors[as.numeric(initialPartonInfo$clus)], pch = markers, xlab="r [fm]", ylab="phi [rad]")
  #scatterplot(r ~ phi | clus, data=initialPartonInfo, xlab="phi", ylab="r")
  
  #Determine, for each cluster
}
