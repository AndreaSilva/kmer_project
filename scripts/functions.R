# methods used for analysis
calc_dissimilarity<-function(data_dir, subset_names, K, method){
  tmp<-vector()
  j=1
  for (l in subset_names){
    tmp[j]<-system.time({
    ptm <- proc.time()
    #l="China"
    filename<-paste(data_dir,l,"_KmerMatr_K",K,".RData", sep="")
    load(file=filename)                       # Read Data from files              
    KmerMatr<-apply(KmerMatr,2,function(x){     # make matrix of normalised frequencies
      x/sum(x)
    }) 
    Dist <- my_dist(KmerMatr, method)                  # Calculating distance for subset not filter  
    filename<-paste("../output/dist_",method,"_",l ,"_K", K,".RData", sep="")
    save(Dist, file=filename)
    if (l=="China"){
      KmerMatr_all<-KmerMatr
    }
    if (l=="USA"){
      KmerMatr_all<-cbind(KmerMatr, KmerMatr_all)
      Dist_all <- my_dist(KmerMatr_all, method) 
      filename<-paste("../output/dist_",method,"_all_K", K,".RData", sep="")
      save(Dist_all, file=filename)
    }
    })[3]
    j=j+1
    #return(proc.time() - ptm)
  }
  return(tmp)
}

# count dispersion
disp <- function(Data){  
  res <- apply(Data, 1, sd)
  res
}

# draw pcoa
draw_pcoa<-function(dist, method){
  library(ape)
  library(ggplot2)
  m.BC<-as.matrix(dist)
  res <- pcoa(m.BC)
  perc_expl <- paste(round(100*(res$values$Relative_eig)[1:10], 1), "%", sep="")
  pn <- paste("PC", 1:length(perc_expl), sep="")
  perc_expl <- paste(pn, perc_expl, sep=": ")
  perc_expl
  cbbPalette <- c("red", "blue", "green", "yellow")
  tmp<-as.data.frame(res$vectors[,1:2])
  tmp$group<-"USA"
  tmp$group[130:281]<-"China"
  #plotPath<-paste("../output/PCoA_",method,".pdf", sep="")
  #pdf(plotPath, height= 8, width = 10)
  p<-ggplot(tmp, aes(x=Axis.1, y=Axis.2, colour=group)) + 
    geom_point(size=4)  + 
    scale_x_continuous(perc_expl[1]) + 
    scale_y_continuous(perc_expl[2]) + 
    guides(colour=FALSE) +
    scale_colour_manual(values=cbbPalette)+theme_bw()
  #print(p)
  return(p)
  #dev.off()
}

#count critical dispersion value for sorting (don't count zeros) treshold from low, high, symm disp
get_crit_filter_value<- function(tmp_disp, fraction, method)
{  
  size=length(tmp_disp)
  j=0
  buf<-vector()
  res<-vector()
  for (i in 1:size){
    if (tmp_disp[[i]]!=0){            # check if Null
      j=j+1
      buf[j]<-tmp_disp[[i]]
    } 
  }
  
  if(method == "low")
  {
    res[1]<-quantile(buf, fraction)
    res[2]<-1000000000
  }
  else
  {
    if(method == "high")
    {
      res[1]<-0
      res[2]<-quantile(buf, 1-fraction)
    }
    else
    {
      if(method == "both")
      {
        res[1]<-quantile(buf, fraction/2)
        res[2]<-quantile(buf, 1-fraction/2)
      }
    }
  }
  return(res)
}

my_dist <- function(x, method) {
  if(method == "COR"){
    new_dist <- 1-cor(x, method="spearman")
    new_dist<-as.dist(new_dist)
  }
  if(method == "BC"){
    new_dist <- bcdist(t(x))
  }
  return(new_dist)
}

m1_vs_m2<-function(d.1.1, d.1.2,d.2.1, d.2.2, method1, method2, xlim){
  cor1<-round(cor(d.1.1, d.1.2, method="spearman"), digits=2)
  cor2<-round(cor(d.2.1, d.2.2, method="spearman"), digits=2)
#   x.ax<-paste("BC ", method1, sep="")
#   y.ax<-paste("BC ", method2, sep="")
#   plotPath<-paste("../output/",method1,"_vs_", method2,".pdf", sep="")
#   pdf(plotPath, height= 5, width = 5)
  plot(d.1.1, d.1.2, col=alpha("blue", 0.1), pch=17,
       xlab=method1, ylab=method2, xlim=c(0,xlim), ylim=c(0,0.7))
  points(d.2.1, d.2.2, col=alpha("red", 0.1), pch=17)
  #   cor_1<-paste("USA (cor = ", round(cor(as.dist(m.1.1), as.dist(m.1.2), method="spearman"), digits=2),")",  sep="")
  #   cor_2<-paste("China (cor = ", round(cor(as.dist(m.2.1), as.dist(m.2.2), method="spearman"), digits=2),")", sep="")
  #   legend("topleft", c(cor_1,cor_2) , pch = 17, pt.bg = "white", 
  #          lty = 0, col = c("blue", red), bty = "n")
  cor_1<-paste("R=",cor1,  sep="")
  cor_2<-paste("R=",cor2,  sep="")
  legend("topleft", c("USA",cor_1,"China",cor_2) , pch = 17, pt.bg = "white", 
         lty = 0, col = c("blue", "white", "red","white"), bty = "n")
  #return(p)
  #dev.off()
  #return(c(cor1, cor2))
}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
