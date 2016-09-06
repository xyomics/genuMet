#' @import grid
multiplot <- function(..., plotlist=NULL, cols) {
  #  require(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # Make the panel
  plotCols = cols                          # Number of columns of plots
  plotRows = ceiling(numPlots/plotCols) # Number of rows needed, calculated from # of cols

  # Set up the page
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(plotRows, plotCols)))
  vplayout <- function(x, y)
    viewport(layout.pos.row = x, layout.pos.col = y)

  # Make each plot, in the correct location
  for (i in 1:numPlots) {
    curRow = ceiling(i/plotCols)
    curCol = (i-1) %% plotCols + 1
    print(plots[[i]], vp = vplayout(curRow, curCol ))
  }

}

#' @title Plot missing rate
#' @import ggplot2
#' @param bds A list returned by \code{\link{falsignal}}.
#' @param method "all" plots the missing rate of all metabolites. "bad1" only plots the missing rate of type 1 bad signal. "bad2" only plots the missing rate of type 2 bad signal. "single" only plots missing rate for a given metabolite.
#' @param mname The name of metabolite if method="single".
#' @return none
#' @author Liu Cao
#' @seealso See \code{\link{falsignal}} for how to generate the input list \code{bds}
#' @export

plot_mrate <- function(bds, type="1", mname){

  if(type == "single"){
    if(sum(bds$mbnames == mname)==0){
      print(paste0("No metabolites called ",mname," in the dataset"))
      return
    }
    mrate = bds$mrate[bds$mbnames == mname,]
    df1 = data.frame(x=1:dim(bds$mrate)[2],y=mrate)
    ggplot(data=df1,aes(x=x,y=y))+geom_line()+ylim(0,1)+labs(title=paste0(mname," ws=",bds$wsize," ss=",bds$ssize),x="#window",y="missing rate")
  }
  else{
    if(type != "1" & type !="2"){
      print("type can only be 1, 2, or single")
      return
    }
    mratebad = NULL
    mrategood = NULL
    if(type == "1"){
      boolbad = bds$boolbad1
    }
    if(type == "2"){
      boolbad = bds$boolbad2
    }
    mratebad = bds$mrate[boolbad,]
    mrategood = bds$mrate[!boolbad,]
    mbnames = bds$mbnames[boolbad]

    ### false
    mrate = mratebad
    pdf(file=paste0("type",type,"_false_mr.pdf"))

    j=1
    while(j <=length(mbnames)){
      df1 = data.frame(x=1:dim(mrate)[2],y=mrate[j,])
      p1 = ggplot(data=df1,aes(x=x,y=y))+geom_line()+ylim(0,1)+labs(title=paste0(mbnames[j]," ws=",bds$wsize," ss=",bds$ssize),x="#window",y="missing rate")
      j=j+1

      if(j<=length(mbnames)){
        df2 = data.frame(x=1:dim(mrate)[2],y=mrate[j,])
        p2 = ggplot(data=df2,aes(x=x,y=y))+geom_line()+ylim(0,1)+labs(title=paste0(mbnames[j]," ws=",bds$wsize," ss=",bds$ssize),x="#window",y="missing rate")
      }
      else{
        multiplot(p1,cols=1)
        break
      }
      j=j+1

      if(j<=length(mbnames)){
        df3 = data.frame(x=1:dim(mrate)[2],y=mrate[j,])
        p3 = ggplot(data=df3,aes(x=x,y=y))+geom_line()+ylim(0,1)+labs(title=paste0(mbnames[j]," ws=",bds$wsize," ss=",bds$ssize),x="#window",y="missing rate")
      }
      else{
        multiplot(p1,p2,cols=1)
        break
      }
      j=j+1

      if(j<=length(mbnames)){
        df4 = data.frame(x=1:dim(mrate)[2],y=mrate[j,])
        p4 = ggplot(data=df4,aes(x=x,y=y))+geom_line()+ylim(0,1)+labs(title=paste0(mbnames[j]," ws=",bds$wsize," ss=",bds$ssize),x="#window",y="missing rate")
        multiplot(p1,p2,p3,p4,cols=2)
      }
      else{
        multiplot(p1,p2,p3,cols=1)
        break
      }
      j=j+1

      if((j-1)%%40 == 0){
        print(paste0(j," metabolites plotted"))
      }
    }
    dev.off()
    print(paste0(dim(mrate)[1], " type ",type,  " false metabolites has been plotted. Saved in type",type,"_false_mr.pdf"  ))

    ### good
    mrate = mrategood
    pdf(file=paste0("type",type,"_good_mr.pdf"))

    j=1
    while(j <=length(mbnames)){
      df1 = data.frame(x=1:dim(mrate)[2],y=mrate[j,])
      p1 = ggplot(data=df1,aes(x=x,y=y))+geom_line()+ylim(0,1)+labs(title=paste0(mbnames[j]," ws=",bds$wsize," ss=",bds$ssize),x="#window",y="missing rate")
      j=j+1

      if(j<=length(mbnames)){
        df2 = data.frame(x=1:dim(mrate)[2],y=mrate[j,])
        p2 = ggplot(data=df2,aes(x=x,y=y))+geom_line()+ylim(0,1)+labs(title=paste0(mbnames[j]," ws=",bds$wsize," ss=",bds$ssize),x="#window",y="missing rate")
      }
      else{
        multiplot(p1,cols=1)
        break
      }
      j=j+1

      if(j<=length(mbnames)){
        df3 = data.frame(x=1:dim(mrate)[2],y=mrate[j,])
        p3 = ggplot(data=df3,aes(x=x,y=y))+geom_line()+ylim(0,1)+labs(title=paste0(mbnames[j]," ws=",bds$wsize," ss=",bds$ssize),x="#window",y="missing rate")
      }
      else{
        multiplot(p1,p2,cols=1)
        break
      }
      j=j+1

      if(j<=length(mbnames)){
        df4 = data.frame(x=1:dim(mrate)[2],y=mrate[j,])
        p4 = ggplot(data=df4,aes(x=x,y=y))+geom_line()+ylim(0,1)+labs(title=paste0(mbnames[j]," ws=",bds$wsize," ss=",bds$ssize),x="#window",y="missing rate")
        multiplot(p1,p2,p3,p4,cols=2)
      }
      else{
        multiplot(p1,p2,p3,cols=1)
        break
      }
      j=j+1

      if((j-1)%%40 == 0){
        print(paste0(j," metabolites plotted"))
      }
    }
    dev.off()
    print(paste0(dim(mrate)[1], " type ", type, " good metabolites has been plotted. Saved in type",type, "_good_mr.pdf"  ))
  }
}

#' @title Plot QC figures
#' @description \code{\link{plot_qc}} plot a set of quality control figures that summarize the distribution of
#'  missing pattern and features of predicted false metabolites and good metabolites.
#' @import ggplot2
#' @param bds A list returned by falsignal().
#' @param method "all" plots the missing rate of all metabolites. "bad1" only plots the missing rate of type 1 bad signal. "bad2" only plots the missing rate of type 2 bad signal. "single" only plots missing rate for a given metabolite.
#' @param mname The name of metabolite if method="single".
#' @return none
#' @author Liu Cao
#' @export
plot_qc <- function(bds){
  pdf(file=paste0("QC.pdf"))
  plot_heatmap(bds)
  plot_variance_distribution(bds)
  plot_meanmrate(bds)
  plot_longestblock_distribution(bds)
  plot_switches_distribution(bds)
  dev.off()
}

plot_heatmap <- function(bds){
    # data_preparation
    bad1 = bds$mrate[bds$boolbad1,]
    bad1good = bds$mrate[(!bds$boolbad1) & bds$meanmrate >0.01,]

    bad2 = bds$mrate[bds$boolbad2,]
    bad2good = bds$mrate[(!bds$boolbad2) & bds$meanmrate >0.01,]

    # heatmap of missing rate
    par(mfrow=c(2,2))

    hrbad1good <- hclust(dist(bad1good),method="complete")
    bad1good=bad1good[rev(hrbad1good$order),]
    image(t(bad1good), xaxt= "n", yaxt= "n", col=gray((0:100)/100),
          xlab=paste0("#windows=",dim(bad1good)[2]),
          ylab=paste0("#matabolites=",dim(bad1good)[1]),
          main=paste0(dim(bad1good)[1]," Type I good metabolic features \n with mean missing rate > 0.01"))

    hrbad1 <- hclust(dist(bad1),method="complete")
    bad1=bad1[rev(hrbad1$order),]
    image(t(bad1), xaxt= "n", yaxt= "n", col=gray((0:100)/100),
          xlab=paste0("#windows=",dim(bad1)[2]),
          ylab=paste0("#matabolites=",dim(bad1)[1]),
          main=paste0(dim(bad1)[1]," Type I Bad metabolic features"))

    hrbad2good <- hclust(dist(bad2good),method="complete")
    bad2good=bad1good[rev(hrbad2good$order),]
    image(t(bad2good), xaxt= "n", yaxt= "n", col=gray((0:100)/100),
          xlab=paste0("#windows=",dim(bad2good)[2]),
          ylab=paste0("#matabolites=",dim(bad2good)[1]),
          main=paste0(dim(bad1good)[1]," Type II good metabolic features \n with mean missing rate > 0.01"))

    hrbad2 <- hclust(dist(bad2),method="complete")
    bad2=bad2[rev(hrbad2$order),]
    image(t(bad2), xaxt= "n", yaxt= "n", col=gray((0:100)/100),
          xlab=paste0("#windows=",dim(bad2)[2]),
          ylab=paste0("#matabolites=",dim(bad2)[1]),
          main=paste0(dim(bad2)[1]," Type II Bad metabolic features"))
}

plot_variance_distribution <- function(bds){
  # variance distribution
  # type I
  df1 = data.frame(tag = bds$boolbad1, mmrate = bds$meanmrate, variance = bds$variance, nswitches = bds$nswitches, longestblock = bds$longestblock)
  p1 = ggplot(subset(df1,tag),aes(x=variance))+geom_histogram(aes(x=variance,y=-..density..,fill="predicted artifact"))+geom_histogram(data=subset(df1,!tag),aes(x=variance,y=..density..,fill="predicted true"))+geom_vline(xintercept=bds$cutoffvar, colour="blue", linetype = "longdash")+
    scale_y_continuous(breaks=NULL)+scale_fill_manual(values=c("#E74C3C", "#1ABC9C"),guide=guide_legend(title=""))+labs(y="Density",title="Distribution of Variance, criterion I",x="") +
    theme_bw()+theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),legend.key = element_blank())

  # type II
  df2 = data.frame(tag = bds$boolbad2, mmrate = bds$meanmrate, variance = bds$variance, nswitches = bds$nswitches, longestblock = bds$longestblock)
  p2 = ggplot(subset(df2,tag),aes(x=variance))+geom_histogram(aes(x=variance,y=-..density..,fill="predicted artifact"))+geom_histogram(data=subset(df2,!tag),aes(x=variance,y=..density..,fill="predicted true"))+geom_vline(xintercept=bds$cutoffvar, colour="blue", linetype = "longdash")+
    scale_y_continuous(breaks=NULL)+scale_fill_manual(values=c("#E74C3C", "#1ABC9C"),guide=guide_legend(title=""))+labs(y="Density",title="Distribution of Variance, criterion II",x="") +
    theme_bw()+theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),legend.key = element_blank())

  # all
  df3 = data.frame(variance=bds$variance)
  p3 = ggplot(df3,aes(x=variance,y=..density..))+geom_histogram()+geom_vline(xintercept=bds$cutoffvar, colour="blue", linetype = "longdash")+
    scale_y_continuous(breaks=NULL)+scale_fill_manual(values=c("#E74C3C", "#1ABC9C"),guide=guide_legend(title=""))+labs(y="Density",title=paste0("Distribution of Variance, all metabolic features"),x=paste0("cutoff =",bds$cutoffvar)) +
    theme_bw()+theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),legend.key = element_blank())

  multiplot(p3, p1,p2,cols=1)
}

plot_longestblock_distribution <- function(bds){
  # longest block distribution
  # bad1
  # type I
  df1 = data.frame(tag = bds$boolbad1, mmrate = bds$meanmrate, variance = bds$variance, nswitches = bds$nswitches, longestblock = bds$longestblock)
  p1 = ggplot(subset(df1,tag),aes(x=longestblock))+geom_histogram(aes(x=longestblock,y=-..density..,fill="predicted artifact"),binwidth=1)+geom_histogram(data=subset(df1,!tag),aes(x=longestblock,y=..density..,fill="predicted true"),binwidth=1)+geom_vline(xintercept=bds$cutofflong, colour="blue", linetype = "longdash")+
    scale_y_continuous(breaks=NULL)+scale_fill_manual(values=c("#E74C3C", "#1ABC9C"),guide=guide_legend(title=""))+labs(y="Density",title="Distribution of Longest Block, criterion I",x="") +
    theme_bw()+theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),legend.key = element_blank())

  # type II
  df2 = data.frame(tag = bds$boolbad2, mmrate = bds$meanmrate, variance = bds$variance, nswitches = bds$nswitches, longestblock = bds$longestblock)
  p2 = ggplot(subset(df2,tag),aes(x=longestblock))+geom_histogram(aes(x=longestblock,y=-..density..,fill="predicted artifact"),binwidth=1)+geom_histogram(data=subset(df2,!tag),aes(x=longestblock,y=..density..,fill="predicted true"),binwidth=1)+geom_vline(xintercept=bds$cutofflong, colour="blue", linetype = "longdash")+
    scale_y_continuous(breaks=NULL)+scale_fill_manual(values=c("#E74C3C", "#1ABC9C"),guide=guide_legend(title=""))+labs(y="Density",title="Distribution of Longest Block, criterion II",x="") +
    theme_bw()+theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),legend.key = element_blank())

  # all
  df3 = data.frame(longestblock=bds$longestblock)
  p3 = ggplot(df3,aes(x=longestblock,y=..density..))+geom_histogram(binwidth=1)+geom_vline(xintercept=bds$cutofflong, colour="blue", linetype = "longdash")+
    scale_y_continuous(breaks=NULL)+scale_fill_manual(values=c("#E74C3C", "#1ABC9C"),guide=guide_legend(title=""))+labs(y="Density",title=paste0("Distribution of Longest Block, all metabolic features"),x=paste0("cutoff =",bds$cutofflong)) +
    theme_bw()+theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),legend.key = element_blank())

  multiplot(p3, p1,p2,cols=1)
}

plot_switches_distribution <- function(bds){
  # #switches distribution
  # type I
  df1 = data.frame(tag = bds$boolbad1, mmrate = bds$meanmrate, variance = bds$variance, nswitches = bds$nswitches, longestblock = bds$longestblock)
  p1 = ggplot(subset(df1,tag),aes(x=nswitches))+geom_histogram(aes(x=nswitches,y=-..density..,fill="predicted artifact"),binwidth=1)+geom_histogram(data=subset(df1,!tag),aes(x=nswitches,y=..density..,fill="predicted true"),binwidth=1)+geom_vline(xintercept=bds$cutoffswt, colour="blue", linetype = "longdash")+
    scale_y_continuous(breaks=NULL)+scale_fill_manual(values=c("#E74C3C", "#1ABC9C"),guide=guide_legend(title=""))+labs(y="Density",title="Distribution of # Switches, criterion I",x="") +
    theme_bw()+theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),legend.key = element_blank())

  # type II
  df2 = data.frame(tag = bds$boolbad2, mmrate = bds$meanmrate, variance = bds$variance, nswitches = bds$nswitches, longestblock = bds$longestblock)
  p2 = ggplot(subset(df2,tag),aes(x=nswitches))+geom_histogram(aes(x=nswitches,y=-..density..,fill="predicted artifact"),binwidth=1)+geom_histogram(data=subset(df2,!tag),aes(x=nswitches,y=..density..,fill="predicted true"),binwidth=1)+geom_vline(xintercept=bds$cutoffswt, colour="blue", linetype = "longdash")+
    scale_y_continuous(breaks=NULL)+scale_fill_manual(values=c("#E74C3C", "#1ABC9C"),guide=guide_legend(title=""))+labs(y="Density",title="Distribution of # Switches, criterion II",x="") +
    theme_bw()+theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),legend.key = element_blank())

  # all
  df3 = data.frame(nswitches=bds$nswitches)
  p3 = ggplot(df3,aes(x=nswitches,y=..density..))+geom_histogram(binwidth=1)+geom_vline(xintercept=bds$cutoffswt, colour="blue", linetype = "longdash")+
    scale_y_continuous(breaks=NULL)+scale_fill_manual(values=c("#E74C3C", "#1ABC9C"),guide=guide_legend(title=""))+labs(y="Density",title=paste0("Distribution of # Switches, all metabolic features"),x=paste0("cutoff =",bds$cutoffswt)) +
    theme_bw()+theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),legend.key = element_blank())

  multiplot(p3, p1,p2,cols=1)
}

plot_meanmrate <- function(bds){
  # mean missing rate distribution
  # type I
  df1 = data.frame(tag = bds$boolbad1, mmrate = bds$meanmrate, variance = bds$variance, nswitches = bds$nswitches, longestblock = bds$longestblock)
  p1 = ggplot(subset(df1,tag),aes(x=mmrate))+geom_histogram(aes(x=mmrate,y=-..density..,fill="predicted artifact"))+geom_histogram(data=subset(df1,!tag),aes(x=mmrate,y=..density..,fill="predicted true"))+geom_vline(xintercept=bds$cutoffmr, colour="blue", linetype = "longdash")+
    scale_y_continuous(breaks=NULL)+scale_x_continuous(breaks=c(0.00,0.5,1))+scale_fill_manual(values=c("#E74C3C", "#1ABC9C"),guide=guide_legend(title=""))+labs(y="Density",title="Distribution of Mean Missing Rate, criterion I",x="") +
    theme_bw()+theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),legend.key = element_blank())

  # type II
  df2 = data.frame(tag = bds$boolbad2, mmrate = bds$meanmrate, variance = bds$variance, nswitches = bds$nswitches, longestblock = bds$longestblock)
  p2 = ggplot(subset(df2,tag),aes(x=mmrate))+geom_histogram(aes(x=mmrate,y=-..density..,fill="predicted artifact"))+geom_histogram(data=subset(df2,!tag),aes(x=mmrate,y=..density..,fill="predicted true"))+geom_vline(xintercept=bds$cutoffmr, colour="blue", linetype = "longdash")+
    scale_y_continuous(breaks=NULL)+scale_x_continuous(breaks=c(0.00,0.5,1))+scale_fill_manual(values=c("#E74C3C", "#1ABC9C"),guide=guide_legend(title=""))+labs(y="Density",title="Distribution of Mean Missing Rate, criterion II",x="") +
    theme_bw()+theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),legend.key = element_blank())

  # all
  df3 = data.frame(mmrate=bds$meanmrate)
  p3 = ggplot(df3,aes(x=mmrate,y=..density..))+geom_histogram()+geom_vline(xintercept=bds$cutoffmr, colour="blue", linetype = "longdash")+
  scale_y_continuous(breaks=NULL)+scale_x_continuous(breaks=c(0.00,0.5,1))+scale_fill_manual(values=c("#E74C3C", "#1ABC9C"),guide=guide_legend(title=""))+labs(y="Density",title=paste0("Distribution of Mean Missing Rate, all metabolic features"),x=paste0("cutoff =",bds$cutoffmr)) +
  theme_bw()+theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),legend.key = element_blank())

  multiplot(p3, p1,p2,cols=1)
}








