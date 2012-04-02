##################################################################################
#' Select plot for a NpdeObject object
#' 
#' Select plot for a NpdeObject object
#' 
#' @usage npde.plot.select(npdeObject,data=FALSE,ecdf=FALSE,qqplot=FALSE, histogram=FALSE,x.scatter=FALSE,pred.scatter=FALSE,x.box=FALSE,pred.box=FALSE, cov.x.scatter=FALSE,cov.pred.scatter=FALSE,cov.x.box=FALSE,cov.pred.box=FALSE, cov.ecdf=FALSE, vpc=FALSE,...)
#' @param npdeObject an object returned by a call to \code{\link{npde}} or \code{\link{autonpde}}
#' @param data boolean, whether to produce a plot of the data
#' @param ecdf boolean, whether to produce a distribution plot of the empirical distribution function
#' @param qqplot boolean, whether to produce a QQ-plot of the empirical distribution function
#' @param histogram boolean, whether to produce a histogram of the metric
#' @param x.scatter boolean, whether to produce a scatterplot of the metric as a function of X
#' @param pred.scatter boolean, whether to produce a scatterplot of the metric as a function of predictions
#' @param x.box boolean, whether to produce whisker plots of the metric as a function of X
#' @param pred.box boolean, whether to produce whisker plots of the metric as a function of predictions
#' @param cov.x.scatter boolean, whether to produce a scatterplot of the metric as a function of X, split by covariate(s)
#' @param cov.pred.scatter boolean, whether to produce a scatterplot of the metric as a function of predictions, split by covariate(s)
#' @param cov.x.box boolean, whether to produce whisker plots of the metric as a function of X, split by covariate(s)
#' @param cov.pred.box boolean, whether to produce whisker plots of the metric as a function of predictions, split by covariate(s)
#' @param cov.ecdf boolean, whether to produce a distribution plot of the empirical distribution function, split by covariate(s)
#' @param vpc boolean, whether to produce a VPC
#' @param \dots additional arguments to be passed on to the function, to control which metric (npde, pd, npd) is used or to override graphical parameters (see the PDF document for details, as well as \code{\link{set.npde.plotoptions}})
#' @author Emmanuelle Comets <emmanuelle.comets@@bichat.inserm.fr>
#' @seealso \code{\link{npde}}, \code{\link{autonpde}}, \code{\link{set.npde.plotoptions}}
#' @references K. Brendel, E. Comets, C. Laffont, C. Laveille, and F.
#' Mentre. Metrics for external model evaluation with an application to the
#' population pharmacokinetics of gliclazide. \emph{Pharmaceutical Research},
#' 23:2036--49, 2006.
#' @keywords plot internal

npde.plot.select<-function(npdeObject,data=FALSE,ecdf=FALSE,qqplot=FALSE, histogram=FALSE,x.scatter=FALSE,pred.scatter=FALSE,x.box=FALSE,pred.box=FALSE, cov.x.scatter=FALSE,cov.pred.scatter=FALSE,cov.x.box=FALSE,cov.pred.box=FALSE, cov.ecdf=FALSE, vpc=FALSE,...) {
# Function selecting which plots are to be drawn
  namObj<-deparse(substitute(npdeObject))
  interactive<-npdeObject["prefs"]$interactive
# ECO TODO: replace with partial matching
  if(data) plot(npdeObject,plot.type="data",...)
  if(ecdf) plot(npdeObject,plot.type="ecdf",...)
  if(qqplot) plot(npdeObject,plot.type="qqplot",...)
  if(histogram) plot(npdeObject,plot.type="histogram",...)
  if(x.scatter) plot(npdeObject,plot.type="x.scatter",...)
  if(pred.box) plot(npdeObject,plot.type="pred.scatter",box=TRUE,...)
  if(x.box) plot(npdeObject,plot.type="x.scatter",box=TRUE,...)
  if(pred.scatter) plot(npdeObject,plot.type="pred.scatter",...)
  if(cov.x.scatter) plot(npdeObject,plot.type="cov.x.scatter",...)
  if(cov.pred.scatter) plot(npdeObject,plot.type="cov.pred.scatter",...)
  if(cov.ecdf) plot(npdeObject,plot.type="cov.ecdf",...)
  if(cov.x.box) plot(npdeObject,plot.type="cov.x.scatter",box=TRUE,...)
  if(cov.pred.box) plot(npdeObject,plot.type="cov.pred.scatter",box=TRUE,...)
  if(vpc) plot(npdeObject,plot.type="vpc",...)
}

#' Default plots for a NpdeObject object
#' 
#' Default plots for a NpdeObject object
#' 
#' @usage default.npde.plots(npdeObject, ...)
#' @param npdeObject an object returned by a call to \code{\link{npde}} or \code{\link{autonpde}}
#' @param \dots additional arguments to be passed on to the function, to control which metric (npde, pd, npd) is used or to override graphical parameters (see the PDF document for details, as well as \code{\link{set.npde.plotoptions}})
#' @keywords plot internal

#### Meta-niveau
default.npde.plots<-function(npdeObject,...) {
# When plot(npdeObject) is called without plot.type  
  par(mfrow=c(2,2),ask=npdeObject["prefs"]$ask)
  npde.plot.select(npdeObject,qqplot=TRUE,histogram=TRUE, x.scatter=TRUE,pred.scatter=TRUE,new=FALSE,...)
}

#' Covariate plots for a NpdeObject object
#' 
#' Covariate plots for a NpdeObject object
#' 
#' @usage npde.plot.covariates(npdeObject, which="x", ...)
#' @param npdeObject an object returned by a call to \code{\link{npde}} or \code{\link{autonpde}}
#' @param which one of "x" (scatterplots of the metric versus X), "pred" (scatterplots of the metric versus predictions) or "ecdf" (empirical distribution function)
#' @param \dots additional arguments to be passed on to the function, to control which metric (npde, pd, npd) is used or to override graphical parameters (see the PDF document for details, as well as \code{\link{set.npde.plotoptions}})
#' @keywords plot internal

npde.plot.covariates<-function(npdeObject,which="x",...) {
# Parameters or random effects versus covariates
  if(which=="x") {
    plot(npdeObject,plot.type="cov.x.scatter",...)
  }
  if(which=="pred") {
    plot(npdeObject,plot.type="cov.pred.scatter",...)
  }
  if(which=="ecdf") {
    plot(npdeObject,plot.type="cov.ecdf",...)
  }
}

#' Plots for pd and npde
#' 
#' Plots for pd and npde
#' 
#' @aliases npde.plot.pd npde.plot.npde
#' @usage npde.plot.pd(npdeObject, ...)
#' @usage npde.plot.npde(npdeObject, ...)
#' @param npdeObject an object returned by a call to \code{\link{npde}} or \code{\link{autonpde}}
#' @param \dots additional arguments to be passed on to the function, to control which metric (npde, pd, npd) is used or to override graphical parameters (see the PDF document for details, as well as \code{\link{set.npde.plotoptions}})
#' @keywords plot internal

npde.plot.pd<-function(npdeObject,...) {
# Advanced goodness of fit plots
  if(npdeObject@options$verbose) cat("Plots for pd\n")
  default.npde.plots(npdeObject,which="pd",...)
}

npde.plot.npde<-function(npdeObject,...) {
# Advanced goodness of fit plots
	if(npdeObject@options$verbose) cat("Plots for npde\n")
  default.npde.plots(npdeObject,...)
}

################################    Data    #####################################


#' Plot functions for NpdeObjects
#' 
#' Workhorse functions for plots (internal)
#' 
#' @aliases npde.plot.data npde.plot.dist npde.plot.scatter npde.plot.vpc npde.plot.loq aux.scatter.box aux.scatter compute.bands.true compute.bands
#' @param npdeObject an object returned by a call to \code{\link{npde}} or \code{\link{autonpde}}
#' @param \dots additional arguments to be passed on to the function, to control which metric (npde, pd, npd) is used or to override graphical parameters (see the PDF document for details, as well as \code{\link{set.npde.plotoptions}})
#' @keywords plot internal
									
npde.plot.data<-function(npdeObject,...) {
    if(length(npdeObject["data"]["icens"])>0) {
      has.cens<-TRUE
      icens<-npdeObject["data"]["icens"]
    } else has.cens<-FALSE
    plot.opt<-npdeObject["prefs"]
    plot.opt$xlab<-paste(npdeObject["data"]["name.predictor"]," (",npdeObject["data"]["units"]$x,")",sep="")
    plot.opt$ylab<-paste(npdeObject["data"]["name.response"]," (",npdeObject["data"]["units"]$y,")",sep="")
    plot.opt<-replace.npde.plotoptions(plot.opt,...)
    if(plot.opt$new) {
      mfrow<-plot.opt$mfrow
      if(length(mfrow)==0) mfrow<-c(1,1)
      par(mfrow=mfrow,ask=plot.opt$ask)
    }
    logtyp<-""
    if(plot.opt$xlog) logtyp<-paste(logtyp,"x",sep="")
    if(plot.opt$ylog) logtyp<-paste(logtyp,"y",sep="")
    plot.type<-plot.opt$type
    x<-npdeObject["data"]["data"][,npdeObject["data"]["name.predictor"]]
    y<-npdeObject["data"]["data"][,npdeObject["data"]["name.response"]]
    if(plot.opt$impute.loq & length(npdeObject["results"]["res"]$ycomp)>0)
    	y<-npdeObject["results"]["res"]$ycomp
    idobs<-npdeObject["data"]["data"][,npdeObject["data"]["name.group"]]
    if(is.null(plot.opt$xlim) & !plot.opt$xlog) plot.opt$xlim<-c(min(x,na.rm=T), max(x,na.rm=T))
    if(is.null(plot.opt$ylim) & !plot.opt$ylog) plot.opt$ylim<-c(min(y,na.rm=T), max(y,na.rm=T))
    if(plot.type=="p" | plot.type=="b") {
    if(has.cens) {
      plot(x,y,xlab=plot.opt$xlab,ylab=plot.opt$ylab,col=plot.opt$col, pch=plot.opt$pch,log=logtyp,xlim=plot.opt$xlim,ylim=plot.opt$ylim, main=plot.opt$main,cex=plot.opt$cex,cex.axis=plot.opt$cex.axis, cex.lab=plot.opt$cex.lab,type="n") 
      points(x[-c(icens)],y[-c(icens)],col=plot.opt$col,pch=plot.opt$pch, cex=plot.opt$cex) 
      if(plot.opt$plot.loq) points(x[icens],y[icens],col=plot.opt$cens.col, pch=plot.opt$cens.pch, cex=plot.opt$cex) 
    } else plot(x,y,xlab=plot.opt$xlab, ylab=plot.opt$ylab,col=plot.opt$col,pch=plot.opt$pch,log=logtyp,xlim=plot.opt$xlim, ylim=plot.opt$ylim,main=plot.opt$main,cex=plot.opt$cex,cex.axis=plot.opt$cex.axis, cex.lab=plot.opt$cex.lab) 
    }
    if(plot.type=="l") {
      plot(x,y,xlab=plot.opt$xlab, ylab=plot.opt$ylab,col=plot.opt$col,lty=plot.opt$lty,lwd=plot.opt$lwd,type="n", log=logtyp,xlim=plot.opt$xlim,ylim=plot.opt$ylim,main=plot.opt$main, cex=plot.opt$cex,cex.axis=plot.opt$cex.axis, cex.lab=plot.opt$cex.lab)
    }
    if(plot.type=="l" | plot.type=="b") {
      for(isuj in unique(idobs)) {
        lines(x[idobs==isuj], y[idobs==isuj],col=plot.opt$col,lty=plot.opt$lty,lwd=plot.opt$lwd)
      }
      if(has.cens & plot.opt$plot.loq) {
        points(x[icens],y[icens],col=plot.opt$cens.col,pch=plot.opt$cens.pch, cex=plot.opt$cex) 
      }    
    }
    if(plot.opt$line.loq & length(npdeObject["data"]["loq"])>0) abline(h=npdeObject["data"]["loq"],lty=plot.opt$ablinelty, col=plot.opt$ablinecol,lwd=plot.opt$ablinelwd)
}

#################    	     Distributions - npde/pd          ###################

npde.plot.dist<-function(npdeObject,which="npde",shape="ecdf",covsplit=FALSE, ...) {
	plot.opt.defaults<-function(plot.opt,plot.type="ecdf",which="npde") {
		plot.opt$new<-TRUE
		plot.opt$main<-""
		if(plot.type=="ecdf") {
			plot.opt$ylab<-"Empirical cdf"
			plot.opt$ylim<-c(0,1)
			plot.opt$xlab<-switch(which,pd="Sample quantiles (pd)",npd="Sample quantiles (npd)",npde="Sample quantiles (npde)")
		}
		if(plot.type=="qqplot") {
			plot.opt$ylab<-"Theoretical Quantiles"
			plot.opt$xlab<-switch(which,pd="Sample quantiles (pd)",npd="Sample quantiles (npd)",npde="Sample quantiles (npde)")
			if(which %in% c("npd","npde")) {
				plot.opt$main<-paste("Q-Q plot versus N(0,1) for ",which)
			} else {
				plot.opt$main<-"Q-Q plot versus U(0,1) for pd"
			}
		}
		if(plot.type=="hist") {
			plot.opt$ylab<-"Frequency"
			plot.opt$xlab<-switch(which,pd="Sample quantiles (pd)",npd="Sample quantiles (npd)",npde="Sample quantiles (npde)")
			plot.opt$xlim<-NULL
			plot.opt$ylim<-NULL
		}	
		return(plot.opt)
	}
	
	if(match(which,c("npde","pd","npd"),nomatch=0)==0) {
		cat("Option which=",which,"not recognised\n")
		return()
	}
	if(match(shape,c("ecdf","qqplot","hist"),nomatch=0)==0) {
		cat("Option shape=",shape,"not recognised\n")
		return()
	}
	if(which=="npde" & length(npdeObject["results"]["res"]$npde)==0)  {
		cat("    Missing npde object to plot.\n")
		return()
	}
	if(which %in% c("pd","npd") & length(npdeObject["results"]["res"]$pd)==0) {
		cat("    Missing pd object to plot.\n")
		return()
	}
	if(which=="pd") distrib<-"unif" else distrib<-"norm"
	if(length(npdeObject["data"]["icens"])>0) has.cens<-TRUE else has.cens<-FALSE
	plot.opt<-npdeObject["prefs"]
	if(covsplit & length(npdeObject["data"]["name.covariates"])==0) {
		if(npdeObject@options$verbose) cat("No covariates in the dataset\n")
		covsplit<-FALSE
	}
	plot.opt<-plot.opt.defaults(plot.opt,plot.type=shape,which=which)
	change.main<-FALSE
	main<-plot.opt$main
	plot.opt<-replace.npde.plotoptions(plot.opt,...)
	if(plot.opt$main!=main) change.main<-TRUE
	if(covsplit) {
		if(plot.opt$which.cov=="") plot.opt$which.cov<-"all"
		if(is.numeric(plot.opt$which.cov)) plot.opt$which.cov<-npdeObject["data"]["name.covariates"][plot.opt$which.cov]
		if(plot.opt$which.cov=="all") lcov<-npdeObject["data"]["name.covariates"] else {
			icov<-match(plot.opt$which.cov,npdeObject["data"]["name.covariates"])
			lcov<-npdeObject["data"]["name.covariates"][icov]
		}
		ncov<-length(lcov)
		if(ncov==0) {
			if(npdeObject@options$verbose) cat("Cannot find covariate(s)", plot.opt$which.cov,"in the dataset\n")
			covsplit<-FALSE
		}
	} else ncov<-0
	if(!covsplit & plot.opt$new) {
		mfrow<-plot.opt$mfrow
		if(length(mfrow)==0) mfrow=c(1,1)
		par(mfrow=mfrow,ask=plot.opt$ask)
	}
	logtyp<-""
	if(plot.opt$xlog) logtyp<-paste(logtyp,"x",sep="")
	if(plot.opt$ylog) logtyp<-paste(logtyp,"y",sep="")
	plot.opt$logtyp<-logtyp
	sim.ypl<-NULL
	if(!plot.opt$approx.pi) {
		if(which %in% c("pd","npd")) {
			if(length(npdeObject["results"]["pd.sim"])==0) {
				cat("You have requested to use simulations to provide the prediction intervals, but the simulated pd are not present.\n")
				plot.opt$approx.pi<-TRUE
			} else sim.ypl<-npdeObject["results"]["pd.sim"]
		}
		if(which=="npde") {
			if(length(npdeObject["results"]["npde.sim"])==0) {
				cat("You have requested to use simulations to provide the prediction intervals, but the simulated npde are not present.\n")
				plot.opt$approx.pi<-TRUE
			} else sim.ypl<-npdeObject["results"]["npde.sim"]
		}
	}
	if(which=="npde") ypl<-npdeObject["results"]["res"]$npde
	if(which %in% c("pd","npd")) ypl<-npdeObject["results"]["res"]$pd    
	if(which=="npd") ypl<-qnorm(ypl)
	if(has.cens) plmat<-data.frame(xpd=ypl, cens=npdeObject["data"]["data"][,npdeObject["data"]["name.cens"]]) else plmat<-data.frame(xpd=ypl,cens=rep(0,length(ypl)))
	keep<-npdeObject["data"]["not.miss"]
	idobs<-npdeObject["data"]["data"][,npdeObject["data"]["name.group"]]
	idobs<-idobs[keep]
	plmat<-plmat[keep,]
	
	if(covsplit) {
		for(icov in 1:ncov) {
			namcov<-lcov[icov]
			zecov<-npdeObject["data"]["data"][keep,namcov]
			idobs<-npdeObject["data"]["data"][keep,npdeObject["data"]["name.group"]]
			ucov<-zecov[match(unique(idobs),idobs)]
			if(!is.numeric(ucov) & length(unique(ucov))<=4) {
				covcont<-FALSE
				ncat<-length(unique(ucov))
				namcat<-paste(namcov,sort(unique(zecov)),sep="=")
			} else {
				covcont<-TRUE
				zecov<-cut(zecov,breaks=quantile(ucov,c(0,1/3,2/3,1)), include.lowest=TRUE, ordered_result=TRUE)
				ncat<-3
				namcat<-paste(namcov,c("<Q1","Q1-Q3",">Q3"),sep=": ")
			}
			if(plot.opt$new) {
				mfrow<-plot.opt$mfrow
				if(length(mfrow)==0) {
					if(covcont) mfrow=c(1,3) else {
						n1<-round(sqrt(ncat))
						n2<-ceiling(ncat/n1)
						mfrow<-c(n1,n2)
					}
				}
				par(mfrow=mfrow,ask=plot.opt$ask)
			}
			zecat<-sort(unique(zecov))
			plmat1<-data.frame(plmat,cov=zecov)
			for(ic in 1:length(zecat)) {
				if(shape=="hist") 
					aux.plot.hist(plmat1[plmat1$cov==zecat[ic],], plot.opt,distrib=distrib, nclass=plot.opt$vpc.bin, sim.ypl=sim.ypl[plmat1$cov==zecat[ic],]) else
						aux.plot.dist(plmat1[plmat1$cov==zecat[ic],],plot.opt, shape=shape,distrib=distrib,nrep=npdeObject["sim.data"]["nrep"], ties=npdeObject["options"]$ties, sim.ypl=sim.ypl[plmat1$cov==zecat[ic],])
			}
		}
	} else {
		if(shape=="hist") aux.plot.hist(plmat,plot.opt,distrib=distrib, nclass=plot.opt$vpc.bin,sim.ypl=sim.ypl) else aux.plot.dist(plmat,plot.opt, shape=shape,distrib=distrib,nrep=npdeObject["sim.data"]["nrep"], ties=npdeObject["options"]$ties, sim.ypl=sim.ypl)
	}
}

#################    	     Scatterplots - npde/pd           ###################

npde.plot.scatter<-function(npdeObject,which="npde",xaxis="x",covsplit=FALSE, ...) {
	if(match(which,c("npde","pd","npd"),nomatch=0)==0) {
		cat("Option which=",which,"not recognised\n")
		return()
	}
	if(xaxis=="cov" & length(npdeObject["data"]["name.covariates"])==0) {
		cat("Option xaxis set to 'cov' but no covariate in dataset\n")
		return()
	}
	if(which=="npde" & length(npdeObject["results"]["res"]$npde)==0) {
		cat("    Missing npde object to plot.\n")
		return()
	}
	if(which %in% c("pd","npd") & length(npdeObject["results"]["res"]$pd)==0)  {
		cat("    Missing pd object to plot.\n")
		return()
	}    
	if(length(npdeObject["data"]["icens"])>0) has.cens<-TRUE else has.cens<-FALSE
	plot.opt<-npdeObject["prefs"]
	if(xaxis=="cov") {
# 		plot.opt$box<-TRUE
 		if(!covsplit & npdeObject@options$verbose) cat("    graph versus covariates requested, setting covsplit to TRUE\n")
 		covsplit<-TRUE
	}
	if(covsplit & length(npdeObject["data"]["name.covariates"])==0) {
		if(npdeObject@options$verbose) cat("No covariates in the dataset\n")
		covsplit<-FALSE
	}
	if(which=="pd") plot.opt$ylim<-c(0,1) else plot.opt$ylim<-NULL
	if(which=="pd") distrib<-"unif" else distrib<-"norm"
	xlab1<-plot.opt$xlab<-switch(xaxis, x=paste(npdeObject["data"]["name.predictor"]," (", npdeObject["data"]["units"]$x,")",sep=""), pred=paste("Predicted ", npdeObject["data"]["name.response"]," (",npdeObject["data"]["units"]$y,")",sep=""), cov="")
	plot.opt$ylab<-switch(which,pd="pd",npd="npd",npde="npde")
	plot.opt$main<-""
	plot.opt<-replace.npde.plotoptions(plot.opt,...)
	change.xlab<-FALSE
	if(xlab1!=plot.opt$xlab) change.xlab<-TRUE
	if(covsplit) {
		if(plot.opt$which.cov=="") plot.opt$which.cov<-"all"
		if(is.numeric(plot.opt$which.cov)) plot.opt$which.cov<-npdeObject["data"]["name.covariates"][plot.opt$which.cov]
		if(plot.opt$which.cov=="all") lcov<-npdeObject["data"]["name.covariates"] else {
			icov<-match(plot.opt$which.cov,npdeObject["data"]["name.covariates"])
			lcov<-npdeObject["data"]["name.covariates"][icov]
		}
		ncov<-length(lcov)
		if(ncov==0) {
			cat("Cannot find covariate",plot.opt$which.cov,"in the dataset\n")
			covsplit<-FALSE
		}
	} else ncov<-0
	if(!covsplit & plot.opt$new) {
		mfrow<-plot.opt$mfrow
		if(length(mfrow)==0) mfrow=c(1,1)
		par(mfrow=mfrow,ask=plot.opt$ask)
	}
	if(xaxis=="cov" & length(plot.opt$mfrow)==0 & plot.opt$new) {
		n1<-round(sqrt(ncov))
		n2<-ceiling(ncov/n1)
		plot.opt$mfrow<-c(n1,n2)
		par(mfrow=plot.opt$mfrow,ask=plot.opt$ask)
	} 
	logtyp<-""
	if(plot.opt$xlog) logtyp<-paste(logtyp,"x",sep="")
	if(plot.opt$ylog) logtyp<-paste(logtyp,"y",sep="")
	plot.opt$logtyp<-logtyp
	sim.ypl<-NULL
	if(!plot.opt$approx.pi) {
		if(which %in% c("pd","npd")) {
			if(length(npdeObject["results"]["pd.sim"])==0) {
				cat("You have requested to use simulations to provide the prediction intervals, but the simulated pd are not present.\n")
				plot.opt$approx.pi<-TRUE
			} else sim.ypl<-npdeObject["results"]["pd.sim"]
		}
		if(which=="npde") {
			if(length(npdeObject["results"]["npde.sim"])==0) {
				cat("You have requested to use simulations to provide the prediction intervals, but the simulated npde are not present.\n")
				plot.opt$approx.pi<-TRUE
			} else sim.ypl<-npdeObject["results"]["npde.sim"]
		}
	}
	if(which=="npde") ypl<-npdeObject["results"]["res"]$npde
	if(which %in% c("pd","npd")) ypl<-npdeObject["results"]["res"]$pd
	if(which=="npd") {
		ypl<-qnorm(ypl)
		if(!plot.opt$approx.pi) sim.ypl<-qnorm(sim.ypl)
	}
	if(xaxis=="x") xpl<-npdeObject["data"]["data"][, npdeObject["data"]["name.predictor"]]
	idobs<-npdeObject["data"]["data"][npdeObject["data"]["not.miss"], npdeObject["data"]["name.group"]]
	if(xaxis=="cov") xpl<-rep(1,length(ypl))
	if(xaxis=="pred") xpl<-npdeObject["results"]["res"]$ypred
	if(has.cens) plmat<-matrix(c(xpl,npdeObject["data"]["data"][, npdeObject["data"]["name.cens"]],ypl),ncol=3) else plmat<-matrix(c(xpl,rep(0,length(xpl)),ypl),ncol=3)
	plmat<-plmat[npdeObject["data"]["not.miss"],]
	if(which %in% c("npd","npde")) {
		x1<-qnorm(0.05)
		dotline<-c(-x1,0,x1)
	}
	if(which=="pd") dotline<-c(0.05,0.5,0.95)
	if(covsplit) {
		# Same limits for the different plots
		if(is.null(plot.opt$xlim) & xaxis!="cov") plot.opt$xlim<-c(min(xpl,na.rm=T), max(xpl,na.rm=T))
		if(is.null(plot.opt$ylim) & xaxis!="cov") plot.opt$ylim<-c(min(ypl,na.rm=T), max(ypl,na.rm=T))
		for(icov in 1:ncov) {
			namcov<-lcov[icov]
			zecov<-npdeObject["data"]["data"][npdeObject["data"]["not.miss"],namcov]
			if(xaxis=="cov") plmat[,1]<-zecov
			ucov<-zecov[match(unique(idobs),idobs)]
			# ECO: SECURISER 
			if(!is.numeric(ucov)) {
				if(length(unique(ucov))>4) cat("Too many categories, the plot will not be informative.\n")
				plot.box<-TRUE
				covcont<-FALSE
				ncat<-length(unique(ucov))
				namcat<-paste(namcov,sort(unique(zecov)),sep="=")
			} else {
				plot.box<-plot.opt$box
				covcont<-TRUE
				zecov<-cut(zecov,breaks=quantile(ucov,c(0,1/3,2/3,1)), include.lowest=TRUE, ordered_result=TRUE)
				ncat<-3
				namcat<-paste(namcov,c("<Q1","Q1-Q3",">Q3"),sep=": ")
			}
			if(plot.opt$new & xaxis!="cov") {
				mfrow<-plot.opt$mfrow
				if(length(mfrow)==0) {
					if(covcont) mfrow=c(1,3) else {
						n1<-round(sqrt(ncat))
						n2<-ceiling(ncat/n1)
						mfrow<-c(n1,n2)
					}
				}
				par(mfrow=mfrow,ask=plot.opt$ask)
			}
			zecat<-sort(unique(zecov))
			plmat1<-data.frame(plmat,zecov)
			if(xaxis=="cov") {
				if(plot.box) {
 				plmat1<-plmat1[,c(4,2,3)]
 				if(!change.xlab) plot.opt$xlab<-namcov
 				aux.scatter.box(plmat1,plot.opt,dotline,main=NULL,namcat=namcat)
				} else {
					aux.scatter(plmat,plot.opt,dotline,distrib=distrib,main=NULL, sim.ypl=sim.ypl)					
				}
			} else {
				for(ic in 1:length(zecat)) {
					if(plot.opt$box) 
						aux.scatter.box(plmat1[plmat1[,4]==zecat[ic],], plot.opt,dotline, main=namcat[ic]) else 
						aux.scatter(plmat1[plmat1[,4]==zecat[ic],],plot.opt,dotline, distrib=distrib,main=namcat[ic], sim.ypl=sim.ypl[plmat1[,4]==zecat[ic],])
				}
			}
		}
	} else {
		if(plot.opt$box) aux.scatter.box(plmat,plot.opt,dotline,main=NULL) else aux.scatter(plmat,plot.opt,dotline,distrib=distrib,main=NULL,sim.ypl=sim.ypl)
	}
}

######################	 Computational functions ################################

# Compute PI bands
compute.bands<-function(nsamp,nseuil=200,quant=c(0.05,0.5,0.95),distrib="norm", alpha=0.95) {
# Compute a prediction interval around selected quantiles of the normal or uniform distribution, for different sizes of samples by randomly sampling from N(0,1) or U(0,1)
### nsamp: a vector giving the sizes of the samples
### size: alpha (defaults to a 95% PI)
### quantile: quant (defaults to 5th, 50th (median) and 95th percentiles)
### distribution: normal (default) or uniform
# When the number of samples isamp is larger than 200, the PI is computed for n=200 and the size of the PI is then adjusted through sqrt(200/isamp)

# Returns
### binf: lower bounds (as many columns as elements in quant) for the PI with level alpha
### bmed: median
### bsup: upper bounds
#  msim<-10000
  msim<-1000 # number of replications used to compute the prediction interval
  idx1<-which(nsamp>nseuil)
  quant.pi<-c((1-alpha)/2,0.5,1-(1-alpha)/2) 
  if(length(idx1)>0) {
     xsamp<-matrix(switch(distrib,norm=rnorm(msim*nseuil),unif=runif(msim*nseuil)), ncol=msim)
     mat<-apply(xsamp,2,quantile,quant)
     xseuil<-apply(mat,1,quantile,quant.pi)
     demi<-(xseuil[3,]-xseuil[1,])/2
  }
  binf<-bsup<-bmed<-matrix(nrow=length(nsamp),ncol=length(quant), dimnames=list(nsamp,quant))
  for(isamp in unique(nsamp)) {
   if(isamp<nseuil) {
     xsamp<-matrix(switch(distrib,norm=rnorm(msim*isamp),unif=runif(msim*isamp)), ncol=msim)
     mat<-apply(xsamp,2,quantile,quant)
     xtab<-apply(mat,1,quantile,quant.pi)
   } else {
     xtab<-matrix(nrow=3,ncol=length(quant))
     xtab[2,]<-switch(distrib,norm=qnorm(quant),unif=quant)
     xtab[1,]<-xtab[2,]-demi*sqrt(nseuil/isamp)
     xtab[3,]<-xtab[2,]+demi*sqrt(nseuil/isamp)
   }
   for(i in which(nsamp==isamp)) {
     binf[i,]<-xtab[1,]
     bmed[i,]<-xtab[2,]
     bsup[i,]<-xtab[3,]
   }
  }
  return(list(binf=binf,bsup=bsup,bmed=bmed))
}

compute.bands.true<-function(sim.ypl,quant=c(0.05,0.5,0.95), distrib="norm", alpha=0.95) {
	# Compute a prediction interval around selected quantiles of the normal or uniform distribution, for different sizes of samples by using the data
	quant.pi<-c((1-alpha)/2,0.5,1-(1-alpha)/2)
	nsim<-dim(sim.ypl)[2]-2
	grp<-unique(sim.ypl[,1])
	binf<-bsup<-bmed<-matrix(nrow=length(grp),ncol=length(quant), dimnames=list(grp,quant))
	for(igrp in grp) {
		mat<-apply(sim.ypl[sim.ypl$grp==igrp,3:dim(sim.ypl)[2]],2,quantile,quant)
		xtab<-apply(mat,1,quantile,quant.pi)
		binf[grp==igrp,]<-xtab[1,]
		bmed[grp==igrp,]<-xtab[2,]
		bsup[grp==igrp,]<-xtab[3,]
	}
	return(list(binf=binf,bsup=bsup,bmed=bmed))
}

# boxplot
aux.scatter.box<-function(plmat,plot.opt,dotline,distrib="norm",main=NULL,namcat=NULL,alpha=0.95) {
	# Auxiliary function used to plot boxplots
	if(!is.null(namcat)) xname<-namcat
	if(is.null(main)) main<-plot.opt$main
	if(is.null(dotline)) dotline<-c(alpha,0.5,1-alpha)
	if(is.numeric(plmat[,1])) {
		plot.opt$xlim<-c(min(plmat[,1]),max(plmat[,1]))
		xbin<-npde.binning(plmat[,1],plot.opt)
		plmat<-data.frame(grp=xbin$xgrp,plmat)
		xcent<-xbin$xat
		if(is.null(namcat)) xname<-format(xbin$xat,digits=2)
	} else {
		xcent<-1:length(unique(plmat[,1]))
		plmat<-data.frame(grp=as.factor(plmat[,1]),plmat)
	}
	boxplot(plmat[,4]~plmat[,1],at=xcent,pch=plot.opt$pch,cex=plot.opt$cex, varwidth=plot.opt$varwidth,boxwex=plot.opt$boxwex, col=plot.opt$fillcol, names=xname,xlab=plot.opt$xlab,ylab=plot.opt$ylab, main=main,cex.lab=plot.opt$cex.lab,xlim=plot.opt$xlim, ylim=plot.opt$ylim)
	if(length(dotline)>0) {
		for(i in dotline) abline(h=i,lty=plot.opt$ablinelty,lwd=plot.opt$ablinelwd, col=plot.opt$ablinecol)
	}
}

# Scatterplot
aux.scatter<-function(plmat,plot.opt,dotline,distrib="norm",main=NULL,namcat=NULL,sim.ypl=NULL) {
	# Auxiliary function used to plot scatterplots
	alpha<-(1-plot.opt$vpc.interval)/2
	if(is.null(main)) main<-plot.opt$main
	if(is.null(dotline)) dotline<-c(alpha,0.5,1-alpha)
	if(plot.opt$bands) {
		if(is.numeric(plmat[,1])) {
			# ECO TODO: should binning be done on full data or after censoring ?
			xbin<-npde.binning(plmat[,1],plot.opt)
			plmat<-data.frame(grp=xbin$xgrp,plmat)
			xcent<-xbin$xat
		} else {
			plmat<-data.frame(grp=plmat[,1],plmat)
			xcent<-unique(plmat[,2])
		}
		nseuil<-200
		alp.pi<-plot.opt$pi.size
		if(alp.pi<0.5) alp.pi<-(1-alp.pi)
		quant<-c(alpha,0.5,1-alpha)
		if(plot.opt$approx.pi) {
			nobs<-tapply(plmat$grp,plmat$grp,length)
			bnds<-compute.bands(nobs,nseuil,quant,distrib,alp.pi)
		} else {
			sim.ypl<-data.frame(grp=plmat$grp,cens=plmat[,3],sim.ypl)
			bnds<-compute.bands.true(sim.ypl,quant,distrib=distrib,alpha=alp.pi)
		}
		binf<-bnds$binf
		bsup<-bnds$bsup
		bmed<-bnds$bmed
		if(is.null(plot.opt$ylim)) plot.opt$ylim<-c(min(plmat[,4],binf,na.rm=T), max(plmat[,4],bsup,na.rm=T))
		plot(plmat[,2],plmat[,4],type="n",xlab=plot.opt$xlab,ylab=plot.opt$ylab, main=main,col=plot.opt$col,pch=plot.opt$pch, cex=plot.opt$cex, cex.lab=plot.opt$cex.lab,ylim=plot.opt$ylim)
		for(icol in 1:3) polygon(c(xcent,rev(xcent)),c(binf[,icol],rev(bsup[,icol])), col=plot.opt$fillcol,lty=plot.opt$ablinelty,border=plot.opt$ablinecol)
		if(length(dotline)==0) {
			for(icol in 1:3) lines(xcent,bmed[,icol],lty=plot.opt$ablinelty, lwd=plot.opt$ablinelwd,col=plot.opt$ablinecol)
		}
		percobs<-matrix(unlist(tapply(plmat[,4],plmat$grp,quantile,quant,na.rm=T)), ncol=3,byrow=T)
		for(icol in 1:3) lines(xcent,percobs[,icol],lty=plot.opt$lty, lwd=plot.opt$lwd,col=plot.opt$col)
	} else {
		plmat<-cbind(grp=plmat[,1],plmat)
		if(is.null(plot.opt$ylim)) plot.opt$ylim<-c(min(plmat[,4],na.rm=T), max(plmat[,4],na.rm=T))
		plot(plmat[,2],plmat[,4],type="n",xlab=plot.opt$xlab,ylab=plot.opt$ylab, main=main,cex=plot.opt$cex,cex.lab=plot.opt$cex.lab,ylim=plot.opt$ylim)
	}  
	points(plmat[plmat[,3]==0,2],plmat[plmat[,3]==0,4],col=plot.opt$obs.col, pch=plot.opt$pch,cex=plot.opt$cex)
	if(plot.opt$plot.loq) points(plmat[plmat[,3]==1,2],plmat[plmat[,3]==1,4], col=plot.opt$cens.col, pch=plot.opt$cens.pch,cex=plot.opt$cex) 
	if(length(dotline)>0) {
		for(i in dotline) abline(h=i,lty=plot.opt$ablinelty,lwd=plot.opt$ablinelwd, col=plot.opt$ablinecol)
	}
}

# histogram
aux.plot.hist<-function(plmat,plot.opt,distrib="norm", nclass=10,sim.ypl=NULL) {
	ndat<-dim(plmat)[1]
	plot.bands<-FALSE
	if(plot.opt$bands) {
		if(!plot.opt$approx.pi & length(sim.ypl)>0) plot.bands<-TRUE
		if(plot.opt$approx.pi) {
			xsamp<-switch(distrib,norm=rnorm(ndat*100),unif=runif(ndat*100))
			sim.ypl<-matrix(xsamp,nrow=ndat)
			plot.bands<-TRUE
		}
	}	
	if(plot.bands) {
		x<-hist(plmat$xpd,breaks=nclass,plot=FALSE)
		alpha<-plot.opt$pi.size
		if(alpha>0.5) alpha<-1-alpha
		tmat<-matrix(nrow=length(x$breaks)-1,ncol=dim(sim.ypl)[2])
		nB<-length(x$breaks)
		for(j in 1:dim(sim.ypl)[2]) {
			xvec<-cut(sim.ypl[,j],breaks=x$breaks,include.lowest=TRUE, ordered_result=TRUE)
			tmat[,j]<-table(xvec)
		}
		row.names(tmat)<-names(table(xvec))
		bnds<-apply(tmat,1,quantile,c(alpha/2,0.5,1-alpha/2))
		if(is.null(plot.opt$xlim)) plot.opt$xlim<-c(min(x$breaks),max(x$breaks))
		if(is.null(plot.opt$ylim)) plot.opt$ylim<-c(0,max(c(bnds,x$counts)))
		plot(x$breaks[-nB],x$counts,xlim=plot.opt$xlim,ylim=plot.opt$ylim, type="n",xlab=plot.opt$xlab,ylab=plot.opt$ylab,main=plot.opt$main, cex.lab=plot.opt$cex.lab)
		# PI on counts
		rect(x$breaks[-nB],bnds[1,],x$breaks[-1],bnds[3,],lwd=plot.opt$lwd.lpi, lty=plot.opt$lty.lpi, border=plot.opt$col.lpi,col=plot.opt$col.fillpi)
		segments(x$breaks[-nB],bnds[2,],x$breaks[-1],bnds[2,], lty=plot.opt$lty.lmed,lwd=plot.opt$lwd.lmed, col=plot.opt$col.lmed)		
		segments(x$breaks[-c(1,nB)],bnds[2,-1],x$breaks[-c(1,nB)], bnds[2,-(nB-1)],lty=plot.opt$lty.lmed,lwd=plot.opt$lwd.lmed, col=plot.opt$col.lmed)
		# Observed data		
		segments(x$breaks[-nB],x$counts,x$breaks[-1],x$counts,lwd=plot.opt$lwd+1, lty=plot.opt$lty,col=plot.opt$col)		
		#		segments(x$breaks[-c(1,nB)],x$counts[-1],x$breaks[-c(1,nB)], x$counts[-(nB-1)],lwd=plot.opt$lwd+1,lty=plot.opt$lty, col=plot.opt$col)
		rect(x$breaks[-nB],0,x$breaks[-1],x$counts,lwd=plot.opt$lwd+1, lty=plot.opt$lty,border=plot.opt$col)
	} else {
		x<-hist(plmat$xpd,breaks=nclass,xlab=plot.opt$xlab,ylab=plot.opt$ylab, main=plot.opt$main,lwd=plot.opt$lwd,border=plot.opt$col,col=plot.opt$fillcol, cex.lab=plot.opt$cex.lab)
		if(distrib=="norm") {
			xpl<-min(plmat$xpd)+c(0:100)/100*(max(plmat$xpd)-min(plmat$xpd))
			ypl<-dnorm(xpl)
			ypl<-ypl/max(ypl)*max(x$counts)
			lines(xpl,ypl,lwd=plot.opt$lwd,lty=plot.opt$lty)  
		}
		if(distrib=="unif") {
			abline(h=ndat/nclass,lty=plot.opt$ablinelty,lwd=plot.opt$lwd)
		}
	}
}

# Empirical cdf/QQplot
aux.plot.dist<-function(plmat,plot.opt,shape="ecdf",distrib="norm",nrep=0, nclass=0,ties=TRUE,sim.ypl=NULL,verbose=FALSE) {
	if(nclass==0) binning<-FALSE else {
		binning<-TRUE
		nbin<-plot.opt$vpc.bin
	}
	ndat<-dim(plmat)[1]
	plot.bands<-FALSE
	if(plot.opt$bands) {
		if(!plot.opt$approx.pi & length(sim.ypl)>0) plot.bands<-TRUE
		if(plot.opt$approx.pi) {
			xsamp<-switch(distrib,norm=rnorm(ndat*100),unif=runif(ndat*100))
			sim.ypl<-matrix(xsamp,nrow=ndat)
			plot.bands<-TRUE
		}
	}
	if(ties & nrep>0) yq<-seq(1/(2*nrep),1-1/(2*nrep),length.out=ndat) else  yq<-seq(0,1,length.out=ndat)
	if(distrib=="norm") yq<-qnorm(yq)		
	if(shape=="ecdf") {
		theo.samp<-yq
		yq<-seq(1/ndat,1,length.out=ndat)
	}
	ymat<-plmat[order(plmat$xpd),]
	ymat<-cbind(ymat,ecdf=yq)
	if(shape=="ecdf") {
		ysh<-data.frame(x=plmat$xpd,y=yq)
		yobs<-data.frame(x=ymat$xpd,y=ymat$ecdf,cens=ymat$cens)
	} else {
		ysh<-data.frame(x=yq,y=plmat$xpd)
		yobs<-data.frame(x=ymat$ecdf,y=ymat$xpd,cens=ymat$cens)
	}
	if(!binning) {
		plot(ysh$x,ysh$y,xlim=plot.opt$xlim,ylim=plot.opt$ylim,type="n", xlab=plot.opt$xlab,ylab=plot.opt$ylab,main=plot.opt$main, cex.lab=plot.opt$cex.lab,log=plot.opt$logtyp)
		if(plot.bands) {
			alpha<-plot.opt$pi.size
			if(alpha>0.5) alpha<-1-alpha
			sim.sort<-colsort(sim.ypl)
			bnds<-apply(sim.sort,1,quantile,c(alpha/2,0.5,1-alpha/2))
			if(shape=="ecdf") {
				rect(bnds[1,-ndat],yq[-ndat],bnds[3,-1],yq[-1], col=plot.opt$fillcol, border=plot.opt$fillcol) 
				lines(bnds[2,],yq,lty=plot.opt$lty.lmed, col=plot.opt$col.lmed, lwd=plot.opt$lwd.lmed)
			} else { 
				rect(yq[-ndat],bnds[1,-ndat],yq[-1],bnds[3,-1],col=plot.opt$fillcol, border=plot.opt$fillcol)
				lines(yq,bnds[2,],lty=plot.opt$lty.lmed, col=plot.opt$col.lmed, lwd=plot.opt$lwd.lmed)
			}
		}
		lines(yobs$x,yobs$y,lty=plot.opt$lty, col=plot.opt$lcol,lwd=plot.opt$lwd)
		if(plot.opt$type=="p" | plot.opt$type=="b") {
			points(yobs$x[yobs$cens==0],yobs$y[yobs$cens==0], col=plot.opt$obs.col, pch=plot.opt$pch,cex=plot.opt$cex)
			if(plot.opt$plot.loq) points(yobs$x[yobs$cens==1],yobs$y[yobs$cens==1], col=plot.opt$cens.col, pch=plot.opt$cens.pch,cex=plot.opt$cex)
		}
		if(!plot.bands) { # theoretical distribution
			if(shape=="qqplot" | distrib=="unif") abline(0,1,lty=plot.opt$ablinelty, lwd=plot.opt$ablinelwd) else lines(theo.samp,c(1:ndat)/ndat, lty=plot.opt$ablinelty, lwd=plot.opt$ablinelwd)
		}
	} else {
		if(verbose) cat("Binning not implemented yet\n")
	}
}

# Binning the X data
npde.binning<-function(xvec,plot.opt,verbose=FALSE) {
  xvec1<-xvec
  xvec<-xvec[!is.na(xvec)]
  if(!is.na(pmatch(plot.opt$vpc.method,"optimal"))) {
  	if(verbose) cat("Optimal binning not yet implemented, reverting to equal binning\n")
    plot.opt$vpc.method<-"equal"
  }
  if(!is.na(pmatch(plot.opt$vpc.method,"user")) & is.null(plot.opt$vpc.breaks)) {
    if(verbose) cat("User-defined method specified, but vpc.breaks is empty; reverting to equal binning\n")
    plot.opt$vpc.method<-"equal"
  }
  if(!is.na(pmatch(plot.opt$vpc.method,c("equal","width"))) & is.null(plot.opt$vpc.bin)) {
    plot.opt$vpc.bin<-10
  }
  nbin<-plot.opt$vpc.bin
  if(!is.na(pmatch(plot.opt$vpc.method,"user"))) {
      bnds<-plot.opt$vpc.breaks
      if(min(bnds)>=min(xvec)) bnds<-c(min(xvec)-1,bnds)
      if(max(bnds)<max(xvec)) bnds<-c(bnds,max(xvec))
  }
  if(!is.na(pmatch(plot.opt$vpc.method,"equal"))) {
      xvec2<-xvec;xvec2[xvec2==min(xvec)]<-min(xvec)-1
      bnds<-unique(quantile(xvec2,(0:nbin)/nbin,type=8))
  }
  if(!is.na(pmatch(plot.opt$vpc.method,"width"))) {
      if(plot.opt$xlog) xvec2<-log(xvec) else xvec2<-xvec
      bnds<-seq(min(xvec2),max(xvec2),length.out=(nbin+1))
      if(plot.opt$xlog) bnds<-exp(bnds)
      bnds[1]<-bnds[1]-1
  }
  if(!is.na(pmatch(plot.opt$vpc.method,c("equal","width","user")))) {
      xgrp<-factor(cut(xvec,bnds,include.lowest=F))
      if(!is.na(pmatch(plot.opt$vpc.method,"equal")) & length(unique(xvec))<=nbin)
        xgrp<-match(xvec,sort(unique(xvec)))
  } else {
      
  }
  nbin<-length(unique(xgrp))
  xpl<-tapply(xvec,xgrp,mean)
  npl<-tapply(xvec,xgrp,length)
  if(!is.na(pmatch(plot.opt$vpc.method,c("equal","width","user")))) {
    tab<-cbind(Interval=names(xpl),Centered.On=format(xpl,digits=2),Nb.obs=npl)
    row.names(tab)<-1:dim(tab)[1]
    if(verbose) {
    	xnam<-switch(EXPR=plot.opt$vpc.method,equal="by quantiles on X", width="equal sized intervals",user="user-defined bins")
    	cat("Method used for binning:",xnam,", dividing into the following intervals\n")
    	print(tab,quote=F)
    }
  }
  xgrp2<-rep(NA,length(xvec1))
  xgrp2[!is.na(xvec1)]<-xgrp
  return(list(xgrp=xgrp2,xat=xpl))
}

# Prediction intervals, vpc
compute.vpc.pi<-function(ysim,xgrp,idrep,nbin,vpc.pi=0.95) {
  nsim<-length(unique(idrep))
  sim.pi.low<-sim.pi.med<-sim.pi.up<-matrix(0,nrow=nbin,ncol=nsim)
  alpha<-(1-vpc.pi)/2
  i0<-1
  for(irep in unique(idrep)) {
    ysim1<-ysim[idrep==irep]
    l1<-unlist(tapply(ysim1,xgrp,function(vec) quantile(vec,c(alpha,0.5,1-alpha))))
    l1<-matrix(l1,ncol=3,byrow=TRUE)
    sim.pi.low[,i0]<-l1[,1]
    sim.pi.med[,i0]<-l1[,2]
    sim.pi.up[,i0]<-l1[,3]
    i0<-i0+1
  }
  return(list(sim.pi.low=sim.pi.low,sim.pi.med=sim.pi.med,sim.pi.up=sim.pi.up))
}

###############################	   VPC	 ########################################

npde.plot.vpc<-function(npdeObject,npc=FALSE,...) {
  has.cens<-FALSE
  if(length(npdeObject["data"]["icens"])>0) {
    has.cens<-TRUE
    is.cens<-(npdeObject["data"]["data"][,npdeObject["data"]["name.cens"]]==1)
  } else is.cens<-rep(FALSE,dim(npdeObject["data"]["data"])[1])
  plot.opt<-npdeObject["prefs"]
  if(npdeObject["options"]["cens.method"]=="omit") plot.opt$impute.loq<-FALSE
  plot.opt$main<-"Visual Predictive Check"
  plot.opt<-replace.npde.plotoptions(plot.opt,...)
  if(plot.opt$new) {
    mfrow<-plot.opt$mfrow
    if(length(mfrow)==0) mfrow<-c(1,1)
    par(mfrow=mfrow,ask=plot.opt$ask)
  }
  logtyp<-""
  if(plot.opt$xlog) logtyp<-paste(logtyp,"x",sep="")
  if(plot.opt$ylog) logtyp<-paste(logtyp,"y",sep="")  

# Binning
  not.miss<-npdeObject["data"]["not.miss"] 
  if(!plot.opt$impute.loq & has.cens) not.miss[npdeObject["data"]["data"][, npdeObject["data"]["name.cens"]]==1]<-FALSE
  xvec<-npdeObject["data"]["data"][,npdeObject["data"]["name.predictor"]]
  if(has.cens & npdeObject["options"]["cens.method"]!="omit" & plot.opt$impute.loq) ydat<-npdeObject["results"]["res"]$ycomp else ydat<-npdeObject["data"]["data"][, npdeObject["data"]["name.response"]]
  xvec<-xvec[not.miss]
  ydat<-ydat[not.miss]
  is.cens<-is.cens[not.miss]
  ysim<-npdeObject["sim.data"]["datsim"]$ysim[not.miss]
  alpha<-(1-plot.opt$vpc.interval)/2
  xbin<-npde.binning(xvec,plot.opt)
  xgrp<-xbin$xgrp
  xpl<-xbin$xat
  nbin<-length(unique(xgrp))
# ECO TODO: implement the optimal binning algorithm of Marc (see library mclust)
# Observed data
    ypl<-tapply(ydat,xgrp,mean)
    obs.bnd<-cbind(tapply(ydat,xgrp,quantile,alpha),tapply(ydat,xgrp,mean), tapply(ydat,xgrp,quantile,1-alpha))
#  }
  if(plot.opt$vpc.pi) {
    idsim<-npdeObject["sim.data"]["datsim"]$idsim[not.miss]
    idrep<-npdeObject["sim.data"]["datsim"]$irsim[not.miss]
    nsim<-npdeObject["prefs"]$bands.nrep
    if(nsim>npdeObject["sim.data"]["nrep"]) nsim<-npdeObject["sim.data"]["nrep"]
    isamp<-sample(1:npdeObject["sim.data"]["nrep"],nsim)
    idx<-match(idrep,isamp,nomatch=0)>0
    sbnd<-compute.vpc.pi(ysim[idx],xgrp,idrep[idx],nbin,0.95)
    pi.low<-apply(sbnd$sim.pi.low,1,quantile,c(0.025,0.5,0.975))
    pi.med<-apply(sbnd$sim.pi.med,1,quantile,c(0.025,0.5,0.975))
    pi.up<-apply(sbnd$sim.pi.up,1,quantile,c(0.025,0.5,0.975))
    vec1<-c(pi.low,obs.bnd[,1])
    vec2<-c(obs.bnd[,3],pi.up)
    if(length(grep("y",logtyp))>0) {
      vec1<-vec1[vec1>0]
      vec2<-vec2[vec2>0]
    }
    limy<-c(min(vec1),max(vec2))
    xvec1<-xvec
    if(length(grep("x",logtyp))>0) xvec1<-xvec1[xvec1>0]
    limx<-c(min(xvec1),max(xvec1))
    plot(xpl,ypl,type="n",xlim=limx,ylim=limy,xlab=plot.opt$xlab, ylab=plot.opt$ylab,main=plot.opt$main,log=logtyp)    
    polygon(c(xpl,rev(xpl)),c(pi.low[1,],rev(pi.low[3,])), col=plot.opt$col.fillpi,lty=plot.opt$lty.lpi, border=plot.opt$col.lpi)
    polygon(c(xpl,rev(xpl)),c(pi.up[1,],rev(pi.up[3,])), col=plot.opt$col.fillpi,lty=plot.opt$lty.lpi, border=plot.opt$col.lpi)
    polygon(c(xpl,rev(xpl)),c(pi.med[1,],rev(pi.med[3,])), col=plot.opt$col.fillmed,lty=plot.opt$lty.lmed, border=plot.opt$col.lmed)
    lines(xpl,pi.low[2,],lty=plot.opt$lty.lpi, col=plot.opt$col.lpi,lwd=plot.opt$lwd.lpi)
    lines(xpl,pi.med[2,],lty=plot.opt$lty.lmed, col=plot.opt$col.lmed,lwd=plot.opt$lwd.lmed)
    lines(xpl,pi.up[2,],lty=plot.opt$lty.lpi, col=plot.opt$col.lpi,lwd=plot.opt$lwd.lpi)
    lines(xpl,obs.bnd[,2],lty=plot.opt$lty.lobs, col=plot.opt$col.lmed,lwd=plot.opt$lwd.lobs)
    for (icol in c(1,3)) lines(xpl,obs.bnd[,icol],lty=plot.opt$lty.lobs, col=plot.opt$col.lobs,lwd=plot.opt$lwd.lobs)
    if(plot.opt$plot.obs) {
    if(has.cens) {
      points(xvec[!is.cens],ydat[!is.cens],pch=plot.opt$pch,col=plot.opt$obs.col, cex=plot.opt$cex)
      if(plot.opt$plot.loq) points(xvec[is.cens],ydat[is.cens], col=plot.opt$cens.col, pch=plot.opt$cens.pch, cex=plot.opt$cex)
      if(plot.opt$line.loq & length(npdeObject["data"]["loq"])>0) abline(h=npdeObject["data"]["loq"], lty=plot.opt$ablinelty, col=plot.opt$ablinecol,lwd=plot.opt$ablinelwd)
    } else points(xvec,ydat,pch=plot.opt$pch,col=plot.opt$obs.col, cex=plot.opt$cex)
    }
  } else {
# Simulated data
    nsim<-npdeObject["sim.data"]["nrep"]
    id.grp<-rep(xgrp,nsim)
    sim.bnd<-cbind(tapply(ysim,id.grp,quantile,alpha),tapply(ysim,id.grp,mean), tapply(ysim,id.grp,quantile,1-alpha))
    vec1<-c(obs.bnd[,1],sim.bnd[,1])
    vec2<-c(obs.bnd[,3],sim.bnd[,3])
    if(length(grep("y",logtyp))>0) {
      vec1<-vec1[vec1>0]
      vec2<-vec2[vec2>0]
    }
    limy<-c(min(vec1),max(vec2))
    xvec1<-xvec
    if(length(grep("x",logtyp))>0) xvec1<-xvec1[xvec1>0]
    limx<-c(min(xvec1),max(xvec1))
    plot(xpl,ypl,type="n",xlim=limx,ylim=limy,xlab=plot.opt$xlab, ylab=plot.opt$ylab,main=plot.opt$main,log=logtyp)
    polygon(c(xpl,rev(xpl)),c(sim.bnd[,3],rev(sim.bnd[,1])), col=plot.opt$fillcol,lty=plot.opt$ablinelty, border=plot.opt$ablinecol)
    lines(xpl,sim.bnd[,2],lty=plot.opt$ablinelty, col=plot.opt$ablinecol,lwd=plot.opt$ablinelwd)
    lines(xpl,obs.bnd[,2],lty=plot.opt$lty, col=plot.opt$lcol,lwd=plot.opt$lwd)
    for (icol in c(1,3)) lines(xpl,obs.bnd[,icol],lty=plot.opt$lty, col=plot.opt$lcol,lwd=plot.opt$lwd)
    if(plot.opt$plot.obs) {
    if(has.cens) {
      points(xvec[!is.cens],ydat[!is.cens],pch=plot.opt$pch,col=plot.opt$obs.col, cex=plot.opt$cex)
      if(plot.opt$plot.loq) points(xvec[is.cens], ydat[is.cens],col=plot.opt$cens.col, pch=plot.opt$cens.pch, cex=plot.opt$cex)
      if(plot.opt$line.loq & length(npdeObject["data"]["loq"])>0) abline(h=npdeObject["data"]["loq"],lty=plot.opt$ablinelty, col=plot.opt$ablinecol,lwd=plot.opt$ablinelwd)
    } else points(xvec,ydat,pch=plot.opt$pch,col=plot.opt$obs.col, cex=plot.opt$cex)
    }
  }
  npc.stat<-c()
  if(npc==TRUE) {
    # ECO TODO: compute NPC - interpolation ? 
  }
  invisible(list(npc=npc.stat))
}

################################  p_LOQ  ########################################

npde.plot.loq<-function(npdeObject,xaxis="x",nsim=200,...) {
# Plot of the probability of an observation being LOQ versus X or predictions, overlaying
### the predicted probability (according to the model)
### the observed probability
### a prediction band obtained using the simulated data
  args1<-match.call(expand.dots=TRUE)
  i1<-match("loq",names(args1))
  if(!is.na(i1)) {
    loq<-as.numeric(args1[[i1]])
  } else {
    if(length(npdeObject["data"]["loq"])==0) {
      ploq<-c()
      yobs<-npdeObject["data"]["data"][,npdeObject["data"]["name.response"]]
      if(length(npdeObject["data"]["loq"])>0) {
        loq<-npdeObject["data"]["loq"]
        if(npdeObject@options$verbose) cat("Computing p(y<LOQ) using LOQ=",loq,"\n")
      } else {
      	yloq<-yobs[npdeObject["data"]["icens"]]
      	if(length(unique(yloq))==1) {
      		if(npdeObject@options$verbose) cat("Same LOQ for all missing data, loq=",loq,"\n")
        	loq<-unique(yloq)
        } else {
        	loq<-min(unique(yloq))
        	if(npdeObject@options$verbose) cat("Computing p(y<LOQ) for the lowest LOQ, loq=",loq,"\n")
      }
      npdeObject["data"]["loq"]<-loq
      }
      if(is.infinite(npdeObject["data"]["loq"])) {
      	if(npdeObject@options$verbose) cat("No loq defined in the data, and no censored data to define it, please call npde.plot.loq with the option loq=XXX where XXX is the value of the LOQ.\n")
        return()
      }
    } else loq<-npdeObject["data"]["loq"]
  }
  has.cens<-FALSE
  if(length(npdeObject["data"]["icens"])>0) {
    has.cens<-TRUE
    icens<-npdeObject["data"]["icens"] #ECO TODO: icens jamais utilise
  }
  plot.opt<-npdeObject["prefs"]
  plot.opt$main<-"Probability of being under the LOQ"
  plot.opt$ylab<-"Pr(Y<LOQ)"
  plot.opt<-replace.npde.plotoptions(plot.opt,...)
  if(plot.opt$new) {
    mfrow<-plot.opt$mfrow
    if(length(mfrow)==0) mfrow<-c(1,1)
    par(mfrow=mfrow,ask=plot.opt$ask)
  }
  logtyp<-""
  if(plot.opt$xlog) logtyp<-paste(logtyp,"x",sep="")
  if(plot.opt$ylog) logtyp<-paste(logtyp,"y",sep="")  
  nsim<-min(nsim,npdeObject["sim.data"]["nrep"])

# Binning
  xvec<-switch(xaxis, x=npdeObject["data"]["data"][,npdeObject["data"]["name.predictor"]], pred=npdeObject["results"]["res"]$ypred, cov="Not implemented yet")
  if(!is.numeric(xvec)) {
  	if(npdeObject@options$verbose) cat(xvec,"\n")
    return()
  }
  if(has.cens) ydat<-npdeObject["data"]["data"][,npdeObject["data"]["name.cens"]] else ydat<-rep(0,length(xvec))
  xbin<-npde.binning(xvec,plot.opt)
  xgrp<-xbin$xgrp
  xpl<-xbin$xat
  nbin<-length(unique(xgrp))
  isamp<-sample(1:npdeObject["sim.data"]["nrep"],nsim)
  ysim<-npdeObject["sim.data"]["datsim"]$ysim
  xtab<-matrix(nrow=nsim,ncol=nbin)
  for(i in 1:nsim) 
    xtab[i,]<-tapply(ysim[npdeObject["sim.data"]["datsim"]$irsim==isamp[i]] < loq,xgrp,mean)
  alpha<-(1-plot.opt$vpc.interval)/2
  quant<-c(alpha,0.5,1-alpha)
  ypl<-apply(xtab,2,quantile,quant)
  xobs<-tapply(ydat,xgrp,mean)
  if(is.null(plot.opt$ylim)) plot.opt$ylim<-c(0,max(c(xtab,xobs),na.rm=T))
  plot(xpl,ypl[1,],type="n",xlim=plot.opt$xlim,ylim=plot.opt$ylim,xlab=plot.opt$xlab, ylab=plot.opt$ylab,main=plot.opt$main,log=logtyp,cex=plot.opt$cex, cex.lab=plot.opt$cex.lab)
  polygon(c(xpl,rev(xpl)),c(ypl[1,],rev(ypl[3,])),col=plot.opt$col.fillpi, lty=plot.opt$lty.lpi, border=plot.opt$col.lpi)
  lines(xpl,ypl[2,],lty=plot.opt$lty.lmed,col=plot.opt$col.lmed, lwd=plot.opt$lwd.lmed)
  lines(xpl,xobs,lty=plot.opt$obs.lty,col=plot.opt$obs.col, lwd=plot.opt$obs.lwd)
}
