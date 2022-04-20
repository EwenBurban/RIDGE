#!/usr/bin/Rscript

#################################################################################################################################
#################################################################################################################################
#####                                                                                                                       #####
#####    This file is part of Demographic Inferences with Linked Selection : DILS.                                          #####
#####                                                                                                                       #####   
#####    DILS is free software: you can redistribute it and/or modify                                                       #####
#####    it under the terms of the GNU General Public License as published by                                               #####
#####    the Free Software Foundation, either version 3 of the License, or                                                  #####
#####    (at your option) any later version.                                                                                #####
#####                                                                                                                       #####    
#####    DILS is distributed in the hope that it will be useful,                                                            #####
#####    but WITHOUT ANY WARRANTY; without even the implied warranty of                                                     #####
#####    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                                      #####
#####    GNU General Public License for more details.                                                                       #####
#####                                                                                                                       #####    
#####    You should have received a copy of the GNU General Public License                                                  #####
#####    along with DILS.  If not, see <https://www.gnu.org/licenses/>.                                                     #####
#####                                                                                                                       #####    
#####    Please send bugreports with examples or suggestions to                                                             #####
#####    camille.roux@univ-lille.fr                                                                                         #####
#####                                                                                                                       #####    
#####    Or write a post on https://groups.google.com/forum/#!forum/dils---demographic-inferences-with-linked-selection     #####
#####                                                                                                                       #####
#################################################################################################################################
#################################################################################################################################

# function to estimate the parameters
abc_nnet_multivar <- function(target,x,sumstat,tol,gwt,rejmethod=F,noweight=F,transf="none",bb=c(0,0),nb.nnet=10,size.nnet=5,trace=T, MaxNWts=10000){
	options(digits=5)
	require(nnet)
	# target is the set of target summary stats
	# x is the parameter vector (long vector of numbers from the simulations) and is the dependent variable for the regression
	# x can also be a matrix for multi-dimensional models. Each column corresponds to a parameter and each row to a simulation.
	# sumstat is an array of simulated summary stats (i.e. independent variables).
	# tol is the required proportion of points nearest the target values
	# gwt is a vector with T/F weights, weighting out any 'bad' values (determined by the simulation program - i.e. nan's etc)
	# if noweight=T, no Epanechnikov weights are calculated
	# if rejmethod=T it doesn't bother with the regression, and just does rejection.
	# nb.nnet>1 is the number of trained neural networks, the more neural nets the more robust is the inference
	# size.nnet is the number of hidden network in the regression. Typically >= the number of parameters in the model
	# transf the vector of transformation for the parameter, ex:transf=c("none","logit","log")
	# bb the vector of bounds for the logit transformation ex:bb=cbind(c(0,0),c(0,2),c(0,0)) (The second column is the only one to be taken into account)
	# If trace=T print messages during the algorithm
	# If rejmethod=F it returns a list with the following components:-
	# $x regression adjusted values
	# $vals - unadjusted values in rejection region (i.e. normal rejection)
	# $wt - the regression weight (i.e. the Epanechnikov weight)
	# $ss - the sumstats corresponding to these points
	# $predmean - estimate of the posterior mean
	if(class(x)=="numeric")
	{
		bb<-cbind(bb)
		x<-cbind(x)
	}

	if(rejmethod)
		transf<-rep("none", dim(x)[2])
	normalise <- function(x,y){
	if(mad(y) == 0)
	return (x)
	else
	return (x/mad(y))
	}

	####Define the weight-decay paramaeter
	repet<-floor(nb.nnet/3)+1
	the_decay<-rep(c(10^(-4),10^(-3),10^(-2)),repet)[1:nb.nnet]
	lt<-dim(x)[2]
	nb_simu<-dim(sumstat)[1]
	for (i in 1:lt)
	{
	if(sum(transf[i] == c("none","log","logit")) == 0){
		stop("transf must be none, log, or logit")
	}
	if(transf[i]=="logit"){
		if(bb[1,i] >= bb[2,i]){
			stop("bounds wrong for logit")
		}
	}
	}

	if(missing(gwt))gwt <- rep(T,length(sumstat[,1]))
	nss <- length(sumstat[1,])
	# scale everything
	    scaled.sumstat <- sumstat
	    for(j in 1:nss){
		scaled.sumstat[,j] <- normalise(sumstat[,j],sumstat[,j][gwt])
	    }
	    target.s.tmp <- target
	    for(j in 1:nss){
		target.s.tmp[,j] <- normalise(target[,j],sumstat[,j][gwt])
	    }

	#déplacé (origine : après calcul de abstol)
		    for (i in 1:lt)
		{
		    if(transf[i] == "log"){
			if(min(x[,i]) <= 0){
				print("log transform: val out of bounds - correcting")
				x.tmp <- ifelse(x[,i] <= 0,max(x[,i]),x[,i])
				x.tmp.min <- min(x.tmp)
				xx[,i] <- ifelse(x[,i] <= 0, x.tmp.min,x[,i])
			}
			x[,i] <- log(x[,i])
		    }
		    else if(transf[i] == "logit"){
			if(min(x[,i]) <= bb[1,i]){
				x.tmp <- ifelse(x[,i] <= bb[1,i],max(x[,i]),x[,i])
				x.tmp.min <- min(x.tmp)
				x[,i] <- ifelse(x[,i] <= bb[1,i], x.tmp.min,x[,i])
			}
			if(max(x[,i]) >= bb[2,i]){
				x.tmp <- ifelse(x[,i] >= bb[2,i],min(x[,i]),x[,i])
				x.tmp.max <- max(x.tmp)
				x[,i] <- ifelse(x[,i] >= bb[2,i], x.tmp.max,x[,i])
			}
			x[,i] <- (x[,i]-bb[1,i])/(bb[2,i]-bb[1,i])
			x[,i] <- log(x[,i]/(1-x[,i]))
		    }
		}

	#boucle le long des targets
	for(L in 1:nrow(target.s.tmp)){
		target.s=as.numeric(target.s.tmp[L,])
		sum1=dst=abstol=wt1=regwt=l1=fit1=ll=predmean=array_pred=mean_averaged=my_residuals=predvar=var_averaged=the_sd=predsd=x.tmp=x.tmp.min=xx=x.tmp.max=NULL
		# calc euclidean distance
		    sum1 <- 0
		    for(j in 1:nss){
			sum1 <- sum1 + (scaled.sumstat[,j]-target.s[j])^2
		   }
		   dst <- sqrt(sum1)
		# includes the effect of gwt in the tolerance
		    dst[!gwt] <- floor(max(dst[gwt])+10)

		# wt1 defines the region we're interested in
		    abstol <- quantile(dst,tol)
		    wt1 <- dst <= abstol
		    if(rejmethod){
			regwt <- 1-dst[wt1]^2/abstol^2
			  l1 <- list(x=cbind(x[wt1,]),wt=regwt,ind=wt1,dst=dst)
		    }
		    else{
			  regwt <- 1-dst[wt1]^2/abstol^2
			if(noweight)
				regwt <- rep(1,length(regwt))
			ll<-NULL
			if(trace==TRUE)
				cat("Regression of the mean ")

			for (i in 1:nb.nnet)
			{
				if(trace==TRUE)
					cat(i," ")
				fit1 <- nnet(scaled.sumstat[wt1,],x[wt1,],weights=regwt,decay=the_decay[i],size=size.nnet,linout=T,maxit=500,trace=F, MaxNWts=MaxNWts)
				ll<-c(ll,list(fit1))
			}
			#Compute the residuals
			predmean<-NULL
			array_pred<-array(dim=c(nb.nnet,sum(wt1),lt))
			for (i in 1:nb.nnet)
			{
				array_pred[i,,]<-ll[[i]]$fitted.values
				predmean<-cbind(predmean,as.numeric(predict(ll[[i]],data.frame(rbind(target.s)))))	
			}
			mean_averaged<-NULL
			for (j in 1:lt)
				mean_averaged<-cbind(mean_averaged,apply(array_pred[,,j],FUN=median,MARGIN=2))
			  predmean<-apply(predmean,FUN=median,MARGIN=1)

			  my_residuals<-(x[wt1,]-mean_averaged)

			#Fit a neural network for predicting the conditional variance
			if(trace==TRUE)
				cat("\nRegression of the variance ")
			ll<-NULL
			for (i in 1:nb.nnet)
			{
				if(trace==TRUE)
					cat(i," ")
				fit2 <- nnet(scaled.sumstat[wt1,],log(my_residuals^2),weights=regwt,decay=the_decay[i],size=size.nnet,linout=T,maxit=500,trace=F, MaxNWts=MaxNWts)
				ll<-c(ll,list(fit2))
			}
			if(trace==TRUE)
				cat("\n")
			predvar<-NULL
			array_pred<-array(dim=c(nb.nnet,sum(wt1),lt))
			for (i in 1:nb.nnet)
			{
				array_pred[i,,]<-ll[[i]]$fitted.values
				predvar<-cbind(predvar,as.numeric(predict(ll[[i]],data.frame(rbind(target.s)))))	
			}
			var_averaged<-NULL
			for (j in 1:lt)
				var_averaged<-cbind(var_averaged,apply(array_pred[,,j],FUN=median,MARGIN=2))
			the_sd<-sqrt(exp(var_averaged))
			predsd<-sqrt(exp(apply(predvar,FUN=median,MARGIN=1)))
			res_correc<-sapply(1:lt,FUN=function(i){predmean[i]+ ((predsd[i]*my_residuals[,i])/the_sd[,i]) })
			l1 <- list(x=cbind(res_correc),vals=cbind(x[wt1,]),wt=regwt,ss=sumstat[wt1,],predmean=predmean)
			}
		    for (i in 1:lt)
		    {
		    if(transf[i] == "log"){
			l1$x[,i] <- exp(l1$x[,i])
			l1$vals[,i] <- exp(l1$vals[,i])
		    }
		    if(transf[i] == "logit"){
		    l1$x[,i] <- exp(l1$x[,i])/(1+exp(l1$x[,i]))
		    l1$x[,i] <- l1$x[,i]*(bb[2,i]-bb[1,i])+bb[1,i]
		    l1$vals[,i] <- exp(l1$vals[,i])/(1+exp(l1$vals[,i]))
		    l1$vals[,i] <- l1$vals[,i]*(bb[2,i]-bb[1,i])+bb[1,i]
		    }
		}
		    return(l1)
	#	    write.table(l1$x, col.names=F, row.names=F, file=paste(output,L,sep=""))
		}
}

# function to plot the prior and posterior
babar<-function(a,b,space=2,breaks="auto",AL=0.5,nameA="A",nameB="B",xl="",yl="",mn="",legx="topright", legende=TRUE){ 
       aprime=a;
       bprime=b;
       if(length(a)>length(b)){ bprime=b; aprime=sample(a,length(b),replace=F) }
       if(length(a)<length(b)){ aprime=a; bprime=sample(b,length(a),replace=F) }

       if(breaks=="auto"){
            bks=hist(c(aprime,bprime),plot=F)$breaks
            bklong=space*length(bks)
            bks=hist(c(aprime,bprime),plot=F,breaks=bklong)$breaks
       }
       else{
            bks=breaks
       }

       h1=hist(a,breaks=bks,plot=F)
       h2=hist(b,breaks=bks,plot=F)
       w1=sum(h1$density)
       w2=sum(h2$density)
       d1=max(h1$density)/w1
       d2=max(h2$density)/w2
       d=max(d1,d2)
       par(lwd=1)
       x=barplot(h1$density/w1,col="white",border=par("fg"),ylim=c(0,d),width=.8,space=.2,ylab=yl,xlab=xl,main=mn, cex.lab=1.2) 
       y=c(x,x[length(x)]+.96)-.5
       axis(side=1,at=y,labels=h1$breaks)
       par(lwd=2,lty=1)
       barplot(h2$density/w2,col=rgb(red=.25,blue=.25,green=.25,alpha=AL),border=NA,ylim=c(0,d/w2),width=.8,space=.2,add=T) 
       if(legende==TRUE){legend(legx,legend=c(nameA,nameB),fill=c("white",rgb(red=.25,blue=.25,green=.25,alpha=.5)),cex=1.5,bty="n")}
}


## get the arguments
#for(i in commandArgs()){
#	tmp = strsplit(i, '=')
#	if(tmp[[1]][1] == 'nameA'){ nameA = tmp[[1]][2] }
#	if(tmp[[1]][1] == 'nameB'){ nameB = tmp[[1]][2] }
#	if(tmp[[1]][1] == 'nCPU'){ nCPU = as.integer(tmp[[1]][2]) }
#	if(tmp[[1]][1] == 'model'){ model = tmp[[1]][2] } # model to simulate
#	if(tmp[[1]][1] == 'nMin'){ nMin = as.integer(tmp[[1]][2]) } # minimal number of sequences
#}


get_posterior<-function(nameA='spA', nameB='spB', nSubdir=10, sub_dir_sim='iteration_2', model='estimation_1', sub_dir_model=1, nPosterior=1000, figure=FALSE, timeStamp='SI_2N', path2observation='SI_2N/iteration_2/pod_0', do_nnet=T, useSFS=0, ncores=ncores){
	library(data.table)
	options(digits=5)
	###################
	# get observed data
	# observed data
	coul = c('#ffffcc', '#c7e9b4', '#7fcdbb', '#41b6c4', '#1d91c0', '#225ea8', '#0c2c84')
	coul = colorRampPalette(coul)

	#obs_ss = read.table(paste(timeStamp, '/ABCstat_global.txt', sep=''), h=T)
	obs_ss_tmp = read.table(paste(path2observation, '/ABCstat_global.txt', sep=''), h=T)

	undesired_stats = c('min', 'max', 'successive', 'pearson', 'ss_sf', 'ss_noSf', 'noSs_sf', 'noSs_noSf','dataset')

	# remove stats
	for(stat in undesired_stats){
		obs_ss_tmp = obs_ss_tmp[, -grep(stat, colnames(obs_ss_tmp))]
	}
	# end of remove stats


	if(useSFS == 1){
		sfs_obs = read.table(paste(path2observation, '/ABCjsfs.txt', sep=''), h=T)
		sfs_obs = sfs_obs[, -c(which(colnames(sfs_obs)=='fA0_fB0'), which(colnames(sfs_obs)=='fA1_fB0'), which(colnames(sfs_obs)=='fA0_fB1'))] # remove the singletons and the (0,0)
		ss_obs = cbind(obs_ss_tmp, sfs_obs)
	}else{
		ss_obs = obs_ss_tmp
	}

	
##### get the number of statistics, number of simulations and number of parameters
####tmp = read.table(paste(timeStamp, '/', sub_dir_sim, '/', model, '_0/ABCstat_global.txt', sep=''), h=T)
####nSimulations = nrow(tmp)
####nStats = ncol(tmp)
####print(nSimulations)
####ss_sim_tmp = matrix(NA, nrow=nSimulations*nSubdir, ncol=nStats)
####colnames(ss_sim_tmp) = colnames(tmp)
####
####tmp = read.table(paste(timeStamp, '/', sub_dir_sim, '/', model, '_0/priorfile.txt', sep=''), h=T)[,-1]
####nParams = ncol(tmp)
####params_sim_tmp = matrix(NA, nrow=nSimulations*nSubdir, ncol=nParams)
####colnames(params_sim_tmp) = colnames(tmp)

####if(useSFS == 1){
####	sfs_sim_tmp = matrix(NA, nrow=nSimulations*nSubdir, ncol=ncol(sfs_obs)); colnames(sfs_sim_tmp) = colnames(sfs_obs)
####}
####	
	# read all these simulated datasets
	ss_sim = list()
	params_sim = list()
    ss_sim_tmp= list()
    sfs_sim_tmp= list()
    params_sim_tmp = list()

	for(rep in seq(0, nSubdir-1, 1)){
        tmp_ss = read.table(file.path(timeStamp,sub_dir_sim,paste0(model,'_',rep),'ABCstat_global.txt'),h=T)
        tmp_ss = tmp_ss[order(tmp_ss[,'dataset']),]
        ss_sim_tmp[[length(ss_sim_tmp)+1]] = tmp_ss

        if(useSFS==1){
            sfs_sim = read.table(file.path(timeStamp,sub_dir_sim,paste0(model,'_',rep),'ABCjsfs.txt'),h=T)
            sfs_sim = sfs_sim[order(sfs_sim[,'dataset']),]
            sfs_sim_tmp[[length(sfs_sim_tmp)+1]] = sfs_sim
        }

        tmp_params = read.table(file.path(timeStamp,sub_dir_sim,paste0(model,'_',rep),'priorfile.txt'),h=T)
        tmp_params = tmp_params[order(tmp_params[,'dataset']),]
        tmp_params = tmp_params[match(tmp_ss$dataset,tmp_params$dataset),]
        params_sim_tmp[[length(params_sim_tmp)+1]] = tmp_params

        stopifnot(nrow(tmp_ss)==nrow(tmp_params))
        if(useSFS==1){stopifnot(nrow(sfs_sim)==nrow(tmp_params))}

	##### statistics
	####tmp_ss = as.matrix(fread(paste(timeStamp, '/', sub_dir_sim, '/', model, '_', rep, '/ABCstat_global.txt', sep=''), h=T))
    ####tmp_ss = tmp_ss[order(tmp_ss[,'dataset']),]
	####ss_sim_tmp[(rep*nSimulations+1):((rep+1)*nSimulations),] = as.matrix(tmp_ss)
	####
	####if(useSFS == 1){
	####	# sfs
	####	sfs_sim = as.matrix(fread(paste(timeStamp, '/', sub_dir_sim, '/', model, '_', rep, '/ABCjsfs.txt', sep=''), h=T))
    ####    sfs_sim = sfs_sim[order(sfs_sim[,'dataset']),]
	####	sfs_sim = sfs_sim[, -c(which(colnames(sfs_sim)=='fA0_fB0'), which(colnames(sfs_sim)=='fA1_fB0'), which(colnames(sfs_sim)=='fA0_fB1'))] # remove the singletons and the (0,0)

	####	sfs_sim_tmp[(rep*nSimulations+1):((rep+1)*nSimulations),] = as.matrix(sfs_sim)
	####}

	##### params
	####tmp_params = as.matrix(fread(paste(timeStamp, '/', sub_dir_sim, '/', model, '_', rep, '/priorfile.txt', sep=''), h=T))
    ####tmp_params = tmp_params[order(tmp_params[,'dataset']),-1]
	####params_sim_tmp[(rep*nSimulations+1):((rep+1)*nSimulations),] = as.matrix(tmp_params)
	}
    ss_sim_tmp = do.call(rbind,ss_sim_tmp)
    if(useSFS==1){    sfs_sim_tmp = do.call(rbind,sfs_sim_tmp)}
    params_sim_tmp = do.call(rbind,params_sim_tmp)
	
	# remove stats
	for(stat in undesired_stats){
		ss_sim_tmp = ss_sim_tmp[, -grep(stat, colnames(ss_sim_tmp))]
	}
	# end of remove stats
	
	# statistics
	if(useSFS == 1){
		ss_sim[[model]] = cbind(ss_sim_tmp, sfs_sim_tmp)  # with SFS
	}else{
		ss_sim[[model]] = ss_sim_tmp  # without SFS
	}

	# remove the uninformative statistics
	toRemove = c(1)
	for(i in 2:ncol(ss_sim[[model]])){
		sd_tmp = sd(as.numeric(c(ss_obs[i], ss_sim[[model]][,i])))
		
		if(sd_tmp<0.00001){
			toRemove = c(toRemove, i)
		}
	}
	
	toRemove = unique(toRemove)
	
	ss_obs = ss_obs[-toRemove]
	ss_sim[[model]] = ss_sim[[model]][, -toRemove]

	# params
	params_sim[[model]] = params_sim_tmp
	nparams = ncol(params_sim[[model]])

	### write the prior
	nPrior = nrow(params_sim[[model]])
	if(nPrior>10000){
		write.table(params_sim[[model]][1:10000,], paste(output_dir, '/priorfile.txt', sep=''), col.names=T, row.names=F, quote=F, sep='\t')
	}else{
		write.table(params_sim[[model]], paste(output_dir, '/priorfile.txt', sep=''), col.names=T, row.names=F, quote=F, sep='\t')
	}
	
	##############
	# inferences
	target_rf = ss_obs 

	# RANDOM FOREST	
        library(abcrf)
        params_model_rf = as.matrix(params_sim[[model]]) # in case of debug
        stats_model_rf = ss_sim[[model]] # in case of debug

        
        
        write("param\tHPD2.5%\tmedian\tHPD%97.5", paste(output_dir, '/posterior_summary_RandomForest_', sub_dir_model, '.txt', sep=''), append=F)
        res_rf = list()
        header_RF = NULL # contains the header of the file posterior_RandomForest
        estimates_RF = NULL # contains the point estimates of the posterior_RandomForest
        for(i in 1:nparams){
        	parameter = params_model_rf[,i]
        	data = data.frame(parameter, stats_model_rf)
        	mod = regAbcrf(parameter~., data, ntree=1000, paral=T, ncores=ncores)
        	estimate = predict(mod, target_rf, data, paral=T, ncores=ncores)

        	param_name = colnames(params_sim[[model]])[i]
        	res_rf[[param_name]] = list()
        	res_rf[[param_name]][['expectation']] = estimate$expectation
        	res_rf[[param_name]][['variance']] = estimate$variance
        	res_rf[[param_name]][['quantile025']] = estimate$quantiles[,1]
        	res_rf[[param_name]][['quantile975']] = estimate$quantiles[,2] 
        	write(paste(c(param_name, round(estimate$quantiles[1],5), round(estimate$expectation, 5), round(estimate$quantiles[2],5)), collapse="\t"), paste(output_dir,'/posterior_summary_RandomForest_', sub_dir_model, '.txt', sep=''), append=T)
        	
        	header_RF = c(header_RF, param_name)
        	estimates_RF = cbind(estimates_RF, estimate$expectation)
        }

        colnames(estimates_RF) = header_RF	
        write.table(estimates_RF, paste(output_dir,'/posterior_RandomForest_', sub_dir_model, '.txt', sep=''), sep="\t", col.names=T, row.names=F, quote=F )
        	
	res_tot = list()
	
	if(do_nnet == TRUE){	
		# NEURAL NETWORK
		library('nnet')
		target = matrix(as.numeric(unlist(ss_obs)), nrow=1)
		x = matrix(as.numeric(unlist(params_sim[[model]])), byrow=F, ncol=ncol(params_sim[[model]]))

		sumstat = matrix(as.numeric(unlist(ss_sim[[model]])), byrow=F, ncol=ncol(ss_sim[[model]]))
		transf_obs = rep("logit", ncol(params_sim[[model]]))
		
		for(param_i in 1:ncol(x)){
			prior_values = sample(x[,param_i], 1000, replace=T)
			if(min(x[,param_i]) <= 0){
				transf_obs[param_i] = 'none'
			}else if(length(table(prior_values)) <= 3){
					transf_obs[param_i] = 'log'
			}
		}
		
		bb = rbind(apply(x, MARGIN=2, FUN="min"), apply(x, MARGIN=2, FUN="max"))
		
		res = abc_nnet_multivar(target=target, x=x, sumstat=sumstat, tol=nPosterior/nrow(x), rejmethod=F, noweight=F, transf=transf_obs, bb=bb, nb.nnet=10, size.nnet=10, trace=T)
	
		posterior = res$x
		colnames(posterior) = colnames(params_sim[[model]])
		write.table(posterior, paste(output_dir,'/posterior_', sub_dir_model, '.txt', sep=''), row.names=F, col.names=T, sep='\t', quote=F)
		res_tot[['neural_network']] = posterior	
	}
	
#	res_tot[['random_forest']] = res_rf

	# plot pdf
	if(figure==T && do_nnet==T){
		library(ggplot2)
		library(ggpubr)
		theme_set(theme_classic())
		figure = list()

		for(i in 1:nparams){
			param_name = colnames(params_sim[[model]])[i]
			
			scale = 1
			if(param_name == 'N1' || param_name == 'N2' || param_name == 'Na'){
				scale = Nref
			}else if(param_name == 'Tsplit' || param_name == 'Tam' || param_name == 'Tsc' || param_name == 'Tdem1' || param_name == 'Tdem2'){ scale = 4*Nref }
			
			prior = x[,i] * scale
			posterior = res$x[,i] * scale
			prior = data.frame(x = prior, label=rep('prior', length(prior)))
			posterior = data.frame(x = posterior, label=rep('posterior', length(posterior)))
			df=rbind(prior, posterior)

			rf_estimate = as.numeric(res_rf[[param_name]][['expectation']]) * scale
			rf_estimate_inf = as.numeric(res_rf[[param_name]][['quantile025']]) * scale
			rf_estimate_sup = as.numeric(res_rf[[param_name]][['quantile975']]) * scale
			
			
			pp = ggdensity(df, x='x', fill='label') +
				theme(axis.text=element_text(size=15), axis.title=element_text(size=15), legend.text=element_text(size=15)) +
#				geom_vline(xintercept = rf_estimate, color = "black", size=1.25) +
#				geom_vline(xintercept = rf_estimate_inf, color = "black", size=0.25, linetype="dashed") +
#				geom_vline(xintercept = rf_estimate_sup, color = "black", size=0.25, linetype="dashed") +
				labs(fill = "") +
				scale_x_continuous(name = param_name) +
				scale_y_continuous(name = 'density') +
				scale_fill_manual(values = c("#f7f7f7", "#f1a340"))
			
			figure[[param_name]] = pp
		}
		ggarrange(plotlist=figure, common.legend = TRUE, labels='AUTO', align='hv')
		ggsave(paste(path2observation,'/posterior_', sub_dir_model, '.pdf', sep=''), bg='white', width=20, height=10)
	}
	
	# retur inferences	
	return(res_tot)
}

