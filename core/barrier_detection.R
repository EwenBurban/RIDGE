#!/usr/bin/env Rscript
library(abcrf)
library(ggpubr)
args = commandArgs(trailingOnly=T)
args = sapply(args,function(x){tmp = unlist(strsplit(x,split='='))
            y = vector()
            y[tmp[1]] = tmp[2]
            return(y)},USE.NAMES=F)

ntree=as.numeric(args['ntree'])
ncores=as.numeric(args['ncores'])
obs_dir=args['obs_dir']
sim_dir=args['sim_dir']
mode=args['mode']
roc_smooth=0.001 # the lower the value is, the smoother the roc will be
n_subset=5
posterior_file=args['posterior']
posterior=read.table(posterior_file,h=T)
pbarrier=sapply(1:nrow(posterior),function(x,...) weighted.mean(posterior[x,c('PbarrierM_ancestral','PbarrierM_current')],w=c(posterior[x,'Tsplit']-posterior[x,'Tam'],posterior[x,'Tsc'])))
mean_prior_ratio=mean((1-pbarrier)/(pbarrier),na.rm=T)

######## defining functions ############
convert2bar = function(mvec){
	bar_vec = ifelse(mvec == 0, 'barrier','non-barrier')
	return(as.factor(bar_vec))
}


param2kp=function(table,cv_threshold=0.1){
	# a function to filter the summary statistic without enough variation
	# if a summary stat have a low variation, abcrf will produce an error and stop
	cv_vec=apply(table,2,sd,na.rm=T)/colMeans(table,na.rm=T)
	if(any(cv_vec<cv_threshold)){
		sel_colnames=c(colnames(table)[-which(cv_vec<cv_threshold)])
	} else {sel_colnames=colnames(table)}
	return(sel_colnames)
}

get_roc = function(prior,est,trsh,...){
	TP = length(which(prior==0 & est>=trsh))
	FP = length(which(prior>0 & est>=trsh))
	TN = length(which(prior>0 & est<trsh))
	FN = length(which(prior==0 & est<trsh))
	TPR = TP / (TP + FN)
	FPR = FP / (FP + TN)
	PPV = TP / (TP + FP)
	nloci= length(which(est>=trsh))
	y=c('TPR'=TPR,'FPR'=FPR,'PPV'=PPV,'nloci'=nloci)
	return(y)
}

auc <-function(x, y, from = min(x, na.rm=TRUE), to = max(x, na.rm=TRUE),absolutearea=F,...)
{
    # Sanity checks
    stopifnot(length(x) == length(y))
    stopifnot(!is.na(from))
    if (length(unique(x)) < 2){
        return(NA)}
	# perform linear interpolation 
	if(absolutearea==F){
	values <- approx(x, y, xout = sort(unique(c(from, to, x[x > from & x < to]))), ...)
	res <- 0.5 * sum(diff(values$x) * (values$y[-1] + values$y[-length(values$y)]))
	}else{
	 o <- order(x)
	ox <- x[o]
	oy <- y[o]

	idx <- which(diff(oy >= 0)!=0)
	newx <- c(x, x[idx] - oy[idx]*(x[idx+1]-x[idx]) / (y[idx+1]-y[idx]))
	newy <- c(y, rep(0, length(idx)))
	values <- approx(newx, newy, xout = sort(unique(c(from, to, newx[newx > from & newx < to]))), ...)
	res <- 0.5 * sum(diff(values$x) * (abs(values$y[-1]) + abs(values$y[-length(values$y)])))
	}
	return(res)

}

#### load observed dataset ###############$
obs_data = read.table(file.path(obs_dir,'ABCstat_locus.txt'),h=T)
if(mode=='test'){
	obs_prior=read.table(file.path(obs_dir,'priorfile_locus.txt'),h=T)
	col=colnames(obs_prior)
	if(any(colnames(obs_prior)=='M_current')){migration='M_current';print('M_current')
	} else if (any(colnames(obs_prior)=='M_ancestral')){migration='M_ancestral';print('M_ancestral')}else{
		migration='null_mig'; obs_prior[,'null_mig']=rep(1000,nrow(obs_prior));print('no_mig')}
	if(length(unique(obs_prior[,migration]))==1){
		migration='null_mig'; obs_prior[,'null_mig']=rep(1000,nrow(obs_prior));print('no_mig')}
	obs_all=merge(obs_data,obs_prior,by='dataset')
}

### load Training set and set up for next usages
train_data =read.table(file.path(sim_dir,'ABCstat_locus.txt'),h=T)
train_prior =read.table(file.path(sim_dir,'priorfile_locus.txt'),h=T)

col_data=colnames(subset(train_data,select=-c(dataset,seed)))
col_prior=colnames(subset(train_prior,select=-c(dataset,seed)))

all = merge(train_data,train_prior,by='seed')
col_data=param2kp(all[,col_data]) ## curate the dataset

obs_roc_all=all

### obs roc & prediction
train_set=obs_roc_all[,col_data]
prior_train_set=obs_roc_all[,col_prior]
if(migration=='null_mig'){prior_train_set[,migration]=rep(1000,nrow(prior_train_set))}
status=	convert2bar(prior_train_set[,migration])
train_set=data.frame(train_set,status)
rf = abcrf(status~.,data=train_set,ncores = ncores,ntree = ntree,paral = T,lda=F)
pred = predict(rf,training = train_set,obs = obs_data[,col_data],paral = T,ncores=ncores)
if(mode=='test'){
	obs_prediction=data.frame('allocation'=pred$allocation,'post.prob'=pred$post.prob,'true_status'=convert2bar(obs_prior[,migration]),'migration'=obs_prior[,migration])
	obs_prediction$post.prob[which(obs_prediction$allocation=='non-barrier')]=1-obs_prediction$post.prob[which(obs_prediction$allocation=='non-barrier')]
	
	threshold_seq=seq(0,1,by=roc_smooth)
	roc_stat=as.data.frame(do.call(rbind,lapply(threshold_seq,function(x,...) get_roc(prior=obs_prediction$migration,est=obs_prediction$post.prob,trsh=x))))
	if(migration=='null_mig'){obs_AUC=0.5}else{
	obs_AUC=auc(roc_stat$FPR,roc_stat$TPR,absolutearea=T)}
	if(is.na(obs_AUC)){obs_AUC=mean(roc_stat$TPR)/(mean(roc_stat$TPR)+mean(roc_stat$FPR))}

	p_barrier_obs_obs=length(which(obs_prediction$true_status=='barrier'))/nrow(obs_prediction)
	#### bayes part ###
	obs_prediction$bayes_factor=mean_prior_ratio * (obs_prediction$post.prob/(1-obs_prediction$post.prob))
	obs_bayes_stat=get_roc(prior=obs_prior[,migration],est=obs_prediction$bayes_factor,trsh=1.000001)
	names(obs_bayes_stat)=paste0('obs_bayes_',names(obs_bayes_stat))
	threshold_seq=seq(0,10,by=roc_smooth)
	obs_roc_bayes_stat=as.data.frame(do.call(rbind,lapply(threshold_seq,function(x,...) get_roc(prior=obs_prior[,migration],est=obs_prediction$bayes_factor,trsh=x))))
	if(migration=='null_mig'){obs_bayes_AUC=0.5}else
	{obs_bayes_AUC=auc(obs_roc_bayes_stat$FPR,obs_roc_bayes_stat$TPR,absolutearea=F)}
	if(is.na(obs_bayes_AUC)){obs_bayes_AUC=mean(obs_roc_bayes_stat$TPR)/(mean(obs_roc_bayes_stat$TPR)+mean(obs_roc_bayes_stat$FPR))}
} else {
	obs_prediction=data.frame('allocation'=pred$allocation,'post.prob'=pred$post.prob)
	obs_prediction$post.prob[which(obs_prediction$allocation=='non-barrier')]=1-obs_prediction$post.prob[which(obs_prediction$allocation=='non-barrier')]
	obs_prediction$bayes_factor=mean_prior_ratio * (obs_prediction$post.prob/(1-obs_prediction$post.prob))
}

p_barrier_obs_est=length(which(obs_prediction$allocation=='barrier'))/nrow(obs_prediction)

####### generate output

write.table(obs_prediction,file=file.path(obs_dir,'Pbarrier.txt'),sep='\t',row.names=F)
### report
report=c('obs_bayes_p'=length(which(obs_prediction$bayes_factor>1))/nrow(obs_data),'p_barrier_obs_est'=p_barrier_obs_est)
if (mode=='test'){
	report=c(report,'obs_AUC'=obs_AUC,'p_barrier_obs_obs'=p_barrier_obs_obs,'obs_bayes_AUC'=obs_bayes_AUC,obs_bayes_stat)
	write.table(roc_stat,file=file.path(obs_dir,'true_roc_table.txt'),sep='\t',row.names=F)

}
write.table(t(report),file.path(obs_dir,'report_barrier_detection.txt'),sep='\t',row.names=F)
	
