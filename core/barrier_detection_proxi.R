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
timeStamp=args['timeStamp']
sim_dir=args['sim_dir']
mode=args['mode']
graph=args['graph'] # expect 'True' to produce graphical output
migration='M_current'
roc_smooth=0.001 # the lower the value is, the smoother the roc will be
n_subset=5
thr=mig_threshold=0.5
######## defining functions ############
convert2bar = function(mvec){
	bar_vec = ifelse(mvec == 0, 'barrier','non-barrier')
	return(as.factor(bar_vec))
}

convert2mig = function(mvec,thr){
	bar_vec=sapply(mvec, function(x) if(x==0){'barrier'} else if(x>0 & x<thr){'pseudo-barrier'} else {'non-barrier'})
	return(as.factor(bar_vec))
}

resample_dataset=function(data,f,prob=NULL){
	rprob=table(f)/length(f)
	prob=prob[names(rprob)]
	fprob=prob/rprob
	print(fprob)	
	m_data=data[sample(1:nrow(data),size=nrow(data),replace=T,prob=as.numeric(fprob[f])),]
	return(m_data)


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
	obs_all=merge(obs_data,obs_prior,by='dataset')
}

### load Training set and set up for next usages
train_data =read.table(file.path(sim_dir,'ABCstat_locus.txt'),h=T)
train_prior =read.table(file.path(sim_dir,'priorfile_locus.txt'),h=T)

col_data=colnames(subset(train_data,select=-c(dataset,seed)))
col_prior=colnames(subset(train_prior,select=-c(dataset,seed)))

all = merge(train_data,train_prior,by='seed')
col_data=param2kp(all[,col_data]) ## curate the dataset

half_sample=sample(1:nrow(all),size=round(nrow(all)*0.5),replace=F)
sim_roc_all=all[half_sample,]
obs_roc_all=all[-half_sample,]

#### sim roc
spl_all=split.data.frame(sim_roc_all,f=sample(1:n_subset,size=nrow(sim_roc_all),replace=T))
sim_prediction_ll=lapply(1:n_subset,function(I,...){
	print('started')	
	test_set = spl_all[[I]][,col_data]
	prior_test_set=spl_all[[I]][,col_prior]
	tmp_spl_bind=do.call(rbind,spl_all[-I])	

	print('started')	
	tmp_spl_bind_table1=resample_dataset(tmp_spl_bind,f=convert2mig(tmp_spl_bind[,migration],thr),prob=c('barrier'=0.25,'pseudo-barrier'=0.25,'non-barrier'=0.5))
	train_set_table1 = tmp_spl_bind_table1[,col_data]
	prior_train_set_table1 = tmp_spl_bind_table1[,col_prior]

	print('started')	
	status=	as.factor(ifelse(prior_train_set_table1>thr,'high_migration','low_migration'))
	train_set=data.frame(train_set_table1,status,row.names=NULL)
	rf = abcrf(status~.,data=train_set,ncores = ncores,ntree = ntree,paral = T,lda=F)
	pred = predict(rf,training = train_set,obs = test_set,paral = T,ncores=ncores)
		
	pred_data=data.frame('low_mig_prob'=pred$post.prob,'true_status'=convert2mig(prior_test_set[,migration],thr),'allocation_low_mig'=pred$allocation)
	pred_data$low_mig_prob[which(pred_data$allocation_low_mig=='high_migration')]=1-pred_data$low_mig_prob[which(pred_data$allocation_low_mig=='high_migration')]
	

	tmp_spl_bind_table2=resample_dataset(tmp_spl_bind,f=convert2mig(tmp_spl_bind[,migration],thr),prob=c('barrier'=0.5,'pseudo-barrier'=0.5,'non-barrier'=0))
	train_set_table2 = tmp_spl_bind_table2[,col_data]
	prior_train_set_table2 = tmp_spl_bind_table2[,col_prior]

	status=	as.factor(ifelse(prior_train_set_table2==0,'barrier','non-barrier'))
	train_set=data.frame(train_set_table2,status,row.names=NULL)
	rf = abcrf(status~.,data=train_set,ncores = ncores,ntree = ntree,paral = T,lda=F)
	pred = predict(rf,training = train_set,obs = test_set,paral = T,ncores=ncores)
	
	pred_data=data.frame(pred_data,'barrier_prob'=pred$post.prob,'allocation'=pred$allocation)
	pred_data$barrier_prob[which(pred_data$allocation=='non-barrier')]=1-pred_data$barrier_prob[which(pred_data$allocation=='non-barrier')]
	pred_data$post.prob=pred_data$low_mig_prob * pred_data$barrier_prob
	print(paste0('subset ',I,' done'))
	return(pred_data)
		})
sim_prediction=do.call(rbind,sim_prediction_ll)

threshold_seq=seq(0,1,by=roc_smooth)
sim_roc_stat=as.data.frame(do.call(rbind,lapply(threshold_seq,function(x,...) get_roc(prior=sim_prediction$migration,est=sim_prediction$post.prob,trsh=x))))
sim_AUC=auc(sim_roc_stat$FPR,sim_roc_stat$TPR,absolutearea=T)
if(is.na(sim_AUC)){sim_AUC=mean(sim_roc_stat$TPR)/(mean(sim_roc_stat$TPR)+mean(sim_roc_stat$FPR))}
p_barrier_sim_obs=length(which(sim_prediction$true_status=='barrier'))/nrow(sim_prediction)
p_barrier_sim_est=length(which(sim_prediction$allocation=='barrier'))/nrow(sim_prediction)


### obs roc & prediction
obs_roc_all_table1=resample_dataset(obs_roc_all,f=convert2mig(obs_roc_all[,migration],thr),prob=c('barrier'=0.25,'pseudo-barrier'=0.25,'non-barrier'=0.5))
obs_roc_all_table2=resample_dataset(obs_roc_all,f=convert2mig(obs_roc_all[,migration],thr),prob=c('barrier'=0.5,'pseudo-barrier'=0.5,'non-barrier'=0))

train_set_table1=obs_roc_all_table1[,col_data]
prior_train_set_table1=obs_roc_all_table1[,col_prior]
status=	as.factor(ifelse(prior_train_set_table1>thr,'high_migration','low_migration'))
train_set=data.frame(train_set_table1,status,row.names=NULL)
rf = abcrf(status~.,data=train_set,ncores = ncores,ntree = ntree,paral = T,lda=F)
pred1 = predict(rf,training = train_set,obs = obs_data[,col_data],paral = T,ncores=ncores)


train_set_table2=obs_roc_all_table2[,col_data]
prior_train_set_table2=obs_roc_all_table2[,col_prior]
status=	as.factor(ifelse(prior_train_set_table2==0,'barrier','non-barrier'))
train_set=data.frame(train_set_table2,status,row.names=NULL)
rf = abcrf(status~.,data=train_set,ncores = ncores,ntree = ntree,paral = T,lda=F)
pred2 = predict(rf,training = train_set,obs = obs_data[,col_data],paral = T,ncores=ncores)

if(mode=='test'){
	obs_prediction=data.frame('allocation'=pred2$allocation,'allocation_low_mig'=pred1$allocation,
							  'low_mig_prob'=pred1$post.prob,'barrier_prob'=pred2$post.prob,
							  'true_status'=convert2mig(obs_prior[,migration],thr),'migration'=obs_prior[,migration])
	obs_prediction$barrier_prob[which(obs_prediction$allocation=='non-barrier')]=1-obs_prediction$barrier_prob[which(obs_prediction$allocation=='non-barrier')]
	obs_prediction$low_mig_prob[which(obs_prediction$allocation_low_mig=='high_migration')]=1-obs_prediction$low_mig_prob[which(obs_prediction$allocation_low_mig=='high_migration')]
	obs_prediction$post.prob=obs_prediction$low_mig_prob * obs_prediction$barrier_prob

	threshold_seq=seq(0,1,by=roc_smooth)
	roc_stat=as.data.frame(do.call(rbind,lapply(threshold_seq,function(x,...) get_roc(prior=obs_prediction$migration,est=obs_prediction$post.prob,trsh=x))))
	obs_AUC=auc(roc_stat$FPR,roc_stat$TPR,absolutearea=T)
	if(is.na(obs_AUC)){obs_AUC=mean(roc_stat$TPR)/(mean(roc_stat$TPR)+mean(roc_stat$FPR))}

	p_barrier_obs_obs=length(which(obs_prediction$true_status=='barrier'))/nrow(obs_prediction)
} else {
	obs_prediction=data.frame('allocation'=pred2$allocation,'allocation_low_mig'=pred1$allocation,
							  'low_mig_prob'=pred1$post.prob,'barrier_prob'=pred2$post.prob)
	obs_prediction$barrier_prob[which(obs_prediction$allocation=='non-barrier')]=1-obs_prediction$barrier_prob[which(obs_prediction$allocation=='non-barrier')]
	obs_prediction$low_mig_prob[which(obs_prediction$allocation_low_mig=='high_migration')]=1-obs_prediction$low_mig_prob[which(obs_prediction$allocation_low_mig=='high_migration')]
	obs_prediction$post.prob=obs_prediction$low_mig_prob * obs_prediction$barrier_prob
}

p_barrier_obs_est=length(which(obs_prediction$allocation=='barrier'))/nrow(obs_prediction)

####### generate output

write.table(obs_prediction,file=file.path(obs_dir,'Pbarrier_proxi.txt'),sep='\t',row.names=F)
write.table(sim_roc_stat,file=file.path(obs_dir,'roc_table_proxi.txt'),sep='\t',row.names=F)
### report
report=c('sim_AUC'=sim_AUC,'p_barrier_sim_obs'=p_barrier_sim_obs,'p_barrier_sim_est'=p_barrier_sim_est,
		 'p_barrier_obs_est'=p_barrier_obs_est)
if (mode=='test'){
	report=c(report,'obs_AUC'=obs_AUC,'p_barrier_obs_obs'=p_barrier_obs_obs)
	write.table(roc_stat,file=file.path(obs_dir,'true_roc_table_proxi.txt'),sep='\t',row.names=F)

}
write.table(t(report),file.path(obs_dir,'report_barrier_detection_proxi.txt'),sep='\t',row.names=F)
	
### graphical output   
if (graph=='True'){

	if(mode=='test'){
		all_roc=merge(sim_roc_stat,roc_stat,by='FPR',all=T)
		p=ggscatter(all_roc,'FPR',c('TPR.x','TPR.y'),merge=T,add='loess',size=0)
		all_roc=merge(sim_roc_stat,roc_stat,by='PPV',all=T)
		ppv=ggscatter(all_roc,'PPV',c('nloci.x','nloci.y'),merge=T,add='loess',size=0,xlim=c(0,1)) + geom_vline(xintercept=0.95)
	} else {
		p=ggscatter(sim_roc_stat,'FPR','TPR',add='loess',size=0)

		ppv_vec=sim_roc_stat$PPV
		test=max(ppv_vec)
		if(test > 0.95){ 
			yy=which(ppv_vec >= 0.95 & ppv_vec < 1)[1] 
			ppv_thr=0.95 
		}else{
				yy=which(ppv_vec==test) 
				ppv_thr=test
		}
		nloci_ppv=sim_roc_stat[yy,'nloci']
		ppv=ggscatter(sim_roc_stat,'PPV','nloci',add='loess',size=0,xlim=c(0,1)) +
			geom_vline(xintercept=0.95) + 
			annotate('text',x=0.7,y=nloci*0.8,paste0('for PPV=',ppv_thr,'nloci=',nloci_ppv)) 

	}

	ggsave(ggarrange(p,ppv),filename=file.path(obs_dir,'barrier_plot.png'),width=14)
}

