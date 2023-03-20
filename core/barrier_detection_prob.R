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
point_posterior_file=args['point_posterior']
point_posterior=read.table(point_posterior_file,h=T)
pbarrier=point_posterior[1,'PbarrierM_current']
print(pbarrier)
######## defining functions ############
convert2bar = function(mvec){
	bar_vec = ifelse(mvec == 0, 'barrier','non-barrier')
	return(as.factor(bar_vec))
}

convert2prob = function(mvec,pbarrier){
	prob_vec = ifelse(mvec == 0, pbarrier,1-pbarrier)
	return(prob_vec)
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

	test_set = spl_all[[I]][,col_data]
	prior_test_set=spl_all[[I]][,col_prior]
	train_set = do.call(rbind,spl_all[-I])[,col_data]
	prior_train_set = do.call(rbind,spl_all[-I])[,col_prior]
	status=	convert2bar(prior_train_set[,migration])
	train_set=data.frame(train_set,status)
	rf = abcrf(status~.,data=train_set,ncores = ncores,ntree = ntree,paral = T,lda=F)
	pred = predict(rf,training = train_set,obs = test_set,paral = T,ncores=ncores)
		
	pred_data=data.frame('allocation'=pred$allocation,'post.prob'=pred$post.prob,'true_status'=convert2bar(prior_test_set[,migration]),'migration'=prior_test_set[,migration])
	pred_data$post.prob[which(pred_data$allocation=='non-barrier')]=1-pred_data$post.prob[which(pred_data$allocation=='non-barrier')]
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
train_set=obs_roc_all[,col_data]
prior_train_set=obs_roc_all[,col_prior]
status=	convert2bar(prior_train_set[,migration])
train_set=data.frame(train_set,status)
###### train set prob #####
# this part generate a prior dataset which respect the expected proportion of barrier 
prob=convert2prob(obs_roc_all[,migration],pbarrier)
prob_obs_roc_all=obs_roc_all[sample(1:nrow(obs_roc_all),size=nrow(obs_roc_all),replace=T,prob=prob),]
status_prob=convert2bar(prob_obs_roc_all[,migration])
prob_train_set=data.frame(prob_obs_roc_all[,col_data],'status'=status_prob)
###### prediction
rf = abcrf(status~.,data=train_set,ncores = ncores,ntree = ntree,paral = T,lda=F)
pred = predict(rf,training = prob_train_set,obs = obs_data[,col_data],paral = T,ncores=ncores)

if(mode=='test'){
	obs_prediction=data.frame('allocation'=pred$allocation,'post.prob'=pred$post.prob,'true_status'=convert2bar(obs_prior[,migration]),'migration'=obs_prior[,migration])
	obs_prediction$post.prob[which(obs_prediction$allocation=='non-barrier')]=1-obs_prediction$post.prob[which(obs_prediction$allocation=='non-barrier')]

	threshold_seq=seq(0,1,by=roc_smooth)
	roc_stat=as.data.frame(do.call(rbind,lapply(threshold_seq,function(x,...) get_roc(prior=obs_prediction$migration,est=obs_prediction$post.prob,trsh=x))))
	obs_AUC=auc(roc_stat$FPR,roc_stat$TPR,absolutearea=T)
	if(is.na(obs_AUC)){obs_AUC=mean(roc_stat$TPR)/(mean(roc_stat$TPR)+mean(roc_stat$FPR))}

	p_barrier_obs_obs=length(which(obs_prediction$true_status=='barrier'))/nrow(obs_prediction)
} else {
	obs_prediction=data.frame('allocation'=pred$allocation,'post.prob'=pred$post.prob)
	obs_prediction$post.prob[which(obs_prediction$allocation=='non-barrier')]=1-obs_prediction$post.prob[which(obs_prediction$allocation=='non-barrier')]
}

p_barrier_obs_est=length(which(obs_prediction$allocation=='barrier'))/nrow(obs_prediction)

####### generate output

write.table(obs_prediction,file=file.path(obs_dir,'Pbarrier_prob.txt'),sep='\t',row.names=F)
write.table(sim_roc_stat,file=file.path(obs_dir,'roc_table_prob.txt'),sep='\t',row.names=F)
### report
report=c('sim_AUC'=sim_AUC,'p_barrier_sim_obs'=p_barrier_sim_obs,'p_barrier_sim_est'=p_barrier_sim_est,
		 'p_barrier_obs_est'=p_barrier_obs_est)
if (mode=='test'){
	report=c(report,'obs_AUC'=obs_AUC,'p_barrier_obs_obs'=p_barrier_obs_obs)
	write.table(roc_stat,file=file.path(obs_dir,'true_roc_table_prob.txt'),sep='\t',row.names=F)

}
write.table(t(report),file.path(obs_dir,'report_barrier_detection_prob.txt'),sep='\t',row.names=F)
	
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

