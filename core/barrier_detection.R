#!/usr/bin/env Rscript
library(abcrf)
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
barrier_type=''
roc_smooth=0.001 # the lower the value is, the smoother the roc will be
n_subset=5
posterior_file=args['posterior']
posterior=read.table(posterior_file,h=T)
pbarrier=sapply(1:nrow(posterior),function(x,...){ 

					# this line weight the barrier proportion estimated for ancestral and current migration by their respective time of effect during divergence
					# TO-DO : directly implement the pbarrier parameter in ABC model, allowing direct estimation of this parameter rather than estimating it at posteriori - E.BURBAN 9/11/2023
				weighted.mean(posterior[x,c('PbarrierM_ancestral','PbarrierM_current')],w=c(posterior[x,'Tsplit']-posterior[x,'Tam'],posterior[x,'Tsc']))})
## in case where there is no barrier, weighted mean produce NaN. So NaN are replace by 0
pbarrier[is.nan(pbarrier)]=0
######## pbarrier ratio computation (pbarrier is called Q in text)
# correct way is to compute the average ratio (1-Q)/Q but this way produce Inf when Q=0. 
# In model averaging, by construction some  posteriors might have a Q=0. 
# Consequently, we remove Inf value from the mean and so increase the mean value. 
# To avoid this pitfall, an aproximate way is to compute the ratio as mean(1-Q)/mean(Q) + var(Q)/mean(Q)³
correct_pbarrier_ratio=(1-pbarrier)/(pbarrier)
correct_mean_pbarrier_ratio=mean(correct_pbarrier_ratio[is.infinite(correct_pbarrier_ratio)==F],na.rm=T)
aproxmate_mean_pbarrier_ratio=(1-mean(pbarrier))/mean(pbarrier) + var(pbarrier)/(mean(pbarrier)^3)

#mean_prior_ratio=

######## defining local functions ############
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
	# AUC function from MASS R package
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

#### load observed dataset ###############
obs_data = read.table(file.path(obs_dir,'ABCstat_locus.txt'),h=T)
if(mode=='test'){
	# when applied to a pseudo-observed dataset where the true status of loci is know, this part define the column where to find migration value
	# and provide a migration column when there is not. 
	obs_prior=read.table(file.path(obs_dir,'priorfile_locus.txt'),h=T)
	col=colnames(obs_prior)
	if(any(colnames(obs_prior)=='M_current')){migration='M_current';print('M_current')
	} else if (any(colnames(obs_prior)=='M_ancestral')){migration='M_ancestral';print('M_ancestral')}else{
		migration='null_mig'; obs_prior[,'null_mig']=rep(0,nrow(obs_prior));print('no_mig')}
	obs_all=merge(obs_data,obs_prior,by='dataset')
}else {
		migration='M_current'
}

### load Training set and set up for next usages
train_data =read.table(file.path(sim_dir,'ABCstat_locus.txt'),h=T)
train_prior =read.table(file.path(sim_dir,'priorfile_locus.txt'),h=T)

col_data=colnames(subset(train_data,select=-c(dataset,seed)))
col_prior=colnames(subset(train_prior,select=-c(dataset,seed)))

all = merge(train_data,train_prior,by='seed')
nb_barr=length(which(all$M_current==0))
nb_nonbarr=length(which(all$M_current!=0))
print(nb_barr/nrow(all))
if (nb_barr > nb_nonbarr){
	all_bar=which(all$M_current==0)
	all_bar=all[all_bar[sample(1:nb_barr,size=nb_nonbarr)],]
	all_nonbar=all[which(all$M_current!=0),]
	all=rbind(all_bar,all_nonbar)
} else if (nb_barr < nb_nonbarr){
	all_nonbar=which(all$M_current!=0)
	all_nonbar=all[all_nonbar[sample(1:nb_nonbarr,size=nb_barr)],]
	all_bar=all[which(all$M_current==0),]
	all=rbind(all_bar,all_nonbar)
}
#col_data=param2kp(all[,col_data]) ## curate the dataset

nb_barr=length(which(all$M_current==0))
print(nb_barr/nrow(all))

obs_roc_all=all

### obs roc & prediction
train_set=obs_roc_all[,col_data]
prior_train_set=obs_roc_all[,col_prior]
#if(migration=='null_mig'){prior_train_set[,migration]=rep(1000,nrow(prior_train_set))}
if(migration=='null_mig'){status=convert2bar(prior_train_set[,'M_current'])}else {status=	convert2bar(prior_train_set[,migration])}
train_set=data.frame(train_set,status)
rf = abcrf(status~.,data=train_set,ncores = ncores,ntree = ntree,paral = T,lda=F)
pred = predict(rf,training = train_set,obs = obs_data[,col_data],paral = T,ncores=ncores)
if(mode=='test'){
	obs_prediction=data.frame('allocation'=pred$allocation,'post.prob'=pred$post.prob,
							  'true_status'=convert2bar(obs_prior[,migration]),'migration'=obs_prior[,migration],
							  'tree_vote'=pred$vote[,'barrier']/rowSums(pred$vote)) ### new feature !!!! to test ! Reminder : line 144 there is another occurance of this code
	obs_prediction$post.prob[which(obs_prediction$allocation=='non-barrier')]=1-obs_prediction$post.prob[which(obs_prediction$allocation=='non-barrier')]
	### roc curve computation	
	threshold_seq=seq(0,1,by=roc_smooth)
	roc_stat=as.data.frame(do.call(rbind,lapply(threshold_seq,function(x,...) get_roc(prior=obs_prediction$migration,est=obs_prediction$post.prob,trsh=x))))
	if(migration=='null_mig'){obs_AUC=NA}else{obs_AUC=auc(roc_stat$FPR,roc_stat$TPR,absolutearea=T)}
	if(is.na(obs_AUC)){obs_AUC=mean(roc_stat$TPR)/(mean(roc_stat$TPR)+mean(roc_stat$FPR))}
	### computing bayes factor ####
	obs_prediction$bayes_factor=correct_mean_pbarrier_ratio * (obs_prediction$post.prob/(1-obs_prediction$post.prob))
	p_barrier_obs_obs=length(which(obs_prediction$true_status=='barrier'))/nrow(obs_prediction)
	#### bayes part ###
} else {
	obs_prediction=data.frame('allocation'=pred$allocation,'post.prob'=pred$post.prob,'tree_vote'=pred$vote[,'barrier']/rowSums(pred$vote))
	obs_prediction$post.prob[which(obs_prediction$allocation=='non-barrier')]=1-obs_prediction$post.prob[which(obs_prediction$allocation=='non-barrier')]

}

### computing bayes factor ####
obs_prediction$bayes_factor=correct_mean_pbarrier_ratio * (obs_prediction$post.prob/(1-obs_prediction$post.prob))
obs_prediction$BF_approxQ=aproxmate_mean_pbarrier_ratio *  (obs_prediction$post.prob/(1-obs_prediction$post.prob))



####### generate output

write.table(cbind(obs_data,obs_prediction),file=file.path(obs_dir,paste0('Pbarrier.txt')),sep='\t',row.names=F)
### report
#report=c('obs_bayes_p'=length(which(obs_prediction$bayes_factor>100))/nrow(obs_data),'p_barrier_obs_est'=p_barrier_obs_est)
report=c('average_pbarrier'=mean(pbarrier),'pbarrier_ratio'=correct_mean_pbarrier_ratio,'aproximate_pbarrier_ratio'=aproxmate_mean_pbarrier_ratio)
if (mode=='test'){
	report=c(report,'obs_AUC'=obs_AUC,'true_pbarrier'=p_barrier_obs_obs)
	write.table(roc_stat,file=file.path(obs_dir,'true_roc_table.txt'),sep='\t',row.names=F)

}
write.table(t(report),file.path(obs_dir,'barrier_proportion_and_ratio.txt'),sep='\t',row.names=F)	
write.table(rf$model.rf$confusion.matrix,sep='\t',row.names=T,col.names=T,quote=F,file=file.path(obs_dir,'confusion_matrix_barrier.txt'))
write.table(rf$model.rf$variable.importance,sep='\t',row.names=T,col.names=T,quote=F,file=file.path(obs_dir,'variable_importance_barrier.txt'))


### gerenate control plot

test_vec=rf$model.rf$predictions == status
test_vec_barr=test_vec[which(status=='barrier')]
test_vec_nonbarr=test_vec[which(status=='non-barrier')]

get_error_class=function(ref_v,s,e,t_v){
	tt=t_v[which(ref_v>=s & ref_v<e)]
	y=length(which(tt==F))/length(t_v)
	return(y)}


h=hist(all$bialsite_avg,breaks=100)
hh=data.frame(s=h$breaks[-length(h$breaks)],e=h$breaks[-1])
error_class=apply(hh,1,function(se,...){get_error_class(all$bialsite_avg,s=se[1],e=se[2],t_v=test_vec)})
error_class_barrier=apply(hh,1,function(se,...){get_error_class(all$bialsite_avg,s=se[1],e=se[2],t_v=test_vec_barr)})
error_class_nonbarrier=apply(hh,1,function(se,...){get_error_class(all$bialsite_avg,s=se[1],e=se[2],t_v=test_vec_nonbarr)})

class_error_data=data.frame(bialsite_breaks=h$mids,error_class,error_class_barrier,error_class_nonbarrier)
write.table(class_error_data,file=file.path(obs_dir,'class_error_data.txt'),sep='\t',row.names=F,col.names=T)
