#!/usr/bin/env Rscript
library(abcrf)
mode = 'class' # class or reg
for(arg in commandArgs()){
  tmp = strsplit(arg, '='); opt = tmp[[1]][1]; val = tmp[[1]][2]
  if(opt == "ntree"){ntree=as.numeric(val)}
  if(opt == "ncores"){ncores=as.numeric(val)}
  if(opt == "obs_dir"){obs_dir=val}
  if(opt == "timeStamp"){timeStamp=val}
  if(opt == "sim_dir"){sim_dir=val}
  if(opt == "model"){m=val}
  if(opt == "param"){p=val}
  if(opt == "output"){output=val}
  
}
## custom func
remove_param = function(table){
	cv = apply(table,2,sd)/colMeans(table)
	y = which(abs(cv) < 0.1 | is.nan(cv)) # TODO check la pertinence de cet modif
	return(y)
}

convert2bar = function(mvec){
	bar_vec = ifelse(mvec == 0, 'barrier','non-barrier')
	return(as.factor(bar_vec))
}

### Param estimation ###
# observed data
obs_ABCstat_file = file.path(obs_dir,'ABCstat_locus.txt')
#obs_ABCjsfs_file = file.path(obs_dir,'ABCjsfs_locus.txt')
#test_data = cbind(read.table(obs_ABCstat_file,h=T),read.table(obs_ABCjsfs_file,h=T))
test_data = read.table(obs_ABCstat_file,h=T)
test_data[is.na(test_data)]=0
# Training set 
train_ABCstat_file_list = list.files(sim_dir,pattern = "ABCstat_locus.txt",recursive = T,full.names = T)
#train_ABCjsfs_file_list = list.files(sim_dir,pattern = "ABCjsfs_locus.txt",recursive = T,full.names = T)
train_priorfile_list = list.files(sim_dir,pattern = "priorfile_locus.txt",recursive = T,full.names = T)

train_targets = do.call(rbind,lapply(train_priorfile_list[grep(m,train_priorfile_list)],read.table,header=T,stringsAsFactors=F))

if( (grepl('SI',m) | grepl('AM',m) ) & grepl('M',p)){model_param_estimation = rep(0,nrow(test_data))} else {

	train_data = do.call(rbind,lapply(train_ABCstat_file_list[grep(m,train_ABCstat_file_list)],read.table,header=T,stringsAsFactors=F))
	train_data = subset(train_data,select=-dataset)
	train_targets = subset(train_targets,select=-dataset)
#	train_data_sfs = do.call(rbind,lapply(train_ABCjsfs_file_list[grep(m,train_ABCjsfs_file_list)],read.table,header=T,stringsAsFactors=F))
#	train_data = cbind(train_data,train_data_sfs)
	all = merge(train_data,train_targets,by='seed')
	train_data = all[,colnames(train_data)]
	train_targets = all[,colnames(train_targets)]
	train_data = subset(train_data,select = -seed)

	## curation
	train_data = as.data.frame(sapply(train_data,as.numeric))
	train_data[is.na(train_data)]=0
	train_targets = as.data.frame(sapply(train_targets,as.numeric))
	train_targets[is.na(train_targets)]=0


#	rejected_statistics = c('minDivAB_avg', 'minDivAB_std', 'maxDivAB_avg', 'maxDivAB_std', 'Gmin_avg', 'Gmin_std', 'Gmax_avg', 'Gmax_std','bialsites_avg','thetaA_avg','thetaB_avg','divAB_avg','ss_sf','noSS_sf','ss_noSf','noSS_noSf','succesive_ss_avg') 
#	rejected_statistics = c('datatset') 
#	rej_stat = na.omit(match(rejected_statistics,colnames(train_data)))
#	param2rm = union(remove_param(train_data),rej_stat)
	param2rm = vector()
	if(length(param2rm)>0){param2kp = colnames(train_data)[-param2rm]} else { param2kp=colnames(train_data)}
#	print(length(param2rm))

	stopifnot(any(grepl(p,colnames(train_targets))))

	### regression / classification ####
	locus = sample(1:nrow(train_data),size=(nrow(train_data)-1)/100,replace=F)
	if(mode == 'reg'){
	rf = regAbcrf(train_targets[locus,p]~.,data=train_data[locus,param2kp],ncores = ncores,ntree = ntree,paral = T)
	pred = predict(rf,training = train_data[locus,param2kp],obs = test_data,paral = T,ncores=ncores)
	model_param_estimation<-pred$expectation
	}
	if(mode == 'class'){
	rf = abcrf(convert2bar(train_targets[locus,p])~.,data=train_data[locus,param2kp],ncores = ncores,ntree = ntree,paral = T,lda=F)
	pred = predict(rf,training = train_data[locus,param2kp],obs = subset(test_data,select=-dataset),paral = T,ncores=ncores)
#	model_param_estimation = ifelse(as.vector(pred$allocation) == 'barrier',0,1)# reconvert into numeric value to be mean with other models 
	model_param_estimation = pred$vote[,'barrier'] / rowSums(pred$vote)
	}
}
model_param_estimation = cbind(model_param_estimation,test_data[,'dataset'])

write.table(model_param_estimation,col.names=c('Pbarrier','dataset'),row.names=F,file = output) 
