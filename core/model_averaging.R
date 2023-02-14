#!/usr/bin/env Rscript
#library(randomForest)
library(abcrf)
# Default values
obs_pattern = "ABCstat.txt"
sim_pattern = "ABCstat.txt"
ntree=1000
ncores=8
test=FALSE
for(arg in commandArgs()){
  tmp = strsplit(arg, '='); opt = tmp[[1]][1]; val = tmp[[1]][2]
  if(opt == "ntree"){ntree=as.numeric(val)}
  if(opt == "obs_dir"){obs_dir=val}
  if(opt == "obs_pattern"){obs_pattern=val}
  if(opt == "sim_pattern"){sim_pattern=val}
  if(opt == 'useSFS'){usesfs=as.numeric(val)} # 0 is F , 1 is T
  if(opt == "sim_dir"){sim_dir=val}
  if(opt == "output_name"){output_name=val}
  if(opt == "timeStamp"){timeStamp=val}
  if(opt == "ncores"){ncores=as.numeric(val)}
  if(opt == "mode"){mode=val}
}
### function
remove_param = function(table){
	cv = apply(table,2,sd)/colMeans(table)
	y = which(abs(cv) < 0.1 | is.nan(cv)) # TODO check la pertinence de cet modif
	return(y)
}

model_target_generator = function(file_dir){
	data = read.table(file_dir,h=T)
	dir = dirname(file_dir)
	model=sub('.*/','',dir)
	model=sub('N_.*','N',model)
	y=rep(model,nrow(data))
	return(y)
}
#### Model weighting #####

# Observed data
if(mode == 'multi'){
list_obs_dir = list.dirs(path=obs_dir,full.names = F,recursive=F)
print(list_obs_dir)
obs_ABCstat_file_list = file.path(obs_dir,list_obs_dir,obs_pattern)
test_data = do.call(rbind,lapply(obs_ABCstat_file_list,read.table,header=T))[,-1]
test=TRUE
}else if (mode == 'single') {
	test_data = read.table(file.path(obs_dir,obs_pattern),h=T)[,-1]
}

# Training set 
train_ABCstat_file_list = list.files(sim_dir,pattern = sim_pattern,recursive = T,full.names = T)
print(train_ABCstat_file_list)
if (usesfs==1){train_ABCsfs_file_list = list.files(sim_dir,pattern='ABCjsfs.txt',recursive=T,full.names=T)}

print(paste('The length of train list file is',length(train_ABCstat_file_list)))
train_models=sub("N_.*","N",list.dirs(path=sim_dir,full.names = F,recursive=F))
train_data = do.call(rbind,lapply(train_ABCstat_file_list,read.table,header=T,fill=T))[,-1]
train_target=factor(unlist(sapply(train_ABCstat_file_list,model_target_generator)))
print(length(train_target))
print(nrow(train_data))
stopifnot(length(train_target)==nrow(train_data))
#RF
train_data_2rf = train_data 
if (usesfs==1){
	train_data_sfs = do.call(rbind,lapply(train_ABCsfs_file_list,read.table,header=T))[,-1]
	train_data_2rf =  data.frame(train_data,train_data_sfs)}

nc_bef = ncol(train_data_2rf)
useless_param = remove_param(train_data_2rf)
if(length(useless_param)!=0){
train_data_2rf = train_data_2rf[,-useless_param] # remove the non-informative parameters
test_data = test_data[,-useless_param]
print(paste('the number of remove parameter is', nc_bef - ncol(train_data_2rf)))
print(paste('the number of remaining parameter is', ncol(train_data_2rf)))
}
test_data[is.na(test_data)] = 0
train_data_2rf[is.na(train_data_2rf)] = 0
train_data_2rf = data.frame(train_data_2rf,'t'=train_target)
rf = abcrf(t~.,train_data_2rf,paral=T,ncores=ncores,lda=F)


## get the vote %
votes = predict(rf,test_data,training=train_data_2rf)$vote
votes=votes/rowSums(votes) # normalize votes

if(test==TRUE){
votes = as.data.frame(votes)
stopifnot(nrow(votes)==length(list_obs_dir))
for(o in list_obs_dir){ # write all the result inside their respective directory (in case of multiple observed data)
	output=file.path(obs_dir,o,output_name)
	index = which(list_obs_dir==o)
	sub_votes = votes[index,,drop=F]
	stopifnot(nrow(sub_votes)==1)
	
	cumsort_v = cumsum(sort(sub_votes))
	usefull_model = names(cumsort_v)[which(cumsort_v >= 0.05)] # remove all the useless models which alltogether doesn’t contribute to result 
	usefull_votes = sub_votes[,usefull_model,drop=F]
	usefull_votes = usefull_votes/sum(usefull_votes) # rescale 
	write.table(usefull_votes, file=output, row.names=F,sep="\t" )
	}
}else{
	output=file.path(obs_dir,output_name)
	tmp = as.vector(votes)
	names(tmp) = colnames(votes)
	votes = tmp
	print(votes)
	cumsort_v = cumsum(sort(votes))
	usefull_model = names(cumsort_v)[which(cumsort_v >= 0.05)] # remove all the useless models which alltogether doesn’t contribute to result 
	usefull_votes = votes[usefull_model]
	usefull_votes = usefull_votes/sum(usefull_votes) # rescale 
	names(usefull_votes) = usefull_model
	print(usefull_votes)
	usefull_votes = as.data.frame(t(usefull_votes))
	write.table(usefull_votes, file=output, row.names=F,sep="\t" )
}
