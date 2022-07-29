#!/usr/bin/env Rscript
library(abcrf)
for(arg in commandArgs()){
  tmp = strsplit(arg, '='); opt = tmp[[1]][1]; val = tmp[[1]][2]
  if(opt == "ntree"){ntree=as.numeric(val)}
  if(opt == "ncores"){ncores=as.numeric(val)}
  if(opt == "timeStamp"){timeStamp=val}
  if(opt == "sim_dir"){sim_dir=val}
  if(opt == "output_dir"){output_dir=val}
  if(opt == "mig_mod"){migration=val}
  
}
## custom func
remove_param = function(table){
	cv = apply(table,2,sd)/colMeans(table)
	y = which(abs(cv) < 0.1 | is.nan(cv)) # TODO check la pertinence de cet modif
	return(y)
}

get_roc_stat = function(prior,est,trsh,...){
	TP = length(which(prior==0 & est>=trsh))
	FP = length(which(prior>0 & est>=trsh))
	TN = length(which(prior>0 & est<trsh))
	FN = length(which(prior==0 & est<trsh))
	TPR = TP / (TP + FN)
	FPR = FP / (FP + TN)
	y=c('TPR'=TPR,'FPR'=FPR)
	return(y)
}
# Training set 
train_ABCstat_file_list = file.path(timeStamp,sim_dir,"ABCstat_locus.txt")
train_priorfile_list = file.path(timeStamp,sim_dir,"priorfile_locus.txt")

train_targets = read.table(train_priorfile_list,header=T,stringsAsFactors=F)
train_targets = train_targets[order(train_targets$dataset),]
train_data = read.table(train_ABCstat_file_list,header=T,stringsAsFactors=F)
train_data = train_data[order(train_data$dataset),]

## curation
dataset_col = which(colnames(train_data)=='dataset')
train_data = as.data.frame(sapply(train_data,as.numeric))[,-dataset_col]
train_targets = as.data.frame(sapply(train_targets,as.numeric))[,-dataset_col]
if(any(is.na(train_data[,1]))){naline= which(is.na(train_data[,1]));  train_data=train_data[-naline,];print(paste(naline,' line removed'))}
if(any(is.na(train_targets[,1]))){naline2= which(is.na(train_targets[,1]));  train_targets=train_targets[-naline2,];print(paste(naline2,' line removed'))}
if((any(is.na(train_data[,1])) & any(is.na(train_targets[,1])))){stopifnot(naline==naline2)}

train_data[is.na(train_data)] = 0
train_targets[is.na(train_targets)] = 0

param2rm = remove_param(train_data)
if(length(param2rm)>0){param2kp = colnames(train_data)[-param2rm]} else { param2kp=colnames(train_data)}
print(length(param2rm))

### classification ####
print(paste('train_data possède',nrow(train_data),'lignes'))
locus_to_train = sample(1:nrow(train_data),size=nrow(train_data)/20,replace=F)
print(length(locus_to_train))
locus_to_eval = sample(1:nrow(train_data)[-locus_to_train],size=length(locus_to_train)/10,replace=F)
tag = as.factor(sapply(train_targets[locus_to_train,migration],function(x) ifelse(x==0,'barrier','non-barrier')))
print(length(tag))
train_set=train_data[locus_to_train,param2kp]
train_set$tag=tag
test_set=train_data[locus_to_eval,param2kp]
rf = abcrf(tag~.,data=train_set,ncores = ncores,ntree = ntree,paral = T,lda=F)
pred = predict(rf,training =train_set ,obs =test_set ,paral = T,ncores=ncores)
model_param_estimation = pred$vote[,'barrier'] / rowSums(pred$vote)
print(train_targets[locus_to_eval,migration])
print(model_param_estimation)
## roc estimation

seq_trsh = seq(0,1,by=0.05)
stats = do.call(rbind,lapply(seq_trsh,function(x) get_roc_stat(train_targets[locus_to_eval,migration],model_param_estimation,x)))
stats = data.frame(cbind(seq_trsh,stats))
colnames(stats)=c('threshold','TPR','FPR')
print(head(stats))
print(class(stats))
write.table(stats,file=file.path(output_dir,'roc_table.txt'),row.names=F)
pdf(file.path(output_dir,'roc.pdf'))
plot(stats[,'FPR'],stats[,'TPR'],type='p')
plot(stats[,'threshold'],stats[,'TPR'],type='l',col='blue')
lines(stats[,'threshold'],stats[,'FPR'],col='red')
legend(x='toprigth',legend=c('TPR','FPR'),col=c('blue','red'))
dev.off()
