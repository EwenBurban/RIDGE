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
  if(opt == "locus_estimation"){locus_estimation=val}
  
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
train_data = read.table(train_ABCstat_file_list,header=T,stringsAsFactors=F)
train_data = subset(train_data,select=-dataset)
train_targets = subset(train_targets,select=-dataset)
all = merge(train_data,train_targets,by='seed')
train_data = all[,colnames(train_data)]
train_targets = all[,colnames(train_targets)]
## curation
train_data[is.na(train_data)] = 0
train_targets[is.na(train_targets)] = 0

param2rm = remove_param(train_data)
if(length(param2rm)>0){param2kp = colnames(train_data)[-param2rm]} else { param2kp=colnames(train_data)}
print(length(param2rm))

### classification ####
print(paste('train_data possède',nrow(train_data),'lignes'))
locus_to_train = sample(1:nrow(train_data),size=round(nrow(train_data)/5),replace=F)
print(length(locus_to_train))
locus_to_eval = sample(1:nrow(train_data)[-locus_to_train],size=round(length(locus_to_train)/8),replace=F)
tag = as.factor(sapply(train_targets[locus_to_train,migration],function(x) ifelse(x==0,'barrier','non-barrier')))
print(length(tag))
train_set=train_data[locus_to_train,param2kp]
train_set$tag=tag
test_set=train_data[locus_to_eval,param2kp]
rf = abcrf(tag~.,data=train_set,ncores = ncores,ntree = ntree,paral = T,lda=F)
pred = predict(rf,training =train_set ,obs =test_set ,paral = T,ncores=ncores)
pred$post.prob[which(pred$allocation=='non-barrier')]=1-pred$post.prob[which(pred$allocation=='non-barrier')]
model_param_estimation=pred$post.prob
## roc estimation
locus_estimation = read.table(locus_estimation,h=T)
seq_trsh = seq(0,1,by=0.005)
stats = do.call(rbind,lapply(seq_trsh,function(x) get_roc_stat(train_targets[locus_to_eval,migration],model_param_estimation,x)))
stats = data.frame(cbind(seq_trsh,stats))
colnames(stats)=c('threshold','TPR','FPR')
stats$nlocus = sapply(stats$threshold,function(x,...) length(which(locus_estimation[,'Pbarrier'] >= x )))
stats[is.na(stats)]=0
print(head(stats))
print(class(stats))
write.table(stats,file=file.path(output_dir,'roc_table.txt'),row.names=F)
pdf(file.path(output_dir,'roc.pdf'))
plot(stats[,'FPR'],stats[,'TPR'],type='p')
plot(stats[,'threshold'],stats[,'TPR'],type='l',col='blue')
lines(stats[,'threshold'],stats[,'FPR'],col='red')
legend(x='topright',legend=c('TPR','FPR'),fill=c('blue','red'))
plot(stats[,'threshold'],stats[,'TPR']/(stats[,'TPR'] + stats[,'FPR']),type='l')
dev.off()
