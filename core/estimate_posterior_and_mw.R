args = commandArgs(trailingOnly=T)
args = sapply(args,function(x){tmp = unlist(strsplit(x,split='='))
            y = vector()
            y[tmp[1]] = tmp[2]
            return(y)},USE.NAMES=F)

library(abcrf)
ntree=as.numeric(args['ntree'])
ncores=as.numeric(args['ncores'])
obs_dir=args['obs_dir']
ref_table_dir=args['sim_dir']
mode=args['mode']
nPosterior=as.numeric(args['nPosterior'])
Pbarrier_max=as.numeric(args['Pbarrier_max'])
if(is.na(args['tag'])==F){tag=args['tag']} else {tag=''}
print(nPosterior)
## zero params ## 
#all_param = c('N1','N2','Na','M_current','M_ancestral','shape_N_a','shape_N_b','shape_M_current_a','shape_M_current_b','shape_M_ancestral_a','shape_M_ancestral_b','Tsc','Tam','Tsplit','PbarrierM_current','PbarrierM_ancestral')
all_param = c('N1','N2','Na','M_current','M_ancestral','shape_N_a','shape_N_b','Tsc','Tam','Tsplit','PbarrierM_current','PbarrierM_ancestral')
all_param = c('M_current','Tsplit','PbarrierM_current')
zero=list()

zero[['SI']] = c('M_current'=0,'M_ancestral'=0,'shape_N_a'=1e4,'shape_N_b'=1e4,'shape_M_current_a'=1e4,
				 'shape_M_current_b'=1e4,'shape_M_ancestral_a'=1e4,'shape_M_ancestral_b'=1e4,'Tsc'=0,
				 'Tam'='x["Tsplit"]','PbarrierM_current'='runif(1)','PbarrierM_ancestral'='runif(1)')
zero[['SC']] = c('M_ancestral'=0,'shape_N_a'=1e4,'shape_N_b'=1e4,'shape_M_current_a'=1e4,
				 'shape_M_current_b'=1e4,'shape_M_ancestral_a'=1e4,'shape_M_ancestral_b'=1e4,
				 'Tam'='x["Tsplit"]','PbarrierM_current'=0,'PbarrierM_ancestral'='runif(1)')
zero[['AM']] = c('M_current'=0,'shape_N_a'=1e4,'shape_N_b'=1e4,'shape_M_current_a'=1e4,
				 'shape_M_current_b'=1e4,'shape_M_ancestral_a'=1e4,'shape_M_ancestral_b'=1e4,'Tsc'=0,
				 'PbarrierM_current'='runif(1)','PbarrierM_ancestral'=0)
zero[['IM']] = c('M_ancestral'=0,'shape_N_a'=1e4,'shape_N_b'=1e4,'shape_M_current_a'=1e4,
				 'shape_M_current_b'=1e4,'shape_M_ancestral_a'=1e4,'shape_M_ancestral_b'=1e4,'Tsc'='as.numeric(x["Tsplit"])*r',
				 'Tam'='as.numeric(x["Tsplit"])*r','PbarrierM_current'=0,'PbarrierM_ancestral'=0)
zero[['SI']] = c('M_current'=0,'M_ancestral'=0,'shape_N_a'=1e4,'shape_N_b'=1e4,
				 'Tsc'=0,
				 'Tam'='x["Tsplit"]','PbarrierM_current'='runif(1,0,Pbarrier_max)','PbarrierM_ancestral'='runif(1,0,Pbarrier_max)')
zero[['SC']] = c('M_ancestral'=0,'shape_N_a'=1e4,'shape_N_b'=1e4,
				 
				 'Tam'='x["Tsplit"]','PbarrierM_current'=0,'PbarrierM_ancestral'='runif(1,0,Pbarrier_max)')
zero[['AM']] = c('M_current'=0,'shape_N_a'=1e4,'shape_N_b'=1e4,
				 'Tsc'=0,
				 'PbarrierM_current'='runif(1,0,Pbarrier_max)','PbarrierM_ancestral'=0)
zero[['IM']] = c('M_ancestral'=0,'shape_N_a'=1e4,'shape_N_b'=1e4,
				 'Tsc'='as.numeric(x["Tsplit"])*r',
				 'Tam'='as.numeric(x["Tsplit"])*r','PbarrierM_current'=0,'PbarrierM_ancestral'=0)
### modif 2 : when no data on Pbarrier => Pbarrier=0
zero[['SI']] = c('M_current'=0,'M_ancestral'=0,'shape_N_a'=1e4,'shape_N_b'=1e4,
				 'Tsc'=0,
				 'Tam'='x["Tsplit"]','PbarrierM_current'=0,'PbarrierM_ancestral'=0)
zero[['SC']] = c('M_ancestral'=0,'shape_N_a'=1e4,'shape_N_b'=1e4,
				 
				 'Tam'='x["Tsplit"]','PbarrierM_current'=0,'PbarrierM_ancestral'=0)
zero[['AM']] = c('M_current'=0,'shape_N_a'=1e4,'shape_N_b'=1e4,
				 'Tsc'=0,
				 'PbarrierM_current'=0,'PbarrierM_ancestral'=0)
zero[['IM']] = c('M_ancestral'=0,'shape_N_a'=1e4,'shape_N_b'=1e4,
				 'Tsc'='as.numeric(x["Tsplit"])*r',
				 'Tam'='as.numeric(x["Tsplit"])*r','PbarrierM_current'=0,'PbarrierM_ancestral'=0)
get_zeros <- function (x,z,...) {# a function to generate the zero values for each row of a model posterior
	r=runif(1)
	m=sub('/','',sub('_.*','',x['model']))
	all_param = c('N1','N2','Na','M_current','M_ancestral','shape_N_a','shape_N_b','shape_M_current_a','shape_M_current_b','shape_M_ancestral_a','shape_M_ancestral_b','Tsc','Tam','Tsplit','PbarrierM_current','PbarrierM_ancestral')
	all_param = c('N1','N2','Na','M_current','M_ancestral','shape_N_a','shape_N_b','Tsc','Tam','Tsplit','PbarrierM_current','PbarrierM_ancestral')
	missing_param = setdiff(all_param,names(x))
	z = zero[[m]][missing_param]
	zz = sapply(z,function(y,...) eval(parse(text=y)))
	zz=as.numeric(zz)
	names(zz)=missing_param
	return(zz)
}
## load of all dataset comings from all models
## each dataset is loaded and at same time, model of origin is extracted from file name
## and zero params are added
ref_table_ld=list.dirs(ref_table_dir,full.names=T,recursive=F)
ref_table_data= lapply(ref_table_ld,function(x,...){
						   ss_data=read.table(file.path(x,'ABCstat_global.txt'),h=T)
						   prior_data=read.table(file.path(x,'priorfile.txt'),h=T)
						   prior_data$model=sub('/','',sub('N_.*','N',sub(ref_table_dir,'',x)))
						   zero_table=t(apply(prior_data,1,function(x,...) get_zeros(x,zero)))
						   prior_data=as.data.frame(cbind(prior_data,zero_table))
						   data=merge(ss_data,prior_data,by='dataset')
						   return(data)
				 })
ref_table_data=do.call(rbind,ref_table_data)
print(unique(ref_table_data$model))
ref_table_ss=ref_table_data[!names(ref_table_data) %in% c(all_param,'model','dataset')]
ref_table_prior=subset(ref_table_data,select=all_param)
print('ref data loaded')
## load observed dataset
if (mode=='single'){
	list_obs_dir=obs_dir
	print(list_obs_dir)
	obs_data=read.table(file.path(obs_dir,'ABCstat_global.txt'),h=T)}
if (mode=='multi'){
	list_obs_dir=obs_dir
	list_obs_dir = list.dirs(path=obs_dir,full.names = T,recursive=F)
#	list_obs_dir = sample(list_obs_dir,size=720,replace=T)
	print('Warning, subsampling has been added line 83 for testrun, remove it for real runs + remove the sf/df_outlier addition line 93')
	print(list_obs_dir)
	obs_data = do.call(rbind,lapply(file.path(list_obs_dir,'ABCstat_global.txt'),read.table,h=T))
}
print('obs data loaded')
## data purification (remove uneccessary summary stats)
cv_vec=apply(ref_table_ss,2,sd,na.rm=T)/colMeans(ref_table_ss,na.rm=T)
if(any(cv_vec<0.01)){
	print(colnames(ref_table_ss)[which(cv_vec<0.01)])
	sel_colnames=c(colnames(ref_table_ss)[-c(which(cv_vec<0.01),which(colnames(ref_table_ss)=='sf_outlier'))])
	print(setdiff(sel_colnames,colnames(obs_data)))
	print(setdiff(colnames(obs_data),sel_colnames))
	ref_table_ss=subset(ref_table_ss,select=sel_colnames)
	obs_data=subset(obs_data,select=sel_colnames)
}
sel_colnames=colnames(ref_table_ss)[-which(colnames(ref_table_ss)=='sf_outlier')]
	print(setdiff(sel_colnames,colnames(obs_data)))
	print(setdiff(colnames(obs_data),sel_colnames))
obs_data=subset(obs_data,select=sel_colnames)
ref_table_ss=subset(ref_table_ss,select=sel_colnames)
print('data filtered')
## generate the regression random forest for each parameter
## predict expected value + generate posterior for the parameter and store rf weigth in a vector
    res=lapply(all_param,function(p,...){
    	target=ref_table_prior[,p]
    	rf_data=data.frame(ref_table_ss,target)
    	rf=regAbcrf(target~.,data=rf_data,ntree=ntree,paral=T,ncores=ncores)
    	pred=predict(rf,training=rf_data,obs=obs_data,paral=T,ncores=ncores,rf.weights=T)
    	## It return a table with O columns (with O the number of observed dataset) and nPosterior rows
    	print(dim(pred$weights))
    	print(nPosterior)
    	posterior=apply(pred$weights,2,function(x,...) sample(target,size=nPosterior-1,replace=T,prob=x))
    	posterior=rbind(posterior,pred$expectation)
    	## It return a table with O columns (with O the number of observed dataset) and nrows(rf_data) rows
    	weights=pred$weights
    	r=list('posterior'=posterior,'weights'=weights)
    	print(paste(p,'done'))
    	return(r)
    })
    names(res)=all_param
## generate posteriors
list_posterior_table=lapply(1:nrow(obs_data),function(O,...){
		   y=do.call(cbind,lapply(res,function(x,...) x[['posterior']][,O,drop=F]))
		   colnames(y)=names(res)
		   return(y)
})

i=1
for ( dir in list_obs_dir){
	output_posterior_file=file.path(dir,paste0('posterior',tag,'.txt'))
	write.table(list_posterior_table[[i]],file=output_posterior_file,sep='\t',quote=F,row.names=F)
	output_point_posterior_file=file.path(dir,paste0('point_posterior',tag,'.txt'))
	write.table(list_posterior_table[[i]][nPosterior,,drop=F],file=output_point_posterior_file,sep='\t',quote=F,row.names=F)
	i=i+1
}
## generate model weight using the sum of all parameter weight 
sum_weights=Reduce('+',lapply(res,function(x) x[['weights']]))
model_weights=lapply(1:nrow(obs_data),function(x){
						 y=tapply(sum_weights[,x],as.factor(ref_table_data$model),sum)/sum(sum_weights[,x])
						 names(y)=levels(as.factor(ref_table_data$model))
						 return(y)
})

i=1
for ( dir in list_obs_dir){
	output_mw_file=file.path(dir,'model_weight.txt')
	write.table(t(model_weights[[i]]),file=output_mw_file,sep='\t',quote=F,row.names=F)
	i=i+1
}
