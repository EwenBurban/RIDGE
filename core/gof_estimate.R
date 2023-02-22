## load argument
args = commandArgs(trailingOnly=T)
args = sapply(args,function(x){tmp = unlist(strsplit(x,split='='))
            y = vector()
            y[tmp[1]] = tmp[2]
            return(y)},USE.NAMES=F)

obs_dir=args['obs_dir']
ref_table_dir=args['sim_dir']
mode=args['mode']
output=args['output']
nb_replicate=as.numeric(args['nb_replicate'])
type=args['type']

## declare function ###################
normalise <- function(x,y){
	if(mad(y) == 0)
	return (x)
	else
	return (x/mad(y))
}


euclidian_dist=function(target,sim){
	d=rowSums((target[rep(1,nrow(sim)),] - sim)^2)
	d=sqrt(d)
	m_d=mean(d)
	return(m_d)
}

gfit=function(target,sumstat,nb.replicate=100){
	### rescale data
	scaled.sumstat=as.data.frame(apply(sumstat,2,function(x) normalise(x,x)))
	scaled.target=as.data.frame(sapply(colnames(target),function(x,...) normalise(target[x],sumstat[,x])))
	## calculate null distance distribution
	print('Calculate null distance distribution')
	pb=txtProgressBar(min=0,max=nb.replicate,style=3)
	eval_row=sample(1:nrow(scaled.sumstat),size=nb.replicate)
	sim.dist=sapply(eval_row,function(x,...) {
						y=euclidian_dist(target=scaled.sumstat[x,],sim=scaled.sumstat[-x,])
						setTxtProgressBar(pb,which(eval_row==x))
						return(y)
			})
	close(pb)
	## calculate the distance of observed data
	if(nrow(target)==1){
	obs.dist=euclidian_dist(target=scaled.target,sim=scaled.sumstat)
	pval=length(which(sim.dist>obs.dist))/length(sim.dist)

	}else if(nrow(target)>1){
	obs.dist=apply(scaled.target,1,function(x,...) euclidian_dist(target=matrix(x,nr=1,nc=length(x)),sim=scaled.sumstat))
	pval=sapply(obs.dist,function(x,...) length(which(sim.dist>x))/length(sim.dist))
	}
	## calculate the pvalue


	res=list(pval=pval,obs.dist=obs.dist,sim.dist=sim.dist)
	return(res)
}

## load ref table and obs dataset

if (type=='prior'){
ref_table_ld=file.path(list.dirs(ref_table_dir,full.names=T,recursive=F),'ABCstat_global.txt')
ref_table_data=do.call(rbind,lapply(ref_table_ld,read.table,h=T))}
if (type=='post'){
	ref_table_data=read.table(file.path(ref_table_dir,'ABCstat_global.txt'),h=T)
}

if (mode=='single'){
	list_obs_dir=obs_dir
	print(list_obs_dir)
	obs_data=read.table(file.path(obs_dir,'ABCstat_global.txt'),h=T)}
if (mode=='multi'){
	list_obs_dir=obs_dir
	list_obs_dir = list.dirs(path=obs_dir,full.names = T,recursive=F)
	print(list_obs_dir)
	obs_data = do.call(rbind,lapply(file.path(list_obs_dir,'ABCstat_global.txt'),read.table,h=T))
}
## calculate gof

gof=gfit(obs_data,ref_table_data,nb.replicate=nb_replicate)

## write gof
i=1
for (d in list_obs_dir){
	write.table(t(c('pval'=gof$pval[i],'dist'=gof$obs.dist[i])),file=file.path(d,output),row.names=F,sep='\t',quote=F)
i=i+1
}
