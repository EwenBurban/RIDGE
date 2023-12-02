library(ggpubr)
args = commandArgs(trailingOnly=T)
args = sapply(args,function(x){tmp = unlist(strsplit(x,split='='))
            y = vector()
            y[tmp[1]] = tmp[2]
            return(y)},USE.NAMES=F)
posterior_data=args['posterior']
prior_dir=args['prior_dir']
Nref=as.numeric(args['Nref'])
output=args['output']
posterior_data='~/remote_serv/ifb_project/crow/poelstra_comp/posterior.txt'
prior_dir='/home/ewenb/remote_serv/ifb_project/crow/poelstra_comp/modelComp/'
Nref=1e5
data=read.table(posterior_data,h=T)
##### define custom function

plot_distr=function(posterior_dist,prior_dist,param=''){
	h=hist(posterior_dist,plot=F,breaks=100)
	p=gghistogram(posterior_dist,bins=100,fill='gray',title=param,xlab=param) +
		geom_vline(xintercept=mean(posterior_dist),color='red') + 
		geom_vline(xintercept=mean(prior_dist),color='black',linetype='dashed') +
		annotate('text',x=1.2* mean(posterior_dist), y=1.2* max(h$counts),label=paste0('mean = ',round(mean(posterior_dist))))
	p = p+ theme_bw() + labs_pubr()
	return(p)}


all_param = c('N1','N2','Na','M_current','M_ancestral','shape_N_a','shape_N_b','Tsc','Tam','Tsplit','PbarrierM_current','PbarrierM_ancestral')
zero=list()

zero[['SI']] = c('M_current'=0,'M_ancestral'=0,'shape_N_a'=1e4,'shape_N_b'=1e4,
				 'Tsc'=0,
				 'Tam'='x["Tsplit"]','PbarrierM_current'=0,'PbarrierM_ancestral'=0)
zero[['SC']] = c('M_ancestral'=0,'shape_N_a'=1e4,'shape_N_b'=1e4,
				 
				 'Tam'='x["Tsplit"]','PbarrierM_current'=0,'PbarrierM_ancestral'=0)
zero[['AM']] = c('M_current'=0,'shape_N_a'=1e4,'shape_N_b'=1e4,
				 'Tsc'=0,
				 'PbarrierM_current'=0,'PbarrierM_ancestral'=0)
zero[['IM']] = c('M_ancestral'='x["M_current"]','shape_N_a'=1e4,'shape_N_b'=1e4,
				 'Tsc'='as.numeric(x["Tsplit"])*r',
				 'Tam'='as.numeric(x["Tsplit"])*r','PbarrierM_current'=0,
				 'PbarrierM_ancestral'='if(is.na(x["PbarrierM_current"])){0}else{x["PbarrierM_current"]}')
get_zeros <- function (x,z,...) {# a function to generate the zero values for each row of a model posterior
	r=runif(1)
	m=sub('/','',sub('_.*','',x['model']))
	all_param = c('N1','N2','Na','M_current','M_ancestral','shape_N_a','shape_N_b','Tsc','Tam','Tsplit','PbarrierM_current','PbarrierM_ancestral')
	missing_param = setdiff(all_param,names(x))
	z = zero[[m]][missing_param]
	zz = sapply(z,function(y,...) eval(parse(text=y)))
	zz=as.numeric(zz)
	names(zz)=missing_param
	return(zz)
}

####### get prior table
prior_ld=list.dirs(prior_dir,full.names=T,recursive=F)
exist=sapply(prior_ld,function(x){file.exists(file.path(x,'ABCstat_global.txt'))})
prior_ld=prior_ld[exist]

prior_data= lapply(prior_ld,function(x){
						   prior_data=read.table(file.path(x,'priorfile.txt'),h=T)
						   prior_data$model=sub('/','',sub('N_.*','N',sub(prior_dir,'',x)))
						   zero_table=t(apply(prior_data,1,function(x,...) get_zeros(x,zero)))
						   prior_data=as.data.frame(cbind(prior_data,zero_table))
						   prior_data=subset(prior_data,select=-model)
						   return(prior_data)
				 })
prior_data=do.call(rbind,prior_data)
##### generate Tsplit & Tsc & Tam plot

Tsplit_plot= plot_distr(data$Tsplit * 4 * Nref,prior_data$Tsplit* 4 * Nref,param='Tsplit')
Tsc_plot= plot_distr(data$Tsc* 4 * Nref,prior_data$Tsc* 4 * Nref,param='Tsc')
Tam_plot= plot_distr(data$Tam* 4 * Nref,prior_data$Tam* 4 * Nref,param='Tam')

### generate N1 & N2 & Na plot
## Because Ne distribution is affected by a beta shape, it recreate first a distribution from mean value and beta parameter
n=100
N1_post_dist=as.vector(apply(data,1,function(vv,...){rbeta(n,vv['shape_N_a'],vv['shape_N_b']) * vv['N1'] /((vv['shape_N_a'])/ (vv['shape_N_a'] + vv['shape_N_b']))}))* Nref
N2_post_dist=as.vector(apply(data,1,function(vv,...){rbeta(n,vv['shape_N_a'],vv['shape_N_b']) * vv['N2'] /((vv['shape_N_a'])/ (vv['shape_N_a'] + vv['shape_N_b']))}))* Nref
Na_post_dist=as.vector(apply(data,1,function(vv,...){rbeta(n,vv['shape_N_a'],vv['shape_N_b']) * vv['Na'] /((vv['shape_N_a'])/ (vv['shape_N_a'] + vv['shape_N_b']))}))* Nref

N1_prior_dist=as.vector(apply(prior_data,1,function(vv,...){rbeta(n,vv['shape_N_a'],vv['shape_N_b']) * vv['N1'] /((vv['shape_N_a'])/ (vv['shape_N_a'] + vv['shape_N_b']))}))* Nref
N2_prior_dist=as.vector(apply(prior_data,1,function(vv,...){rbeta(n,vv['shape_N_a'],vv['shape_N_b']) * vv['N2'] /((vv['shape_N_a'])/ (vv['shape_N_a'] + vv['shape_N_b']))}))* Nref
Na_prior_dist=as.vector(apply(prior_data,1,function(vv,...){rbeta(n,vv['shape_N_a'],vv['shape_N_b']) * vv['Na'] /((vv['shape_N_a'])/ (vv['shape_N_a'] + vv['shape_N_b']))}))* Nref

N1_plot=plot_distr(N1_post_dist,N1_prior_dist,param='N1')
N2_plot=plot_distr(N2_post_dist,N2_prior_dist,param='N2')
Na_plot=plot_distr(Na_post_dist,Na_prior_dist,param='Na')

### generate M_curr and M_anc plot
### Because m distribution is affected by a Pbarrier, it reacreate first a distribution from m value and Pbarrier 
n=100
M_cur_post_dist=as.vector(apply(data,1,function(vv,...){rbinom(n,size=1,prob=vv['PbarrierM_current']) * vv['M_current'] * 4* Nref}))
M_anc_post_dist=as.vector(apply(data,1,function(vv,...){rbinom(n,size=1,prob=vv['PbarrierM_ancestral']) * vv['M_ancestral'] * 4* Nref}))

M_cur_prior_dist=as.vector(apply(prior_data,1,function(vv,...){rbinom(n,size=1,prob=vv['PbarrierM_current']) * vv['M_current'] * 4* Nref}))
M_anc_prior_dist=as.vector(apply(prior_data,1,function(vv,...){rbinom(n,size=1,prob=vv['PbarrierM_ancestral']) * vv['M_ancestral'] * 4* Nref}))

M_cur_plot=plot_distr(M_cur_post_dist,M_cur_prior_dist,param='M_cur')
M_anc_plot=plot_distr(M_anc_post_dist,M_anc_prior_dist,param='M_anc')

all_plot=ggarrange(N1_plot,N2_plot,Na_plot,Tsplit_plot,Tsc_plot,Tam_plot,M_cur_plot,M_anc_plot,ncol=3,nrow=3)
ggsave(all_plot,filename=output,width=14,height=14)
