library(ggpubr)
args = commandArgs(trailingOnly=T)
args = sapply(args,function(x){tmp = unlist(strsplit(x,split='='))
            y = vector()
            y[tmp[1]] = tmp[2]
            return(y)},USE.NAMES=F)
data_locus=args['locus_data']
data_global=args['global_data']
output=args['output']
data_locus=read.table(data_locus,h=T)
data_global=read.table(data_global,h=T)
ss_2compare=c('pi','theta','Dtaj','sx')
ss_2view=c('FST','divAB','netDivAB','sf','ss')
### comparition
compare_plot=lapply(ss_2compare,function(param,...){
						param_vec=paste0(param,LETTERS[1:2],'_avg')
						p=ggboxplot(data_locus,NULL,c(param_vec),merge=T,xlab=param) + 
							geom_segment(x=0,xend=1,y=data_global[1,param_vec[1]],yend=data_global[1,param_vec[1]],linetype='dashed',color='red')+
							geom_segment(x=1,xend=2,y=data_global[1,param_vec[2]],yend=data_global[1,param_vec[2]],linetype='dashed',color='red')
						return(p)})
density_plot=lapply(ss_2view,function(param,...){
						 param=paste0(param,'_avg')
						 p=ggdensity(data_locus,param,title=param)+geom_vline(xintercept=data_global[1,param],color='red')
						 return(p)})
all_plot=c(compare_plot,density_plot)
ggexport(plotlist=all_plot,ncol=2,nrow=3,filename=output)




