library(ggpubr)
library(FactoMineR)
library(factoextra)
args = commandArgs(trailingOnly=T)
args = sapply(args,function(x){tmp = unlist(strsplit(x,split='='))
            y = vector()
            y[tmp[1]] = tmp[2]
            return(y)},USE.NAMES=F)
dir=args['dir']
output_dir=args['output_dir']
ddir=file.path(dir,'modelComp')
prior_data=do.call(rbind,lapply(list.files(ddir,recursive=T,full.names=T,pattern='ABCstat_global.txt'),read.table,h=T))
prior_data$tag='prior'

obs_data=read.table(file.path(dir,'ABCstat_global.txt'),h=T)
obs_data$tag='obs'

post_data=read.table(file.path(dir,'sim_posterior/ABCstat_global.txt'),h=T)
post_data$tag='post'

all_data=rbind(prior_data,obs_data,post_data)
all_data_pca=subset(all_data,select=-c(dataset,tag))
pca=PCA(all_data_pca,graph=F)
ind_p=fviz_pca_ind(pca,habillage=as.factor(all_data$tag),label='n',addEllipse=T)
var_p=fviz_pca_var(pca)
ggsave(ggarrange(ind_p,var_p),filename=file.path(output_dir,'QC_posterior_acp.pdf'),width=14)


param_2view=paste0(c('piA','piB','thetaA','thetaB','FST','divAB','netDivAB','sf','ss'),'_avg')
print(param_2view)

plot_list=lapply(param_2view,function(param,...){
					 print(param)
					 p=ggdensity(all_data,param,color='tag',fill='tag',title=param)+geom_vline(xintercept=obs_data[1,param],color='red')
					 return(p)
			})
ggexport(ggarrange(plotlist=plot_list,ncol=2,nrow=3),filename=file.path(output_dir,'QC_posterior_density.pdf'))


