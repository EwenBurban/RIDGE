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


all_data=rbind(prior_data,obs_data)
all_data_pca=subset(all_data,select=-c(dataset,tag))
all_data_pca=all_data_pca[,colnames(all_data_pca)[-grep(colnames(all_data_pca),pattern='std|median')]]
pca=PCA(all_data_pca,graph=F)
ind_p=fviz_pca_ind(pca,habillage=as.factor(all_data$tag),label='n',addEllipse=T)
var_p=fviz_pca_var(pca,repel=T)
ggsave(ggarrange(ind_p,var_p),filename=file.path(output_dir,'QC_prior_acp.pdf'),width=14)


param_2view=paste0(c('piA','piB','thetaA','thetaB','FST','divAB','netDivAB','sf','ss'),'_avg')
param_2view=c(param_2view,'sf_outlier','fst_outlier','divAB_outlier','netDivAB_outlier','piA_outlier','piB_outlier')

plot_list=lapply(param_2view,function(param,...){
					 p=ggdensity(all_data,param,color='tag',fill='tag',title=param)+geom_vline(xintercept=obs_data[1,param],color='red')
					 return(p)
			})
ggexport(ggarrange(plotlist=plot_list,ncol=2,nrow=3),filename=file.path(output_dir,'QC_prior_density.pdf'))


