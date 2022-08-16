#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=T)
args = sapply(args,function(x){tmp = unlist(strsplit(x,split='='))
		y = vector()
		y[tmp[1]] = tmp[2]
		return(y)},USE.NAMES=F)
obs_dataset = read.table(args['obs_data'],h=T)
first_round_dir = args['first_round_dir']
est_model_dir = first_round_dir
output = args['output']
output_list = list()
library(ggpubr)
library(FactoMineR)

sim_pattern = 'ABCstat_global.txt'
sim_files = list.files(first_round_dir,pattern = sim_pattern,recursive = T,full.names = T)

models=sub("N_.*","N",list.dirs(path=first_round_dir,full.names = F,recursive=F))
sim_data = lapply(sim_files,read.table,h=T)
stopifnot(length(sim_data)==length(models))
model_tag = vector()
for(i in 1:length(sim_files)){
	model_tag = c(model_tag,rep(models[i],nrow(sim_data[[i]])))
}

acp_data = do.call(rbind,sim_data)[,-1]
### add obs dataset ###
acp_data = rbind(acp_data,obs_dataset[,-1])
col_vec = c(rep('grey',nrow(acp_data)-1),'red')
### acp
acp = PCA(acp_data,graph=F)
acp_plot = plot(acp,choix='ind',label='none',habillage='ind',col.hab=col_vec,cex=0.8)

ggsave(acp_plot,filename=output)
#### density plots ### 
sim_res_data= do.call(rbind,sim_data)
sim_res_data$model = model_tag
all_var = c('FST_avg','divAB_avg','netDivAB_avg','sf_avg','piA_avg','piB_avg','thetaA_avg','thetaB_avg','DtajA_avg','DtajB_avg')
ll = list()
for (var in all_var){
ll[[var]] = ggdensity(sim_res_data,var,color='model') + geom_vline(xintercept=obs_dataset[,var],color='red')
if(var!=all_var[1]){ll[[var]]=ggpar(ll[[var]],legend='none',font.legend = c(6, "plain", "black"))}else{
	ll[[var]] = ggpar(ll[[var]],font.legend = c(6, "plain", "black"),legend.title=list(color='model',shape=19))
}
}
density_plot = ggarrange(plotlist = ll,ncol=2,nrow=2)
print('density_plot done')
#### posterior plots ###
posterior_pattern = 'posterior'

posterior_files = list.files(est_model_dir,pattern = posterior_pattern,recursive = F,full.names = T)
posterior_files  = grep(posterior_files,pattern='RandomForest',value=T,invert=T)
posterior_namefiles = list.files(est_model_dir,pattern = posterior_pattern,recursive = F,full.names = F)
posterior_namefiles  = grep(posterior_namefiles,pattern='RandomForest',value=T,invert=T)
models=sub("posterior_","",posterior_namefiles)
models = sub('.txt','',models)
print(models)
list_posterior_data = lapply(posterior_files,read.table,h=T)
print(length(list_posterior_data))
all_param = unique(unlist(lapply(list_posterior_data,colnames)))
if(any(all_param=='dataset')){print('dataset found');all_param = all_param[-which(all_param=='dataset')]}
for(i in 1:length(models)){
	list_posterior_data[[i]]$tag = rep(models[i],nrow(list_posterior_data[[i]]))
}
est_density_plotlist = list()
colnames_posterior = lapply(list_posterior_data,colnames)
##### Tsplit
for (param in c('Tsplit','Tam','Tsc')){
tmp_list= list_posterior_data[grep(x=colnames_posterior,pattern=param)]
if(length(tmp_list)==0){next()}	
for(i in 1:length(tmp_list)){tmp_list[[i]] = subset(tmp_list[[i]],select=c(param,'tag'))}
tmp_data = do.call(rbind,tmp_list)
if(param=='Tsplit'){save_Tsplit_data = tmp_data}
est_density_plotlist[[param]] = ggdensity(tmp_data,param,color='tag')
est_density_plotlist[[param]] = ggpar(est_density_plotlist[[param]],font.legend = c(6, "plain", "black")) + rremove('legend.title')
}


##### N plot
for (param in c('Na','N1','N2')){

	tmp_list= list_posterior_data[grep(x=colnames_posterior,pattern=param)]
	if(length(tmp_list)==0){next()}	
	for(i in 1:length(tmp_list)){
		if(any(grepl(colnames(tmp_list[[i]]),pattern='shape'))){			
		tmp_df = tmp_list[[i]]
		tmp_list[[i]] = subset(tmp_df,select=c(param,'tag'))
		} else {
		tmp_list[[i]] = subset(tmp_list[[i]],select=c(param,'tag')) }
	}
	tmp_data = do.call(rbind,tmp_list)
	if(param=='Na'){save_na_data = tmp_data}
	est_density_plotlist[[param]] = ggdensity(tmp_data,param,color='tag')
	est_density_plotlist[[param]] = ggpar(est_density_plotlist[[param]],font.legend = c(6, "plain", "black")) + rremove('legend.title')

}
#### Migration plot 
for (param in c('M_current','M_ancestral')){

	tmp_list= list_posterior_data[grep(x=colnames_posterior,pattern=param)]
	if(length(tmp_list)==0){next()}	

	for(i in 1:length(tmp_list)){
	 	if(any(grepl(colnames(tmp_list[[i]]),pattern='shape_M'))){			
		tmp_df = tmp_list[[i]]
		tmp_df = apply(tmp_df,1,function(x,...){
			size = 1e2
			a = as.numeric(x[paste0('shape_',param,'_a')])
			b = as.numeric(x[paste0('shape_',param,'_b')])
			barrier = as.numeric(x[paste0('Pbarrier',param)])
			M = as.numeric(rep(x[param],size)) * (rbeta(size,a,b) / (a/(a + b)))# * sample(c(0,1),size,replace=T,prob=c(barrier,1-barrier))
			y = as.data.frame(cbind(M,rep(x['tag'],length(M))))
			colnames(y) = c(param,'tag')
			return(y)})
		tmp_list[[i]] = do.call(rbind,tmp_df)
		} else {
		tmp_list[[i]] = subset(tmp_list[[i]],select=c(param,'tag')) }
	}
	tmp_data = do.call(rbind,tmp_list)
	tmp_data = as.data.frame(tmp_data)
	rownames(tmp_data)= NULL
	tmp_data[,param] = as.numeric(tmp_data[,param])
	est_density_plotlist[[param]] = ggdensity(tmp_data,param,color='tag')
	est_density_plotlist[[param]] = ggpar(est_density_plotlist[[param]],xlim=c(0,150),font.legend = c(6, "plain", "black")) + rremove('legend.title')
}
#### 2D plot Tspilt vs Na

tmp_data = as.data.frame(cbind(save_na_data,save_Tsplit_data))
print(head(tmp_data))
tmp_data = tmp_data[,-2]
tmp_data[,c(1,2)] = apply(tmp_data[,c(1,2)],2,as.numeric)
est_density_plotlist[['2Dplot']] = ggline(tmp_data,'Na','Tsplit',color='tag',numeric.x.axis=T)
est_density_plotlist[['2Dplot']] = ggpar(est_density_plotlist[['2Dplot']],font.legend = c(6, "plain", "black")) + rremove('legend.title')
print('2d plot done')

est_density_plot = ggarrange(plotlist=est_density_plotlist,ncol=2,nrow=2)
ggexport(plotlist=list(acp_plot,density_plot,est_density_plot),filename=output)
