#!/usr/bin/env Rscript

for(arg in commandArgs()){
	tmp= unlist(strsplit(arg,split="="))
	opt = tmp[1];val = tmp[2]
	if(opt=='master_post_dir'){m_post_dir=val}
	if(opt=='nPosteriors_locus'){nPost=as.numeric(val)}
	if(opt=='weight_file'){weight_file=val}
	if(opt=='output'){output=val}

}
zero_vec=c()

model_list = c('SC_1M_1N', 'SC_1M_2N', 'SC_2M_1N', 'SC_2M_2N', 'AM_1M_1N', 'AM_1M_2N', 'AM_2M_1N', 'AM_2M_2N', 'IM_1M_1N', 'IM_1M_2N', 'IM_2M_1N', 'IM_2M_2N', 'SI_1N', 'SI_2N')
all_param = c('N1','N2','Na','founders1','founders2','Tdem1','Tdem2','M12','M21','shape_N_a','shape_N_b','shape_M12_a','shape_M12_b','shape_M21_a','shape_M21_b','Tsc','Tam','Tsplit','nBarriersM12','nBarriersM21')

## zero params ## 
zero=list()

zero[['SI']] = c('M12'=0,'M21'=0,'shape_N_a'=9e6,'shape_N_b'=1,'shape_M12_a'=9e6,'shape_M12_b'=1,'shape_M21_a'=9e6,'shape_M21_b'=1,'Tsc'=0,'Tam'=1e6,'nBarriersM12'=1e3,'nBarriersM21'=1e3)
zero[['SC']] = c('M12'=0,'M21'=0,'shape_N_a'=9e6,'shape_N_b'=1,'shape_M12_a'=9e6,'shape_M12_b'=1,'shape_M21_a'=9e6,'shape_M21_b'=1,'Tsc'=1e6,'Tam'=1e6,'nBarriersM12'=0,'nBarriersM21'=0)
zero[['AM']] = c('M12'=0,'M21'=0,'shape_N_a'=9e6,'shape_N_b'=1,'shape_M12_a'=9e6,'shape_M12_b'=1,'shape_M21_a'=9e6,'shape_M21_b'=1,'Tsc'=0,'Tam'=1e6,'nBarriersM12'=0,'nBarriersM21'=0)
zero[['IM']] = c('M12'=0,'M21'=0,'shape_N_a'=9e6,'shape_N_b'=1,'shape_M12_a'=9e6,'shape_M12_b'=1,'shape_M21_a'=9e6,'shape_M21_b'=1,'Tsc'=1e6,'Tam'=1e6,'nBarriersM12'=0,'nBarriersM21'=0)
zero = do.call(rbind,zero) 
weight_data = read.table(weight_file,h=T)
selected_model = colnames(weight_data)
print('weight data loaded')
list_post_files = file.path(m_post_dir,paste0('posterior_',selected_model,'.txt'))
stopifnot(length(list_post_files)>0)
list_post_data=list()
for (m in selected_model){
	dem_model = substring(m,1,2)
	m_post_file = grep(list_post_files,pattern=paste0(m,'.txt'),value=T)
	print(m_post_file)
      	stopifnot(length(m_post_file)==1)
	m_post = read.table(m_post_file,h=T)
	missing_param = setdiff(all_param,colnames(m_post))
	zero_data = zero[rep(dem_model,nrow(m_post)),missing_param,drop=F]
	m_post = cbind(m_post,zero_data)
	print(ncol(m_post))
	list_post_data[[m]] = m_post	
}

#list_post_data = lapply(list_post_files,read.table,h=T)

sample_size_model = table(sample(colnames(weight_data),nPost,prob=weight_data[1,],replace=T))
tmp_list = list()
for(model in names(sample_size_model)){
	sample_size = sample_size_model[model]
	sel_lines = sample(1:nPost,size=as.numeric(sample_size))
	print(model)
	tmp_list[[length(tmp_list)+1]] = data.frame(list_post_data[[model]][sel_lines,,drop=F])
	print( colnames(data.frame(list_post_data[[model]][sel_lines,,drop=F])))

}

tmp_data = do.call(rbind,tmp_list)

write.table(tmp_data,file=output,row.names=F,quote=F,sep='\t')
