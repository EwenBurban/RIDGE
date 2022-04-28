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
all_param = c('N1','N2','Na','M_current','M_ancestral','shape_N_a','shape_N_b','shape_M_current_a','shape_M_current_b','shape_M_ancestral_a','shape_M_ancestral_b','Tsc','Tam','Tsplit','PbarrierM_current','PbarrierM_ancestral')

## zero params ## 
zero=list()

zero[['SI']] = c('M_current'=0,'M_ancestral'=0,'shape_N_a'=90,'shape_N_b'=1,'shape_M_current_a'=90,
				 'shape_M_current_b'=1,'shape_M_ancestral_a'=90,'shape_M_ancestral_b'=1,'Tsc'=0,
				 'Tam'='x["Tsplit"]','PbarrierM_current'='runif(1)','PbarrierM_ancestral'='runif(1)')
zero[['SC']] = c('M_ancestral'=0,'shape_N_a'=90,'shape_N_b'=1,'shape_M_current_a'=90,
				 'shape_M_current_b'=1,'shape_M_ancestral_a'=90,'shape_M_ancestral_b'=1,
				 'Tam'='x["Tsplit"]','PbarrierM_current'='runif(1)','PbarrierM_ancestral'='runif(1)')
zero[['AM']] = c('M_current'=0,'shape_N_a'=90,'shape_N_b'=1,'shape_M_current_a'=90,
				 'shape_M_current_b'=1,'shape_M_ancestral_a'=90,'shape_M_ancestral_b'=1,'Tsc'=0,
				 'PbarrierM_current'='runif(1)','PbarrierM_ancestral'='runif(1)')
zero[['IM']] = c('M_ancestral'=0,'shape_N_a'=90,'shape_N_b'=1,'shape_M_current_a'=90,
				 'shape_M_current_b'=1,'shape_M_ancestral_a'=90,'shape_M_ancestral_b'=1,'Tsc'='x["Tsplit"]',
				 'Tam'=0,'PbarrierM_current'='runif(1)','PbarrierM_ancestral'='runif(1)')
print(zero)
get_zeros <- function (x,z,m,missing_param) {# a function to generate the zero values for each row of a model posterior
	z = zero[[m]][missing_param]
	zz = sapply(z,function(y) eval(parse(text=y)))
	return(as.vector(zz))
}
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
	zero_data = t(apply(m_post,1,function(x) get_zeros(x,zero,dem_model,missing_param)))
	colnames(zero_data) = missing_param
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
	tmp_list[[length(tmp_list)+1]] = data.frame(list_post_data[[model]][sel_lines,,drop=F])

}

tmp_data = do.call(rbind,tmp_list)

write.table(tmp_data,file=output,row.names=F,quote=F,sep='\t')
