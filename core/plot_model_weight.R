
library(ggpubr)
args = commandArgs(trailingOnly=T)
args = sapply(args,function(x){tmp = unlist(strsplit(x,split='='))
            y = vector()
            y[tmp[1]] = tmp[2]
            return(y)},USE.NAMES=F)
model_weight=args['model_weight']
prior_dir=args['prior_dir']
output=args['output']


prior_ld=list.dirs(prior_dir,full.names=T,recursive=F)
model_list=sub('/','',sub('N_.*','N',sub(prior_dir,'',prior_ld))) 
neutral_weight=table(model_list)/length(model_list)
weight=read.table(model_weight,h=T)
diff_weight=weight-neutral_weight
dd=data.frame('weight'=as.vector(t(diff_weight)),'model'=names(diff_weight))
plot_weigth_change=ggbarplot(dd,'model','weight',fill='gray',ylab='w - neutral_w') + theme_bw() + labs_pubr()
ggsave(plot_weigth_change,filename=output,width=14)
