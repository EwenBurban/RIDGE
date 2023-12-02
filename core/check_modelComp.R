args = commandArgs(trailingOnly=T)
args = sapply(args,function(x){tmp = unlist(strsplit(x,split='='))
            y = vector()
            y[tmp[1]] = tmp[2]
            return(y)},USE.NAMES=F)

ref_table_dir=args['dir']

ref_table_ld=list.dirs(ref_table_dir,full.names=T,recursive=F)
exist=sapply(ref_table_ld,function(x){file.exists(file.path(x,'ABCstat_global.txt'))})
ref_table_ld=ref_table_ld[exist]

ref_table_lf=file.path(ref_table_ld,'ABCstat_global.txt')
ref_table_loaded=lapply(ref_table_lf,read.table,h=T)

#### check the presence of header

status=sapply(ref_table_loaded,function(x){if(colnames(x)[1]=='dataset'){return(T)}else{return(F)}})

### correct corrupted files

header=colnames(ref_table_loaded[status][[1]])
print(header)
print(ref_table_lf[status==F])
for (f in ref_table_lf[status==F]){
	print(f)
	d=read.table(f)
	colnames(d)=header
	write.table(d,file=f,row.names=F,sep='\t',col.names=T,quote=F)
}
