require(stringr)

args=commandArgs(T)
# read file
z=read.csv(args[1],header=F)

# select ones with RMSD < 2.0 to the canonical P-loop fragment
x=z[which(z$V10<=2.0),]

# concatenate name and sheetid "NAME SHEETID" is taken as a single string
y=paste(x$V1,x$V15)
r=rle(sort(y))

# find the entry which appears twice (meaning having two P-loop on the same beta-sheet)
have2=r$value[which(r$lengths==2)]

# split the concatenated lines
spr=str_split(have2," ")

# make vector
uspr=unlist(spr)

# find the length of vector. This makes vector of "Name1 sheetid1 name2 sheetid2...." 
len=length(uspr)

# extract the name of entires (odd numbered elements are names)
name=uspr[seq(1,len,2)]

# get sheetid. (even numbered elements are sheetid)
sheetid=uspr[seq(1+1,len,2)]

# extract entries of interest from original data and save it
ind=is.element(x$V1,name)
xx=x[ind,]
write.csv(file="DualPloop_RMSDlt2p0.csv",xx,quote = F,row.names = F)

# remove "sheetid = -1" which means "No-sheet" and save it
mo=name[which(sheetid!=-1)]
ind2=is.element(x$V1,mo)
xxx=x[ind2,]
write.csv(file="DualPloop_RMSDlt2p0_samesheet.csv",xxx,quote = F,row.names = F)
