# Mapping of Function and Comments
df <- read.csv(file = "faprotax_report.csv", header = FALSE)
res <- list()
col_name <- c()
ncol_name <- 0
res_col <- 0

for (i in 1:nrow(df)) {
  if (grepl("#", df[i, 1])) {
    res_row <- 1
    res_col <- res_col + 1
    res[[res_col]] <- list(df[i, 1])
    ncol_name <- ncol_name + 1
    col_name[ncol_name] = df[i, 1]
  } else {
    res[[res_col]][[res_row]] <- df[i, 1]
    res_row <- res_row + 1
  }
}

# Standardize the length of the list
library(purrr)
res <- map(res, ~ {
  len <- max(lengths(res))
  if (length(.x) < len) {
    c(.x, rep("nothing", len - length(.x)))
  } else {
    .x
  }
})

df_output <- as.data.frame(Reduce(cbind, res))
colnames(df_output) <- col_name

# Convert all columns to character type
df_output[] <- lapply(df_output, as.character)

# Write the results to a file
write.table(df_output, "FuncAnnoMap.csv", sep = ",", row.names = FALSE, quote = FALSE)

# Remove duplicates in each column and fill with "nothing" to match the length of the longest column
df_output_clean <- as.data.frame(lapply(df_output, function(x) {c(unique(x), rep("nothing", max(lengths(df_output)) - length(unique(x))))}))

# View the processed dataframe
head(df_output_clean)

# Save the result
write.csv(df_output_clean, file = "df_output_clean.csv", row.names = FALSE)

# Map Functions to otutab
uniqList = read.csv("df_output_clean.csv",header = T)
library(tidyr)
uniqMatrix = gather(uniqList,func,anno)
res = uniqMatrix[which(uniqMatrix$anno!=""),]
otutab = read.csv("16S_otutab_anno.csv")
otuFuncRes = merge(otutab,res,by.x = "taxonomy","anno")
write.csv(otuFuncRes,"otuFuncRes.csv",row.names = F)

# An OTU can correspond to multiple functions, handle duplicate functions
# Extract all OTU-function pairs
df = subset(otuFuncRes,select = c(OUT_ID,func))
colnames(df) <- c("OTU_ID","func")

# Get a list of all unique functions
func_list <- unique(df$func)

# Calculate the maximum length for the OTU ID column list
max_len <- max(table(df$func))

# Define an empty list
res_list <- list()
funcCount = 1
for (func in func_list) {
  # Get all OTU IDs for the current function
  otuset <- df[df$func == func, "OTU_ID"]
  # Add all OTU IDs for the current function to res_list, using func as the column name
  res_list[[func]] <- c(otuset, rep("nothing", max_len-length(otuset)))
  res_list[[paste("func",funcCount,sep = "")]] <- c(rep(func,length(otuset)),rep("nothing", max_len-length(otuset)))
  funcCount = funcCount + 1
}

# Convert the result list to a dataframe, while converting func to the corresponding column names
res_df <- data.frame(res_list)
# Save the result dataframe to a file
write.csv(res_df, file = "FinalFuncRes.csv", row.names = FALSE)

# Quantifying the process of community assembly
# the OTU table file (Tab delimited txt file)
com.file="otutab_16S_agriculture_rarefy_A.csv"
# the phylogenetic tree file
tree.file="root_tree_16s(1).nwk"
# the function information
clas.file="FinalFuncRes.csv"
# the environmental varialbes
env.file="allenvfactor_16SA.csv"
# the treatment informaiton table
treat.file="treat16SA.csv"

# key parameter setting
prefix="LKF_ITSA"  # prefix of the output file names. usually use a project ID.
rand.time=1000 # randomization time, 1000 is usually enough. 
nworker=102 # nworker is thread number for parallel computing
memory.G=400 #to set the memory size as you need (but should be less than the available space in your hard disk

# load R packages and data
library(iCAMP)
library(ape)
# the OTU table file (Tab delimited txt file)
comm=read.csv(com.file, header = TRUE,row.names = 1) 
# the phylogenetic tree file
tree=read.tree(file = tree.file)
# the classification (taxonomy) information
clas=read.csv(clas.file, header = TRUE)
# the environmental varialbes
env=read.csv(env.file, header = TRUE, row.names = 1)
treat=read.csv(treat.file, header = TRUE,row.names = 1)

# match sample IDs in OTU table and treatment information table
sampid.check=match.name(rn.list=list(comm=comm,env=env,treat=treat))
comm=sampid.check$comm
comm=comm[,colSums(comm)>0,drop=FALSE] # if some unmatched samples were removed, some OTUs may become ghosts, then you may use this line to remove them if necessary.
env=sampid.check$env # skip this if you do not have env.file
treat=sampid.check$treat

# Load calculation results
load("Fan16SA.iCAMP.Summary.rda")
# Initialization result
res = data.frame(HeS="",HoS="",DL="",HD="",DR="",Stochasticity="")
i = 1
flag_count = 1
while(i<ncol(clas)){
  ini_comm = comm
  ini_clas = clas
  ini_tree = tree
  # Obtain the i-th functional group
  guild = data.frame()
  guild <- clas[,c(i,i+1)]
  colnames(guild) <- c("otuid","func")
  guild = guild[which(guild[,1]!=""),]
  rownames(guild) <- guild[,1]
  guild = subset(guild,select = func)
  ini_clas = guild
  
  #  match OTU IDs in OTU table and tree file
  spid.check=match.name(cn.list=list(comm=ini_comm),rn.list=list(clas=ini_clas),tree.list=list(tree=ini_tree))
  ini_comm=spid.check$comm
  ini_clas=spid.check$clas
  ini_tree=spid.check$tree
  flag = ini_clas[1,1]
  
  Guild=icamp.cate(icamp.bins.result = icbin, comm = ini_comm, cate =ini_clas,
                   treat=treat, silent = FALSE, between.group = FALSE) #treat is necessary
  if (length(Guild$Ptx[,which(grepl(flag,colnames(Guild$Ptx)))])==0) {
    i = i + 2
    next
  }
  
  res[flag_count,] = as.character(Guild$Ptx[,which(grepl(flag,colnames(Guild$Ptx)))])
  rownames(res)[flag_count] = flag
  print(res[flag_count,])
  i = i + 2
  flag_count = flag_count + 1
}
write.csv(res,"16SAresult.csv")