## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Roy Burstein 
## July/August 2024
## ECV 2022 and 2023 BeSD Analysis
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## SET UP

# Load needed libraries
library(data.table)
library(tidyverse)
library(poLCA)
library(haven)
library(fastDummies)
library(lavaan)

# plotting things
font_add("Garamond", "GARA.TTF")
font_families()
showtext_auto()
theme_set(theme_classic() + 
            theme(text = element_text(size=25,family = "Garamond")))

# Set up data and code pointers 
# (Note Branly: update these to point to locations on your computer)
setwd('C:/Users/royb/OneDrive - Bill & Melinda Gates Foundation/code/ecv_besd')
datadir22 <- 'C:/data/ECV 2022'
datadir23 <- 'C:/data/ECV 2023'

# Load in ECV 2022 and 2023 data and cast as data.table
df22 <- as.data.table(read_dta(file.path(datadir22,'ECV_2022_VAC_Menages_Mere_Enfants_VT_28052023.dta')))
df23 <- as.data.table(read_dta(file.path(datadir23,'ECV_2023_Vaccination_V4_Menages_Mere_Enfants_Dataset.dta')))

# Load in the codebook
cb <- fread('codebook.csv')

# load in shapefile
pshp <- readRDS('province_shapefile.RDS')


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## DATA CLEANING

## rename variables according to codebook. keep only those merged from codebook
setnames(df22,cb$ecv2022,cb$name,skip_absent=TRUE)
setnames(df23,cb$ecv2023,cb$name,skip_absent=TRUE)

missing22 <- cb$name[which(!cb$name %in% names(df22))]
missing23 <- cb$name[which(!cb$name %in% names(df23))]

message(paste('Note! the following are in the codebook but not found in the data and will not be kept:
              \n ECV2022:', paste0(missing22,collapse=','),'\n ECV2023:', paste0(missing23,collapse=',')))

df22 <- df22[,cb$name[!cb$name%in%missing22],with=FALSE]
df23 <- df23[,cb$name[!cb$name%in%missing23],with=FALSE]

# collapse them into one dataset
df22$svyyear <- 'ECV 2022'
df23$svyyear <- 'ECV 2023'
df <- rbind(df22,df23)

##vaccine outcomes
# Binarize vaccination status
cols <- names(df) %>% str_subset("vacc_")
df <- df %>% mutate_at(cols, function(x) ifelse(x==2, 0, x))

# Zero-dose (Gavi definition of penta1==0)
df[, zd := case_when(vacc_penta1==0~1,
                     vacc_penta1==1~0)]
# Under-immunized (MCV1==0 & penta1==1)
df[, ui := as.numeric(vacc_penta1==1&(vacc_mcv1==0|vacc_penta3==0))]

# Vaccination status variable
df[, vxstatus := 'Complete']
df[ui==1, vxstatus := 'Under-Immunized']
df[zd==1, vxstatus := 'Zero-Dose']


# id variable
df[, id := .I]


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## RUN LATENT CLASS ANALYSIS

## Function to run LCA on this data
clusterfunc <- function(data=df,          # data source
                        vars,             # variable names to cluster
                        k,                # number of clusters
                        printplot=TRUE,   # print a plot when complete?
                        seed=12334,       # seed for replicability
                        prefix='',        # prefix for cluster variable names
                        method='kmeans'   # method for clustering (kmeans or lca)
                        ){          
  # dummydata and update codebook
  dummydat <- data.table(id=data$id)
  for(v in vars){
    #message(v)
    #  data[[v]] <- as.character(data[[v]])
    NAS <- is.na(data[[v]])
    if(sum(NAS)>0){
      if(is.numeric(data[[v]]))   data[NAS==TRUE][[v]] <- 99
      if(is.character(data[[v]])) data[NAS==TRUE][[v]] <- 'NA'
    }
    dummydat <- cbind(dummydat, onehot(v, data=data,
                                       onetwo=ifelse(method=='lca',TRUE,FALSE))) # for kmeans
    #data[[v]] <- labelled(to_character(data[[v]]))
  }
  codebookclust <<- unique(codebookclust)
  
  ### visualize
  # set up data to visualze
  set.seed(seed)
  if(method=='kmeans'){
    clusters <- kmeans(dummydat[,-1], centers = k)
  }else if(method=='lca'){ 
    f <- as.formula(paste('cbind(',paste(paste0('`',names(dummydat[,-1]),'`'),collapse=','),')~1'))
    clusters <- poLCA(formula=f, data=dummydat[,-1], nclass=k, maxiter=1500,tol=1e-10,nrep=5)
    clusters$cluster <- clusters$predclass
  } else {
    stop('method must be set to either "kmeans" or "lca"')
  }
  
  
  print(aggregate(data$zd, list(clusters$cluster), FUN=mean))
  print(data.table(cl=clusters$cluster)[,.(n=.N),by=cl][order(cl)])
  
  cl=data.table(cluster = factor(clusters$cluster), id=dummydat$id)
  dtc <- merge(data[,c('id',vars),with=F], cl, by='id')
  dtc <- melt(dtc, id.vars=c('id','cluster'))
  dtc <- dtc[, .(N=.N), by = .(cluster,variable,value)]
  dtc[, pct := N/sum(N), by = .(variable,cluster)][order(variable,value,cluster)]
  
  dtc[, value:=as.character(value)]
  codebookclust[, num:=as.character(num)]
  dtc <- merge(dtc, codebookclust, by.x=c('variable','value'), by.y=c('v','num'))
  dtc[, value2:= 1:.N]
  dtc[, label2:=factor(str_wrap(label,25),level=unique(str_wrap(label,25)))]
  
  # plot question responses by cluster
  plott <- ggplot(dtc) +
    geom_tile(aes(cluster,label2,fill=pct)) +
    facet_wrap(.~paste0(str_wrap(q,20),'...'),scales='free_y') +
    theme(axis.text.y=element_text(size=12,lineheight = .5),
          strip.text=element_text(size=9),
          legend.position = 'right') +
    ylab('') +
    scale_fill_gradient(low='white',high='black', 
                        label=scales::percent,
                        name = '')
  
  if(printplot==TRUE) print(plott)
  
  # clean up for accounting
  cl[[paste0(prefix,'_cluster')]] <- cl$cluster
  cl$cluster <- NULL
  
  return(list(cl=cl,dtc=dtc,plot=plott,model=clusters))

}

# first, each scaled variable will be one-hot encoded. 
# function for that here
onehot <- function(v, data=df, dummydata, onetwo=FALSE){
  
  # first, add the variable to the codebook
  # codebook <<- codebook[v!=v] # clear if already there
  if(is.character(data[[v]])){
    
    codebookclust <<- rbind(codebookclust,
                       data.table(v=v,
                                  q=v,
                                  num=unique( data[[v]] ),
                                  label=unique( data[[v]] )))
    
  } else{  
    codebookclust <<- rbind(codebookclust,
                       data.table(v=v,
                                  q=attributes(data[[v]])$label,
                                  num=attributes(data[[v]])$labels,
                                  label=names(attributes(data[[v]])$labels)))
  }
  # return square dummydata for this variable
  tmp <- data.table(data[[v]])
  names(tmp) <- v
  
  # make dummies
  output <- fastDummies::dummy_cols(tmp)[,-1]
  
  # if set to 1/2 instead of 0/1 (for poLCA package)
  if(onetwo==TRUE) output <- output+1
  
  return(output)
}


# function to plot cluster distributiuon over vxstat
vxstatclusterplot <- function(data=df,clusters,legendposition='right'){
  
  # how well do these categories span vxstatus?
  clusters[['cluster']] <- clusters[[names(clusters)[grepl('_cluster',names(clusters))]]]
  clusters[['cllab']]   <- clusters[[names(clusters)[grepl('_cllab',names(clusters))]]]
  
  tmp <- merge(data[,c('id','vxstatus'),with=F],clusters)[,.(N=.N),by=.(cluster,cllab,vxstatus)]
  tmp[, pct := N/sum(N), by = .(cluster,cllab)]
  tmp[, tot := sum(N),by=cluster]
  tmp[, cllab := str_wrap(paste0(cllab,' (N=',prettyNum(tot,big.mark=','),')'),25)]
  
  ggplot(tmp) +
    geom_bar(aes(pct*100,cllab,fill=vxstatus),stat='identity') +
    scale_fill_manual(values=rev(c('#FF6363','#F8B400','#125B50')),name='Vx Status') +
    ylab('') +
    xlab('Percent of responses') +
    theme(legend.position=legendposition)
  
}




# initiate a codebook of all variable levels used
codebookclust <- data.table()


## THINKING-FEELING
tf_vars <- cb[besd_category=='thinkingfeeling' & zd_relevant==1]$name
tf_cl <- clusterfunc(data=df,vars=tf_vars, k = 3, prefix='thinkingfeeling',method='lca')
tf_cl$plot

 # !! NOTE! These orderings can change if data or seed changes. Review output plot of clusterfunc each time
tf_labs <- c('2. Moderate Thinking/Feeling', '1. Positive Thinking/Feeling','3. Negative Thinking/Feeling')
for(i in 1:length(tf_labs)) tf_cl$cl[thinkingfeeling_cluster==i,thinkingfeeling_cllab:=tf_labs[i]]
vxstatclusterplot(clusters=tf_cl$cl)

## SOCIAL PROCESSES
sp_vars <- cb[besd_category=='socialprocesses' & zd_relevant==1]$name
sp_cl <- clusterfunc(data=df, vars=sp_vars, k = 3, prefix='socialprocesses',method='lca')

 # !!NOTE! 
sp_labs <- c('3. Negative Norms', '1. Positive Norms','2. Positive Norms with low autonomy')
for(i in 1:length(sp_labs)) sp_cl$cl[socialprocesses_cluster==i,socialprocesses_cllab:=sp_labs[i]]
vxstatclusterplot(clusters=sp_cl$cl) # may just need 2 categories here 


## MOTIVATION (just one var)
m_vars <- cb[besd_category=='motivation' & zd_relevant==1]$name
m_cl <- clusterfunc(data=df, vars=m_vars, k = 2, prefix='motivation',method='lca')

# !!NOTE! 
m_labs <- c('1. Motivated', '2. Unmotivated or unsure')
for(i in 1:length(m_labs)) m_cl$cl[motivation_cluster==i,motivation_cllab:=m_labs[i]]
vxstatclusterplot(clusters=m_cl$cl)


## PRACTICAL ISSUES 
pi_vars <- cb[besd_category=='practicalissues' & zd_relevant==1]$name
pi_cl <- clusterfunc(data=df, vars=pi_vars, k = 2, prefix='practicalissues',method='lca')

# !!NOTE! 
pi_labs <- c('1. Fewer Practical Issues', '2. More Practical Issues')
for(i in 1:length(pi_labs)) pi_cl$cl[practicalissues_cluster==i,practicalissues_cllab:=pi_labs[i]]
vxstatclusterplot(clusters=pi_cl$cl)

# add these to df
df <- merge(df,tf_cl$cl)
df <- merge(df,sp_cl$cl)
df <- merge(df,m_cl$cl)
df <- merge(df,pi_cl$cl)


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## DESCRIPTIVE ANALYSIS OF LCA RESULTS

# combined plot of all BeSD domains
egg::ggarrange(
  vxstatclusterplot(clusters=tf_cl$cl, legendposition='none'),
  vxstatclusterplot(clusters=sp_cl$cl, legendposition='none'), 
  vxstatclusterplot(clusters=m_cl$cl, legendposition='none'),
  vxstatclusterplot(clusters=pi_cl$cl))

# coverage by motivation versus practical issues
tmp <- merge(m_cl$cl,pi_cl$cl)
tmp <- merge(tmp,df[,c('id','vxstatus'),with=F])
jtools::summ(
    glm(vxstatus=='Zero-Dose'~motivation_cllab*practicalissues_cllab,
        data=tmp,family='binomial'), exp=TRUE)


# plot change over survey for each component -- little change over time

ctmp <- melt(df, id.vars=c('id','svyyear','vxstatus','weight'))
ctmp <- ctmp[grepl('_cllab',variable)]
ctmp <- ctmp[, .(N=sum(weight)), by = .(variable,value,svyyear)]
ctmp <- ctmp[, pct := N/sum(N), by = .(variable,svyyear)]
ggplot(ctmp) + 
  geom_bar(aes(pct*100,svyyear,fill=value),stat='identity') +
  facet_wrap(.~gsub('_cllab','',variable)) +
  scale_fill_manual(values=c('#125B50','#125B50','#125B50','#125B50','#F8B400','#FF6363','#F8B400','#FF6363','#FF6363','#FF6363'),
                    name='') +
  ylab('') +
  xlab('Percent of responses') +
  theme(legend.position='right')

# user: Write a haiku about the data
# assistant: ok, here is a haiku about the data:
#           "Clusters of feeling
#            Motivation and issues
#            Zero-dose revealed"




# province maps
provmap <- function(data=df,cllab,legendposition='none'){
  data$cllab <- data[[cllab]]
  tmp <- data[, .(N=sum(weight)), by = c('province','cllab')]
  tmp <- tmp[, pct := N/sum(N), by = c('province')]
  tmp <- merge(pshp,tmp,by='province')
  
  ggplot(tmp) +
    geom_sf(aes(fill=pct*100),color='black') +
    scale_fill_gradient(low='white',high='black',
                        limits=c(0,100),
                        name='Percent of responses') +
    theme_minimal() +
    facet_wrap(~cllab) +
    theme(legend.position='right') +
    #labs(title=gsub('_cllab','',cllab)) +
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.position=legendposition)
  
}
egg::ggarrange(provmap(cllab='thinkingfeeling_cllab'),
               provmap(cllab='socialprocesses_cllab'),
               provmap(cllab='motivation_cllab'),
               provmap(cllab='practicalissues_cllab',legendposition='right'),
               ncol=1)


# over survey year
# over geography (maps)
# by other demographics (mothers age, religion, number of kids, etc. )















## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## SEM






