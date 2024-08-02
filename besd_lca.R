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
library(sysfonts)
library(showtext)

# plottishowtext# plotting things
font_add("Garamond", "GARA.TTF")
font_families()
showtext_auto()
theme_set(theme_classic() + 
            theme(text = element_text(size=25,family = "Garamond")))
cols <- c('#478CCF','#FF4E88','#88D66C','#FF8225')
bcnsord <- c("Thinking/Feeling", "Social Processes", "Motivation", "Practical Issues")

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
besd_cb <- fread('besd_codebook.csv')

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
vxstatclusterplot <- function(data=df,clusters,legendposition='right',title=''){
  
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
    ggtitle(title) +
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

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# combined plot of all BeSD domains and outcomes
pdf('figure2 -- LCA and vx outcomes.pdf',width=27,height=18,onefile=FALSE)
egg::ggarrange(
  vxstatclusterplot(clusters=tf_cl$cl, legendposition='none',title='Thinking/Feeling'),
  vxstatclusterplot(clusters=sp_cl$cl, legendposition='none',title='Social Processes'), 
  vxstatclusterplot(clusters=m_cl$cl,  legendposition='none',title='Motivation'),
  vxstatclusterplot(clusters=pi_cl$cl, title='Practical Issues'))
dev.off()


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


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# province maps
provmap <- function(data=df,cllab,legendposition='right',maxcol='black',title=''){
  data$cllab <- data[[cllab]]
  tmp <- data[, .(N=sum(weight)), by = c('province','cllab')]
  tmp <- tmp[, pct := N/sum(N), by = c('province')]
  tmp <- merge(pshp,tmp,by='province')
  
  # get centroid
  ggplot(tmp) +
    geom_sf(aes(fill=pct*100),color='black') +
    scale_fill_gradient(low='white',high=maxcol,
                       # limits=c(0,100),
                        name='Percent of\nresponses') +
    theme_minimal() +
    facet_wrap(~cllab) +
    theme(legend.position='right') +
    ggtitle(title) +
    #labs(title=gsub('_cllab','',cllab)) +
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.position=legendposition,
          text=element_text(size=20,family = "Garamond"))

  
}
pdf('figure3 -- LCA and province.pdf',width=16,height=27,onefile=FALSE)
egg::ggarrange(provmap(cllab='thinkingfeeling_cllab',maxcol=cols[1],title='Thinking/Feeling'),
               provmap(cllab='socialprocesses_cllab',maxcol=cols[2],title='Social Processes'),
               provmap(cllab='motivation_cllab'     ,maxcol=cols[3],title='Motivation'),
               provmap(cllab='practicalissues_cllab',maxcol=cols[4],title='Practical Issues'),
               ncol=1)
dev.off()

# over survey year
# over geography (maps)
# by other demographics (mothers age, religion, number of kids, etc. )




## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## FIGURE 1 - Question breakdown across each BeSD domain/component (one large figure)


# read in the besd_codebook and format it
besd_cbl <- melt(besd_cb, id.vars=c('name',	'description',	'besd_category','besd_category_nice',	'question'),
                 variable.name = 'value', value.name = 'value_label')
besd_cbl[,value := gsub('response_','',value)]
besd_cbl <- besd_cbl[value_label!='']
besd_cbl[, value_label:=paste0(value,'. ',value_label)]

# collapse
cllabvars <- names(df)[grepl('_cllab',names(df))]
tmp <- df[,c('id',pi_vars,tf_vars, sp_vars ,m_vars,cllabvars),with=F]
tmp <- melt(tmp, id.vars=c('id',cllabvars))
tmp[, value:=as.character(value)]

# merge codebook
tmp <- merge(tmp, besd_cbl, by.x=c('variable','value'), by.y=c('name','value'))

# call out specific LCA class for each row (since they are question associated)
tmp <- tmp %>%
  mutate(cllab_associated = case_when(
    besd_category == "motivation"      ~ motivation_cllab,
    besd_category == "thinkingfeeling" ~ thinkingfeeling_cllab,
    besd_category == "socialprocesses" ~ socialprocesses_cllab,
    besd_category == "practicalissues" ~ practicalissues_cllab,
    TRUE ~ NA_character_
  ))

# get pct
tmp <- tmp[, .(N=.N), by = c( 'question','value_label','cllab_associated','besd_category','besd_category_nice')]
tmp <- tmp[, pct := N/sum(N), by = c('question','cllab_associated','besd_category','besd_category_nice')]

# add zeros
for(q in unique(tmp$question)){
  for(v in unique(tmp[question==q]$value_label)){
    for(c in unique(tmp[question==q]$cllab_associated)){
      if(nrow(tmp[question==q & value_label==v & cllab_associated==c])==0){
        tmp <- rbind(tmp,data.table(question=q,value_label=v,
                                    cllab_associated=c,N=0,pct=0,
                                    besd_category=unique(tmp[question==q]$besd_category),
                                    besd_category_nice=unique(tmp[question==q]$besd_category_nice)))
      }
    }
  }
}


# plot a facetted tile plot

listofplots <- list()
for(i in 1:4){
  col <- cols[i]
  bcn <- bcnsord[i]
  
  listofplots[[i]] <-
  ggplot(tmp[besd_category_nice==bcn]) +
    geom_tile(aes(value_label,cllab_associated,fill=pct),color=col) +
    facet_wrap(~str_wrap(question, width = 40),scales='free_x') +
   theme(axis.text.y=element_text(size=18,lineheight = .5),
         axis.text.x=element_text(size=18,lineheight = .5,angle=45,hjust=1),
        strip.text=element_text(size=18),
        legend.position = 'none') +
    # add percent text
    geom_text(aes(value_label,cllab_associated,label=paste0(round(pct*100,0),'%')),
              size=6, color='black') +
    ylab('') + #Latent Class Assigned') +
    xlab('') + #Question Response') +
    ggtitle(paste0(bcn)) +
    scale_fill_gradient(low='white',high=col, 
                        label=scales::percent, 
                        limits=c(0,1),
                        name = '') 
} 
pdf('figure1 -- LCA and questions.pdf',width=20,height=35,onefile=FALSE)
egg::ggarrange(listofplots[[1]],listofplots[[2]],listofplots[[3]],listofplots[[4]],ncol=1)
dev.off()


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Correlation of BeSD components with BeSD components



## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## SEM

















## ________________
# quick one, card retention over time for carolina
df[, dayssincestart := as.numeric(int_date-min(int_date)), by = .(svyyear,interviewer)]
df[,interviewer2:=factor(paste0(svyyear,'_',interviewer))]
df[, Ninterviewer := .N, by = .(interviewer2)]

ggplot(df[Ninterviewer>10]) +
  geom_histogram(aes(dayssincestart))

jtools::summ(glm(opv1_card==1 ~ dayssincestart + svyyear + province ,
                 data=df[Ninterviewer>10 & dayssincestart<40],
                 family=binomial))

library(mgcv)
m <- mgcv::bam(opv1_card==1 ~ s(dayssincestart) ,
               data=df[Ninterviewer>30 & dayssincestart<30])
plot(m)

tmp <- df[Ninterviewer>10 & dayssincestart<30,.(N=.N),by=.(dayssincestart,opv1_card)]
tmp[, pct := N/sum(N), by = .(dayssincestart)]
ggplot(tmp[opv1_card==1 & dayssincestart>0]) +
  geom_point(aes(dayssincestart,pct*100,size=N)) +
  ylab('Proportion Reported Card Retention') +
  xlab('Days since interviewer started') +
  theme(legend.position='right')

# int hour
df[, hour := as.numeric(substr(int_time,1,2))]
df[,hrsincestart :=hour-min(hour), by = .(int_date,interviewer2)]
ggplot(df[Ninterviewer>10 & svyyear=='ECV 2022' & hrsincestart<18]) +
  geom_histogram(aes(hrsincestart))
jtools::summ(glm(opv1_card==1 ~ hour, data=df[svyyear=='ECV 2022'], family='binomial'),exp=T)

tmp <- df[svyyear=='ECV 2022' & hrsincestart<20,.(N=.N),by=.(hrsincestart,opv1_card)]
tmp[, pct := N/sum(N), by = .(hrsincestart)]

ggplot(tmp[opv1_card==1 &N>1000]) +
  geom_point(aes(hrsincestart,pct*100,size=N)) +
  ylab('Proportion Reported Card Retention') +
  xlab('Hrs since interviewer started working that day') +
  scale_x_continuous(breaks=0:10) +
  theme(legend.position='right')



# per interviewer average retention by hour since started working
tmp <- df[svyyear=='ECV 2022' & Ninterviewer>30]
tmp[,m := mean(opv1_card==1), by = .(interviewer2)]
tmp[,hrsincestart :=hour-min(hour), by = .(int_date,interviewer2)]
tmp <- tmp[,.(ratio=100*mean(opv1_card==1)/m), by = .(interviewer2,hrsincestart)]
tmp <- tmp[,.(ratio=median(ratio,na.rm=T)), by = .(hrsincestart)]
ggplot(tmp[hrsincestart<=11]) +
  geom_point(aes(hrsincestart,ratio*1)) +
  geom_hline(yintercept=100,color='pink') +
  ylab('Avg. ratio of card retention reported') +
  xlab('Hrs since interviewer started working that day') +
  scale_x_continuous(breaks=0:10) +
  theme(legend.position='right')

# 
# library(lme4)
# m<-glmer(opv1_card==1 ~ dayssincestart + svyyear + province + (1|interviewer2), data=df, family='binomial')
# jtools::summ(m,exp=TRUE)
# 


