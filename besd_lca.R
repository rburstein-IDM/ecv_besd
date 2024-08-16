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
library(egg)
library(patchwork)
library(sjPlot)
library(broom)
library(jtools)
library(Hmisc)
library(forcats)
library(tableone)
library(survey)


# plottishowtext# plotting things
font_add("Garamond", "GARA.TTF")
font_families()
showtext_auto()
theme_set(theme_classic() + 
            theme(text = element_text(size=25,family = "Garamond")))
colz <- c('#478CCF','#FF4E88','#88D66C','#FF8225')
bcnsord <- c("Thinking/Feeling", "Social Processes", "Motivation", "Practical Issues")

# Set up data and code pointers 
# (Note Branly: update these to point to locations on your computer)
setwd('C:/Users/royb/OneDrive - Bill & Melinda Gates Foundation/code/ecv_besd')
datadir22 <- 'C:/data/ECV 2022'
datadir23 <- 'C:/data/ECV 2023'


skipprep <- TRUE

  if(skipprep==FALSE){
    
  # Load in ECV 2022 and 2023 data and cast as data.table
  df22 <- as.data.table(read_dta(file.path(datadir22,'ECV_2022_VAC_Menages_Mere_Enfants_VT_28052023.dta')))
  df23 <- as.data.table(read_dta(file.path(datadir23,'ECV_2023_Vaccination_V4_Menages_Mere_Enfants_Dataset.dta')))
  
  # Load in the codebook
  cb <- fread('codebook_bm.csv')
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
  df[, fvxstatus := factor(vxstatus, levels=c('Zero-Dose','Under-Immunized','Complete'))]
  
  # id variable
  df[, id := .I]
  
  # haven_labelled variables to factors using their haven labels, starting with education_level_hhlead
  labelledvars <- c('strate', 'sex_hhlead', 'etat_civil_hhlead', 'education_level_hhlead','religion_hhlead',
                    'ethnie_hhlead','etat_civil_caregiver','age_caregiver','education_level_caregiver',
                    'religion_caregiver','wealth_quintile')
  for(v in labelledvars) df[[v]] <- as_factor(df[[v]])
  
  
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
  # test between 2 and 6 latent classes for each domain
  # this is slow, so only run if needed
  
  if(TRUE == FALSE) {
    
    # run the models
    tflcalist <- splcalist <- mlcalist <- pilcalist <- list()
    for(k in 2:6){
      message(k)
      tflcalist[[k]] <- clusterfunc(data=df,vars=tf_vars, k = k, prefix='thinkingfeeling',method='lca',printplot=FALSE)
      splcalist[[k]] <- clusterfunc(data=df,vars=sp_vars, k = k, prefix='socialprocesses',method='lca',printplot=FALSE)
      if(k<4) mlcalist[[k]]  <- clusterfunc(data=df,vars=m_vars, k = k, prefix='motivation',method='lca',printplot=FALSE)
      pilcalist[[k]] <- clusterfunc(data=df,vars=pi_vars, k = k, prefix='practicalissues',method='lca',printplot=FALSE)
    }
    
    # likelihood ratio tests to compare models (get pval)
    lrt <- function(mod1, mod2){
      pchisq(-2 * (mod1$llik - mod2$llik),
             mod2$npar - mod1$npar,
             lower.tail = FALSE)
    }
    entropyfunc <- function(model){ 
      if(is.null(model)){
        return(NA)
      } else {
        posterior_probs <- model$posterior
        entropy_individual <- function(p) { p <- p[p > 0]; -sum(p * log(p)) }
        entropy <- mean(apply(posterior_probs, 1, entropy_individual))
        return(entropy)
      }
    }
    
    # extract goodness of fit stats (AIC, BIC, G^2, X^2) for each from the $model object
    grabfitstats <- function(modlist,ks,domainname){
      message(domainname)
      out <-
      data.table(
        domain  = domainname,
        k       = ks,
        loglik  = unlist(lapply(modlist[ks],  function(x) x$model$llik)),
        AIC     = unlist(lapply(modlist[ks],  function(x) x$model$aic)),
        BIC     = unlist(lapply(modlist[ks],  function(x) x$model$bic)),
        G2      = unlist(lapply(modlist[ks],  function(x) x$model$Gsq)),
        X2      = unlist(lapply(modlist[ks],  function(x) x$model$Chisq)),
        entropy = lapply(lapply(modlist[ks],  function(x) x$model),entropyfunc),
        npar    = unlist(lapply(modlist[ks],  function(x) x$model$npar))
      )
     for(K in ks[2]:ks[length(ks)]){ # llik ratio test
      out[k==K, lrt_pval := lrt(modlist[[K-1]]$model,modlist[[K]]$model)]
     }
      return(out)
    }
    # make a table across all indicators
    fitstats <- rbindlist(
      list(grabfitstats(tflcalist,2:6,'Thinking/Feeling'),
           grabfitstats(splcalist,2:6,'Social Processes'),
           grabfitstats(mlcalist,2:3,'Motivation'),
           grabfitstats(pilcalist,2:6,'Practical Issues'))
    )
    fitstats
    fwrite(fitstats,'fitstats.csv')
    
    # Note that we get statistically significant improvements when add classes
    #  though, given the large SS, this is expected. Practical considerations are
    #  important here too. 
    
    # TODO. could bootstrap to test stability
    
  }
  
  
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #  how does each question do in explaining the classes (were the BeSD core questions the right ones?)
  #  TODO
  
  
  
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
  summ(
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
  egg::ggarrange(provmap(cllab='thinkingfeeling_cllab',maxcol=colz[1],title='Thinking/Feeling'),
                 provmap(cllab='socialprocesses_cllab',maxcol=colz[2],title='Social Processes'),
                 provmap(cllab='motivation_cllab'     ,maxcol=colz[3],title='Motivation'),
                 provmap(cllab='practicalissues_cllab',maxcol=colz[4],title='Practical Issues'),
                 ncol=1)
  dev.off()
  
  # over survey year
  # over geography (maps)
  # by other demographics (mothers age, religion, number of kids, etc. )
  
  
  
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## FIGURE 1 - Question breakdown across each BeSD domain/component (one large figure)
  ###           Item-Response Probabilities (Conditional Probabilities):
  ###           The probability of endorsing each item given class membership, which helps in understanding the characteristics of each class.
  
  
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
    col <- colz[i]
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
  ## Demographic Variables plots
  
  # make age of caregiver into a group of 5 year bands, as a factor
  df[, age_caregiver_grp := cut(as.numeric(as.character(age_caregiver)), 
                                breaks = seq(15, 65, by = 5), 
                                labels = c('18-20','21-25','26-30','31-35','36-40','41-45','46-50','51-55','56-60','60+'))]
  df[as.numeric(as.character(age_caregiver))>50, age_caregiver_grp := '> 50']
  df[as.numeric(as.character(age_caregiver))<18, age_caregiver_grp := '< 18']
  df[, age_caregiver_grp := factor(age_caregiver_grp, levels=rev(c('< 18','18-20','21-25','26-30','31-35','36-40',
                                                               '41-45','46-50','> 50')))]
  # AS PLOT
  demdistplot <- function(dv, weighted = TRUE) {
    message(paste(dv,'----------'))
    tmp <- copy(df)
    
    if(dv=='age_caregiver_grp') 
      tmp <- tmp[as.numeric(as.character(age_caregiver)) < 46]
    
    if(weighted==FALSE) tmp$weight <- 1
    tmp <- tmp[,.(N=sum(weight)),by=c(dv,cllabvars)]
    tmp <- melt(tmp, id.vars=c(dv,'N'), variable.name = 'besd_category', value.name = 'value_label')
    tmp <- tmp[, NN := sum(N), by = c(dv,'besd_category')]
    tmp[, pct := N/NN]
    tmp[, dv := get(dv)]
    tmp[, besd_category := gsub('_cllab','',besd_category)]
    tmp <- merge(tmp, unique(besd_cbl[,c('besd_category','besd_category_nice')]), by='besd_category')
    
      
    plotlist <- list()
    for(c in 1:4){
      cat <- c("Thinking/Feeling","Social Processes","Motivation" , "Practical Issues")[c]    
      message(cat)
      tmpp <- tmp[besd_category_nice==cat]
      nvals <- length(unique(tmpp$value_label))
      CLZ <- c('#F9DBBA','#5B99C2','#1A4870')
      if(nvals==2) CLZ <- CLZ[c(1,3)]
      
      plotlist[[c]] <-
        ggplot(tmpp) +
          geom_bar(aes(pct,dv,fill=value_label),position='stack',stat='identity') +
          ggtitle(cat) +
          theme(axis.text.y=element_text(angle=0,hjust=1),
                legend.position='top')+
          scale_x_continuous(labels=scales::percent, name = '') +
          scale_fill_manual(values=rev(CLZ), name='') + ylab('') +
         guides(fill = guide_legend(ncol = 1))
    }
    
    plotlist[[5]] <-
      ggplot(tmpp) +
        geom_bar(aes(dv,N),fill='darkgrey',stat='identity') +
        theme(axis.text.y=element_blank(),
            axis.line.y =element_blank(),
            axis.ticks.y=element_blank(),
            legend.position='right')+
        ggtitle(dv) + 
        scale_fill_manual(values=CLZ, name='') + xlab('') + ylab('') 
      
    # patchwork plot them
    layout <- plotlist[[5]] /
      (plotlist[[1]]+plotlist[[2]]) /
      (plotlist[[3]]+plotlist[[4]]) 
      plot_layout(heights = c(0.3, 1, 1))
    
    # Display the combined plot
    return(layout)
  }
  
  
  demvars <- 
    c('strate', 'sex_hhlead', 'etat_civil_hhlead', 'education_level_hhlead','religion_hhlead',
      'ethnie_hhlead','etat_civil_caregiver','age_caregiver_grp','education_level_caregiver',
      'religion_caregiver','wealth_quintile')
  
  pdf('figure4 -- LCA and demographics.pdf',width=15,height=17)
  for(dvv in demvars) plot(demdistplot(dvv))
  dev.off()
  
  # AS TABLE (TODO)
  


# set up variables to include in regressions
# set up as factors with hypothesized 'worser' levels higher up (will be >1 OR in regression)
demvarsreg <- 
  c('strate', 'sex_hhlead', 'ethnie_hhlead','etat_civil_caregiver',
    'age_caregiver_grp','education_level_caregiver',
    'religion_caregiver','wealth_quintile')

               
                         
# relabel factors
label(df$strate) <- 'Setting'
levels(df$strate) <- c('Urban','Rural')

levels(df$sex_hhlead) <- c('male','female')
label(df$sex_hhlead)  <- 'Sex of HH Head'

label(df$ethnie_hhlead) <- 'Ethnicity of HH Head'

label(df$etat_civil_caregiver)  <- 'Marital Status of Caregiver'
levels(df$etat_civil_caregiver) <- c('Married','Common-law','Separated','Single','Divorced','Widowed')
df$etat_civil_caregiver <- relevel(df$etat_civil_caregiver, ref = "Married")

label(df$age_caregiver_grp) <- 'Age of Caregiver'

label(df$education_level_caregiver) <- 'Education Level of Caregiver'
levels(df$education_level_caregiver) <- c('Never been to school','Primary','Secondary',
                                          'Higher','Do not know','No response')
df$education_level_caregiver <- relevel(df$education_level_caregiver, ref = "Higher")

label(df$religion_caregiver) <- 'Religion of Caregiver'
levels(df$religion_caregiver) <- c('No religion','Catholic','Protestant',
                                    'Kimbanguist','Muslim','Revival/Independent Church',
                                    'Other','Dont know','No response ')
df$religion_caregiver <- relevel(df$religion_caregiver, ref = "Catholic")


label(df$wealth_quintile) <- 'Wealth Quintile'
levels(df$wealth_quintile) <- c('Poorest','Poorer','Middle','Richer','Richest')
# order all levels richest to poorest 
df$wealth_quintile <- fct_relevel(df$wealth_quintile, 'Richest','Richer','Middle','Poorer','Poorest')

label(df$thinkingfeeling_cllab) <- 'Thinking/Feeling'
label(df$socialprocesses_cllab) <- 'Social Processes'
label(df$motivation_cllab) <- 'Motivation'
label(df$practicalissues_cllab) <- 'Practical Issues'

label(df$province) <- 'Province'
label(df$age_caregiver_grp) <- 'Age of Caregiver'





## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## SAVE ENVIRONMENT HERE 
#save.image('besd_lca.RData')
## LOAD HERE
  } # end skipprep

load('besd_lca.RData')




## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Table 1

table1 <- CreateTableOne(vars = c(cllabvars,demvarsreg,'province','svyyear'),
                         strata = 'vxstatus', data = df[c_age_m>11])
print(table1, showAllLevels = TRUE, test = TRUE)
table1_df <- as.data.frame()
write.csv(print(table1, showAllLevels = TRUE, test=FALSE), file = "table1.csv", row.names = TRUE)





## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## simple regression




## ~~~~
## How much do tf and sp explain motivation? 
forml1 <- as.formula(paste0('motivation_cllab=="1. Motivated" ~',paste0(tf_vars,collapse='+'),
                            '+',paste0(sp_vars,collapse='+')))
forml2 <- as.formula(paste0('motivation_cllab=="1. Motivated" ~ thinkingfeeling_cllab + socialprocesses_cllab'))

# with the full model, only about 45% of the variance in motivation is explained, meaning we should keep it in a full model
summ(glm(forml1,data=df, family='binomial'),exp=TRUE) # R2 ~ 0.45
summ(glm(forml2,data=df, family='binomial'),exp=TRUE) # R2 ~ 0.39



## ~~~~
## How much do the besd categories explain ZD (versus both UI and FIC)
forml1 <- as.formula(paste0('zd ~',paste0(tf_vars,collapse='+'),'+',paste0(sp_vars,collapse='+'),
                            '+',paste0(pi_vars,collapse='+'),'+',paste0(m_vars,collapse='+')))
forml2 <- as.formula(paste0('zd ~ thinkingfeeling_cllab + socialprocesses_cllab + practicalissues_cllab + motivation_cllab'))

# 
summ(glm(forml1,data=df, family='binomial'),exp=TRUE) # R2 ~ 0.40
summ(glm(forml2,data=df, family='binomial'),exp=TRUE) # R2 ~ 0.36


## ~~~~
## Add in demographic and geog variables to the above 
forml1 <- as.formula(paste0('zd ~',paste0(tf_vars,collapse='+'),'+',
                                   paste0(sp_vars,collapse='+'),'+',
                                   paste0(pi_vars,collapse='+'),'+',
                                   paste0(m_vars,collapse='+'),
                                   '+ province +',
                                   paste0(demvarsreg,collapse='+')))
forml2 <- as.formula(paste0('zd ~ thinkingfeeling_cllab + socialprocesses_cllab + practicalissues_cllab + motivation_cllab + ',
                            paste0(demvarsreg,collapse='+'), '+ province'))

# The lack of diff in R2 may indicate that the explanatory value differntial in 
#  full Qs vs BeSD summary Qs is explained by the demographic variables
summ(glm(forml1,data=df, family='binomial'),exp=TRUE) # R2 ~ 0.50
modzd <- glm(forml2,data=df, family='binomial')
summ(modzd,exp=TRUE) # R2 ~ 0.50
# demographic/geog only: R2 ~ 0.29 (so a big diff coming in BeSD only) 
#  <<GOOD STORY FOR PAPER, LEAVE WITH THIS SIMPLE MODEL AND LEAVE MORE COMPLEX MODELLING FOR NEXT PAPER>>





## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## FULL MODELS FOR PAPER
## TODO: DO USING SVY PACKAGE
# (ZD, >12m, ZD versus any vx)
zddf <- df[c_age_m>=12]
formlzd <- as.formula(paste0('zd ~ thinkingfeeling_cllab + socialprocesses_cllab + practicalissues_cllab + motivation_cllab + ',
                            paste0(demvarsreg,collapse='+'), '+ province'))
modzd <- glm(formlzd,data=zddf, family='binomial') #, weights = zddf$weight)

# (UI, >12m, UI versus FIC
uidf <- df[c_age_m>=12 & zd==0]
formlui <- as.formula(paste0('ui ~ thinkingfeeling_cllab + socialprocesses_cllab + practicalissues_cllab + motivation_cllab + ',
                            paste0(demvarsreg,collapse='+'), '+ province'))
modui <- glm(formlui,data=uidf, family='binomial')



### TABLE
# export the summary table in a nice format
cleanregressiontable <- function(modobj,datlab=NULL){
  if(is.null(modobj$survey.design)) {
    tmp <- as.data.table(summ(modobj,exp=T)$coeftable)
    names(tmp) <- c('estimate','conf.low','conf.high','z','p')
  }else {
    tmp <-  as.data.table(summ(modobj,exp=T)$coeftable)
    tmp[, conf.low := exp(confint(modobj)[,1])]
    tmp[, conf.high := exp(confint(modobj)[,2])]
    tmp <- tmp[, c('exp(Est.)','conf.low','conf.high','t val.','p')]
    names(tmp) <- c('estimate','conf.low','conf.high','z','p')
  }
  
  tmp$level <- rownames(summ(modobj,exp=T)$coeftable)
  # sep out variable name from varlevl variable
  for(v in all.vars(modobj$formula)[-1]){
    tmp[grepl(v,level), variable := v]
    tmp[grepl(v,level), level := gsub(v,'',level)]
  }
  tmp$variable[1] <- '_Intercept'
  # make order_id variable by variable
  tmp[, order_id := 1:.N, by = variable]
  # add the ref category from each
  reftemplate<-copy(tmp)[1,]
  reftemplate[, c('estimate','conf.low','conf.high') := 1]
  reftemplate[, c('z','p') := NA]
  reftemplate[, order_id := 0]
  reftab <- copy(tmp)[0,]
  for(v in all.vars(modobj$formula)[-1]){
    lvlsdat <- unique(modobj$data[[v]])
    lvlstab <- unique(tmp[variable==v]$level)
    ref     <- lvlsdat[!lvlsdat%in%lvlstab]
    if(length(ref)==1){
      reftemplate[, variable := v]
      reftemplate[, level := paste0(ref,' (ref)')]
      reftab <- rbind(reftab,reftemplate)
    } else {
      message(paste('No reference found for var',v))
    }
  }
  tmp[, variable := factor(variable, levels=c('_Intercept',all.vars(modobj$formula)[-1]))]
  tmp <- rbind(tmp,reftab)[order(variable,order_id)]
  
  # nice variable names and orderings
  for(v in unique(tmp$variable)){
    if(label(modobj$data[[v]])=='') tmp[variable==v,varnice := v]
    else tmp[variable==v,varnice := label(modobj$data[[v]])]
  }
  tmp[,varnice := factor(varnice, levels=unique(tmp[order(variable)]$varnice))]
  
  tmp[, level := factor(level, levels=tmp$level)]
  
  #clean up table
  tmp <- tmp[, c('varnice','level','estimate','conf.low','conf.high','z','p')]
  tmp[, c('estimate','conf.low','conf.high','z','p') := lapply(.SD, function(x) round(x,2)), .SDcols = c('estimate','conf.low','conf.high','z','p')]
  names(tmp) <- c('Variable','Level','Odds Ratio','95% CI Lower','95% CI Upper','Z','p')
  if(!is.null(datlab)) tmp$datlab <- datlab
  return(tmp)
}
zdtab <- cleanregressiontable(modzd,datlab='Zero-Dose versus any vaccination')
uitab <- cleanregressiontable(modui,datlab='Under-Immunized versus Fully-Immunized')

# export the summary tables 
as.data.table(zdtab) %>% fwrite('ZD model coefficients.csv')
as.data.table(uitab) %>%  fwrite('UI model coefficients.csv')


## PLOT
pdf('figure6 -- all regression coefficients.pdf',width=25,height=20)
ggplot(rbind(zdtab,uitab)) +
  geom_vline(xintercept=1,linetype='dashed',color='grey') +
  
  geom_pointrange(aes(y=Level,x=`Odds Ratio`,xmin=`95% CI Lower`,xmax=`95% CI Upper`,color=datlab),
                  position=position_dodge(width=0.3)) +
  facet_wrap(~Variable,scales='free_y') +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        axis.text.y=element_text(size=15)) +
  scale_color_manual(values=colz, name='') +
  scale_x_continuous(limits=c(0,NA), name='Odds Ratio')+
  ylab('') +
  ggtitle('Odds Ratios of BeSD Categories in Predicting Zero-Dose and Under-Immunized Status') +
  theme(legend.position='bottom')
dev.off()


pdf('figure5 -- BeSD regression coefficients.pdf',width=18,height=12)
ggplot(rbind(zdtab,uitab)[Variable %in% c('Thinking/Feeling','Social Processes','Motivation','Practical Issues')]) +
  geom_vline(xintercept=1,linetype='dashed',color='grey') +
  
  geom_pointrange(aes(y=Level,x=`Odds Ratio`,xmin=`95% CI Lower`,xmax=`95% CI Upper`,color=datlab),
                  position=position_dodge(width=0.3)) +
  facet_wrap(~Variable,scales='free_y') +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        axis.text.y=element_text(size=15)) +
  scale_color_manual(values=colz, name='') +
  scale_x_continuous(limits=c(0,NA), name='Odds Ratio')+
  ylab('') +
  ggtitle('Odds Ratios of BeSD Categories in Predicting Zero-Dose and Under-Immunized Status') +
  theme(legend.position='bottom')
dev.off()







## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Regression using survey weights/complex design

# TODO: HOW TO DEAL WITH DIFF SURVEY YEARS? JUST STRATA? ASK DALE
survey_design_zd <- svydesign(
  id = ~zone + area,   # ZS as strata, AS as primary sampling units (clusters)
  strata = ~interaction(province, svyyear),  # ?
  weights = ~weight,   # The weights calculated based on the sampling probabilities
  data = zddf,          # The data frame
  nest = TRUE          # Nesting within each higher-level sampling unit
)

modzd_svy <- svyglm(
  formlzd,
  design = survey_design_zd,
  family = binomial  
)

ggplot(
  rbind(cleanregressiontable(modzd,datlab='Zero-Dose (unweighted)'),
        cleanregressiontable(modzd_svy,datlab='Zero-Dose (weighted)')
        )) + #[Variable %in% c('Thinking/Feeling','Social Processes','Motivation','Practical Issues')]) +
  geom_vline(xintercept=1,linetype='dashed',color='grey') +
  
  geom_pointrange(aes(y=Level,x=`Odds Ratio`,xmin=`95% CI Lower`,xmax=`95% CI Upper`,color=datlab),
                  position=position_dodge(width=0.3)) +
  facet_wrap(~Variable,scales='free_y') +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        axis.text.y=element_text(size=15)) +
  scale_color_manual(values=colz, name='') +
  scale_x_continuous(limits=c(0,NA), name='Odds Ratio')+
  ylab('') +
  theme(legend.position='bottom')

# same for UI
survey_design_ui <- svydesign(
  id = ~zone + area,   
  strata = ~interaction(province, svyyear),  # ?
  weights = ~weight,
  data = uidf,         
  nest = TRUE        
)
modui_svy <- svyglm(
  formlui,
  design = survey_design_ui,
  family = binomial  
)
ggplot(
  rbind(
        cleanregressiontable(modzd,datlab='Zero-Dose (unweighted)'),
        cleanregressiontable(modzd_svy,datlab='Zero-Dose (weighted)')
  )) +
  geom_vline(xintercept=1,linetype='dashed',color='grey') +
  
  geom_pointrange(aes(y=Level,x=`Odds Ratio`,xmin=`95% CI Lower`,xmax=`95% CI Upper`,color=datlab),
                  position=position_dodge(width=0.3)) +
  facet_wrap(~Variable,scales='free_y') +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        axis.text.y=element_text(size=15)) +
  scale_color_manual(values=colz, name='') +
  scale_x_continuous(limits=c(0,NA), name='Odds Ratio')+
  ylab('') +
  theme(legend.position='bottom')



## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# demographic/geog only: R2 ~ 0.29 (so a big diff coming in BeSD only) 
#  <<GOOD STORY FOR PAPER, LEAVE WITH THIS SIMPLE MODEL AND LEAVE MORE COMPLEX MODELLING FOR NEXT PAPER>>





## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## SEM


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
## ordinal regression
#library(MASS)
#m <- polr(fvxstatus ~ thinkingfeeling_cllab + socialprocesses_cllab + motivation_cllab + practicalissues_cllab, data=data.frame(df))
#summary(m)


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
## multinomial regression
## Spell out the different assumptions from each model
# m <- multinom(vxstatus ~ thinkingfeeling_cllab + socialprocesses_cllab + motivation_cllab + practicalissues_cllab, data=df)
# summary(m)



















# 
# 
# 
# ## ________________
# # quick one UNRELATED, card retention over time for carolina
# df[, dayssincestart := as.numeric(int_date-min(int_date)), by = .(svyyear,interviewer)]
# df[,interviewer2:=factor(paste0(svyyear,'_',interviewer))]
# df[, Ninterviewer := .N, by = .(interviewer2)]
# 
# ggplot(df[Ninterviewer>10]) +
#   geom_histogram(aes(dayssincestart))
# 
# summ(glm(opv1_card==1 ~ dayssincestart + svyyear + province ,
#                  data=df[Ninterviewer>10 & dayssincestart<40],
#                  family=binomial))
# 
# library(mgcv)
# m <- mgcv::bam(opv1_card==1 ~ s(dayssincestart) ,
#                data=df[Ninterviewer>30 & dayssincestart<30])
# plot(m)
# 
# tmp <- df[Ninterviewer>10 & dayssincestart<30,.(N=.N),by=.(dayssincestart,opv1_card)]
# tmp[, pct := N/sum(N), by = .(dayssincestart)]
# ggplot(tmp[opv1_card==1 & dayssincestart>0]) +
#   geom_point(aes(dayssincestart,pct*100,size=N)) +
#   ylab('Proportion Reported Card Retention') +
#   xlab('Days since interviewer started') +
#   theme(legend.position='right')
# 
# # int hour
# df[, hour := as.numeric(substr(int_time,1,2))]
# df[,hrsincestart :=hour-min(hour), by = .(int_date,interviewer2)]
# ggplot(df[Ninterviewer>10 & svyyear=='ECV 2022' & hrsincestart<18]) +
#   geom_histogram(aes(hrsincestart))
# summ(glm(opv1_card==1 ~ hour, data=df[svyyear=='ECV 2022'], family='binomial'),exp=T)
# 
# tmp <- df[svyyear=='ECV 2022' & hrsincestart<20,.(N=.N),by=.(hrsincestart,opv1_card)]
# tmp[, pct := N/sum(N), by = .(hrsincestart)]
# 
# ggplot(tmp[opv1_card==1 &N>1000]) +
#   geom_point(aes(hrsincestart,pct*100,size=N)) +
#   ylab('Proportion Reported Card Retention') +
#   xlab('Hrs since interviewer started working that day') +
#   scale_x_continuous(breaks=0:10) +
#   theme(legend.position='right')
# 
# 
# 
# # per interviewer average retention by hour since started working
# tmp <- df[svyyear=='ECV 2022' & Ninterviewer>30]
# tmp[,m := mean(opv1_card==1), by = .(interviewer2)]
# tmp[,hrsincestart :=hour-min(hour), by = .(int_date,interviewer2)]
# tmp <- tmp[,.(ratio=100*mean(opv1_card==1)/m), by = .(interviewer2,hrsincestart)]
# tmp <- tmp[,.(ratio=median(ratio,na.rm=T)), by = .(hrsincestart)]
# ggplot(tmp[hrsincestart<=11]) +
#   geom_point(aes(hrsincestart,ratio*1)) +
#   geom_hline(yintercept=100,color='pink') +
#   ylab('Avg. ratio of card retention reported') +
#   xlab('Hrs since interviewer started working that day') +
#   scale_x_continuous(breaks=0:10) +
#   theme(legend.position='right')
# 
# # 
# # library(lme4)
# # m<-glmer(opv1_card==1 ~ dayssincestart + svyyear + province + (1|interviewer2), data=df, family='binomial')
# # summ(m,exp=TRUE)
# # 
# 
# 
