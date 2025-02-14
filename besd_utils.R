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

# set up for plotting
colz <- c('#4F8FD0','#FF4E88','#88D66C','#FF8225')
goodbadcolz <- c('#FF6363','#F8B400','#125B50')
font_add("Garamond", "GARA.TTF")
font_families()
showtext_auto()
theme_set(theme_classic() + 
            theme(text = element_text(size=31,family = "Garamond")))


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


# function to plot cluster distributiuon over vxstat
vxstatclusterplot2 <- function(data=df,clusters,legendposition='right',title=''){
  
  # how well do these categories span vxstatus?
  clusters[['cluster']] <- clusters[[names(clusters)[grepl('_cluster',names(clusters))]]]
  clusters[['cllab']]   <- clusters[[names(clusters)[grepl('_cllab',names(clusters))]]]
  
  tmp <- merge(data[,c('id','vxstatus'),with=F],clusters)[,.(N=.N),by=.(cluster,cllab,vxstatus)]
  tmp[, pct := N/sum(N), by = .(cluster,cllab)]
  tmp[, tot := sum(N),by=cluster]
  tmp[, pctoftot := round((tot/sum(unique(tot)))*100)]
  tmp[, cllab := paste0(cllab,'\n(N=',prettyNum(tot,big.mark=','),
                        ' [',pctoftot,'%])')]
  
  ggplot(tmp) +
    geom_bar(aes(pct,cllab,fill=vxstatus),stat='identity') +
    scale_fill_manual(values=rev(c('#FF6363','#F8B400','#125B50')),name='Vx Status') +
    ylab('') +
    ggtitle(title) +
    theme(axis.text.x = element_text(lineheight = .1)) +
    scale_x_continuous(labels=scales::percent, name='Percent of Responses') +
    theme(legend.position=legendposition)
  
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



## Demographic variable plot
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


# extract r-squared
calculate_r2 <- function(modobj) {
  # Check the family of the model
  family_type <- family(modobj)$family
  # Calculate log-likelihood of the full model
  logLik_full <- logLik(modobj)
  # Fit the null model (intercept only)
  null_model <- update(modobj, . ~ 1)
  logLik_null <- logLik(null_model)
  # Calculate R-squared or pseudo R-squared based on the model type
  if (family_type == "gaussian") {
    # Linear regression: R-squared
    r_squared <- summary(modobj)$r.squared
  } else {
    # For other GLM families: Pseudo R-squared (McFadden)
    r_squared <- 1 - (as.numeric(logLik_full) / as.numeric(logLik_null))
  }
  return(r_squared)
}






# Province maps
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






# Regression result table
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
    lvlstab <- gsub('`','',unique(tmp[variable==v]$level))  # updated 16-JAN-2025
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
  
  # added 16-JAN-2025
  tmp[, level := gsub('`','',level)]
  tmp[, level := factor(level, levels=unique(tmp$level))]
  
  #clean up table
  tmp <- tmp[, c('varnice','level','estimate','conf.low','conf.high','z','p')]
  tmp[, c('estimate','conf.low','conf.high','z','p') := lapply(.SD, function(x) round(x,2)), .SDcols = c('estimate','conf.low','conf.high','z','p')]
  names(tmp) <- c('Variable','Level','Odds Ratio','95% CI Lower','95% CI Upper','Z','p')
  if(!is.null(datlab)) tmp$datlab <- datlab
  return(tmp)
}

