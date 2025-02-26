## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Roy Burstein 
## Script for analysis, figures, and tables for the manuscript:
##   Measuring behavioral and social drivers (BeSD) and their association 
##   with zero-dose and under-vaccinated status in children in the Democratic 
##   Republic of Congo. 
## February 2025
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Set up data and code pointers 
setwd('C:/Users/royb/OneDrive - Bill & Melinda Gates Foundation/code/ecv_besd')
datadir22 <- 'C:/data/ECV 2022'
datadir23 <- 'C:/data/ECV 2023'

# script to load libraries, bespoke functions, and plot setup
source('./besd_utils.R')

# if we skip prep (to save time), load prepped data instead
skipprep <- TRUE

## Data prep.. or skip to load
if(skipprep==FALSE){
    
  # Load in ECV 2022 and 2023 data and cast as data.table
  df22o <- as.data.table(read_dta(file.path(datadir22,'ECV_2022_VAC_Menages_Mere_Enfants_VT_28052023.dta')))
  df23o <- as.data.table(read_dta(file.path(datadir23,'ECV_2023_Vaccination_V4_Menages_Mere_Enfants_Dataset.dta')))
  
  # active copy for analysis (to avoid long load times in interactive coding)
  df22 <- copy(df22o)
  df23 <- copy(df23o)
  
  # Load in the codebook
  cb      <- fread('codebook_bm.csv')
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
  
  ## vaccine outcomes
  # Binarize vaccination status
  cols <- names(df) %>% str_subset("vacc_")
  df <- df %>% mutate_at(cols, function(x) ifelse(x==2, 0, x))
  
  
  ## ORIGINAL WAY WE DEFINED OUTCOMES 
  # Zero-dose (Gavi definition of penta1==0)
  df[, zd_orig := case_when(vacc_penta1==0~1,
                            vacc_penta1==1~0)]
  # Under-immunized (MCV1==0 & penta1==1)
  df[, ui_orig := as.numeric(vacc_penta1==1&(vacc_mcv1==0|vacc_penta3==0))]
  
  # Vaccination status variable
  df[, vxstatus_orig := 'Complete']
  df[ui_orig==1, vxstatus_orig := 'Under-Immunized']
  df[zd_orig==1, vxstatus_orig := 'Zero-Dose']
  df[, fvxstatus_orig := factor(vxstatus_orig, levels=c('Zero-Dose','Under-Immunized','Complete'))]
 
  ## JAN 2025: Update definition to DPT0 and eight basic vaccines 
  ## Make an alternative definition of UI based on a set of vx: (BCG+Polio x 3 + Penta x3 + MCV x 1).
  df[, numcorevx := as.numeric(vacc_bcg    +  vacc_opv1   + vacc_opv2   + 
                               vacc_opv3   +  vacc_penta1 + vacc_penta2 + 
                               vacc_penta3 +  vacc_mcv1)]
  
  df[, zd_hard := as.numeric(numcorevx == 0)] # hard definition (True ZD)
  df[, zd := zd_orig] # Using the penta0 definition of ZD
  df[, ui  := as.numeric(numcorevx > 0 & numcorevx < 8)]
  
  # this is the outcome variable we will use
  df[, vxstatus := 'Fully-Vaccinated']
  df[ui==1, vxstatus := 'Under-Vaccinated']
  df[zd==1, vxstatus := 'Zero-Dose']
  df[, fvxstatus := factor(vxstatus, levels=c('Zero-Dose','Under-Vaccinated','Fully-Vaccinated'))]

  # id variable
  df[, id := .I]
  
  # haven_labelled variables (from original stata file) to factors using their haven labels
  labelledvars <- c('strate', 'sex_hhlead', 'etat_civil_hhlead', 'education_level_hhlead','religion_hhlead',
                    'ethnie_hhlead','etat_civil_caregiver','age_caregiver','education_level_caregiver',
                    'religion_caregiver','wealth_quintile')
  for(v in labelledvars) df[[v]] <- as_factor(df[[v]])
  
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## RUN LATENT CLASS ANALYSIS
  
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
  
  sp_labs <- c('3. Negative Norms', '1. Positive Norms','2. Positive Norms with low autonomy')
  for(i in 1:length(sp_labs)) sp_cl$cl[socialprocesses_cluster==i,socialprocesses_cllab:=sp_labs[i]]
  vxstatclusterplot(clusters=sp_cl$cl) # may just need 2 categories here 
  
  
  ## MOTIVATION (just one var)
  m_vars <- cb[besd_category=='motivation' & zd_relevant==1]$name
  m_cl <- clusterfunc(data=df, vars=m_vars, k = 2, prefix='motivation',method='lca')
  
  m_labs <- c('1. Higher motivation', '2. Lower motivation')
  for(i in 1:length(m_labs)) m_cl$cl[motivation_cluster==i,motivation_cllab:=m_labs[i]]
  vxstatclusterplot(clusters=m_cl$cl)
  
  
  ## PRACTICAL ISSUES 
  pi_vars <- cb[besd_category=='practicalissues' & zd_relevant==1]$name
  pi_cl <- clusterfunc(data=df, vars=pi_vars, k = 2, prefix='practicalissues',method='lca')
  
  pi_labs <- c('1. Fewer Practical Issues', '2. More Practical Issues')
  for(i in 1:length(pi_labs)) pi_cl$cl[practicalissues_cluster==i,practicalissues_cllab:=pi_labs[i]]
  vxstatclusterplot(clusters=pi_cl$cl)
  
  # add these to df
  df <- merge(df,tf_cl$cl)
  df <- merge(df,sp_cl$cl)
  df <- merge(df,m_cl$cl)
  df <- merge(df,pi_cl$cl)
  
  
  
  ## Add Demographics
  # set up variables to include in regressions
  # set up as factors with hypothesized 'worser' levels higher up (will be >1 OR in regression)
  demvarsreg <- 
    c('strate', 'sex_hhlead', 'ethnie_hhlead','etat_civil_caregiver',
      'age_caregiver_grp','education_level_caregiver',
      'religion_caregiver','wealth_quintile','occupation_hhlead')
  
  # make age of caregiver into a group of 5 year bands, as a factor
  df[, age_caregiver_grp := cut(as.numeric(as.character(age_caregiver)), 
                                breaks = seq(15, 65, by = 5), 
                                labels = c('18-20','21-25','26-30','31-35','36-40','41-45','46-50','51-55','56-60','60+'))]
  df[as.numeric(as.character(age_caregiver))>50, age_caregiver_grp := '> 50']
  df[as.numeric(as.character(age_caregiver))<18, age_caregiver_grp := '< 18']
  df[, age_caregiver_grp := factor(age_caregiver_grp, levels=rev(c('< 18','18-20','21-25','26-30','31-35','36-40',
                                                                   '41-45','46-50','> 50')))]
  
  
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
  # order from richest to poorest 
  df$wealth_quintile <- fct_relevel(df$wealth_quintile, 'Richest','Richer','Middle','Poorer','Poorest')
  
  label(df$thinkingfeeling_cllab) <- 'Thinking/Feeling'
  label(df$socialprocesses_cllab) <- 'Social Processes'
  label(df$motivation_cllab) <- 'Motivation'
  label(df$practicalissues_cllab) <- 'Practical Issues'
  
  label(df$province) <- 'Province'
  label(df$age_caregiver_grp) <- 'Age of Caregiver'
  
  df$occupation_hhlead <- as.factor(df$occupation_hhlead)
  label(df$occupation_hhlead) <- 'Occupation of HH Head'
  levels(df$occupation_hhlead) <- c('Without profession','Teacher','Civil servant','Farmer/herder',
                                    'Fisherman','Trader','Laborer','Other Profession','Student')
  

  # expanded varlist for plotting later on
  demvars <- 
    c('strate', 'sex_hhlead', 'etat_civil_hhlead', 'education_level_hhlead','religion_hhlead',
      'ethnie_hhlead','etat_civil_caregiver','age_caregiver_grp','education_level_caregiver',
      'religion_caregiver','wealth_quintile','occupation_hhlead')


  ## individual BeSD questions to  nice format
  # read in the besd_codebook and format it
  besd_cbl <- melt(besd_cb, id.vars=c('name',	'description',	'besd_category','besd_category_nice',	'question'),
                   variable.name = 'value', value.name = 'value_label')
  besd_cbl[,value := gsub('response_','',value)]
  besd_cbl <- besd_cbl[value_label!='']

  nicelistofBQs <- c()
  for(BQ in unique(cb[zd_relevant==1]$name)){
    message(BQ)
    newname <- besd_cbl[name==BQ]$description[1]
    df[[newname]] <-  factor(df[[BQ]], 
                             levels=as.numeric(besd_cbl[name==BQ]$value), 
                             labels=besd_cbl[name==BQ]$value_label)
    nicelistofBQs <- c(nicelistofBQs,newname)
  }
  
  # hotfix for some level orderings to be best to worst
  df[, `Service satisfaction` := factor(`Service satisfaction`, 
        levels=c('Very satisfied','Moderately satisfied','A little satisfied','Not at all satisfied '))]
  
  df[, `Affordability` := factor(`Affordability`,
         levels=c('Very easy','Moderately easy','Quite difficult','Not at all easy'))]
  
  df[,`Mother's travel autonomy` := factor(`Mother's travel autonomy`,
         levels=c('No','Yes'))]
  
  df[, `Ease of access` := factor(`Ease of access`,
         levels=c('Very easy','Moderately easy','Quite difficult','Not at all easy'))]
  
  df[, `Confidence in  healthworkers` := factor(`Confidence in  healthworkers`,
         levels=c('Very much','Moderately','A little','Not at all'))]
  
  df[, `Confidence in vaccine safety` := factor(`Confidence in vaccine safety`,
         levels=c('Very safe','Moderately safe','A little safe','Not at all safe'))]
  
  df[, `Confidence in vaccine benefits` := factor(`Confidence in vaccine benefits`,
         levels=c('Very important','Moderately important','A little important','Not at all important'))]
  
  

  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## SAVE ENVIRONMENT HERE 
  # remove unneeded large memory objects
  rm(list=c('df22','df23','df22o','df23o'))
  save.image('besd_lca4.RData')
  
} # end skip prep

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## LOAD ENVIRONMENT HERE
load('besd_lca4.RData')




## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Table 1  

table1 <- CreateTableOne(vars   = c(cllabvars,
                                    nicelistofBQs,
                                    demvarsreg,
                                    'province',
                                    'svyyear'),
                         strata = 'vxstatus', 
                         data   = df[c_age_m>11])
print(table1, showAllLevels = TRUE, test = TRUE)
write.csv(print(table1, showAllLevels = TRUE, test=FALSE), 
          file = "manu_output/table1.csv", row.names = TRUE)


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Figure 1: Each BeSD Variable and Outcome. 


# read in the besd_codebook and format it
besd_cbl <- melt(besd_cb, id.vars=c('name',	'description',	'besd_category','besd_category_nice',	'question'),
                 variable.name = 'value', value.name = 'value_label')
besd_cbl[,value := gsub('response_','',value)]
besd_cbl <- besd_cbl[value_label!='']
besd_cbl[, value_label:=paste0(value,'. ',value_label)]

# collapse
cllabvars <- names(df)[grepl('_cllab',names(df))]
tmp <- df[c_age_m>=12,c('id',pi_vars,tf_vars, sp_vars ,m_vars,cllabvars,'vxstatus','province'),with=F]
tmp <- melt(tmp, id.vars=c('id',cllabvars,'vxstatus','province'))
tmp[, value:=as.character(value)]

# merge codebook
tmp <- merge(tmp, besd_cbl, by.x=c('variable','value'), by.y=c('name','value'))

## FIGURE
tmptab <- tmp[, .(N=.N), by = c('question','value_label','vxstatus','besd_category_nice')]
tmptab[, pct := N/sum(N), by = c('question','value_label')]
tmptab[, labN := sum(N), by = c('question','value_label')]
tmptab[, vallab := paste0(value_label,' (N=',prettyNum(labN,big.mark=','),')')]

# plot a facetted tile plot
bcnsord <- c('Thinking/Feeling','Social Processes','Motivation','Practical Issues')
listofplots <- list()
for(i in 1:4){
  bcn <- bcnsord[i]
  
  listofplots[[i]] <-
    
    ggplot(tmptab[besd_category_nice==bcn]) +
    geom_bar(aes(vallab,pct,fill=vxstatus),stat='identity') +
    facet_grid(besd_category_nice~str_wrap(question, width = 17),scales='free_x') +
    scale_y_continuous(labels=scales::percent) +
    theme(axis.text.y=element_text(size=18,lineheight = .5),
          axis.text.x=element_text(size=18,lineheight = 1,angle=90,hjust=1),
          strip.text=element_text(size=18),
          legend.position = 'right') +
    ylab('') + #Latent Class Assigned') +
    xlab('') + #Question Response') +
    ggtitle('') +
    scale_fill_manual(values=rev(c('#FF6363','#F8B400','#125B50')),name='')
  
}


pdf('manu_output/figure 1 -- besd questions and vxstat.pdf',width=20,height=40,onefile=FALSE)
egg::ggarrange(listofplots[[1]],listofplots[[2]],listofplots[[3]],listofplots[[4]],
               ncol=1,
               heights=c(5,5,5,5))
dev.off()



## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  SI Figure 1: specific predictive validity of survey questions

output <- data.table()
for(q in  cb[zd_relevant==1]$name){
  message(q)
  zdmod<-glm(as.formula(paste0('zd ~ as.factor(',q,')')),data=df[c_age_m>=12], family='binomial')
  uimod<-glm(as.formula(paste0('ui ~ as.factor(',q,')')),data=df[c_age_m>=12 & zd==0], family='binomial')
  output <- rbind(output,
                  data.table(name=q,
                             `ZD versus any vaccination (UV+FV)`  = calculate_r2(zdmod),
                             `UV versus FV` = calculate_r2(uimod)))
}
output <- merge(output, besd_cb[,c('name','besd_category_nice','description')], by='name')
output$besd_category_nice <- factor(output$besd_category_nice, 
                                    levels=c("Thinking/Feeling", "Social Processes","Motivation","Practical Issues" ))
output[name %in% c('tf_benefits','m_intent','sp_familynorms','pi_knowwhere','pi_affordability'),
       priority := 'Priority Question']
output[is.na(priority),priority := '']

pdf('manu_output/figure SI -- question specific predictive validity.pdf',width=18,height=13)
ggplot(melt(output,id=c('name','besd_category_nice','priority','description'))) +
  geom_point(aes(value,description,color=priority),size=6) +
  facet_grid(besd_category_nice~variable,scales='free') +
  xlab('Pseudo R-Squared') +
  ylab('') +
  scale_color_manual(values=c('black','red'),name='') +
  theme(panel.grid.major.x = element_line(color = "grey"),
        panel.grid.major.y = element_line(color = "grey"),
        panel.grid.minor.x = element_line(color = "grey"),
        panel.grid.minor.y = element_line(color = "grey"),
        legend.position = 'top',
        # Adjust text and spacing
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.spacing = unit(1, "lines"))  # Add space between facets
dev.off()





## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  Figure 2: Item - repsonse for LCA


# call out specific LCA class for each row (since they are question associated)
tmp <- tmp %>%
  mutate(cllab_associated = case_when(
    besd_category == "motivation"      ~ as.character(motivation_cllab),
    besd_category == "thinkingfeeling" ~ as.character(thinkingfeeling_cllab),
    besd_category == "socialprocesses" ~ as.character(socialprocesses_cllab),
    besd_category == "practicalissues" ~ as.character(practicalissues_cllab),
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
bcnsord <- c('Thinking/Feeling','Social Processes','Motivation','Practical Issues')
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
pdf('manu_output/figure 2 -- LCA and questions.pdf',width=20,height=35,onefile=FALSE)
egg::ggarrange(listofplots[[1]],listofplots[[2]],listofplots[[3]],listofplots[[4]],ncol=1)
dev.off()




## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## SI Figure 3: Corelation between demographics and LCA classes

pdf('manu_output/figure SI -- LCA and demographics.pdf',width=20,height=20)
for(dvv in demvarsreg) plot(demdistplot(dvv))
dev.off()



## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  Figure 3: LCA and vaccination outcomes


pdf('manu_output/figure 3 -- LCA and vx outcomes.pdf',width=27,height=18,onefile=FALSE)
egg::ggarrange(
  vxstatclusterplot2(data=df[c_age_m>=12], clusters=tf_cl$cl, legendposition='none',title='Thinking/Feeling'),
  vxstatclusterplot2(data=df[c_age_m>=12], clusters=sp_cl$cl, legendposition='none',title='Social Processes'), 
  vxstatclusterplot2(data=df[c_age_m>=12], clusters=m_cl$cl,  legendposition='none',title='Motivation'),
  vxstatclusterplot2(data=df[c_age_m>=12], clusters=pi_cl$cl, title='Practical Issues'))
dev.off()




## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  Table 2: Logistic Regression Results with LCA
##  Table 3: Logistic Regression Results with all Qs
##  Figure 4: Forest plot of logistic regression results
##  Supplementary table X: Question-specific regression


## FULL MODELS 
# (ZD, >12m, ZD versus any vx)

# set up survey data and formulas for regression
zddf <- df[c_age_m>=12]
formlzd <- as.formula(paste0('zd ~ thinkingfeeling_cllab + socialprocesses_cllab + practicalissues_cllab + motivation_cllab + ',
                             paste0(demvarsreg,collapse='+'), '+ province + svyyear'))
# (UI, >12m, UI versus FIC
uidf <- df[c_age_m>=12 & zd==0]
formlui <- as.formula(paste0('ui ~ thinkingfeeling_cllab + socialprocesses_cllab + practicalissues_cllab + motivation_cllab + ',
                             paste0(demvarsreg,collapse='+'), '+ province + svyyear'))

# patchfix missing vars as levels
zddf[is.na(`Intention to get the child vaccinated`),`Intention to get the child vaccinated`:='No Response']
zddf[is.na(`Community leader norms`),`Community leader norms`:='No Response']
zddf[is.na(`Vaccination availability`),`Vaccination availability`:='No Response']
uidf[is.na(`Intention to get the child vaccinated`),`Intention to get the child vaccinated`:='No Response']
uidf[is.na(`Community leader norms`),`Community leader norms`:='No Response']
uidf[is.na(`Vaccination availability`),`Vaccination availability`:='No Response']

# Regression set up using survey weights and design variables
survey_design_zd <- svydesign(
  id = ~id, 
  strata = ~province,  
  weights = ~weight,   
  data = zddf,         
  nest = TRUE         
)

survey_design_ui <- svydesign(
  id = ~id, 
  strata = ~province,  
  weights = ~weight,   
  data = uidf,         
  nest = TRUE         
)


# run the two main logistics regression models
modzd_svy <- svyglm(
  formlzd,
  design = survey_design_zd,
  family = binomial  
)

modui_svy <- svyglm(
  formlui,
  design = survey_design_ui,
  family = binomial  
)


## TABLE
# export the summary table in a nice forma
zdtab <- cleanregressiontable(modzd_svy,datlab='Zero-Dose versus any vaccination\n(ZD vs. UV+FV)')
uitab <- cleanregressiontable(modui_svy,datlab='Under-Vaccinated versus Fully-Vaccinated\n(UV vs. FV)')

# export the summary tables 
as.data.table(zdtab) %>%  fwrite('manu_output/table2 -- ZD model coefficients.csv')
as.data.table(uitab) %>%  fwrite('manu_output/table2 -- UI model coefficients.csv')



# Main figure of forest plot

# clean up data for plotting (new svy year convention as year of survey)
plottmp <- rbind(zdtab,uitab)
plottmp[Variable=='svyyear',Variable:='Survey Year']
plottmp[Variable=='Survey Year'&Level=='ECV 2023',Level:='VCS 2024']
plottmp[Variable=='Survey Year'&Level=='ECV 2022 (ref)',Level:='VCS 2023 (ref)']
plottmp[Variable=='Province',Level:=gsub(' Province','',substr(Level,4,10000))]

# one plot thats just the LCA classes


pdf('manu_output/figure 4 -- regression forest plot.pdf',width=27,height=27,onefile=FALSE)
egg::ggarrange(
  ggplot(plottmp[Variable %in% bcnsord]) +
    geom_vline(xintercept=1,linetype='dashed',color='black') +
    geom_pointrange(aes(y=Level,x=`Odds Ratio`,xmin=`95% CI Lower`,xmax=`95% CI Upper`,color=datlab),
                    position=position_dodge(width=0.3),size=1.5,lwd=1.5) +
    facet_wrap(~Variable,scales='free_y',ncol=1) +
    theme(axis.text.y=element_text(size=15),
          panel.grid.major.x = element_line(color = "grey"),
    ) +
    scale_color_manual(values=colz, name='') +
    scale_x_continuous(limits=c(0,9), name='Odds Ratio',
                       breaks = 0:9)+
    ylab('') +
    theme(legend.position='bottom'),
  ggplot(plottmp[!Variable %in% bcnsord & ! Variable %in% c('_Intercept','Province')]) + # save space exclude intercept and prov intercepts
    geom_vline(xintercept=1,linetype='dashed',color='black') +
    geom_pointrange(aes(y=Level,x=`Odds Ratio`,xmin=`95% CI Lower`,xmax=`95% CI Upper`,color=datlab),
                    position=position_dodge(width=0.3),size=1.5,lwd=1.5) +
    facet_wrap(~Variable,scales='free_y',ncol=1) +
    theme(axis.text.y=element_text(size=15),
          # add x gridlines
          panel.grid.major.x = element_line(color = "grey"),
    ) +
    scale_color_manual(values=colz, name='') +
    scale_x_continuous(limits=c(0,9), name='Odds Ratio',
                       breaks = 0:9)+
    ylab('') +
    theme(legend.position='none'),
  ncol = 2)

dev.off()


## Table 3

modzd_svy2 <- svyglm(
  as.formula(paste0('zd ~ ', paste0(c(paste0(paste0("`", nicelistofBQs, "`")),demvarsreg),collapse='+'), '+ province + svyyear')),
  design = survey_design_zd,
  family = binomial  
)


modui_svy2 <- svyglm(
  as.formula(paste0('ui ~ ',  paste0(c(paste0(paste0("`", nicelistofBQs, "`")),demvarsreg),collapse='+'), '+ province + svyyear')),
  design = survey_design_ui,
  family = binomial  
)


zdtab2 <- cleanregressiontable(modzd_svy2,datlab='Zero-Dose versus any vaccination\n(ZD vs. UV+FV)')
uitab2 <- cleanregressiontable(modui_svy2,datlab='Under-Vaccinated versus Fully-Vaccinated\n(UV vs. FV)')

as.data.table(zdtab2) %>%  fwrite('manu_output/table3 -- ZD model coefficients FULL BESD Qs.csv')
as.data.table(uitab2) %>%  fwrite('manu_output/table3 -- UI model coefficients FULL BESD Qs.csv')



## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Supplementary Table 3: VIF
vifzd <- data.table(cbind(gsub('`','',rownames(car::vif(modzd_svy2))),car::vif(modzd_svy2)))
vifui <- data.table(cbind(gsub('`','',rownames(car::vif(modui_svy2))),car::vif(modui_svy2)))
fwrite(vifzd, 'manu_output/supp table X -- ZD model VIF.csv')
fwrite(vifui, 'manu_output/supp table X -- UI model VIF.csv')




## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##  Table 4: Predictive power by specification

# write formulae 
forml_zd_allBeSDQs <- as.formula(paste0('zd ~',paste0(tf_vars,collapse='+'),'+',paste0(sp_vars,collapse='+'),
                                        '+',paste0(pi_vars,collapse='+'),'+',paste0(m_vars,collapse='+')))
forml_zd_prioty    <- as.formula('zd ~ tf_benefits + m_intent + sp_familynorms + pi_knowwhere + pi_affordability')
forml_zd_lca       <- as.formula('zd ~ thinkingfeeling_cllab + socialprocesses_cllab + practicalissues_cllab + motivation_cllab')

forml_ui_allBeSDQs <- as.formula(paste0('ui ~',paste0(tf_vars,collapse='+'),'+',paste0(sp_vars,collapse='+'),
                                        '+',paste0(pi_vars,collapse='+'),'+',paste0(m_vars,collapse='+')))
forml_ui_prioty    <- as.formula('ui ~ tf_benefits + m_intent + sp_familynorms + pi_knowwhere + pi_affordability')
forml_ui_lca       <- as.formula('ui ~ thinkingfeeling_cllab + socialprocesses_cllab + practicalissues_cllab + motivation_cllab')

foml_zd_allBeSDQs_dem <- as.formula(paste0('zd ~',paste0(tf_vars,collapse='+'),'+',paste0(sp_vars,collapse='+'),
                                           '+',paste0(pi_vars,collapse='+'),'+',paste0(m_vars,collapse='+'),'+',
                                           paste0(demvarsreg,collapse='+'),'+ province'))
foml_zd_prioty_dem    <- as.formula(paste0('zd ~ tf_benefits + m_intent + sp_familynorms + pi_knowwhere + pi_affordability +',
                                           paste0(demvarsreg,collapse='+'),'+ province'))
foml_zd_lca_dem       <- as.formula(paste0('zd ~ thinkingfeeling_cllab + socialprocesses_cllab + practicalissues_cllab + motivation_cllab +',
                                           paste0(demvarsreg,collapse='+'),'+ province'))
foml_zd_dem <- as.formula(paste0('zd ~',paste0(demvarsreg,collapse='+'),'+ province'))

foml_ui_allBeSDQs_dem <- as.formula(paste0('ui ~',paste0(tf_vars,collapse='+'),'+',paste0(sp_vars,collapse='+'),
                                           '+',paste0(pi_vars,collapse='+'),'+',paste0(m_vars,collapse='+'),'+',
                                           paste0(demvarsreg,collapse='+'),'+ province'))
foml_ui_prioty_dem    <- as.formula(paste0('ui ~ tf_benefits + m_intent + sp_familynorms + pi_knowwhere + pi_affordability +',
                                           paste0(demvarsreg,collapse='+'),'+ province'))
foml_ui_lca_dem       <- as.formula(paste0('ui ~ thinkingfeeling_cllab + socialprocesses_cllab + practicalissues_cllab + motivation_cllab +',
                                           paste0(demvarsreg,collapse='+'),'+ province'))
foml_ui_dem <- as.formula(paste0('ui ~',paste0(demvarsreg,collapse='+'),'+ province'))

formlist <- list(forml_zd_allBeSDQs,forml_zd_prioty,forml_zd_lca,foml_zd_dem,
                 foml_zd_allBeSDQs_dem,foml_zd_prioty_dem,foml_zd_lca_dem,
                 forml_ui_allBeSDQs,forml_ui_prioty,forml_ui_lca,foml_ui_dem,
                 foml_ui_allBeSDQs_dem,foml_ui_prioty_dem,foml_ui_lca_dem)

# remove ZD and UV from the above
modnamesnice <- c('All 16 BeSD Qs only','5 Priority BeSD Qs only','4 LCA-derived BeSD domains only','Demographics only',
                  'All 16 BeSD Qs + Demographics','5 Priority BeSD Qs + Demographics','4 LCA-derived BeSD domains + Demographics')

# note: not using svyglm because it doesnt use MLE and cant get Psuedo-RS that way
modlist <- list()
for(f in 1:length(formlist)){
  message(f)
  if(as.character(formlist[[f]][2])=='zd'){
    modlist[[f]] <- glm(formlist[[f]],data=zddf, family='binomial')
  } else if(as.character(formlist[[f]][2])=='ui'){
    modlist[[f]] <- glm(formlist[[f]],data=uidf, family='binomial')
  }
}

rsqs <- sapply(modlist,calculate_r2)
tab4 <- data.table(outcome=c(rep('Zero-Dose',7),rep('Under-Vaccinated',7)),model=modnamesnice,rsquared=round(rsqs,2))
tab4 <- reshape(tab4, direction = "wide", idvar = "model", timevar = "outcome")

write.csv(tab4, file = "manu_output/table4 -- predictive power by specification.csv", row.names = FALSE)



## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Supplementary Figure 2: Province maps

pdf('./manu_output/figure SI -- LCA and province.pdf',width=16,height=27,onefile=FALSE)
egg::ggarrange(provmap(cllab='thinkingfeeling_cllab',maxcol=colz[1],title='Thinking/Feeling'),
               provmap(cllab='socialprocesses_cllab',maxcol=colz[2],title='Social Processes'),
               provmap(cllab='motivation_cllab'     ,maxcol=colz[3],title='Motivation'),
               provmap(cllab='practicalissues_cllab',maxcol=colz[4],title='Practical Issues'),
               ncol=1)
dev.off()





## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Supplementary Figure 4: Pairwise correlation heatmaps of cllab variables


# first get all unique pairwise combinations of cllab vars
cllabpairs <- combn(cllabvars,2,simplify=FALSE)

groblist <- list()
for(i in 1:length(cllabpairs)){
  
  tmp <- df[,.(N=sum(weight)),by=c(cllabpairs[[i]],'vxstatus')]
  tmp[,pct:=N/sum(N),by=vxstatus]
  tmp2 <- df[,.(N=sum(weight)),by=c(cllabpairs[[i]])]
  tmp2[,pct:=N/sum(N)]
  tmp2[,vxstatus:='All']
  tmp <- rbind(tmp,tmp2)
  nice1 <- besd_cb[besd_category==gsub('_cllab','',cllabpairs[[i]][1])]$besd_category_nice[1]
  nice2 <- besd_cb[besd_category==gsub('_cllab','',cllabpairs[[i]][2])]$besd_category_nice[1]
  tmp$x<-tmp[[cllabpairs[[i]][1]]]
  tmp$y<-tmp[[cllabpairs[[i]][2]]]
  groblist[[i]] <-
    ggplot(tmp) +
    geom_tile(aes(x,y,fill=pct)) +
    scale_fill_gradient(low='white',high='black', 
                        label=scales::percent,
                        name = 'Percent of responses') +
    geom_text(aes(x,y,label=paste0(round(pct*100,0),'%')),
              size=8, color='#800080') +
    facet_wrap(~vxstatus, nrow=1) + 
    xlab('') + ylab('') +
    ggtitle(paste(nice1,'and',nice2)) +
    theme(axis.text.x=element_text(angle=45,hjust=1),
          axis.text.y=element_text(angle=0,hjust=1),
          legend.position='none')
  
}

pdf('manu_output/figure SI -- LCA and pairwise correlations with vxstatus.pdf',width=22,height=37,onefile=FALSE)
egg::ggarrange(plots=groblist,ncol=1)
dev.off()


























## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Supp Figure (Not used in Manu.)
##  plot change over survey for each component -- little change over time
ctmp <- melt(df, id.vars=c('id','svyyear','vxstatus','weight'))
ctmp <- ctmp[grepl('_cllab',variable)]
ctmp <- ctmp[, .(N=sum(weight)), by = .(variable,value,svyyear)]
ctmp <- ctmp[, pct := N/sum(N), by = .(variable,svyyear)]
png('./manu_output/figure SI -- LCA and survey year.png',width=1200,height=800)
ggplot(ctmp) + 
  geom_bar(aes(pct*100,svyyear,fill=value),stat='identity') +
  facet_wrap(.~gsub('_cllab','',variable)) +
  scale_fill_manual(values=c('#125B50','#125B50','#125B50','#125B50','#F8B400','#FF6363','#F8B400','#FF6363','#FF6363','#FF6363'),
                    name='') +
  ylab('') +
  xlab('Percent of responses') +
  theme(legend.position='right')
dev.off()





## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Suppl. Figure (Unused): distribution of motivation in ZD pop across BESD variables

besdqs <- names(df)[grep('tf_|sp_|pi_',names(df))]
tmp <- copy(df)
for(q in c('m_intent',besdqs)){
  # make haven_labelled characters using their labels as chars
  tmp[[q]] <- haven::as_factor(tmp[[q]])
  tmp[[q]] <- as.character(tmp[[q]])
}
tmp <- tmp[zd==1]
tmp <- melt(tmp[,c('id','m_intent',besdqs),with=F], id.vars=c('id','m_intent'))

pdf('./manu_output/figure SI -- LCA and motivation by BeSD question.pdf',width=15,height=25)
ggplot(tmp[m_intent!='4'&!value %in% c('8','NA') & !is.na(value)]) +
  geom_bar(aes(value,fill=m_intent),position='fill',stat='count') +
  facet_wrap(~variable,scales='free_x') +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        axis.text.y=element_text(size=12,lineheight = .5),
        strip.text=element_text(size=9),
        legend.position = 'right') +
  ylab('') +
  xlab('') +
  ggtitle('Motivation by BeSD Question') +
  scale_fill_manual(values=c('#F8B400','#FF6363','#125B50'),name='Motivation') +
  theme(legend.position='right')
dev.off()




### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Suppl. Figure (Unused): Look at info source wrt to besd categories

sourcevars <- names(df)[grep('infosource',names(df))]

tmp <- melt(df[,c('id','vxstatus','weight','svyyear',cllabvars,sourcevars),with=F], 
            id.vars=c('id','vxstatus','weight','svyyear',cllabvars),
            variable.name = 'infosource', value.name = 'value')
tmp[, infosource := gsub('infosource_','',infosource)]
tmp[, value := as.numeric(value==1)]

# loop over all besd vars
groblist <- list()
for(v in cllabvars){
  tmp$v <- tmp[[v]]
  bnice <- besd_cb[besd_category==gsub('_cllab','',v)]$besd_category_nice[1]
  
  # save groblist

  groblist[[v]] <-
  ggplot(tmp[value==1]) +
    geom_bar(aes(y=infosource,fill=vxstatus)) +
    facet_wrap(~v) +
    theme(axis.text.x=element_text(angle=45,hjust=1)) +
    scale_fill_manual(values=rev(goodbadcolz),name='') +
    ggtitle(bnice) +
    ylab('') +
    theme(legend.position=ifelse(bnice=='Practical Issues','bottom','none'))
} 

pdf('./manu_output/figure SI -- info source and besd categories.pdf',width=20,height=30,onefile=FALSE)
ggarrange(plots=groblist,ncol=1)
dev.off()











## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Supplementary Table 2: Fit statistics for LCA models
#   test between 2 and 6 latent classes for each domain
#   this is slow, so only run if needed

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

}






## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Suppl. Experiments: Start testing SEM approach (out of scope for current paper)
# Load required libraries
library(lavaan)
library(semPlot)
library(data.table)


# Assuming zddf is a data.table, create one-hot encoded (dummy) variables for thinkingfeeling_cllab, socialprocesses_cllab, and practicalissues_cllab
zddf[, `:=`(
  tf_pos = as.numeric(thinkingfeeling_cllab == "1. Positive Thinking/Feeling"),
  tf_mod = as.numeric(thinkingfeeling_cllab == "2. Moderate Thinking/Feeling"),
  
  sp_pos = as.numeric(socialprocesses_cllab == "1. Positive Norms"),
  sp_pos_low_auto = as.numeric(socialprocesses_cllab == "2. Positive Norms with low autonomy"),
  
  pi_more = as.numeric(practicalissues_cllab == "2. More Practical Issues")  # Binary variable for practical issues
)]

# Motivation as a binary variable
zddf[, motiv := as.numeric(motivation_cllab == "1. Motivated")]

# Define the SEM model, adding practical issues (pi) as a downstream variable of motivation and a predictor of zd
sem_model <- '
  # Direct effects
  zd ~ c1 * tf_pos + c2 * tf_mod +
        c3 * sp_pos + c4 * sp_pos_low_auto + 
        b * motiv + pi * pi_more
  
  # Mediator model: Motivation explained by one-hot encoded TF and SP variables
  motiv ~ a1 * tf_pos + a2 * tf_mod + 
           a3 * sp_pos + a4 * sp_pos_low_auto
  
  # Practical issues explained by motivation
  pi_more ~ d * motiv

  # Indirect effects via motivation
  indirect_tf_pos := a1 * b
  indirect_tf_mod := a2 * b
  indirect_sp_pos := a3 * b
  indirect_sp_pos_low_auto := a4 * b

  # Indirect effects via practical issues
  indirect_pi := d * pi

  # Total effects
  total_tf_pos := c1 + (a1 * b)
  total_tf_mod := c2 + (a2 * b)
  total_sp_pos := c3 + (a3 * b)
  total_sp_pos_low_auto := c4 + (a4 * b)
  total_pi := pi + (d * b)
'

# Fit the SEM model with WLSMV estimator
fit <- sem(sem_model, 
           data = zddf, 
           ordered = c("zd", "motiv"),  # Binary variables declared as ordered
           estimator = "WLSMV")

# Summary of the SEM model
summary(fit, fit.measures = TRUE, standardize = TRUE, rsquare = TRUE)

# Plot the SEM model
semPaths(fit, "std", what = "paths", layout = "tree", edge.label.cex = 1.2)




