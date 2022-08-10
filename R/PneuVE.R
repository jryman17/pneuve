#' @title PneuVE
#' 
#' @description This package predicts pneumococcal conjugate vaccine (PCV) effectiveness for shared serotypes with current PCV using summary-level (mean and 95% CI) current PCV serotype specific effectiveness data, summary-level (geometric mean and 95% CI) current PCV serotype specific IgG data, and summary-level (geometric mean and 95% CI) new PCV serotype specific IgG data.
#' 
#' @param nsim the number of simulations
#' @param seed seed=realization
#' @param obs.VEs observed vaccine effectiveness (starting with the mean, followed by the lower, and then the upper 95% CI)
#' @param nSubPlacebo number of subjects in the placebo treated arm
#' @param nSubCurrentPCV number of subjects in the current PCV treated arm
#' @param nSubNewPCV number of subjects in the next PCV treated arm
#' @param serotype serotype being analyzed (in "")
#' @param GMC.Placebo geometric mean concentration in the placebo treated arm
#' @param upperCI.Placebo upper 95% CI bound in the placebo treated arm
#' @param GMC.CurrentPCV geometric mean concentration in the current PCV treated arm
#' @param upperCI.CurrentPCV upper 95% CI bound in the current PCV treated arm
#' @param GMC.NewPCV geometric mean concentration in the next PCV treated arm
#' @param upperCI.NewPCV upper 95% CI bound in the new PCV treated arm
#' 
#' @return NULL
#' 
#' @examples 
#' ser4sim <- VEpred(nsim=10,seed=realization,obs.VEs=c(97, 65, 100),nSubPlacebo=189,nSubCurrentPCV=123,nSubNewPCV=235,serotype="pn4",GMC.Placebo=0.03,upperCI.Placebo=0.04,GMC.CurrentPCV=1.61,upperCI.CurrentPCV=1.84,GMC.NewPCV=1.55,upperCI.NewPCV=1.70)
#' @export 


VEpred <- function(nsim,seed,obs.VEs,nSubPlacebo,nSubCurrentPCV,nSubNewPCV,serotype,
                   GMC.Placebo,upperCI.Placebo,GMC.CurrentPCV,upperCI.CurrentPCV,
                   GMC.NewPCV,upperCI.NewPCV) {
  
  #Setting up Summary level DFs
  placeboDF <- data.frame(serotype,GMC.Placebo,nSubPlacebo,upperCI.Placebo) %>% 
    mutate(logMean=log(GMC.Placebo),
           logCI95=log(upperCI.Placebo),
           logSEM=(logCI95-logMean)/1.96, 
           logSD=logSEM*sqrt(nSubPlacebo),
           logVAR=logSD^2) 
  
  currentPCV.DF <- data.frame(serotype,GMC.CurrentPCV,nSubCurrentPCV,upperCI.CurrentPCV) %>% 
    mutate(logMean=log(GMC.CurrentPCV),
           logCI95=log(upperCI.CurrentPCV),
           logSEM=(logCI95-logMean)/1.96, 
           logSD=logSEM*sqrt(nSubCurrentPCV),
           logVAR=logSD^2)
  
  newPCV.DF <- data.frame(serotype,GMC.NewPCV,nSubNewPCV,upperCI.NewPCV) %>% 
    mutate(logMean=log(GMC.NewPCV),
           logCI95=log(upperCI.NewPCV),
           logSEM=(logCI95-logMean)/1.96, 
           logSD=logSEM*sqrt(nSubNewPCV),
           logVAR=logSD^2) 
  
  
  logMean_Generator <- function(x,seed,vaccine,numberOfSubjects,pcvAggregate,serotype) {
    set.seed(seed)
    
    n <- 1 
    
    Mlog <- as.numeric(pcvAggregate %>% filter(serotype==serotype) %>% 
                         dplyr::select(logMean))
    
    SElog <- as.numeric(pcvAggregate %>% filter(serotype==serotype) %>% 
                          dplyr::select(logSEM)) #use SEM instead
    
    
    samples <- rnorm(n,Mlog,SElog)
    
    sample.avgs <- data.frame(logMean=samples)
    
    
    
    return(sample.avgs)
    
  }
  
  
  
  #function to simulate new SDs
  logSD_Generator <- function(x,seed,vaccine,numberOfSubjects,pcvAggregate,serotype) {
    set.seed(seed)
    
    n <- as.numeric(numberOfSubjects)
    
    logVar <- as.numeric(pcvAggregate %>% filter(serotype==serotype) %>%
                           dplyr::select(logVAR))
    
    
    iter.logSD <- data.frame(logSD=sqrt((logVar*rchisq(n=1, df=n-1))/(n-1))) 
    
    
    
    return(iter.logSD)
  }
  
  
  #VE Generator 
  VE_dist_func <- function(x, seed, obs.VEs , nsim=1) {
    set.seed(seed)
    
    ve.all <- c(obs.VEs[1], obs.VEs[2], obs.VEs[3]) #[1]=mean VE, [2]=lower 95%CI VE, [3]=upper 95%CI VE
    
    ve.all[ve.all == 100] <- 99.99 # Prevents an infinite log_OR
    
    log_OR <- log(1 - ve.all / 100)
    
    names(log_OR) <- c('mean', 'ucl', 'lcl')
    
    log_or_sd <-  (log_OR['ucl'] - log_OR['mean'])/1.96
    
    sim.log.or <- rnorm(n = nsim, mean = log_OR['mean'], sd = log_or_sd)
    
    
    
    sim.VE <- (1 - exp(sim.log.or))
    
    return(sim.VE)
  }
  
  
  #function to create a placebo, PCV13 and PCV7 titer dataset
  SimData <- function(x,seed,pcvAggregate,serotype,logMean,logSD,numberOfSubjects) { 
    set.seed(seed)
    
    
    
    
    
    logSimData <- rnorm(n = numberOfSubjects, 
                        mean = as.numeric(logMean),
                        sd = as.numeric(logSD)) 
    simTiter <- data.frame(ELISA=exp(logSimData))
    
    
    return(simTiter)
  }
  
  
  
  #calculate ELISA RCDC
  ELISARCDC <- function(x,seroDataset,vaccine) { 
    #Sequence for rcdc to evaluate - smoothes the RCDC 
    z = seq(0, 10, by=0.01)
    
    PCVecdf <- ecdf(seroDataset) #ecdf is a cumulative distribution function
    GoBetween = PCVecdf(z) #runs the ecdf function at the specified intervals (z)
    PCVrcdc<-data.frame(z,GoBetween) %>% 
      mutate(percent=100*(1-GoBetween)) %>% # becasue I want the reverse cumulative distribution I use 1 - the ecdf output
      dplyr::rename(ELISA=z) %>% 
      mutate(vax=vaccine)
    
    return(PCVrcdc)
    
  }
  
  
  #This function is looking for the fraction of interest for each assay (pv/pc)
  CoPSolver <- function(realization,ELISAplaceborcdc,ELISAPCV7rcdc,ELISAmedianVE) { 
    #ELISA normalized COPs
    ELISAmedianPnVE=1-ELISAmedianVE #median fraction of interest
    
    #joining The PCV13 and placebo rcdc's and then calculating (pv/pc)
    ELISAcalc <- right_join(ELISAplaceborcdc,ELISAPCV7rcdc, by="ELISA") %>% 
      mutate(percentLessThenUnvax=100-percent.x,
             percentLessThenVax=100-percent.y,
             ratioVaxUnvax=percentLessThenVax/percentLessThenUnvax) %>% 
      filter(!is.na(ratioVaxUnvax), !ratioVaxUnvax=="Inf") %>% 
      arrange(desc(ELISA)) #evaluates the largest titer associated with a ratioVaxUnvax - this is done because there can be multiple correct ratios associated with an ELISA range - I wanted to be conservative and chose the highest ELISA threshold that is associated with the correct ratio
    #finds the row that contains the fraction closest to the fraction of interest 
    ELISAmedianCOPvector <- ELISAcalc[which.min(abs(ELISAmedianPnVE-ELISAcalc$ratioVaxUnvax)),]
    
    #Pulls out the titer associated with the fraction of interest
    ELISAmedianCOP = as.numeric(as.character(ELISAmedianCOPvector$ELISA))
    
    #Put together all of the COPs for Heatmap
    COPs <- data.frame(ELISAmedianCOP)
    
    return(COPs)
  }
  
  
  #This function is looking at the whole ELISAcalc dataset - I used this to QC the solution from COPsolver function above
  CoPSolverDataset <- function(realization,ELISAplaceborcdc,ELISAPCV7rcdc,ELISAmedianVE) { 
    #ELISA normalized COPs
    ELISAmedianPnVE=1-ELISAmedianVE #median fraction of interest
    
    #joining The PCV13 and placebo rcdc's and then calculating (pv/pc)
    ELISAcalc <- right_join(ELISAplaceborcdc,ELISAPCV7rcdc, by="ELISA") %>% 
      mutate(percentLessThenUnvax=100-percent.x,
             percentLessThenVax=100-percent.y,
             ratioVaxUnvax=percentLessThenVax/percentLessThenUnvax) %>% 
      filter(!is.na(ratioVaxUnvax), !ratioVaxUnvax=="Inf") %>% 
      arrange(desc(ELISA)) #evaluates the largest titer associated with a ratioVaxUnvax
    
    return(ELISAcalc)
  }
  
  
  
  #Function that allows you to look at the whole RCDC dataset and ratios used to calculate Efficacy 
  #after we have CoP
  #this function is used within the CoPVEtable function below - the newVEvector contains the ratioVaxUnvax which is used to back calculate vaccine efficacy
  CoPdatasetVE <- function(placeborcdc,rcdc,COP,tolAdjust) {
    calc <- merge(placeborcdc,rcdc, by="ELISA", all = T) %>% 
      mutate(percentLessThenUnvax=100-percent.x,
             percentLessThenVax=100-percent.y,
             ratioVaxUnvax=percentLessThenVax/percentLessThenUnvax) %>% 
      filter(!is.na(ratioVaxUnvax), !ratioVaxUnvax=="Inf") %>% 
      arrange(desc(ELISA)) #evaluates the largest titer associated with a ratioVaxUnvax
    
    newVEvector <- data.frame(calc %>% filter(near(ELISA,COP,tol = tolAdjust))) #adjustment of tolerance may be needed for a single solution, but really has not needed to be adjusted from 0.005 in my experience
    
    return(newVEvector)
  }
  
  
  #function for the Efficacy table for each vaccine
  #function leverages CoPdatasetVE function above
  CoPVEtable <- function(realization,PCV7ELISArcdc,placeboELISArcdc,PCV13ELISArcdc,COPE) {
    #PCV13
    #ELISA
    PCV13ELISAefficacyV <- CoPdatasetVE(placeborcdc = placeboELISArcdc,rcdc = PCV13ELISArcdc,COP = COPE,
                                        tolAdjust = 0.005) %>% 
      mutate(efficacy=(1-ratioVaxUnvax)*100) #This is Sibers Equation VE=1-(pv/pc)
    PCV13ELISAefficacy <- as.numeric(PCV13ELISAefficacyV$efficacy) #pulls out the efficacy from dataset
    
    #PCV7
    #ELISA
    PCV7ELISAefficacyV <- CoPdatasetVE(placeborcdc = placeboELISArcdc,rcdc = PCV7ELISArcdc,
                                       COP = COPE,tolAdjust = 0.005) %>% 
      mutate(efficacy=(1-ratioVaxUnvax)*100) #This is Sibers Equation VE=1-(pv/pc)
    PCV7ELISAefficacy <- as.numeric(PCV7ELISAefficacyV$efficacy)
    
    
    #Table Output with vaccine specific efficacy associated with each assay 
    TableEffSero <- data.frame(bind_rows(vaccine=c("PCV7","PCV13"),
                                         ELISAefficacy=c(PCV7ELISAefficacy,PCV13ELISAefficacy))) 
    
    return(TableEffSero)
  }
  
  ###Serotype Summary####
  sersim <- as_tibble(expand.grid(realization=seq(1,nsim))) %>% 
    rowwise() %>% 
    mutate(vaccineEfficacy=purrr::map(realization,.f=VE_dist_func,seed=realization,
                                      obs.VEs=obs.VEs,nsim=1)) %>% 
    mutate(VE=as.numeric(vaccineEfficacy)) %>% 
    rowwise() %>% 
    mutate(placebologMean=purrr::map(realization,.f=logMean_Generator,seed=realization,vaccine="placebo",
                                     numberOfSubjects=nSubPlacebo,
                                     serotype=serotype,pcvAggregate=placeboDF)) %>% 
    mutate(placebologSD=purrr::map(realization,.f=logSD_Generator,seed=realization,vaccine="placebo",
                                   numberOfSubjects=nSubPlacebo,
                                   serotype=serotype,pcvAggregate=placeboDF)) %>% 
    mutate(pcv7logMean=purrr::map(realization,.f=logMean_Generator,seed=realization,vaccine="PCV7",
                                  numberOfSubjects=nSubCurrentPCV,
                                  serotype=serotype,pcvAggregate=currentPCV.DF)) %>% 
    mutate(pcv7logSD=purrr::map(realization,.f=logSD_Generator,seed=realization,vaccine="PCV7",
                                numberOfSubjects=nSubCurrentPCV,
                                serotype=serotype,pcvAggregate=currentPCV.DF)) %>% 
    mutate(pcv13logMean=purrr::map(realization,.f=logMean_Generator,seed=realization,vaccine="PCV13",
                                   numberOfSubjects=nSubNewPCV,
                                   serotype=serotype,pcvAggregate=newPCV.DF)) %>% 
    mutate(pcv13logSD=purrr::map(realization,.f=logSD_Generator,seed=realization,vaccine="PCV13",
                                 numberOfSubjects=nSubNewPCV,
                                 serotype=serotype,pcvAggregate=newPCV.DF)) %>% 
    rowwise() %>% 
    mutate(placebologMeanNum=list(as.data.frame(placebologMean)$logMean)) %>% 
    mutate(placebologSDNum=list(as.data.frame(placebologSD)$logSD)) %>% 
    mutate(pcv7logMeanNum=list(as.data.frame(pcv7logMean)$logMean)) %>% 
    mutate(pcv7logSDNum=list(as.data.frame(pcv7logSD)$logSD)) %>% 
    mutate(pcv13logMeanNum=list(as.data.frame(pcv13logMean)$logMean)) %>% 
    mutate(pcv13logSDNum=list(as.data.frame(pcv13logSD)$logSD)) %>% 
    rowwise() %>% 
    mutate(placeboVector=purrr::map(realization,.f=SimData,seed=realization,pcvAggregate=placeboDF,
                                    serotype=serotype,logMean=placebologMeanNum,
                                    logSD=placebologSDNum,numberOfSubjects=nSubPlacebo)) %>% 
    mutate(pcv7vector=purrr::map(realization,.f=SimData,seed=realization,pcvAggregate=currentPCV.DF,
                                 serotype=serotype,logMean=pcv7logMeanNum,logSD=pcv7logSDNum, 
                                 numberOfSubjects=nSubCurrentPCV)) %>% 
    mutate(pcv13vector=purrr::map(realization,.f=SimData,seed=realization,pcvAggregate=newPCV.DF,
                                  serotype=serotype,logMean=pcv13logMeanNum,logSD=pcv13logSDNum, 
                                  numberOfSubjects=nSubNewPCV)) %>% 
    rowwise() %>% 
    mutate(placeboTiters=list(as.data.frame(placeboVector)$ELISA),
           pcv7Titers=list(as.data.frame(pcv7vector)$ELISA),
           pcv13Titers=list(as.data.frame(pcv13vector)$ELISA)) %>% 
    ungroup() %>% 
    mutate(placeboELISArcdc=purrr::pmap(.l=list(realization,placeboTiters),.f=ELISARCDC,vaccine="placebo")) %>% 
    mutate(pcv7ELISArcdc=purrr::pmap(.l=list(realization,pcv7Titers),.f=ELISARCDC,vaccine="PCV7")) %>% 
    mutate(pcv13ELISArcdc=purrr::pmap(.l=list(realization,pcv13Titers),.f=ELISARCDC,vaccine="PCV13")) %>% 
    
    
    mutate(COPsolve=purrr::pmap(.l=list(realization,placeboELISArcdc,pcv7ELISArcdc,VE),
                                .f=CoPSolver)) %>% 
    rowwise() %>% 
    mutate(COP=COPsolve$"ELISAmedianCOP") %>% 
    ungroup() %>% 
    mutate(Eff=purrr::pmap(.l=list(realization,pcv7ELISArcdc,placeboELISArcdc,pcv13ELISArcdc,COP),
                           .f=CoPVEtable)) %>% 
    rowwise() %>% 
    mutate(pcv7Eff=Eff$'ELISAefficacy'[[1]],
           pcv13Eff=Eff$'ELISAefficacy'[[2]])
  
  
  summarizeEff <- sersim %>% 
    dplyr::select(pcv13Eff) %>% 
    group_by() %>% 
    summarise(medianEff=median(pcv13Eff),
              upper95CiEff=quantile(pcv13Eff, 0.975),
              lower95CiEff=quantile(pcv13Eff, 0.025)) 
  
  
  sumarizeCOP <- sersim %>% 
    dplyr::select(COP) %>% 
    group_by() %>% 
    summarise(medianCOP=median(COP),
              upper95CiCOP=quantile(COP, 0.975),
              lower95CiCOP=quantile(COP, 0.025))
  
  
  
  return(list(summarizeEff,sumarizeCOP,sersim))
  
  
}