#' @title backboinetr
#'
#' @description Obtain the operating characteristics of the backfill bayesian
#' optimal interval design using efficacy and toxicity outcomes for dose
#' optimization within random scenarios
#'
#' @param target_T Target toxicity probability. The default value is
#' \code{target_T=0.3}. When observing 1 DLT out of 3 patients and the target
#' DLT rate is between 0.25 and 0.279, the decision is to stay at the current
#' dose due to a widely accepted practice.
#' @param target_Tr The upper boundary for the toxicity when generating the
#' random scenarios. The default value is \code{target_Tr=0.359}.
#' @param target_E The minimum required efficacy probability. The default value
#' is \code{target_E=0.25}.
#' @param target_Er The lower boundary for the efficacy when generating the
#' random scenarios. The default value is \code{target_Er=0.197}.
#' @param n.dose Number of dose.
#' @param startdose Starting dose. The lowest dose is generally recommended.
#' @param ncohort Number of cohort.
#' @param cohortsize Cohort size.
#' @param pT.saf Highest toxicity probability that is deemed sub-therapeutic
#' such that dose-escalation should be pursued. The default value is
#' \code{pT.saf=target_T*0.6}.
#' @param pT.tox Lowest toxicity probability that is deemed overly toxic such
#' that dose de-escalation is needed. The default value is
#' \code{pT.tox=target_T*1.4}.
#' @param pE.saf Minimum probability deemed efficacious such that the dose
#' levels with less than delta1 are considered sub-therapeutic.
#' The default value is \code{pE.saf=target_E*0.6}.
#' @param alpha.T1 Probability that toxicity event occurs in the late half of
#' toxicity assessment window. The default value is \code{alpha.T1=0.5}.
#' @param alpha.E1 Probability that efficacy event occurs in the late half of
#' assessment window. The default value is \code{alpha.E1=0.5}.
#' @param tau.T Toxicity assessment windows (months).
#' @param tau.E Efficacy assessment windows (months).
#' @param te.corr Correlation between toxicity and efficacy probability,
#' specified as Gaussian copula parameter. The default value is
#' \code{te.corr=0.2}.
#' @param gen.event.time Method to generate the time to first toxicity and
#' efficacy outcome. Weibull distribution is used when
#' \code{gen.event.time ="weibull"}. Uniform distribution is used when
#' \code{gen.event.time="uniform"}.
#' The default value is \code{gen.event.time="weibull"}.
#' @param accrual Accrual rate (months) (patient accrual rate per month).
#' @param gen.enroll.time Method to generate enrollment time. Uniform
#' distribution is used when
#' \code{gen.enroll.time="uniform"}. Exponential distribution is used when
#' \code{gen.enroll.time="exponential"}. The default
#' value is \code{gen.enroll.time="uniform"}.
#' @param n.elimination a minimum sample size for dose elimination. If the number
#' of patients treated at the current dose reaches \code{n.elimination} and meet
#' elimination dose level criteria, eliminate current dose level and higher doses
#' when meet toxicity criteria and eliminate current dose level when meet efficacy
#' criteria. The default value is \code{n.elimination=6}.
#' @param stopping.npts Early study termination criteria for the number of
#' patients in the dose-escalation and backfill cohorts. If the number of
#' patients at the current dose reaches this criteria and the same dose level
#' is recommended as the next dose level, the study is terminated.
#' The default value is \code{stopping.npts=12}.
#' @param suspend the suspension rule that holds off the decision on dose
#' allocation for the dose-escalation cohort until sufficient toxicity
#' information is available. For example, setting as 0.33 which means one-third
#' of the patients had not completed the toxicity evaluation at the current dose
#' level in the dose escalation cohort. The default value \code{suspend=0}
#' essentially turns off this type of suspending rule, that is all patients
#' should complete the toxicity evaluation at the current dose level in the dose
#' escalation cohort
#' @param stopping.prob.T Early study termination criteria for toxicity,
#' taking a value between 0 and 1. If the posterior probability that toxicity
#' outcome is less than the target toxicity probability (\code{target_T}) is
#' larger than this criteria, the dose levels are eliminated from the study.
#' The default value is \code{stopping.prob.T=0.95}.
#' @param stopping.prob.E Early study termination criteria for efficacy,
#' taking a value between 0 and 1. If the posterior probability that efficacy
#' outcome is less than the minimum efficacy probability (\code{target_E}) is
#' larger than this criteria, the dose levels are eliminated from the study.
#' The default value is \code{stopping.prob.E=0.90}.
#' @param ppsi01 Score for toxicity=yes and efficacy=no in utility defined by
#' scoring.The default value is \code{psi01=0}.
#' @param ppsi00 Score for toxicity=no and efficacy=no in utility defined by
#' scoring. The default value is \code{psi00=40}.
#' @param ppsi11 Score for toxicity=yes and efficacy=yes in utility defined by
#' scoring. The default value is \code{psi11=60}.
#' @param ppsi10 Score for toxicity=no and efficacy=yes in utility defined by
#' scoring. The default value is \code{psi10=100}.
#' @param n.sim Number of simulated trial. The default value is
#' \code{n.sim=10000}.
#' @param seed.sim Seed for random number generator. The default value is
#' \code{seed.sim=30}.
#' @details The \code{backboinetr} is a function which generates the operating
#' characteristics of the backfill bayesian optimal interval design using
#' efficacy and toxicity outcomes for dose optimization by a simulation study.
#' Users can specify a variety of study settings to simulate studies. The
#' operating characteristics of the design are summarized by the percentage of
#' times that each dose level was selected as optimal biological dose and the
#' average number of patients who were treated at each dose level. The
#' percentage of times that the study was terminated and the expected study
#' duration are also provided.
#' @return
#' The \code{backboinetr} returns a list containing the following components:
#' \item{toxprob}{The random true toxicity probability.}
#' \item{effprob}{The random true efficacy probability.}
#' \item{phi}{Target toxicity probability.}
#' \item{delta}{Target efficacy probability.}
#' \item{target_Tr}{The upper boundary for the toxicity when generating the
#' random scenarios.}
#' \item{target_Er}{The lower boundary for the efficacy when generating the
#' random scenarios.}
#' \item{bd.true}{The target optimal dose (OD) level when generating the random
#' scenarios.}
#' \item{mtd.true}{The maximum tolerated dose (MTD) level when generating the
#' random scenarios.}
#' \item{lambda1}{Lower toxicity boundary in dose escalation/de-escalation.}
#' \item{lambda2}{Upper toxicity boundary in dose escalation/de-escalation.}
#' \item{eta1}{Lower efficacy boundary in dose escalation/de-escalation.}
#' \item{tau.T}{Toxicity assessment windows (months).}
#' \item{tau.E}{Efficacy assessment windows (months).}
#' \item{suspend}{The suspension rule that holds off the decision on dose
#' allocation for the dose-escalation cohort until sufficient toxicity
#' information is available.}
#' \item{accrual}{Accrual rate (months) (patient accrual rate per month).}
#' \item{n.patient}{Average number of patients who were treated at each dose
#' level in dose-esclation and backfill cohorts}
#' \item{n.bpatient}{Average number of back filled patients who were treated
#' at each dose level}
#' \item{n.tox.patient}{Average number of patients who experienced toxicity at
#' each dose level in dose-esclation and backfill cohorts}
#' \item{n.eff.patient}{Average number of patients who experienced efficacy at
#' each dose level in dose-esclation and backfill cohorts}
#' \item{n.tox.bpatient}{Average number of patients who experienced toxicity at
#' each dose level in backfill cohort}
#' \item{n.eff.bpatient}{Average number of patients who experienced efficacy at
#' each dose level in backfill cohort}
#' \item{prop.select}{Percentage of times that each dose level was selected as
#' optimal biological dose.}
#' \item{prop.stop}{Percentage of times that the study was terminated.}
#' \item{duration}{Expected study duration (months)}
#' \item{totaln}{Total patients}
#' \item{data.obs.n}{Record the number of patients in each dose level within the
#' simulations during the trial}
#' \item{obd}{Record the optimal dose in each simulation during the trial}
#' \item{backfilltimes}{Record how may times we back-filled during the trial}
#' \item{backfillcount}{Record the number of back-filled patients in dose level
#' within the simulations during the trial}
#' \item{PCS}{The percentage of trials that the optimal dose was correctly
#' selected.}
#' \item{PCA}{The percentage of patients that were correctly allocated to the
#' optimal dose.}
#' \item{PTS}{The percentage of toxic doses selection.}
#' \item{PTA}{The percentage of patients who were allocated to toxic doses.}
#' @references
#' Takeda, K., Zhu, J. and Hirakawa, A. (2025), BF-BOIN-ET: A Backfill Bayesian
#' Optimal Interval Design Using Efficacy and Toxicity Outcomes for Dose
#' Optimization. Pharmaceutical Statistics, 24: e2470. https://doi.org/10.1002/pst.2470
#' @examples
#'
#' target_T=0.3
#' target_E=0.25
#'\dontrun{
#' get.oc.backboinetr(target_T=target_T,target_Tr=0.359,target_E=target_E,
#' target_Er=0.197,n.dose=5,startdose=1,ncohort=10,cohortsize=3,
#' pT.saf=0.6 * target_T,pT.tox = 1.4 * target_T,pE.saf = 0.6 * target_E,
#' alpha.T1=0.5,alpha.E1=0.5,tau.T=1,tau.E=1,te.corr=0.2,
#' gen.event.time="weibull",accrual=3,gen.enroll.time="uniform",n.elimination=6,
#' stopping.npts=12,suspend=0,stopping.prob.T=0.95,stopping.prob.E=0.90,
#' ppsi01=0,ppsi00=40,ppsi11=60,ppsi10=100,n.sim=2,seed.sim=30)
#'}
#' @import Iso copula dplyr tidyselect magrittr
#' @importFrom stats binomial dbinom pbeta pbinom rmultinom runif rexp rbeta
#' @export


get.oc.backboinetr <- function (target_T=0.3,target_Tr=0.359,target_E=0.25,target_Er=0.197,n.dose,startdose,ncohort,cohortsize,
                               pT.saf=0.6 * target_T,pT.tox = 1.4 * target_T,pE.saf = 0.6 * target_E,
                               alpha.T1=0.5,alpha.E1=0.5,tau.T,tau.E,te.corr=0.2,gen.event.time="weibull",
                               accrual,gen.enroll.time="uniform",n.elimination=6,stopping.npts=12,
                               suspend=0,stopping.prob.T=0.95,stopping.prob.E=0.90,ppsi01=0,ppsi00=40,ppsi11=60,
                               ppsi10=100,n.sim=10000,seed.sim=30){


  if(!((pT.saf<target_T)&(target_T<pT.tox))){
    stop("Design parameters must satisfy a condition of pT.saf < target_T < pT.tox.")

  }else if(!(pE.saf<target_E)){
    stop("Design parameters must satisfy a condition of pE.saf < target_E.")

  }

  #####Get utility score by measuring the efficacy-toxicity trade off####;
  futility_score <- function(pe,pt,psi00,psi01,psi10,psi11)
  {
    psi.e0t1 <- psi01
    psi.e0t0 <- psi00
    psi.e1t1 <- psi11
    psi.e1t0 <- psi10

    ut = (  psi.e0t1*(1-pe)*pt
            + psi.e0t0*(1-pe)*(1-pt)
            + psi.e1t1*pe    *pt
            + psi.e1t0*pe    *(1-pt))

    return(ut)
  }

  #####Get optimal dose by measuring the efficacy-toxicity trade off####;
  obd.select <- function(probt, probe,tterm, eterm, ph,stopT, stopE,psi00, psi01,psi10,psi11)
  {

    candose <- which((tterm>=(1-stopT))&(eterm>=(1-stopE)))

    if( (length(probt[candose])>0) & (length(candose)>0)  ){
      mdif <- min(abs(probt[candose]-ph))
      mtd1  <- max(which(abs(probt-ph)==mdif))

      effdose <- intersect(candose,1:mtd1)
    }else{
      effdose <- 0
    }
    fu <- futility_score(pt=probt,pe=probe,
                         psi00=psi00,psi01=psi01,psi10=psi10,psi11=psi11)

    if(length(fu[effdose])>0){
      fu.max <- which(fu==max(fu[effdose]))
      if(length(intersect(fu.max,effdose))>0){
        re     <- min(intersect(fu.max,effdose))
      }else{
        re     <- 0
      }
    }else{
      re       <-0
    }


    return(re)
  }

  ####calculate toxicity escalation and de-escalation boundaries#####;
  lambda1 = log((1 - pT.saf)/(1 - target_T))/log(target_T *
                                                   (1 - pT.saf)/(pT.saf * (1 - target_T)))

  lambda2 = log((1 - target_T)/(1 - pT.tox))/log(pT.tox * (1 -
                                                             target_T)/(target_T * (1 - pT.tox)))

  ####calculate efficacy boundaries#####;

  eta1 = log((1-pE.saf)/(1-target_E))/log(target_E*(1-pE.saf)/(pE.saf* (1-target_E)))

######generate random scenarios#####;
  ndose<-n.dose # number of doses
  targetT<-target_T # highest acceptable toxicity rate
  targetE<-target_Er # lowest acceptable efficacy rate
  u1<-ppsi11
  u2<-ppsi00
  ntrial<-n.sim
  obd.vec<-mtd.vec<-OBD.vec<-c()
  pE.mat<-pT.mat<-u.mat<-c()
  obdi<-0

  set.seed(seed.sim)

  for (iii in startdose:ndose){
    obdi <- iii

    for(s in 1:ntrial){

      obd<-obdi
      mtd<-obd # initialization
      prob_plateau<-0.5
      obd.temp<-0
      jj<-kk<-rep(0,ndose)
      uu<-runif(1)

      if(uu<prob_plateau & obd<ndose){
        # generate plateau cases
        while(obd.temp!=obd | kk[obd]>target_Tr | jj[obd]<targetE+0.05){
          mtd<-(obd-1)+sample(ndose-obd+1,1)
          D=1:ndose
          if(mtd<ndose){temp<-sort(runif(length(D)-mtd,targetT,1));
          bornesup<-max(targetT+0.1,temp[length(D)-mtd])}
          if(mtd==ndose){bornesup<-targetT +(1-targetT)*rbeta(1,0.5,1)}
          mtd.temp<-0

          while(mtd.temp!=mtd){
            kk<-sort(runif(ndose,0,bornesup))
            mtd.temp<-which.min(abs(kk-targetT))
          }

          bornesup<-runif(1,targetE+0.1,1)
          jj[1:obd]<-sort(runif(obd,0,bornesup))
          jj[obd:ndose]<-jj[obd]

          utility<-u1*jj+(1-kk)*u2

          obd.temp<-which.max(utility[1:mtd])

        }
      }


      if(uu>=prob_plateau & obd<ndose){
        # generate non-plateau cases
        while(obd.temp!=obd | kk[obd]>target_Tr | jj[obd]<targetE+0.05){
          mtd<-(obd-1)+sample(ndose-obd+1,1)
          D=1:ndose
          if(mtd<ndose){temp<-sort(runif(length(D)-mtd,targetT,1));
          bornesup<-max(targetT+0.1,temp[length(D)-mtd])}
          if(mtd==ndose){bornesup<-targetT +(1-targetT)*rbeta(1,0.5,1)}
          mtd.temp<-0

          while(mtd.temp!=mtd ){
            kk<-sort(runif(ndose,0,bornesup))
            mtd.temp<-which.min(abs(kk-targetT))
          }

          med<-(obd-1)+sample(ndose-obd+1,1)
          bornesup<-runif(1,targetE+0.1,1)
          jj[med]<-max(runif(ndose,0,bornesup))
          if(med>1){
            jj[1:(med-1)]<-sort(runif(med-1,0,jj[med]))
          }
          if(med<ndose){
            jj[(med+1):ndose]<- sort(runif(ndose-med,0,jj[med]),decreasing=TRUE)
          }
          utility<-u1*jj+(1-kk)*u2
          obd.temp<-which.max(utility[1:mtd])

        }
      }

      if(obd==ndose){
        while(obd.temp!=obd | kk[obd]>target_Tr | jj[obd]<targetE+0.05){
          mtd<-(obd-1)+sample(ndose-obd+1,1)
          D=1:ndose
          if(mtd<ndose){temp<-sort(runif(length(D)-mtd,targetT,1));
          bornesup<-max(targetT+0.1,temp[length(D)-mtd])}
          if(mtd==ndose){bornesup<-targetT +(1-targetT)*rbeta(1,0.5,1)}
          mtd.temp<-0

          while(mtd.temp!=mtd){
            kk<-sort(runif(ndose,0,bornesup))
            mtd.temp<-which.min(abs(kk-targetT))
          }
          bornesup<-runif(1,targetE+0.1,1)
          jj<-sort(runif(ndose,0,bornesup))
          utility<-u1*jj+(1-kk)*u2
          obd.temp<-which.max(utility[1:mtd])
        }
      }

      pE.mat<-rbind(pE.mat,jj)
      pT.mat<-rbind(pT.mat,kk)
      u.mat<-rbind(u.mat,utility)
      obd.vec<-c(obd.vec,obd)
      mtd.vec<-c(mtd.vec,mtd)

    }
    OBD.vec<-c(OBD.vec,rep(iii,ntrial))
    all<-list(OBD=OBD.vec,obd.true=obd.vec,mtd.true=mtd.vec,Pe=pE.mat,Pt=pT.mat)
    alld<-as.data.frame(all)
  }
  eff<-alld %>% select(starts_with("Pe."))
  tox<-alld %>% select(starts_with("Pt."))


####end#####;
  an.sim<-n.sim*(n.dose-startdose+1)

  data.obs.n <- array(0,dim=c(an.sim,n.dose))
  toxprobr   <- array(0,dim=c(an.sim,n.dose))
  effprobr   <- array(0,dim=c(an.sim,n.dose))
  data.dur   <- array(0,dim=c(an.sim))

  obd <- array(0,dim=c(an.sim))
  mtd.true <- array(0,dim=c(an.sim))
  bd <- array(0,dim=c(an.sim))

  toxicity=matrix(nrow=an.sim,ncol=n.dose)
  efficacy=matrix(nrow=an.sim,ncol=n.dose)

  btoxicity=matrix(nrow=an.sim,ncol=n.dose)
  befficacy=matrix(nrow=an.sim,ncol=n.dose)


  prop.select <- array(0,dim=c(an.sim,n.dose))
  prop.stop <- array(0,dim=c(an.sim))

  dose.curr=startdose

  backfilltimes=rep(0,an.sim) ## record how may times we back-filled during the trial
  backfillcount=matrix(nrow=an.sim,ncol=n.dose) ## record the location of backfill

  nmax=ncohort*cohortsize  ####max number of patients for dose-escalation#####;

  dosen <- 1:n.dose
  dose  <- paste("Dose",dosen,sep="")
  pr.alpha <- 1
  pr.beta  <- 1

  alpha.T1 <- alpha.T1
  alpha.T2 <- 0.5
  alpha.E1 <- alpha.E1
  alpha.E2 <- 0.5


  for(simu in 1:an.sim){
    set.seed(seed.sim)
    effprob = eff[simu,]
    toxprob = tox[simu,]
    bd[simu]=alld[simu,2]
    mtd.true[simu]=alld[simu,3]

    toxp <- data.frame(toxprob)
    colnames(toxp) <- dose

    effp <- data.frame(effprob)
    colnames(effp) <- dose

    efftoxp <- list(toxp=toxp,effp=effp)

    ncop    <- copula::normalCopula(te.corr,dim=2,dispstr="ex")
    mv.ncop <- NULL

    if(gen.event.time=="weibull"){

      for(i in 1:n.dose){
        psi.T    <- efftoxp$toxp[i][[1]]
        zetta.T1 <- log(log(1-psi.T)/log(1-psi.T+alpha.T1*psi.T))/log(1/(1-alpha.T2))
        zetta.T2 <- tau.T/(-log(1-psi.T))^(1/zetta.T1)

        psi.E    <- efftoxp$effp[i][[1]]
        zetta.E1 <- log(log(1-psi.E)/log(1-psi.E+alpha.E1*psi.E))/log(1/(1-alpha.E2))
        zetta.E2 <- tau.E/(-log(1-psi.E))^(1/zetta.E1)

        mv.ncop <- append(mv.ncop,copula::mvdc(copula       = ncop,
                                               margins      = c("weibull","weibull"),
                                               paramMargins = list(list(shape=zetta.T1,scale=zetta.T2),
                                                                   list(shape=zetta.E1,scale=zetta.E2))))
      }

    }else if(gen.event.time=="uniform"){

      for(i in 1:n.dose){
        psi.T <- efftoxp$toxp[i][[1]]
        psi.E <- efftoxp$effp[i][[1]]

        mv.ncop <- append(mv.ncop,copula::mvdc(copula       = ncop,
                                               margins      = c("unif","unif"),
                                               paramMargins = list(list(min=0,max=tau.T*(1/psi.T)),
                                                                   list(min=0,max=tau.E*(1/psi.E)))))
      }

    }

    #####end of generating event time for efficacy and toxicity######;
    yE=NULL
    yT=NULL
    localt.yT=NULL
    localt.yE=NULL

    obs.n     <- numeric(n.dose)
    obs.tox   <- numeric(n.dose)
    obs.tox.n <- numeric(n.dose)
    obs.eff   <- numeric(n.dose)
    obs.eff.n <- numeric(n.dose)
    pe        <- numeric(n.dose)
    pt        <- numeric(n.dose)

    obs.btox   <- numeric(n.dose)
    obs.beff   <- numeric(n.dose)

    dose.curr <- startdose

    bgamma  <- numeric(n.dose)

    d=NULL ####record dose allocation####
    t.enter=NULL
    localt.enter=NULL
    localt.finish=NULL

    t.eventT=NULL
    t.eventE=NULL
    t.onset=NULL
    t.finish=NULL
    t.curr=0
    tite.df    <- NULL
    nvector=rep(0,n.dose)

    bdosekeep=NULL
    backfill=0 # how many times of backfill.
    backfillvector=rep(0,length=n.dose)  ## which dose is backfilled.
    elimi_tox<-rep(0,n.dose)
    elimi_eff<-rep(0,n.dose)
    elimi<-rep(0,n.dose)
    earlystop<-0
    type=NULL

    for(i in 1:ncohort){
      ### assign three maxpatients to dose.curr
      # t.curr; t.finish;t.onset

      localt.enter=NULL

      for(j in 1:cohortsize){
        ## local t.enter
        if(j==1){
          localt.enter=c(localt.enter,t.curr)

        }else{
          if(gen.enroll.time=="uniform"){ localt.enter =c(localt.enter, localt.enter[length(localt.enter)]+ runif(1, 0, 2/accrual))}
          if(gen.enroll.time=="exponential"){ localt.enter = c(localt.enter,  localt.enter[length(localt.enter)]+ rexp(1, rate=accrual))}

        }
      }

      localt.finish=localt.enter+tau.T

      t.enter=c(t.enter,localt.enter)

      t.finish=c(t.finish,localt.finish)

      time.te <- copula::rMvdc(cohortsize,mv.ncop[[dose.curr]]) ####event time#####

      localt.onset=localt.enter+time.te[,1]

      localt.yT=as.numeric(time.te[,1]<=tau.T)
      localt.yE=as.numeric(time.te[,2]<=tau.E)

      yT=c(yT,as.numeric(time.te[,1]<=tau.T))
      yE=c(yE,as.numeric(time.te[,2]<=tau.E))

      d=c(d,rep(dose.curr,cohortsize))
      type=c(type,rep(1,cohortsize))

      t.eventT=c(t.eventT,time.te[,1])
      t.eventE=c(t.eventE,time.te[,2])

      t.onset=t.enter+apply(data.frame(
        t.eventT=t.eventT,
        t.eventE=t.eventE
      ), 1, max, na.rm=TRUE)

      endt=apply(data.frame(
        endt1=t.enter+tau.T,
        endt2=t.enter+t.eventT
      ), 1, min, na.rm=TRUE)

      #####added for more than or equal 1-suspend patients completed the DLT assessment####;
      if (suspend!=0){
        quantile=quantile(endt[which(d==dose.curr)], probs = (1-suspend))
        quantile50=max(endt[which(d==dose.curr)][endt[which(d==dose.curr)]<=quantile])
      }
      #####end######;

      ende=apply(data.frame(
        ende1=t.enter+tau.E,
        ende2=t.enter+t.eventE
      ), 1, min, na.rm=TRUE)

      nvector[dose.curr]=nvector[dose.curr]+cohortsize


      ## renew t.curr until there is open doses
      flag=0
      t.curr=max(t.enter)
      if (suspend!=0){
        t.bench=min(max(localt.finish),max(localt.onset),quantile50)
      }else{
        t.bench=min(max(localt.finish),max(localt.onset))
      }
      if(i==ncohort){t.bench=t.bench+100}
      ###bacause in the last cohort, we need to enroll enough patients to the backfilling doses, setting as 100 months#######;
      queue=NULL

      # t.enter;t.onset;t.finish

      while(flag==0){

        flag=1

        if(gen.enroll.time=="uniform"){
          t.curr=t.curr+runif(1, 0, 2/accrual)
          queue=c(queue,t.curr)

        }else if(gen.enroll.time=="exponential"){
          t.curr=t.curr+rexp(1, rate=accrual)
          queue=c(queue,t.curr)
        }

        if(t.bench>t.curr){

          flag=0
        }

      }

      queue=queue[which(queue!=t.curr)]

      db=NULL ##record dose allocation for backfilling
      localt.enterb=NULL
      localt.finishb=NULL
      localt.yTb=NULL
      localt.yEb=NULL
      time.teba=NULL
      titeb.df=NULL
      s_back=NULL

      # t.curr;queue
      # d
      ## check open doses
      ## be adviced, we do not eliminate a dose from bf set if it turns of be overtoxic during the trial.
      ## assume response is immediately available after dose assignment.
      ## this simplifies the situation, during backfilling, no lower doses watch response
      ## if not, change the code easily, decide the bfset everytime we queueing
      flag2=0

      while(flag2==0){

        flag2=1

        if(dose.curr!=1){ # dynnamically change the bfset as each person arrives.
          for (q in queue){
            ####select max utility score for backfilling dose####;
            if(!is.null(bdosekeep)){
              if (q!=min(queue)){
                #####using the latest information to get the max utility score#####
                titeb.df <- data.frame(dose   = d,
                                       enter  = t.enter,
                                       endtox = endt,
                                       dlt    = yT,
                                       endeff = ende,
                                       orr    = yE)

                for(ds in 1:n.dose){
                  if(sum(titeb.df$dose==ds)>0){

                    bcompsub.T <- titeb.df[(titeb.df$endtox<=q)&(titeb.df$dose==ds),]
                    bpendsub.T <- titeb.df[(titeb.df$endtox >q)&(titeb.df$dose==ds),]
                    bcompsub.E <- titeb.df[(titeb.df$endeff<=q)&(titeb.df$dose==ds),]
                    bpendsub.E <- titeb.df[(titeb.df$endeff >q)&(titeb.df$dose==ds),]

                    bx.DLT  <- sum(bcompsub.T$dlt)
                    bn.DLT  <- bx.DLT+sum(1-bcompsub.T$dlt)+sum(q-bpendsub.T$enter)/tau.T

                    bx.ORR  <- sum(bcompsub.E$orr)
                    bn.ORR  <- bx.ORR+sum(1-bcompsub.E$orr)+sum(q-bpendsub.E$enter)/tau.E

                    obs.tox[ds]   <- bx.DLT
                    obs.tox.n[ds] <- bn.DLT
                    pt[ds]        <- bx.DLT/bn.DLT

                    obs.eff[ds]   <- bx.ORR
                    obs.eff.n[ds] <- bn.ORR
                    pe[ds]        <- bx.ORR/bn.ORR

                    titeb.curdose <- titeb.df[titeb.df$dose==ds,]
                    bgamma.T <- as.numeric(titeb.curdose$endtox<=q)
                    bgamma.E <- as.numeric(titeb.curdose$endeff<=q)
                    bgamma[ds]<- mean(bgamma.T*bgamma.E)

                  }}
              }


              ###select the backfilling dose using utility score#####;
              ####add utility score#####;
              evadose <- intersect(intersect(dosen[nvector!=0],bdosekeep),dosen[nvector<stopping.npts])

              if(!identical(evadose, integer(0))) {
                estpt <- Iso::pava(pt[evadose])


                utility_b<-futility_score(pe=pe[evadose],pt=estpt,psi00=ppsi00,psi01=ppsi01,psi10=ppsi10,psi11=ppsi11)
                s_back=max(evadose[utility_b==max(utility_b)]) ###backfilling dose#####;
              }else{
                s_back=NULL
              }
            }


            ## check the status of current dos, if it all finishes and reaches n.stop and decision is stay, terminate the trial
            if(!is.null(s_back)){
              s=s_back
              if(nvector[s]<stopping.npts){
                ## note, I do not care if toxicity of this backfill dose is currently over toxic or not.
                db=c(db,s)
                d=c(d,s)
                type=c(type,0)

                localt.enterb=c(localt.enterb,q)
                localt.finishb=c(localt.finishb,q+tau.T)

                t.enter=c(t.enter,q)
                t.finish=c(t.finish, q+tau.T)

                time.teb<-copula::rMvdc(1,mv.ncop[[s]])

                time.teba<-rbind(time.teba,time.teb)

                localt.yTb=c(localt.yTb,as.numeric(time.teb[,1]<=tau.T))
                localt.yEb=c(localt.yEb,as.numeric(time.teb[,2]<=tau.E))

                yT=c(yT,as.numeric(time.teb[,1]<=tau.T))
                yE=c(yE,as.numeric(time.teb[,2]<=tau.E))

                t.eventT=c(t.eventT,time.teb[,1])
                t.eventE=c(t.eventE,time.teb[,2])

                t.onset=t.enter+apply(data.frame(
                  t.eventT=t.eventT,
                  t.eventE=t.eventE
                ), 1, max, na.rm=TRUE)

                endt=apply(data.frame(
                  endt1=t.enter+tau.T,
                  endt2=t.enter+t.eventT
                ), 1, min, na.rm=TRUE)

                ende=apply(data.frame(
                  endeb1=t.enter+tau.E,
                  endeb2=t.enter+t.eventE
                ), 1, min, na.rm=TRUE)

                nvector[s]=nvector[s]+1
                backfillvector[s]=backfillvector[s]+1
                backfill=backfill+1
              }

            }



          } # end of queue
          # queue;t.enter;t.event;t.onset;t.finish;d;yT;yE;nvector; backfillvector; d



        } # dose.curr !=1
      }


      ####decide next dose#####;

      tite.df <- rbind(tite.df,
                       data.frame(dose   = dose.curr,
                                  enter  = localt.enter,
                                  endtox = localt.enter+apply(as.matrix(1:cohortsize),1,function(x){min(time.te[x,1],tau.T)}),
                                  dlt    = localt.yT,
                                  endeff = localt.enter+apply(as.matrix(1:cohortsize),1,function(x){min(time.te[x,2],tau.E)}),
                                  orr    = localt.yE,
                                  backfill = rep(0,cohortsize)))

      if(!is.null(localt.enterb)){
        tite.df <- rbind(tite.df,
                         data.frame(dose   = db,
                                    enter  = localt.enterb,
                                    endtox = localt.enterb+apply(as.matrix(1:length(localt.enterb)),1,function(x){min(time.teba[x,1],tau.T)}),
                                    dlt    = localt.yTb,
                                    endeff = localt.enterb+apply(as.matrix(1:length(localt.enterb)),1,function(x){min(time.teba[x,2],tau.E)}),
                                    orr    = localt.yEb,
                                    backfill = rep(1,length(localt.enterb))))
      }

      tite.curdose <- tite.df[tite.df$dose==dose.curr,]
      gamma.T <- as.numeric(tite.curdose$endtox<=t.curr)
      gamma.E <- as.numeric(tite.curdose$endeff<=t.curr)

      gamma.all.T <- as.numeric(tite.df$endtox<=t.curr)
      gamma.all.E <- as.numeric(tite.df$endeff<=t.curr)

      for(ds in 1:n.dose){
        if(sum(tite.df$dose==ds)>0){

          compsub.T <- tite.df[(tite.df$endtox<=t.curr)&(tite.df$dose==ds),]
          pendsub.T <- tite.df[(tite.df$endtox >t.curr)&(tite.df$dose==ds),]
          compsub.E <- tite.df[(tite.df$endeff<=t.curr)&(tite.df$dose==ds),]
          pendsub.E <- tite.df[(tite.df$endeff >t.curr)&(tite.df$dose==ds),]

          x.DLT  <- sum(compsub.T$dlt)
          n.DLT  <- x.DLT+sum(1-compsub.T$dlt)+sum(t.curr-pendsub.T$enter)/tau.T

          x.ORR  <- sum(compsub.E$orr)
          n.ORR  <- x.ORR+sum(1-compsub.E$orr)+sum(t.curr-pendsub.E$enter)/tau.E

          obs.tox[ds]   <- x.DLT
          obs.tox.n[ds] <- n.DLT
          pt[ds]        <- x.DLT/n.DLT

          obs.eff[ds]   <- x.ORR
          obs.eff.n[ds] <- n.ORR
          pe[ds]        <- x.ORR/n.ORR

          #####create two paramters for -	Average number of patients who experienced toxicity/efficacy at each dose level in backfill cohort####;
          x.bDLT  <- sum(compsub.T[compsub.T$backfill==1,]$dlt)
          x.bORR  <- sum(compsub.E[compsub.E$backfill==1,]$orr)
          obs.btox[ds]   <- x.bDLT
          obs.beff[ds]   <- x.bORR

        }}

      #####dose allocation#####;
      if(pt[dose.curr]<=lambda1){
        nxtdose <- dose.curr+1

        ####backfill doses#####;
        if(length(dosen[which(pe>=eta1)])>0){
          jmin<-min(dosen[which(pe>=eta1)])
        }else{
          jmin<-NULL
        }
        if(!is.null(jmin)){
          if (jmin<=dose.curr){
            backdose<-c(jmin:dose.curr)
          }else{
            backdose<-NULL
          }
        }else{
          backdose<-NULL
        }
      }else if( ((pt[dose.curr]>lambda1) & (pt[dose.curr]<=lambda2)) | ( (obs.tox[dose.curr]==1) & (nvector[dose.curr]==3))  ){
        nxtdose <- dose.curr

        ####backfill doses#####;
        if(length(dosen[which(pe>=eta1)])>0){
          jmin<-min(dosen[which(pe>=eta1)])
        }else{
          jmin<-NULL
        }
        if(!is.null(jmin)){
          if (jmin<=dose.curr-1){
            backdose<-c(jmin:(dose.curr-1))
          }else{
            backdose<-NULL
          }
        }else{
          backdose<-NULL
        }

      }else if(pt[dose.curr]>lambda2){
        nxtdose <- dose.curr-1

        ####backfill doses#####;
        if(length(dosen[which(pe>=eta1)])>0){
          jmin<-min(dosen[which(pe>=eta1)])
        }else{
          jmin<-NULL
        }
        if(!is.null(jmin)){
          if (jmin<=dose.curr-2){
            backdose<-c(jmin:(dose.curr-2))
          }else{
            backdose<-NULL
          }
        }else{
          backdose<-NULL
        }
      }


      ####eliminate dose#####;
      po.shape1 <- pr.alpha + obs.tox
      po.shape2 <- pr.beta  + (nvector-obs.tox)
      tterm     <- pbeta(target_T,po.shape1,po.shape2)

      po.shape1 <- pr.alpha + obs.eff
      po.shape2 <- pr.beta  + (nvector-obs.eff)
      eterm     <- 1-pbeta(target_E,po.shape1,po.shape2)

      ###determine which dose level should be eliminated

      if(length(which((tterm<(1-stopping.prob.T))&(nvector>=n.elimination)))>0){
        elimi_tox[(min(which((tterm<(1-stopping.prob.T))&(nvector>=n.elimination)))):n.dose]<-1
      }
      if(length(which((eterm<(1-stopping.prob.E))&(nvector>=n.elimination)))>0){
        elimi_eff[which((eterm<(1-stopping.prob.E))&(nvector>=n.elimination))]<-1
      }

      elimi<-elimi_tox+elimi_eff


      admflg  <- (elimi==0)
      admdose <- dosen[admflg]

      if(sum(admflg)==0){
        earlystop=1
        break
      }else if((length(d[which(type==1)]))==nmax){ ## dose-escalation patients are consumed! ##enrollment is stopped when 1)
        break
      }else if((nvector[dose.curr]>=stopping.npts) & ( ##enrollment is stopped when 2)
        ((pt[dose.curr]>lambda1) & (pt[dose.curr]<=lambda2)) | (dose.curr==1 & pt[dose.curr]>=lambda2) | ((dose.curr==n.dose||admflg[dose.curr+1]==0) & pt[dose.curr]<=lambda1)
      )){
        break
      }else{
        if(nxtdose==0){
          if(admflg[1]){
            dose.curr <- 1
          }else{
            break
          }
        }else if(nxtdose==(n.dose+1)){
          if(admflg[n.dose]){
            dose.curr <- n.dose
          }else{
            break
          }
        }else if(is.element(nxtdose,admdose)){
          dose.curr <- nxtdose
        }else if(dose.curr<nxtdose){
          if(sum(admdose>=nxtdose)!=0){
            dose.curr <- min(admdose[admdose>=nxtdose])
          }
        }else if(dose.curr>=nxtdose){
          if(sum(admdose<=nxtdose)!=0){
            dose.curr <- max(admdose[admdose<=nxtdose])
          }else{
            break
          }
        }

        ####for backfilling doses#####;
        if(!is.null(backdose)){
          if(length(intersect(backdose,admdose))>0){
            bdosekeep<-intersect(backdose,admdose)
          }else{
            bdosekeep<-NULL
          }
        }else{
          bdosekeep<-NULL
        }

        ####end;
      }
    }

    t.curr=max(endt,ende)
    data.obs.n[simu,] <- nvector
    data.dur[simu]    <- t.curr

    toxprobr[simu,] <- as.vector(t(toxprob))
    effprobr[simu,] <- as.vector(t(effprob))

    efficacy[simu,]=obs.eff
    toxicity[simu,]=obs.tox

    befficacy[simu,]=obs.beff
    btoxicity[simu,]=obs.btox


    backfilltimes[simu]=backfill
    backfillcount[simu,]=backfillvector

    evadose <- intersect(dosen[nvector!=0],dosen[])
    obspt <- obs.tox[evadose]/nvector[evadose]
    obspe <- obs.eff[evadose]/nvector[evadose]

    tterm.obd <- numeric(n.dose)
    eterm.obd <- numeric(n.dose)

    for(i in evadose){
      po.shape1    <- pr.alpha + obs.tox[i]
      po.shape2    <- pr.beta  + (nvector[i]-obs.tox[i])
      tterm.obd[i] <- pbeta(target_T,po.shape1,po.shape2)

      po.shape1    <- pr.alpha + obs.eff[i]
      po.shape2    <- pr.beta  + (nvector[i]-obs.eff[i])
      eterm.obd[i] <- 1-pbeta(target_E,po.shape1,po.shape2)
    }

    if(length(d)<nmax & earlystop==1){

      obd[simu] <- 0

    }else if(length(evadose)==1){

      if((tterm.obd[evadose]>=(1-stopping.prob.T))&(eterm.obd[evadose]>=(1-stopping.prob.E))){
        obd[simu] <- evadose
      }

    }else if(sum((tterm.obd[evadose]>=(1-stopping.prob.T))&(eterm.obd[evadose]>=(1-stopping.prob.E)))>=1){

      estpt <- Iso::pava(obspt)

      estpe <- obspe

      obd[simu] <- obd.select(probt=estpt, probe=estpe,tterm=tterm.obd[evadose], eterm=eterm.obd[evadose],
                              stopT=stopping.prob.T, stopE=stopping.prob.E , ph=target_T ,psi00=ppsi00, psi01=ppsi01,psi10=ppsi10,psi11=ppsi11)
    }

    ####output results####;

    for(i in 1:n.dose){
      prop.select[simu,i] <- round(mean(obd[simu]==i)*100,digits=1)
    }
    colnames(prop.select) <- dose


    prop.stop[simu] <- round(mean(obd[simu]==0)*100,digits=1)

  }



  n.patient <- round(apply(data.obs.n,2,mean),digits=2)
  names(n.patient) <- dose

  n.bpatient <- round(apply(backfillcount,2,mean),digits=2)
  names(n.bpatient) <- dose

  n.tox.patient <- round(apply(toxicity,2,mean),digits=2)
  names(n.tox.patient) <- dose

  n.eff.patient <- round(apply(efficacy,2,mean),digits=2)
  names(n.eff.patient) <- dose

  n.tox.bpatient <- round(apply(btoxicity,2,mean),digits=2)
  names(n.tox.bpatient) <- dose

  n.eff.bpatient <- round(apply(befficacy,2,mean),digits=2)
  names(n.eff.bpatient) <- dose

  names(toxprob)      <- dose
  names(effprob)      <- dose

  duration  <- round(mean(data.dur),digits=1)
  names(duration) <- "Trial duration"

  totaln = sum(data.obs.n)/an.sim
  names(totaln)<-"Total patients"

  names(target_T)          <- "Target toxicity prob."
  names(target_E)        <- "Target efficacy prob."
  names(lambda1)      <- "Lower toxicity boundary"
  names(lambda2)      <- "Upper toxicity boundary"
  names(eta1)         <- "Lower efficacy boundary"
  names(tau.T)        <- "Tox. assessment window (months)"
  names(tau.E)        <- "Eff. assessment window (months)"
  names(accrual)      <- "Accrual rate (months)"
  names(suspend)      <- "Suspension rule"
  names(target_Tr)    <- "The upper boundary for the toxicity"
  names(target_Er)    <- "The lower boundary for the efficacy"

  colnames(toxprobr) <- dose
  colnames(effprobr) <- dose
  colnames(data.obs.n) <- dose
  colnames(backfillcount) <- dose

  PCSA   <- array(0,dim=c(an.sim))
  PCAA   <- array(0,dim=c(an.sim))
  POAA   <- array(0,dim=c(an.sim))

  POS=sum(obd>mtd.true)*100/an.sim

  for (k in 1:an.sim){

    if(bd[k]==0){
      PCSA[k]=prop.stop[k]
      PCAA[k]=0
    }else{
      PCSA[k]=prop.select[k,bd[k]]
      PCAA[k]=(data.obs.n[k,]*100/sum(data.obs.n[k,]))[bd[k]]
    }

    if (mtd.true[k]==n.dose){POAA[k]=0}else{
      POAA[k]= sum(data.obs.n[k,(mtd.true[k]+1):n.dose])*100/sum(data.obs.n[k,])# percent of overdose allocation
    }
  }

  POA=round(mean(POAA),digits=1)
  PCS=round(mean(PCSA),digits=1)
  PCA=round(mean(PCAA),digits=1)

  names(PCS) <-"The percentage of trials that the optimal dose was correctly selected"
  names(PCA) <-"The percentage of patients that were correctly allocated to the optimal dose"

  names(POS) <-"The percentage of toxic doses selection"
  names(POA) <-"The percentage of patients who were allocated to toxic doses"

  result <- list(toxprob      = toxprobr,
                 effprob      = effprobr,
                 phi          = target_T,
                 delta        = target_E,
                 target_Tr    = target_Tr,
                 target_Er    = target_Er,
                 bd.true      = bd,
                 mtd.true     = mtd.true,
                 lambda1      = lambda1,
                 lambda2      = lambda2,
                 eta1         = eta1,
                 tau.T        = tau.T,
                 tau.E        = tau.E,
                 suspend      = suspend,
                 accrual      = accrual,
                 n.patient    = n.patient,
                 n.bpatient   = n.bpatient,
                 n.tox.patient= n.tox.patient,
                 n.eff.patient= n.eff.patient,
                 n.tox.bpatient= n.tox.bpatient,
                 n.eff.bpatient= n.eff.bpatient,
                 prop.select  = prop.select,
                 prop.stop    = prop.stop,
                 duration     = duration,
                 totaln       = totaln,
                 data.obs.n   = data.obs.n,
                 obd          = obd,
                 backfilltimes= backfilltimes,
                 backfillcount= backfillcount,
                 PCS=PCS,
                 PCA=PCA,
                 PTS=POS,
                 PTA=POA
  )

  class(result) <- "BF-BOIN-ET-R"
  return(result)

}
