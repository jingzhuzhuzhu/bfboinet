#' @title bfboinet_rp2
#'
#' @description Obtain the operating characteristics of a seamless two-stage phase
#' I/II trial design with backfill and joint monitoring for dose optimization
#' within fixed scenarios
#'
#' @param target_T Target toxicity probability. The default value is
#' \code{target_T=0.3}. When observing 1 DLT out of 3 patients and the target
#' DLT rate is between 0.25 and 0.279, the decision is to stay at the current
#' dose due to a widely accepted practice.
#' @param toxprob Vector of true toxicity probability.
#' @param target_E The minimum required efficacy probability. The default value
#' is \code{target_E=0.25}.
#' @param effprob Vector of true efficacy probability.
#' @param n.dose Number of dose for stage 1.
#' @param startdose Starting dose. The lowest dose is generally recommended.
#' @param ncohort Number of cohort for stage 1.
#' @param cohortsize Cohort size for stage 1.
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
#' @param n.elimination to avoid allocating patients to severely toxic doses,
#' dose elimination criteria are applied before the dose allocation decision when
#' n treated at the current dose reaches \code{n.elimination},it is a minimum sample
#' size for dose elimination.The default value is \code{n.elimination=6}.
#' @param stopping.npts Early study termination criteria for the number of
#' patients in the dose-escalation and backfill cohorts. If the number of
#' patients at the current dose reaches this criteria and the same dose level
#' is recommended as the next dose level, the study is terminated.
#' The default value is \code{stopping.npts=12}.
#' @param suspend The suspension rule that holds off the decision on dose
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
#' @param Nesc the total number of patients () in the dose-escalation cohort (Stage 1) reaches
#' the maximum total number of patients in the dose-escalation cohort (\code{Nesc}).
#' The default value is \code{Nesc=36}.
#' @param boundMTD set \code{boundMTD=TRUE} to impose the condition: the isotonic
#' estimate of toxicity probability for the selected MTD must be less than de-escalation
#' boundary. The default value is \code{boundMTD=FALSE}.
#' @param estpt.method Method to estimate the efficacy probability. Fractional
#' polynomial logistic regression is used when \code{estpt.method="fp.logistic"}.
#' Model averaging of multiple unimodal isotopic regression is used when
#' \code{estpt.method="multi.iso"}. Observed efficacy probability is used when
#' \code{estpt.method="obs.prob"}.
#' @param obd.method Method to select the optimal biological dose. Utility
#' defined by weighted function is used when \code{obd.method="utility.weighted"}.
#' Utility defined by truncated linear function is used when
#' \code{obd.method="utility.truncated.linear"}. Utility defined by scoring is
#' used when \code{obd.method="utility.scoring"}. Highest estimated efficacy
#' probability is used when \code{obd.method="max.effprob"}.
#' @param w1 Weight for toxicity-efficacy trade-off in utility defined by
#' weighted function. This must be specified when using
#' \code{obd.method="utility.weighted"}. The default value is \code{w1=0.33}.
#' @param w2 Weight for penalty imposed on toxic doses in utility defined by
#' weighted function. This must be specified when using
#' \code{obd.method="utility.weighted"}. The default value is \code{w2=1.09}.
#' @param plow.ast Lower threshold of toxicity linear truncated function. This
#' must be specified when using \code{obd.method="utility.truncated.linear"}.
#' The default value is \code{plow.ast=pT.saf}.
#' @param pupp.ast Upper threshold of toxicity linear truncated function. This
#' must be specified when using \code{obd.method="utility.truncated.linear"}.
#' The default value is \code{pupp.ast=pT.tox}.
#' @param qlow.ast Lower threshold of efficacy linear truncated function. This
#' must be specified when using \code{obd.method="utility.truncated.linear"}.
#' The default value is \code{qlow.ast=pE.saf/2}.
#' @param qupp.ast Upper threshold of efficacy linear truncated function. This
#' must be specified when using \code{obd.method="utility.truncated.linear"}.
#' The default value is \code{qupp.ast=target_E}.
#' @param stage1.method Method to patient assignment for backfilling patients.
#' Pick-the-winner is used when \code{stage1.method="PW"}. Equal randomization is
#' used when \code{stage1.method="ER"}.
#' @param H0 Stage 2: A numeric value for the response rate under the null hypothesis
#' (toxicity - OR, no toxicity - OR, toxicity - no OR, no toxicity - No OR).
#' @param H1 Stage 2: A numeric value for the response rate under the alternative hypothesis.
#' @param nIA.sample Stage 2: A numeric vector representing the additional patients enrolled at each interim analysis.
#' The value at index 'i' indicates the number of patients at interim analysis 'i'.
#' For example, for four interim analyses with total sample sizes of 10, 15, 20, and 30,
#' the vector would be represented as 'n = c(10, 15, 20, 30)'.
#' @param nIA Stage 2: A numeric value for the number of interim analysis. The default value is \code{nIA=length(nIA.sample)}.
#' @param method Stage 2: A character string specifying the method to use for calculating cutoff values for the efficacy stopping.
#' Options are "power" (default) or "OF" for "O'Brien-Fleming".
#' @param t1e_optimal_pars Stage 2: Desired Type - I error rate. If specified it will only return results with type I error rate less the specified value.
#' @param lambda1_optimal_pars Stage 2: Starting value for 'lambda' values to search.
#' @param lambda2_optimal_pars Stage 2: Ending value for 'lambda' values to search.
#' @param grid1_optimal_pars Stage 2: Number of 'lambda' values to consider between lambda1 and lambda2. A fine grid by 0.01 is recommended.
#' @param gamma1_optimal_pars Stage 2: Starting value for 'gamma' values to search.
#' @param gamma2_optimal_pars Stage 2: Ending value for 'gamma' values to search.
#' @param grid2_optimal_pars Stage 2: Number of 'gamma' values to consider between gamma1 and gamma2. A fine grid by 0.01 is recommended.
#' @param eta1_optimal_pars Stage 2: Starting value for 'eta' values to search.
#' @param eta2_optimal_pars Stage 2: Ending value for 'eta' values to search.
#' @param grid3_optimal_pars Stage 2: Number of eta values to consider between eta1 and eta2. A fine grid by 0.01 is recommended.
#' @param ppsi01 Score for toxicity=yes and efficacy=no in utility defined by
#' scoring.The default value is \code{psi01=0}.
#' @param ppsi00 Score for toxicity=no and efficacy=no in utility defined by
#' scoring. The default value is \code{psi00=40}.
#' @param ppsi11 Score for toxicity=yes and efficacy=yes in utility defined by
#' scoring. The default value is \code{psi11=60}.
#' @param ppsi10 Score for toxicity=no and efficacy=yes in utility defined by
#' scoring. The default value is \code{psi10=100}.
#' @param n.sim Number of simulated trial. The default value is
#' \code{n.sim=1000}.
#' @param seed.sim Seed for random number generator. The default value is
#' \code{seed.sim=100}.
#' @details The \code{bfboinet_rp2} is a function which generates the operating
#' characteristics of the seamless two-stage Phase I/II trial design integrating
#' dose optimization with efficacy evaluation by a simulation study.
#' Users can specify a variety of study settings to simulate studies. The
#' operating characteristics of the design are summarized by the percentage of
#' times that each dose level was selected as optimal biological dose and the
#' average number of patients who were treated at each dose level. The
#' percentage of times that the study was terminated and the expected study
#' duration are also provided.
#' @return
#' The \code{get.oc.backboinet_rp2} returns a list containing the following components:
#' \item{toxprob}{True toxicity probability.}
#' \item{effprob}{True efficacy probability.}
#' \item{phi}{Target toxicity probability.}
#' \item{delta}{Target efficacy probability.}
#' \item{lambda1}{Lower toxicity boundary in dose escalation/de-escalation.}
#' \item{lambda2}{Upper toxicity boundary in dose escalation/de-escalation.}
#' \item{eta1}{Lower efficacy boundary in dose escalation/de-escalation.}
#' \item{tau.T}{Toxicity assessment windows (months).}
#' \item{tau.E}{Efficacy assessment windows (months).}
#' \item{suspend}{The suspension rule that holds off the decision on dose
#' allocation for the dose-escalation cohort until sufficient toxicity
#' information is available.}
#' \item{accrual}{Accrual rate (months) (patient accrual rate per month).}
#' \item{n.patient.all}{Average number of patients who were treated at each dose
#' level at stage 1 and stage 2.}
#' \item{nptsdosepct.all}{The percentage of patients who were treated at each dose
#' level at stage 1 and stage 2.}
#' \item{n.tox.patient.all}{Average number of patients who experienced toxicity at each dose
#' level at stage 1 and stage 2.}
#' \item{n.eff.patient.all}{Average number of patients who experienced efficacy at each dose
#' level at stage 1 and stage 2.}
#' \item{n.patient.stage2}{Average number of patients who were treated at each dose
#' level at stage 2.}
#' \item{n.bpatient}{Average number of back filled patients who were treated
#' at each dose level at stage 1.}
#' \item{prop.select}{Percentage of times that each dose level was selected as
#' optimal biological dose at stage 1 and stage 2.}
#' \item{MTD.select}{Percentage of times that each dose level was selected as
#'  maximum tolerated dose at stage 1 and stage 2.}
#' \item{claim.select}{Percentage of times that each dose level was claimed efficacy
#' at stage 1 and stage 2.}
#' \item{prop.stop}{Percentage of times that the study was terminated at stage 1.}
#' \item{duration}{Expected study duration (months) at stage 1 and stage 2.}
#' \item{duration1}{Expected study duration (months) at stage 1.}
#' \item{duration2}{Expected study duration (months) at stage 2.}
#' \item{totaln}{Total patients at stage 1 and stage 2.}
#' \item{data.obs.n}{Record the number of patients in each dose level within the
#' simulations during the trial at stage 1 and stage 2.}
#' \item{data.obs.n.stage2}{Record the number of patients in each dose level within the
#' simulations during the trial at stage 2.}
#' \item{data.obs.n.stage1}{Record the number of patients in each dose level within the
#' simulations during the trial at stage 1.}
#' \item{obd}{Record the optimal dose in each simulation during the trial at stage 1 and stage 2.}
#' \item{claim}{Record the optimal dose with efficacy in each simulation during the trial at stage 1 and stage 2.}
#' \item{PCS}{The percentage of trials that the optimal dose was correctly
#' selected at stage 1 and stage 2.}
#' \item{PCC}{The percentage of patients that efficacy were correctly allocated to the
#' optimal dose at stage 1 and stage 2.}
#' \item{PTS}{The percentage of toxic doses selection at stage 1 and stage 2.}
#' \item{PTA}{The percentage of patients who were allocated to toxic doses at stage 1 and stage 2.}
#' \item{boundary_tab.all.out}{Futility and efficacy stopping boundaries at stage 2.}
#' \item{lambda}{Lambda values for cut-off probabilitys at stage 2.}
#' \item{gamma}{Gamma values for cut-off probability at stage 2.}
#' \item{eta}{Eta values for cut-off probability at stage 2.}
#' @references
#' A seamless two-stage phase I/II trial design with backfill
#' and joint efficacy and toxicity monitoring as described in (Takeda et al (2026) <doi:10.1002/pst.70092>.
#' @examples
#'
#' target_T=0.3
#' target_E=0.25
#' pttt=c(0.02, 0.05, 0.07, 0.10, 0.15)
#' pee=c(0.05, 0.08, 0.15, 0.30, 0.45)
#' #####Stage 2####;
#' H0=c(0.10, 0.15, 0.30, 0.45)
#' H1=c(0.05, 0.45, 0.15, 0.35)
#' nIA.sample=c(24,30,36,42,48)
#' nIA=length(nIA.sample)
#' t1e_optimal_pars=0.1
#' lambda1_optimal_pars=0
#' lambda2_optimal_pars=1
#' grid1_optimal_pars=101
#' gamma1_optimal_pars=0
#' gamma2_optimal_pars=1
#' grid2_optimal_pars=101
#' eta1_optimal_pars=0
#' eta2_optimal_pars=3
#' grid3_optimal_pars=301
#' \dontrun{
#'get.oc.backboinet_rp2(target_T=target_T, toxprob=pttt,target_E=target_E,effprob=pee,n.dose=5,
#'startdose=1,ncohort=40,cohortsize=3,pT.saf=0.6 * target_T,pT.tox = 1.4 * target_T,
#'pE.saf = 0.6 * target_E,alpha.T1=0.5,alpha.E1=0.5,tau.T=1,tau.E=1,te.corr=0.2,
#'gen.event.time="weibull",accrual=3,gen.enroll.time="uniform",n.elimination=6,
#'stopping.npts=12,suspend=0,stopping.prob.T=0.95,stopping.prob.E=0.90,Nesc=36,
#'boundMTD = FALSE,estpt.method="obs.prob", obd.method="utility.scoring",
#'w1= 0.33, w2=1.09,plow.ast=pT.saf, pupp.ast=pT.tox, qlow.ast=pE.saf/2, qupp.ast=target_E,
#'stage1.method="ER",H0=H0,H1=H1,nIA.sample=nIA.sample,nIA=length(nIA.sample),
#'t1e_optimal_pars=t1e_optimal_pars,lambda1_optimal_pars=lambda1_optimal_pars,
#'lambda2_optimal_pars=lambda2_optimal_pars,grid1_optimal_pars=grid1_optimal_pars,
#'gamma1_optimal_pars=gamma1_optimal_pars,gamma2_optimal_pars=gamma2_optimal_pars,
#'grid2_optimal_pars=grid2_optimal_pars,eta1_optimal_pars=eta1_optimal_pars,
#'eta2_optimal_pars=eta2_optimal_pars,grid3_optimal_pars=grid3_optimal_pars,
#'ppsi01=0,ppsi00=40,ppsi11=60, ppsi10=100,n.sim=1,seed.sim=100)
#'}
#'
#' @import Iso copula BOP2FE dplyr magrittr boinet BOIN
#' @importFrom stats binomial dbinom pbeta pbinom rmultinom runif rexp
#' @export

get.oc.backboinet_rp2 <- function (target_T=0.3, toxprob,target_E=0.25,effprob,n.dose,startdose,ncohort,cohortsize,
                                pT.saf=0.6 * target_T,pT.tox = 1.4 * target_T,pE.saf = 0.6 * target_E,
                                alpha.T1=0.5,alpha.E1=0.5,tau.T,tau.E,te.corr=0.2,gen.event.time="weibull",
                                accrual,gen.enroll.time="uniform",n.elimination=6,stopping.npts=12,
                                suspend=0,stopping.prob.T=0.95,stopping.prob.E=0.90,Nesc=36,boundMTD = FALSE,
                                estpt.method, obd.method,
                                w1= 0.33, w2=1.09,
                                plow.ast=pT.saf, pupp.ast=pT.tox, qlow.ast=pE.saf/2, qupp.ast=target_E,stage1.method,
                                H0,H1,nIA.sample,nIA=length(nIA.sample),method="power",t1e_optimal_pars,
                                lambda1_optimal_pars,lambda2_optimal_pars,grid1_optimal_pars,
                                gamma1_optimal_pars,gamma2_optimal_pars,grid2_optimal_pars,eta1_optimal_pars,
                                eta2_optimal_pars,grid3_optimal_pars,
                                ppsi01=0,ppsi00=40,ppsi11=60,
                                ppsi10=100,n.sim=1000,seed.sim=100){

  if(length(toxprob)!=n.dose){
    stop("Number of dose must be the same as the length of true toxicity probability.")

  }else if(length(effprob)!=n.dose){
    stop("Number of dose must be the same as the length of true efficacy probability.")

  }else if(!((pT.saf<target_T)&(target_T<pT.tox))){
    stop("Design parameters must satisfy a condition of pT.saf < target_T < pT.tox.")

  }else if(!(pE.saf<target_E)){
    stop("Design parameters must satisfy a condition of pE.saf < target_E.")

  }

  ####calculate toxicity escalation and de-escalation boundaries#####;
  lambda1 = log((1 - pT.saf)/(1 - target_T))/log(target_T *
                                                   (1 - pT.saf)/(pT.saf * (1 - target_T)))

  lambda2 = log((1 - target_T)/(1 - pT.tox))/log(pT.tox * (1 -
                                                             target_T)/(target_T * (1 - pT.tox)))

  ####calculate efficacy boundaries#####;

  eta1 = log((1-pE.saf)/(1-target_E))/log(target_E*(1-pE.saf)/(pE.saf* (1-target_E)))

  ####generate event time for efficacy and toxicity######;
  dosen <- 1:n.dose
  dose  <- paste("Dose",dosen,sep="")

  toxp <- data.frame(t(toxprob))
  colnames(toxp) <- dose

  effp <- data.frame(t(effprob))
  colnames(effp) <- dose

  pr.alpha <- 1
  pr.beta  <- 1

  alpha.T1 <- alpha.T1
  alpha.T2 <- 0.5
  alpha.E1 <- alpha.E1
  alpha.E2 <- 0.5

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

  data.obs.n <- array(0,dim=c(n.sim,n.dose))
  data.obs.dose.stage2 <- array(0,dim=c(n.sim,n.dose))

  data.obs.n.stage2 <- array(0,dim=c(n.sim,n.dose))
  data.obs.n.stage1 <- array(0,dim=c(n.sim,n.dose))
  data.obs.earlystopfuti <- array(0,dim=c(n.sim,n.dose))
  data.obs.earlystopsupe <- array(0,dim=c(n.sim,n.dose))

  data.dur   <- array(0,dim=c(n.sim))
  data.dur1   <- array(0,dim=c(n.sim))
  data.dur2   <- array(0,dim=c(n.sim))

  obd <- array(0,dim=c(n.sim))
  claim<- array(0,dim=c(n.sim))
  MTD<- array(0,dim=c(n.sim))

  stage2.earlystop.f<-array(0,dim=c(n.sim))
  stage2.earlystop.e<-array(0,dim=c(n.sim))

  toxicity=matrix(nrow=n.sim,ncol=n.dose)
  efficacy=matrix(nrow=n.sim,ncol=n.dose)

  dose.curr=startdose

  backfilltimes=rep(0,n.sim) ## record how may times we back-filled during the trial
  backfillcount=matrix(nrow=n.sim,ncol=n.dose) ## record the location of backfill

  nmax=Nesc  ####max number of patients for dose-escalation#####;

  set.seed(seed.sim)


  #####apply futility and efficacy stopping rule for stage 2#######;
  nIA.sample.b<-numeric(length(nIA.sample))
  for (i in 1: length(nIA.sample)){
    if (i==1){
      nIA.sample.b[i]<-nIA.sample[i]
    }else{
      nIA.sample.b[i]<-nIA.sample[i]-nIA.sample[i-1]
    }

  }

  ####Please update nsim to 5000 for robust result#####;
  Optimal<- BOP2FE_jointefftox(
    H0=H0,
    H1= H1,
    n = nIA.sample.b,
    nsim = 5000, t1e = t1e_optimal_pars, method = method,
    lambda1 = lambda1_optimal_pars, lambda2 = lambda2_optimal_pars, grid1 = grid1_optimal_pars,
    gamma1 = gamma1_optimal_pars, gamma2 = gamma2_optimal_pars, grid2 = grid2_optimal_pars,
    eta1 = eta1_optimal_pars, eta2 = eta2_optimal_pars, grid3 = grid3_optimal_pars,
    seed = seed.sim
  )

  Optimal2<-summary(Optimal)

  lambda=Optimal2$opt_pars$lambda
  gamma=Optimal2$opt_pars$gamma
  eta=Optimal2$opt_pars$eta

  cn11f_max<-Optimal2$boundary["Futility boundary (OR)",]

  cn12f_max<-Optimal2$boundary["Futility boundary (Tox)",]

  cn11s_min<-Optimal2$boundary["Efficacy boundary (OR)",]

  cn12s_min<-Optimal2$boundary["Efficacy boundary (Tox)",]

  Number_patients<-nIA.sample

  boundary_tab.all<-as.matrix(cbind(Number_patients,cn11f_max,cn12f_max,cn11s_min,cn12s_min))

  boundary_tab.all.out<-boundary_tab.all

  colnames(boundary_tab.all.out) <- c("N","Futility-# of responses<=","Futility-# of toxicities>=","Efficacy-# of responses>=","Efficacy-# of toxicities<=")

  for(simu in 1:n.sim){
    set.seed(seed.sim+simu)

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
    nvector.stage2=rep(0,n.dose)

    bdosekeep=NULL
    backfill=0 # how many times of backfill.
    backfillvector=rep(0,length=n.dose)  ## which dose is backfilled.
    elimi_tox<-rep(0,n.dose)
    elimi_eff<-rep(0,n.dose)
    elimi<-rep(0,n.dose)
    earlystop<-0
    type=NULL
    stage=NULL
    stage.IA=NULL
    earlystopfuti=rep(0,n.dose)
    earlystopsupe=rep(0,n.dose)

    ######Stage 1: ncohort#####;
    for(i in 1:ncohort){

      ### assign three maxpatients to dose.curr
      # t.curr; t.finish;t.onset

      localt.enter=NULL
      remains=cohortsize

      if((length(d[d==dose.curr])+cohortsize)>stopping.npts & (length(d[d==dose.curr])<stopping.npts) ){ ## do not have a full cohort remains;
        remains=stopping.npts-(length(d[d==dose.curr]))
      }

      if(remains==0){
        break
      }



      for(j in 1:remains){
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

      time.te <- copula::rMvdc(remains,mv.ncop[[dose.curr]]) ####event time#####

      localt.onset=localt.enter+time.te[,1]

      localt.yT=as.numeric(time.te[,1]<=tau.T)
      localt.yE=as.numeric(time.te[,2]<=tau.E)

      yT=c(yT,as.numeric(time.te[,1]<=tau.T))
      yE=c(yE,as.numeric(time.te[,2]<=tau.E))

      d=c(d,rep(dose.curr,remains))
      type=c(type,rep(1,remains))

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



      ende=apply(data.frame(
        ende1=t.enter+tau.E,
        ende2=t.enter+t.eventE
      ), 1, min, na.rm=TRUE)

      nvector[dose.curr]=nvector[dose.curr]+remains

      stage=c(stage,rep(1,remains))
      stage.IA=c(stage.IA,rep(0,remains))
      #####added for more than or equal 1-suspend patients completed the DLT assessment####;
      if (suspend!=0){
        quantile=quantile(endt[which(d==dose.curr)], probs = (1-suspend))
        quantile50=max(endt[which(d==dose.curr)][endt[which(d==dose.curr)]<=quantile])
      }
      #####end######;

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
                titeb.df <- data.frame(dose   = d,
                                       enter  = t.enter,
                                       endtox = endt,
                                       dlt    = yT,
                                       endeff = ende,
                                       orr    = yE)

                if(stage1.method=="PW"){
                  #####using the latest information to get the max utility score#####

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
              }


              evadose <- intersect(intersect(dosen[nvector!=0],bdosekeep),dosen[nvector<stopping.npts])

            if(stage1.method=="PW"){
              ###select the backfilling dose using utility score#####;
              if(!identical(evadose, integer(0))) {
                estpt <- Iso::pava(pt[evadose])

                if(estpt.method=="multi.iso"){
                  estpe <- multi.iso(obs=obs.eff[evadose],n=obs.eff.n[evadose])

                }else if(estpt.method=="fp.logistic"){
                  estpe <- fp.logit(obs=obs.eff[evadose],n=obs.eff.n[evadose],dose=evadose)

                }else if(estpt.method=="obs.prob"){
                  estpe <- pe[evadose]
                }

                if(obd.method=="utility.weighted"){
                  utility_b <- utility.weighted(probt=estpt,probe=estpe,
                                          w1=w1,w2=w2,tox.upper=pT.tox)

                  utility_b.max <- evadose[utility_b==max(utility_b)]

                  s_back     <- min(intersect(utility_b.max,evadose)) ###backfilling dose#####;


                }else if(obd.method=="utility.truncated.linear"){
                  utility_b <- utility.truncated.linear(probt=estpt,probe=estpe,
                                                  tlow=plow.ast,tupp=pupp.ast,
                                                  elow=qlow.ast,eupp=qupp.ast)

                  utility_b.max <- evadose[utility_b==max(utility_b)]

                  s_back     <- min(intersect(utility_b.max,evadose)) ###backfilling dose#####;

                }else if(obd.method=="utility.scoring"){
                  utility_b <- utility.scoring(probt=estpt,probe=estpe,
                                         psi00=ppsi00,psi11=ppsi11)

                  utility_b.max <- evadose[utility_b==max(utility_b)]

                  s_back     <- min(intersect(utility_b.max,evadose)) ###backfilling dose#####;

                }

              }else{
                s_back=NULL
              }
            }else if(stage1.method=="ER"){
              ###select the backfilling dose using equal randomization approach#####;
              if(!identical(evadose, integer(0))) {
                s_back=sample(evadose,1)
              }else{
                s_back=NULL
              }
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
                stage=c(stage,1)
                stage.IA=c(stage.IA,0)

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
                                  endtox = localt.enter+apply(as.matrix(1:remains),1,function(x){min(time.te[x,1],tau.T)}),
                                  dlt    = localt.yT,
                                  endeff = localt.enter+apply(as.matrix(1:remains),1,function(x){min(time.te[x,2],tau.E)}),
                                  orr    = localt.yE))

      if(!is.null(localt.enterb)){
        tite.df <- rbind(tite.df,
                         data.frame(dose   = db,
                                    enter  = localt.enterb,
                                    endtox = localt.enterb+apply(as.matrix(1:length(localt.enterb)),1,function(x){min(time.teba[x,1],tau.T)}),
                                    dlt    = localt.yTb,
                                    endeff = localt.enterb+apply(as.matrix(1:length(localt.enterb)),1,function(x){min(time.teba[x,2],tau.E)}),
                                    orr    = localt.yEb))
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
      #if(length(which((eterm<(1-stopping.prob.E))&(nvector>=n.earlystop)))>0){
      #  elimi_eff[which((eterm<(1-stopping.prob.E))&(nvector>=n.earlystop))]<-1
      #}

      #elimi<-elimi_tox+elimi_eff
      elimi<-elimi_tox


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

        if(nvector[dose.curr]>=stopping.npts){
          break
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
    data.dur1[simu]    <- t.curr
    ####start to stage 2######;
    evadose_stage1=NULL  ###inital set up#####;
    admose_stage1<-NULL
    obd_stage1 <- 0
    ######1, select admissible doses for stage II#####;
    ####1), MTD is selected #####;
    if (earlystop==1) {
      MTD_stage1 = 0
    }else {
      MTD_stage1 = select_mtd(target=target_T, npts=nvector, ntox=obs.tox, cutoff.eli = stopping.prob.T, extrasafe = FALSE,
                              offset = 0, p.tox=pT.tox, boundMTD = boundMTD, n.elimination=n.elimination)$MTD
    }
    #####end#####;

    #####2), select the dose that maximizes the utility other than MTD#####;

    if(MTD_stage1 != 0){
      evadose_stage1 <- intersect(intersect(dosen[nvector!=0],dosen[]),dosen[dosen<MTD_stage1])
      obspt_stage1 <- obs.tox[evadose_stage1]/nvector[evadose_stage1]
      obspe_stage1 <- obs.eff[evadose_stage1]/nvector[evadose_stage1]

      tterm.obd_stage1 <- numeric(n.dose)
      eterm.obd_stage1 <- numeric(n.dose)

      for(i in evadose_stage1){
        po.shape1_stage1    <- pr.alpha + obs.tox[i]
        po.shape2_stage1    <- pr.beta  + (nvector[i]-obs.tox[i])
        tterm.obd_stage1[i] <- pbeta(target_T,po.shape1_stage1,po.shape2_stage1)

        po.shape1_stage1    <- pr.alpha + obs.eff[i]
        po.shape2_stage1    <- pr.beta  + (nvector[i]-obs.eff[i])
        eterm.obd_stage1[i] <- 1-pbeta(target_E,po.shape1_stage1,po.shape2_stage1)
      }

      evadose_stage1 <- which((tterm.obd_stage1[evadose_stage1]>=(1-stopping.prob.T)))

      if(length(evadose_stage1)==1){

        obd_stage1 <- evadose_stage1

      }else if(length(evadose_stage1)>=1){

        estpt_stage1 <- Iso::pava(obspt_stage1)

        if(estpt.method=="multi.iso"){
          estpe_stage1 <- multi.iso(obs=obs.eff[evadose_stage1],n=nvector[evadose_stage1])

        }else if(estpt.method=="fp.logistic"){
          estpe_stage1 <- fp.logit(obs=obs.eff[evadose_stage1],n=nvector[evadose_stage1],dose=evadose_stage1)

        }else if(estpt.method=="obs.prob"){
          estpe_stage1 <- obspe_stage1
        }

        if(obd.method=="utility.weighted"){
          utility_stage1 <- utility.weighted(probt=estpt_stage1,probe=estpe_stage1,
                                  w1=w1,w2=w2,tox.upper=pT.tox)

        }else if(obd.method=="utility.truncated.linear"){
          utility_stage1 <- utility.truncated.linear(probt=estpt_stage1,probe=estpe_stage1,
                                          tlow=plow.ast,tupp=pupp.ast,
                                          elow=qlow.ast,eupp=qupp.ast)

        }else if(obd.method=="utility.scoring"){
          utility_stage1 <- utility.scoring(probt=estpt_stage1,probe=estpe_stage1,
                                 psi00=ppsi00,psi11=ppsi11)

        }

        obd_stage1=max(evadose_stage1[utility_stage1==max(utility_stage1)])
      }

    }

    #####admissible dose sets for stage 2#####;

    if(MTD_stage1!=0 & obd_stage1!=0){
      admose_stage1<-c(obd_stage1,MTD_stage1)
    }else if(MTD_stage1!=0 & obd_stage1==0){
      admose_stage1=MTD_stage1
    }else if(MTD_stage1==0 & obd_stage1!=0){
      admose_stage1=obd_stage1
    }

    ###end#####;
    #####end######;


    for (ia in 1: nIA){

      if (is.null(admose_stage1) |   length(admose_stage1)==0){
        break
      }

      admose_stage1_futility<-admose_stage1

      for(dose_stage2 in 1:length(admose_stage1)){
        dose.curr=admose_stage1[dose_stage2]

        data.obs.dose.stage2[simu,dose.curr]=1

        ### assign remaining maxpatients to dose.curr
        # t.curr; t.finish;t.onset

        localt.enter=NULL
        remains=nIA.sample[ia]-nvector[dose.curr]

        for(j2 in 1:remains){
          ## local t.enter
          if(j2==1){
            localt.enter=c(localt.enter,t.curr)

          }else{
            if(gen.enroll.time=="uniform"){ localt.enter =c(localt.enter, localt.enter[length(localt.enter)]+ runif(1, 0, 2/accrual))}
            if(gen.enroll.time=="exponential"){ localt.enter = c(localt.enter,  localt.enter[length(localt.enter)]+ rexp(1, rate=accrual))}

          }
        }

        localt.finish=localt.enter+tau.T

        t.enter=c(t.enter,localt.enter)

        t.finish=c(t.finish,localt.finish)

        time.te <- copula::rMvdc(remains,mv.ncop[[dose.curr]]) ####event time#####

        localt.onset=localt.enter+time.te[,1]

        localt.yT=as.numeric(time.te[,1]<=tau.T)
        localt.yE=as.numeric(time.te[,2]<=tau.E)

        yT=c(yT,as.numeric(time.te[,1]<=tau.T))
        yE=c(yE,as.numeric(time.te[,2]<=tau.E))

        d=c(d,rep(dose.curr,remains))
        type=c(type,rep(1,remains))

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

        ende=apply(data.frame(
          ende1=t.enter+tau.E,
          ende2=t.enter+t.eventE
        ), 1, min, na.rm=TRUE)


        nvector[dose.curr]=nvector[dose.curr]+remains

        nvector.stage2[dose.curr]=nvector.stage2[dose.curr]+remains

        type=c(type,rep(2,remains))


        stage=c(stage,rep(2,remains))
        stage.IA=c(stage.IA,rep(ia,remains))


        tite.df <- rbind(tite.df,
                         data.frame(dose   = dose.curr,
                                    enter  = localt.enter,
                                    endtox = localt.enter+apply(as.matrix(1:remains),1,function(x){min(time.te[x,1],tau.T)}),
                                    dlt    = localt.yT,
                                    endeff = localt.enter+apply(as.matrix(1:remains),1,function(x){min(time.te[x,2],tau.E)}),
                                    orr    = localt.yE))

        t.curr <- t.enter[length(t.enter)]

      }

      t.curr<- max(endt,ende)

      gamma.all.T <- as.numeric(tite.df$endtox<=t.curr)
      gamma.all.E <- as.numeric(tite.df$endeff<=t.curr)

      for(ds in 1:n.dose){
        if(sum(tite.df$dose==ds)>0){

          x.DLT  <- sum(tite.df[(tite.df$dose==ds),]$dlt)
          n.DLT  <- nrow(tite.df[(tite.df$dose==ds),])

          x.ORR  <- sum(tite.df[(tite.df$dose==ds),]$orr)
          n.ORR  <- nrow(tite.df[(tite.df$dose==ds),])

          obs.tox[ds]   <- x.DLT
          obs.tox.n[ds] <- n.DLT
          pt[ds]        <- x.DLT/n.DLT

          obs.eff[ds]   <- x.ORR
          obs.eff.n[ds] <- n.ORR
          pe[ds]        <- x.ORR/n.ORR

        }}

      estpt_stage2 <- Iso::pava(pt)

      if(estpt.method=="multi.iso"){
        estpe_stage2 <- multi.iso(obs=obs.eff,n=obs.eff.n)

      }else if(estpt.method=="fp.logistic"){
        estpe_stage2 <- fp.logit(obs=obs.eff,n=obs.eff.n,dose=dosen[])

      }else if(estpt.method=="obs.prob"){
        estpe_stage2 <- pe
      }


      if(obd.method=="utility.weighted"){
        utility_stage2 <- utility.weighted(probt=estpt_stage2,probe=estpe_stage2,
                                w1=w1,w2=w2,tox.upper=pT.tox)

      }else if(obd.method=="utility.truncated.linear"){
        utility_stage2 <- utility.truncated.linear(probt=estpt_stage2,probe=estpe_stage2,
                                        tlow=plow.ast,tupp=pupp.ast,
                                        elow=qlow.ast,eupp=qupp.ast)

      }else if(obd.method=="utility.scoring"){
        utility_stage2 <- utility.scoring(probt=estpt_stage2,probe=estpe_stage2,
                               psi00=ppsi00,psi11=ppsi11)

      }


      #####apply futility and efficacy stopping rule for stage 2#######;
      boundary_tab<-boundary_tab.all[ia,]

      for(dose_stage2 in 1:length(admose_stage1)){
        dose.curr=admose_stage1[dose_stage2]
        if(obs.eff[dose.curr]<=boundary_tab["cn11f_max"] | obs.tox[dose.curr]>=boundary_tab["cn12f_max"] ){
          earlystopfuti[dose.curr]<-1
        }

        if(obs.eff[dose.curr]>=boundary_tab["cn11s_min"] & obs.tox[dose.curr]<=boundary_tab["cn12s_min"] ){
          if(ia<nIA){
            if(max(utility_stage2[admose_stage1])==utility_stage2[dose.curr]) {earlystopsupe[dose.curr]<-1}
          }else{
            earlystopsupe[dose.curr]<-1
          }
        }

        if(ia<nIA ){
          if(earlystopsupe[dose.curr]==1){admose_stage1=NULL; stage2.earlystop.e[simu]<-1; break;}
          if(all(earlystopfuti[admose_stage1_futility] == 1)){ admose_stage1=NULL; stage2.earlystop.f[simu]<-1; break;}
        }

      }

      if(is.null(admose_stage1)){
        break
      }else{
        admose_stage1=intersect(intersect(dosen[earlystopfuti==0],admose_stage1),dosen[earlystopsupe==0])
      }





      ######end#######;


    }



    data.obs.n[simu,] <- nvector
    data.obs.n.stage2[simu,]<- nvector.stage2
    data.obs.n.stage1[simu,]<- nvector-nvector.stage2
    data.obs.earlystopfuti[simu,] <- earlystopfuti
    data.obs.earlystopsupe[simu,] <- earlystopsupe

    data.dur[simu]    <- t.curr
    data.dur2[simu]    <- t.curr-data.dur1[simu]

    efficacy[simu,]=obs.eff
    toxicity[simu,]=obs.tox

    backfilltimes[simu]=backfill
    backfillcount[simu,]=backfillvector

    evadose <- intersect(dosen[nvector>stopping.npts],dosen[])
    obspt <- obs.tox[evadose]/nvector[evadose]
    obspe <- obs.eff[evadose]/nvector[evadose]



    if(earlystop==1 | length(stage[stage==2])==0 | (max(nvector)<max(nIA.sample)  &  stage2.earlystop.e[simu]==0 ) ){

      obd[simu] <- 0

    }else if(stage2.earlystop.e[simu]==1 & sum(earlystopsupe)==1){

      if(length(intersect(dosen[earlystopsupe==1],evadose))==0) {
        obd[simu] <- 0
      }else{
        obd[simu] <- intersect(dosen[earlystopsupe==1],evadose)
      }

    }else if(length(evadose)>0){

      estpt <- Iso::pava(obspt)

      if(estpt.method=="multi.iso"){
        estpe <- multi.iso(obs=obs.eff[evadose],n=nvector[evadose])

      }else if(estpt.method=="fp.logistic"){
        estpe <- fp.logit(obs=obs.eff[evadose],n=nvector[evadose],dose=evadose)

      }else if(estpt.method=="obs.prob"){
        estpe <- obspe
      }

      if(obd.method=="utility.weighted"){
        utility <- utility.weighted(probt=estpt,probe=estpe,
                                w1=w1,w2=w2,tox.upper=pT.tox)

      }else if(obd.method=="utility.truncated.linear"){
        utility <- utility.truncated.linear(probt=estpt,probe=estpe,
                                        tlow=plow.ast,tupp=pupp.ast,
                                        elow=qlow.ast,eupp=qupp.ast)

      }else if(obd.method=="utility.scoring"){
        utility <- utility.scoring(probt=estpt,probe=estpe,
                               psi00=ppsi00,psi11=ppsi11)

      }

      obd[simu]=min(evadose[utility==max(utility)])

    }

    if(length(intersect(obd[simu],dosen[earlystopsupe==1]))>0){
      claim[simu]<-intersect(obd[simu],dosen[earlystopsupe==1])
    }

    ####MTD based on stage1 and stage 2 #####;
    if (earlystop==1) {
      MTD[simu] = 0
    }else {
      MTD[simu] = select_mtd(target=target_T, npts=nvector, ntox=obs.tox, cutoff.eli = stopping.prob.T, extrasafe = FALSE,
                             offset = 0, p.tox=pT.tox, boundMTD = boundMTD, n.elimination=n.elimination)$MTD
    }

  }



  ####output results####;
  prop.select <- array(0,dim=c(n.dose))
  for(i in 1:n.dose){
    prop.select[i] <- round(mean(obd==i)*100,digits=1)
  }
  names(prop.select) <- dose

  claim.select <- array(0,dim=c(n.dose))
  for(i in 1:n.dose){
    claim.select[i] <- round(mean(claim==i)*100,digits=1)
  }
  names(claim.select) <- dose

  MTD.select <- array(0,dim=c(n.dose))
  for(i in 1:n.dose){
    MTD.select[i] <- round(mean(MTD==i)*100,digits=1)
  }
  names(MTD.select) <- dose


  prop.stop <- round(mean(obd==0)*100,digits=1)
  names(prop.stop) <- "Stop %"

  prop.stop.futility <- round(mean(stage2.earlystop.f==1)*100,digits=1)
  names(prop.stop.futility) <- "Stop % for futility"

  prop.stop.efficacy <- round(mean(stage2.earlystop.e==1)*100,digits=1)
  names(prop.stop.efficacy) <- "Stop % for efficacy"

  n.patient <- round(apply(data.obs.n,2,mean),digits=2)
  names(n.patient) <- dose

  n.patient.stage2 <- round(apply(data.obs.n.stage2,2,mean),digits=2)
  names(n.patient.stage2) <- dose

  n.bpatient <- round(apply(backfillcount,2,mean),digits=2)
  names(n.bpatient) <- dose

  doselevel.stage2 <- round(apply(data.obs.dose.stage2,2,mean),digits=2)
  names(doselevel.stage2) <- dose

  duration  <- round(mean(data.dur),digits=1)
  names(duration) <- "Trial duration"

  duration1  <- round(mean(data.dur1),digits=1)
  names(duration1) <- "Trial duration 1"

  duration2  <- round(mean(data.dur2),digits=1)
  names(duration2) <- "Trial duration 2"

  totaln = round(sum(data.obs.n)/n.sim,1)
  names(totaln)<-"Total patients"

  nptsdose = apply(data.obs.n,2,mean); # number of selecting each dose
  nptsdosepct=round(nptsdose*100/sum(nptsdose),1)

  n.tox.patient <- round(apply(toxicity,2,mean),digits=2)
  names(n.tox.patient) <- dose

  n.eff.patient <- round(apply(efficacy,2,mean),digits=2)
  names(n.eff.patient) <- dose


  names(toxprob)      <- dose
  names(effprob)      <- dose
  names(target_T)     <- "Target toxicity prob."
  names(target_E)     <- "Target efficacy prob."
  names(lambda1)      <- "Lower toxicity boundary"
  names(lambda2)      <- "Upper toxicity boundary"
  names(eta1)         <- "Lower efficacy boundary"
  names(tau.T)        <- "Tox. assessment window (months)"
  names(tau.E)        <- "Eff. assessment window (months)"
  names(accrual)      <- "Accrual rate (months)"
  names(suspend)      <- "Suspension rule"


  colnames(data.obs.n) <- dose
  colnames(backfillcount) <- dose



  if(boundMTD){
    phat<-Iso::pava(toxprob)
    phat<-phat[phat<=lambda2]
    mdif <- min(abs(phat-target_T))
    mtd1m  <- intersect(which(abs(phat-target_T)==mdif),which(phat<=lambda2))
    if(length(mtd1m)>0){
      mtd1<-max(mtd1m)
    }else{
      mtd1<-NULL
    }
  }else{
    phat<-Iso::pava(toxprob)
    phat<-phat[phat<=target_T]
    mdif <- min(abs(phat-target_T))
    mtd1m  <- intersect(which(abs(phat-target_T)==mdif),which(phat<=target_T))
    if(length(mtd1m)>0){
      mtd1<-max(mtd1m)
    }else{
      mtd1<-NULL
    }
  }


  eff1 <-which(effprob>=eta1)
  if(length(eff1)>0){
    eff1min<-min(eff1)
  }else{
    eff1min<-NULL
  }

  if(length(mtd1)>0 & length(eff1min)>0){
    evadoseu <- intersect(1:mtd1,eff1min:max(dosen))
  }else{
    evadoseu<-NULL
  }


  if (length(evadoseu)>0){

    if(obd.method=="utility.weighted"){
      utility.score <- utility.weighted(probt=Iso::pava(toxprob),probe=effprob,
                              w1=w1,w2=w2,tox.upper=pT.tox)

    }else if(obd.method=="utility.truncated.linear"){
      utility.score <- utility.truncated.linear(probt=Iso::pava(toxprob),probe=effprob,
                                      tlow=plow.ast,tupp=pupp.ast,
                                      elow=qlow.ast,eupp=qupp.ast)

    }else if(obd.method=="utility.scoring"){
      utility.score <- utility.scoring(probt=Iso::pava(toxprob),probe=effprob,
                             psi00=ppsi00,psi11=ppsi11)

    }

    obd.true=min(intersect(dosen[utility.score==max(utility.score[evadoseu])],evadoseu))
    mtd.true=max(evadoseu)
  }else if(length(evadoseu)==0 & length(mtd1)!=0){
    obd.true<-0
    mtd.true<-mtd1
  }else{
    obd.true<-0
    mtd.true<-0
  }

  if(obd.true==0){
    PCS=prop.stop
    PCC=NA_integer_
    PCA=0
  }else{
    PCS=prop.select[obd.true]
    PCC=claim.select[obd.true]
    PCA=round((n.patient*100/sum(n.patient))[obd.true],1)
  }

  POS=sum(obd>mtd.true)*100/n.sim


  if (mtd.true==n.dose){POA=0}else{
    POA= round(sum(data.obs.n[,(mtd.true+1):n.dose])*100/sum(data.obs.n),digits=1)# percent of overdose allocation
  }

  names(PCS) <-"The percentage of trials that the optimal dose was correctly selected"
  names(PCC) <-"The percentage of trials that the correct claim was at a correct optimal dose"
  names(PCA) <-"The percentage of patients that were correctly allocated to the optimal dose"

  names(POS) <-"The percentage of toxic doses selection"
  names(POA) <-"The percentage of patients who were allocated to toxic doses"


  result <- list(toxprob      = toxprob,
                 effprob      = effprob,
                 phi          = target_T,
                 delta        = target_E,
                 lambda1      = lambda1,
                 lambda2      = lambda2,
                 eta1         = eta1,
                 tau.T        = tau.T,
                 tau.E        = tau.E,
                 suspend      = suspend,
                 accrual      = accrual,
                 n.patient.all = n.patient,
                 nptsdosepct.all   = nptsdosepct,
                 n.tox.patient.all= n.tox.patient,
                 n.eff.patient.all= n.eff.patient,
                 n.patient.stage2 =n.patient.stage2,
                 n.bpatient   = n.bpatient,
                 prop.select  = prop.select,
                 MTD.select   = MTD.select,
                 claim.select = claim.select,
                 prop.stop    = prop.stop,
                 #prop.stop.futility=prop.stop.futility,
                 #prop.stop.efficacy=prop.stop.efficacy,
                 duration     = duration,
                 duration1     = duration1,
                 duration2     = duration2,
                 totaln       = totaln,
                 #data.obs.earlystopfuti=data.obs.earlystopfuti,
                 #data.obs.earlystopsupe=data.obs.earlystopsupe,
                 data.obs.n   = data.obs.n,
                 data.obs.n.stage2=data.obs.n.stage2,
                 data.obs.n.stage1=data.obs.n.stage1,
                 #doselevel.stage2=doselevel.stage2,
                 obd          = obd,
                 claim=claim,
                 #backfilltimes= backfilltimes,
                 #backfillcount= backfillcount,
                 PCS=PCS,
                 PCC=PCC,
                 #PCA=PCA,
                 PTS=POS,
                 PTA=POA,
                 boundary_tab.all.out=boundary_tab.all.out,
                 lambda=lambda,
                 gamma=gamma,
                 eta=eta
  )

  if(stage1.method=="PW"){
    class(result) <- "BF-BOIN-ET-U"
  }else if(stage1.method=="ER"){
    class(result) <- "BF-BOIN-ET-E"
  }
  return(result)

}



