#' @title select_mtd
#'
#' @description Obtain the maximum tolerated dose (MTD) of the backfill 
#' bayesian optimal interval design using efficacy and toxicity outcomes.
#'
#' @param target Target toxicity probability. The default value is \code{target_T=0.3}.
#' @param npts The number of patients enrolled at each dose level.
#' @param ntox Number of patients with dose limiting toxicity (DLT).
#' @param cutoff.eli The cutoff to eliminate an overly toxic dose for safety.
#' We recommend the default value of (\code{cutoff.eli=0.95}) for general use.
#' @param extrasafe Set \code{extrasafe=TRUE} to impose a more stringent stopping rule.
#' The default value is \code{extrasafe=FALSE}.
#' @param offset A small positive number (between 0 and 0.5) to control how strict the
#' stopping rule is when \code{extrasafe=TRUE}. A larger value leads to a more
#' strict stopping rule. The default value \code{offset=0.05} generally works well.
#' @param p.tox the lowest toxicity probability that is deemed overly toxic such
#' that deescalation is required. The default value is \code{p.tox=1.4*target}.
#' @param boundMTD set \code{boundMTD=TRUE} to impose the condition: the isotonic 
#' estimate of toxicity probability for the selected MTD must be less than 
#' de-escalation boundary. The default value is \code{boundMTD=FALSE}.
#' @param n.elimination The sample size cutoff for elimination. 
#' The default is \code{n.elimination=3}.
#' 
#' @return \code{select_mtd()} returns the selected dose.
#' 
#' @references
#'1. Liu S. and Yuan, Y. (2015). Bayesian optimal interval designs for phase I clinical trials, Journal of the Royal Statistical Society: Series C , 64, 507-523.
#'2. Yuan, Y., Hess, K. R., Hilsenbeck, S. G., & Gilbert, M. R. (2016). Bayesian optimal interval design: a simple and well-performing design for phase I oncology trials. Clinical Cancer Research, 22(17), 4291-4301.
#'3. Zhou, H., Yuan, Y., & Nie, L. (2018). Accuracy, safety, and reliability of novel phase I trial designs. Clinical Cancer Research, 24(18), 4357-4364.
#'4. Zhou, Y., Lin, R., Kuo, Y. W., Lee, J. J., & Yuan, Y. (2021). BOIN Suite: A Software Platform to Design and Implement Novel Early-Phase Clinical Trials. JCO Clinical Cancer Informatics, 5, 91-101.
#'5. Takeda K, Xia Q, Liu S, Rong A. TITE-gBOIN: Time-to-event Bayesian optimal interval design to accelerate dose-finding accounting for toxicity grades. Pharm Stat. 2022 Mar;21(2):496-506. doi: 10.1002/pst.2182. Epub 2021 Dec 3. PMID: 34862715.
#'6. Yuan, Y., Lin, R., Li, D., Nie, L. and Warren, K.E. (2018). Time-to-event Bayesian Optimal Interval Design to Accelerate Phase I Trials. Clinical Cancer Research, 24(20): 4921-4930.
#'7. Rongji Mu, Ying Yuan, Jin Xu, Sumithra J. Mandrekar, Jun Yin, gBOIN: A Unified Model-Assisted Phase I Trial Design Accounting for Toxicity Grades, and Binary or Continuous End Points, Journal of the Royal Statistical Society Series C: Applied Statistics, Volume 68, Issue 2, February 2019, Pages 289–308, https://doi.org/10.1111/rssc.12263.
#'8. Lin R, Yuan Y. Time-to-event model-assisted designs for dose-finding trials with delayed toxicity. Biostatistics. 2020 Oct 1;21(4):807-824. doi: 10.1093/biostatistics/kxz007. PMID: 30984972; PMCID: PMC8559898.
#'9. Hsu C, Pan H, Mu R (2022). _UnifiedDoseFinding: Dose-Finding Methods for Non-Binary Outcomes_. R package version 0.1.9, <https://CRAN.R-project.org/package=UnifiedDoseFinding>.
#' @examples
#' target<-0.3
#' y<-c(0,0,1,2,3,0)
#' n<-c(3,3,6,9,9,0)
#' select_mtd(target=target,npts=n,ntox=y)
#' 
#' @importFrom stats pbeta qbeta rexp rmultinom runif var
#' @export


select_mtd <- function (target=0.3, npts, ntox, cutoff.eli = 0.95, extrasafe = FALSE,
                        offset = 0.05, p.tox=1.4*target, boundMTD = FALSE, n.elimination=3)
{
  pava <- function(x, wt = rep(1, length(x))) {
    n <- length(x)
    if (n <= 1)
      return(x)
    if (any(is.na(x)) || any(is.na(wt))) {
      stop("Missing values in 'x' or 'wt' not allowed")
    }
    lvlsets <- (1:n)
    repeat {
      viol <- (as.vector(diff(x)) < 0)
      if (!(any(viol)))
        break
      i <- min((1:(n - 1))[viol])
      lvl1 <- lvlsets[i]
      lvl2 <- lvlsets[i + 1]
      ilvl <- (lvlsets == lvl1 | lvlsets == lvl2)
      x[ilvl] <- sum(x[ilvl] * wt[ilvl])/sum(wt[ilvl])
      lvlsets[ilvl] <- lvl1
    }
    x
  }
  
  lambda_d = log((1 - target)/(1 - p.tox))/log(p.tox * (1 -target)/(target * (1 - p.tox)))
  
  y = ntox
  n = npts
  ndose = length(n)
  elimi = rep(0, ndose)
  for (i in 1:ndose) {
    if (n[i] >= n.elimination) {
      if (1 - pbeta(target, y[i] + 1, n[i] - y[i] + 1) >
          cutoff.eli) {
        elimi[i:ndose] = 1
        break
      }
    }
  }
  if (extrasafe) {
    if (n[1] >= n.elimination) {
      if (1 - pbeta(target, y[1] + 1, n[1] - y[1] + 1) >
          cutoff.eli - offset) {
        elimi[1:ndose] = 1
      }
    }
  }
  if (elimi[1] == 1 || sum(n[elimi == 0]) == 0) {
    selectdose = 0
  }
  else {
    adm.set = (n != 0) & (elimi == 0)
    adm.index = which(adm.set == T)
    y.adm = y[adm.set]
    n.adm = n[adm.set]
    phat = (y.adm + 0.05)/(n.adm + 0.1)
    phat.var = (y.adm + 0.05) * (n.adm - y.adm + 0.05)/((n.adm +
                                                           0.1)^2 * (n.adm + 0.1 + 1))
    phat = pava(phat, wt = 1/phat.var)
    phat = phat + (1:length(phat)) * 1e-10
    if(boundMTD){
      if(all(phat>lambda_d)){selectdose=0}else{
        phat=phat[phat<=lambda_d]
        selectd = sort(abs(phat - target), index.return = T)$ix[1]
        selectdose = adm.index[selectd]
      }
      
    }else{
      selectd = sort(abs(phat - target), index.return = T)$ix[1]
      selectdose = adm.index[selectd]
    }
    
    
  }
  
  trtd = (n != 0)
  poverdose = pava(1 - pbeta(target, y[trtd] + 0.05, n[trtd] -
                               y[trtd] + 0.05))
  phat.all = pava((y[trtd] + 0.05)/(n[trtd] + 0.1), wt = 1/((y[trtd] +
                                                               0.05) * (n[trtd] - y[trtd] + 0.05)/((n[trtd] + 0.1)^2 *
                                                                                                     (n[trtd] + 0.1 + 1))))
  lowerCIs=pava(qbeta(0.025, y[trtd] + 0.05,n[trtd] - y[trtd] + 0.05),wt = 1/((y[trtd] +
                                                                                 0.05) * (n[trtd] - y[trtd] + 0.05)/((n[trtd] + 0.1)^2 *
                                                                                                                       (n[trtd] + 0.1 + 1))))
  upperCIs=pava(qbeta(0.975, y[trtd] + 0.05,n[trtd] - y[trtd] + 0.05),wt = 1/((y[trtd] +
                                                                                 0.05) * (n[trtd] - y[trtd] + 0.05)/((n[trtd] + 0.1)^2 *
                                                                                                                       (n[trtd] + 0.1 + 1))))
  
  A1 = A2 = A3 = A4 = NULL
  k = 1
  for (i in 1:ndose) {
    if (n[i] > 0) {
      A1 = append(A1, formatC(phat.all[k], digits = 2,
                              format = "f"))
      A2 = append(A2, formatC(lowerCIs[i], digits = 2, format = "f"))
      A3 = append(A3, formatC(upperCIs[i], digits = 2, format = "f"))
      A4 = append(A4, formatC(poverdose[k], digits = 2,
                              format = "f"))
      k = k + 1
    }
    else {
      A1 = append(A1, "----")
      A2 = append(A2, "----")
      A3 = append(A3, "----")
      A4 = append(A4, "----")
    }
  }
  p_est = data.frame(cbind(dose = 1:length(npts), phat = A1,
                           CI = paste("(", A2, ",", A3, ")", sep = "")))
  out = list(target = target, MTD = selectdose, p_est = p_est,
             p_overdose = A4)
  return(out)
}