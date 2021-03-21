
# USAGE NOTES ----------------------

# Algorithm structure:
# IC_evalue > RD_distance > RDt_bound < RDt_var

# If you eventually move these fns to package, note that:
#  - Will need to think through fns' assumptions about signs 
#    (e.g., RDc < 0 case)
#  - Have not dealt with weird cases like if IC^c is already less than IC^t 
#    IC_evalue will probably complain about sqrt(RR * (RR - 1)) being undef'nd in that case
#  - Search "#@" for notes on other assumptions in fns that would need to be gnlz'd  

# BIG STATS FUNCTIONS ----------------------

# bias-corrected variance for one stratum
# Ding & VanderWeele, eAppendix 5 (delta method)
# handles either positive or negative bias
#@think about the fact that recoding depends on confounding rather than true p1, p0
RDt_var = function(f, p1, p0, n1, n0, .maxB) {
  
  #@assumes that we always consider bias away from null 
  # if risk difference < 0, reverse strata coding in order to
  #  apply same bound
  if ( p1 - p0 < 0 ){
    .p0 = p1
    .p1 = p0
    .n0 = n1
    .n1 = n0
  } else {
    .p0 = p0
    .p1 = p1
    .n0 = n0
    .n1 = n1
  }
  
  f.var = var_prop(p = f, n = .n1 + .n0)
  p1.var = var_prop(p = .p1, n = .n1)
  p0.var = var_prop(p = .p0, n = .n0)
  
  term1 = p1.var + p0.var * .maxB^2
  term2 = ( f + (1 - f) / .maxB )^2
  term3 = ( .p1 - .p0 * .maxB )^2 * ( 1 - 1/.maxB )^2 * f.var
  
  return( term1 * term2 + term3 )
}


# corrected RD for a given amount of bias
# @needs to accept argument "true"
# .maxB: the bias factor, NOT the Evalue!
#@always shifts stratum W up and stratum M down, so assumes they are ordered
# and shifts each stratum by same bias factor .maxB
RDt_bound = function( pw_1,
                      pw_0,
                      nw_1,
                      nw_0,
                      fw,
                      maxB_w,
                      biasDir_w,
                      
                      pm_1,
                      pm_0,
                      nm_1,
                      nm_0,
                      fm,
                      maxB_m = NA,
                      biasDir_m,
                      
                      alpha = 0.05 ) {
  
  if ( is.na(maxB_m) ) {
    maxB_m = maxB_w
    message("Assuming same bias in each stratum because you didn't provide maxB_m")
  }
  
  # if ( is.na(biasDir_m) ) {
  #   biasDir_m = "negative"
  #   message("Assuming same bias in each stratum because you didn't provide maxB_m")
  # }
  
  ### Corrected point estimate
  # corrected RD for each stratum - pg 376
  if ( biasDir_w == "positive" ) RDtW = ( pw_1 - pw_0 * maxB_w ) * ( fw + ( 1 - fw ) / maxB_w )
  if ( biasDir_w == "negative" ) RDtW = ( pw_1 * maxB_w - pw_0 ) * ( fw + ( 1 - fw ) / maxB_w )
  
  if ( biasDir_m == "positive" ) RDtM = ( pm_1 - pm_0 * maxB_m ) * ( fm + ( 1 - fm ) / maxB_m )
  if ( biasDir_m == "negative" ) RDtM = ( pm_1 * maxB_m - pm_0 ) * ( fm + ( 1 - fm ) / maxB_m )
  
  # # old version (agrees numerically with the above)
  # RDtW = ( pw_1 - pw_0 * .maxB ) * ( fw + ( 1 - fw ) / .maxB )
  # # corrected RD for X=0 (men) stratum (shift downward) - pg 376
  # RDtM = ( pm_1 * .maxB - pm_0 ) * ( fm + ( 1 - fm )/.maxB )  # without recoding f
  
  RDt = RDtW - RDtM
  
  # sanity check
  # should recover uncorrected RDs when there is no bias
  if ( maxB_w == 1 ) expect_equal(RDtW, pw_1 - pw_0)
  if ( maxB_m == 1 ) expect_equal(RDtM, pm_1 - pm_0)

  
  ### Corrected confidence interval
  # calculate var for each stratum (W and M) separately
  # as in Ding & VanderWeele, eAppendix 5 (delta method)
  VarRDtW = RDt_var(f = fw,
                    p1 = pw_1,
                    p0 = pw_0,
                    n1 = nw_1,
                    n0 = nw_0,
                    .maxB = maxB_w )
  
  VarRDtM = RDt_var(f = fm,
                    p1 = pm_1,
                    p0 = pm_0,
                    n1 = nm_1,
                    n0 = nm_0,
                    .maxB = maxB_m )
  
  VarRDt = VarRDtW + VarRDtM
  
  
  ### Organize and return results
  res = data.frame( stratum = c("1", "0", "effectMod"),
                    
                    RD = c( RDtW, RDtM, RDt),
                    
                    se = c( sqrt(VarRDtW), sqrt(VarRDtM), sqrt(VarRDt) ) )
  
  crit = qnorm( 1 - alpha/2 )
  res$lo = res$RD - crit * res$se
  res$hi = res$RD + crit * res$se
  res$pval = 2 * ( 1 - pnorm( abs( res$RD / res$se ) ) )
  
  return(res)
}




# returns distance of a statistic for a given bias factor from desired true value ("true")
# varName: the statistic (of those returned by RDt_bound: "RD", "lo", "hi") whose distance from "true"
#  should be measured
# ...: args to be passed to RDt_bound

RD_distance = function(stratum,
                       varName,
                       true,
                       ...) {
  .RDs = RDt_bound(...)
  
  abs( .RDs[[varName]][ RDs$stratum == stratum ] - true )
}


# RD_distance( stratum = "effectMod",
#              varName = "RD",
#              true = 0.01,
#              
#              # everything to follow is passed to RD_distance
#              pw_1 = pw_1,
#              pw_0 = pw_0,
#              nw_1 = nw_1,
#              nw_0 = nw_0,
#              fw = fw,
#              biasDir_w = "positive",
#              maxB_w = 1,
# 
#              pm_1 = pm_1,
#              pm_0 = pm_0,
#              nm_1 = nm_1,
#              nm_0 = nm_0,
#              fm = fm,
#              biasDir_m = "positive",
#              maxB_m = 1,
# 
#              alpha = 0.05 )




# varName: as in RD_distance, lets you choose point estimate or CI limit for E-value
# monotonicBias: "no" (non-monotonic), "positive", "negative"
# does NOT consider monontonic bias in arbitrary direction; for that, need to call 
#  the wrapper IC_evalue_outer, which calls IC_evalue twice for each candidate bias direction
# also gives E-values for each stratum separately if wanted (based on varName)
IC_evalue = function( stratum,
                      varName,
                      true = 0,
                      monotonicBias = "no",
                      
                      pw_1,
                      pw_0,
                      nw_1,
                      nw_0,
                      fw,
                      
                      pm_1,
                      pm_0,
                      nm_1,
                      nm_0,
                      fm,
                      
                      alpha = 0.05 ) {
  

  # prepare to pass all arguments to another fn
  # https://stackoverflow.com/questions/29327423/use-match-call-to-pass-all-arguments-to-other-function
  # "-1" removes the name of the fn that was called ("IC_evalue")
  .args = as.list(match.call()[-1])
  # remove other args that are not to be passed to RD_distance
  .args = .args[ !names(.args) %in% c("monotonicBias") ]
  
  # # test match.call situation
  # .args$.maxB = 5
  # do.call(RD_distance, .args)
  
  ### Set up the bounding factor fn to be maximized to get the E-value
  # depends on biasDir assumptions
  #@assumes W stratum > 0 and M is < 0
  #so basically need to warn user to recode exposure if IC^c < 0 
  if ( monotonicBias == "no" ) {
     boundfn = function(x){
       .args$maxB_w = x
       .args$biasDir_w = "positive"
       .args$maxB_m = x
       .args$biasDir_m = "negative"
        do.call( RD_distance, .args )
     }
  }
  
  if ( monotonicBias == "positive" ) {
    boundfn = function(x){
      .args$maxB_w = x
      .args$biasDir_w = "positive"
      .args$maxB_m = 1  # no bias in this stratum
      .args$biasDir_m = "positive"
      do.call( RD_distance, .args )
    }
  }
  
  if ( monotonicBias == "negative" ) {
    boundfn = function(x){
      .args$maxB_w = 1 # no bias in this stratum
      .args$biasDir_w = "negative"
      .args$maxB_m = x
      .args$biasDir_m = "negative"
      do.call( RD_distance, .args )
    }
  }
  
  #@ revisit upper bound of search space (500)
  opt = optimize( f = boundfn,
                  interval = c(0, 500),
                  maximum = FALSE )
                  
  if ( abs( opt$objective - true ) > 0.001 ) warning("E-value didn't move estimate close enough to true value; look into optimize() call")

  return( data.frame( evalue = g(opt$minimum),
                      biasFactor = opt$minimum,  # not the bias factor, but the regular bias
                      bound = opt$objective ) )  # should be equal to true
  
}


# looks at monotonic bias without assuming direction by calling IC_evalue twice
# NOT in shape for package
# lots of dataset-specific things in here
#bm
IC_evalue_outer = function(varName,
                           
                           stratum = "effectMod",
                           true = 0,
                           monotonicBias = "positive",
                           
                           pw_1 = pw_1,
                           pw_0 = pw_0,
                           nw_1 = nw_1,
                           nw_0 = nw_0,
                           fw = fw,
                           
                           pm_1 = pm_1,
                           pm_0 = pm_0,
                           nm_1 = nm_1,
                           nm_0 = nm_0,
                           fm = fm,
                           
                           alpha = 0.05
                           ) {
  
  
  # prepare to pass all arguments to another fn
  # https://stackoverflow.com/questions/29327423/use-match-call-to-pass-all-arguments-to-other-function
  # "-1" removes the name of the fn that was called ("IC_evalue")
  .args = as.list(match.call()[-1])
  
  .args$stratum = 
  
  # remove other args that are not to be passed to RD_distance
  #.args = .args[ !names(.args) %in% c("monotonicBias") ]
  
  # # test match.call situation
  # .args$.maxB = 5
  # do.call(RD_distance, .args)
  
  
  # E-value candidate 1: Shift stratum W down to match stratum M
  ( cand1 = IC_evalue( stratum = "effectMod",
                       varName = varName,
                       true = 0,
                       monotonicBias = "positive",
                       
                       pw_1 = pw_1,
                       pw_0 = pw_0,
                       nw_1 = nw_1,
                       nw_0 = nw_0,
                       fw = fw,
                       
                       pm_1 = pm_1,
                       pm_0 = pm_0,
                       nm_1 = nm_1,
                       nm_0 = nm_0,
                       fm = fm,
                       
                       alpha = 0.05 ) )
  
  # E-value candidate 2: Shift stratum M up to match stratum W
  ( cand2 = IC_evalue( stratum = "effectMod",
                       varName = varName,
                       true = 0,
                       monotonicBias = "negative",
                       
                       pw_1 = pw_1,
                       pw_0 = pw_0,
                       nw_1 = nw_1,
                       nw_0 = nw_0,
                       fw = fw,
                       
                       pm_1 = pm_1,
                       pm_0 = pm_0,
                       nm_1 = nm_1,
                       nm_0 = nm_0,
                       fm = fm,
                       
                       alpha = 0.05 ) )
  
  # Choose candidate E-value that is smaller
  winner = min(cand1$evalue, cand2$evalue)
  
  return( list( evalue = winner,
                evalueBiasDir = ifelse( winner == cand1$evalue, "positive", "negative"),
                candidates = data.frame( biasDir = c("positive", "negative"),
                                         evalue = c(cand1$evalue, cand2$evalue),
                                         biasFactor = c(cand1$biasFactor, cand2$biasFactor),
                                         isMin = c(cand1$evalue == winner, cand2$evalue == winner) ) ) )
}



# SMALL STATS FUNCTIONS ----------------------

# quick variance of proportion
var_prop = function(p, n) ( p * (1 - p) ) / n


quick_ci = function( est, var ) {
  c( est - qnorm(.975) * sqrt(var),
     est + qnorm(.975) * sqrt(var) )
}

quick_pval = function( est, var ) {
  2 * ( 1 - pnorm( abs( est / sqrt(var) ) ) )
}



var_RD = function(p1, p0, n1, n0) {
  num = ( p1 * ( 1 - p1 ) ) / n1
  denom = ( p0 * ( 1 - p0 ) ) / n0
  num + denom
}

# basic E-value transformation
g = function(RR) {
  RR + sqrt( RR * (RR - 1) )
}



# FORMATTING ----------------------

# writes estimate, CI, and pval to results file
# and returns them as a df
#**if using RRs, must pass log-RR, not RR itself
write_est_inf = function( est, var, prefix, takeExp ) {
  
  CIs = quick_ci( est = est, var = var )
  pval = quick_pval( est = est, 
                     var = var)
  
  transf = function(x) x
  if (takeExp == TRUE ) transf = function(x) exp(x)
  
  # write results
  update_result_csv( name = paste( prefix, "est" ),
                     value = round( transf(est), 2) )
  
  update_result_csv( name = paste( prefix, "lo" ),
                     value = round( transf(CIs[1]), 2) )
  
  update_result_csv( name = paste( prefix, "hi" ),
                     value = round( transf(CIs[2]), 2) )
  
  update_result_csv( name = paste( prefix, "pval" ),
                     value = format.pval( pval, eps=0.0001) )
  
  # also return results as (unrounded) df
  res = data.frame( est = transf(est),
                    lo = transf(CIs[1]),
                    hi = transf(CIs[2]),
                    pval = pval )
}

# for reproducible manuscript-writing
# adds a row to the file "stats_for_paper" with a new statistic or value for the manuscript
# optionally, "section" describes the section of code producing a given result
update_result_csv = function( name,
                              section = NA,
                              value = NA,
                              print = FALSE ) {
  setwd(results.dir)
  
  new.rows = data.frame( name,
                         value = as.character(value),
                         section = as.character(section) )
  
  # to avoid issues with variable types when overwriting
  new.rows$name = as.character(new.rows$name)
  new.rows$value = as.character(new.rows$value)
  new.rows$section = as.character(new.rows$section)
  
  
  if ( "stats_for_paper.csv" %in% list.files() ) {
    .res = read.csv( "stats_for_paper.csv",
                     stringsAsFactors = FALSE,
                     colClasses = rep("character", 3 ) )
    
    # if this entry is already in the results file, overwrite the
    #  old one
    if ( all(name %in% .res$name) ) .res[ .res$name %in% name, ] = new.rows
    else .res = rbind(.res, new.rows)
  }
  
  if ( !"stats_for_paper.csv" %in% list.files() ) {
    .res = new.rows
  }
  
  write.csv( .res, 
             "stats_for_paper.csv",
             row.names = FALSE,
             quote = FALSE )
  
  # also write to Overleaf
  setwd(overleaf.dir)
  write.csv( .res, 
             "stats_for_paper.csv",
             row.names = FALSE,
             quote = FALSE )
  
  if ( print == TRUE ) {
    View(.res)
  }
}


# stands for "wipe results"
wr = function(){
  setwd(results.dir)
  if( "stats_for_paper.csv" %in% list.files() ) system("rm stats_for_paper.csv")
  setwd(overleaf.dir)
  if( "stats_for_paper.csv" %in% list.files() ) system("rm stats_for_paper.csv")
}

# stands for "view results"
vr = function(){
  setwd(results.dir)
  View( read.csv("stats_for_paper.csv") )
}



