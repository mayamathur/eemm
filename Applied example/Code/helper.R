# SMALL STATS FUNCTIONS ----------------------

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


# BIG STATS FUNCTIONS ----------------------


# quick variance of proportion
var_prop = function(p, n) ( p * (1 - p) ) / n

#@same for lower and upper?
# bias-corrected variance for one stratum
# Ding & VanderWeele, eAppendix 5 (delta method)
RDt_var = function(f, p1, p0, n1, n0, .maxB) {
  
  f.var = var_prop(p = f, n = n1 + n0)
  p1.var = var_prop(p = p1, n = n1)
  p0.var = var_prop(p = p0, n = n0)
  
  term1 = p1.var + p0.var * .maxB^2
  term2 = ( f + (1 - f) / .maxB )^2
  term3 = ( p1 - p0 * .maxB )^2 * ( 1 - 1/.maxB )^2 * f.var
  
  return( term1 * term2 + term3 )
}


# corrected RD for a given amount of bias
# @needs to accept argument "true"
# .maxB: the bias factor, NOT the Evalue!
RDt_bound = function( pw_1,
                      pw_0,
                      nw_1,
                      nw_0,
                      fw,
                      
                      pm_1,
                      pm_0,
                      nm_1,
                      nm_0,
                      fm,
                      
                      alpha = 0.05,
                      
                      .maxB ) {
  
  #browser()
  # TEST ONLY
  #.maxB = 1
  

  # c.f. R package:
  # est.BF = (sqrt((true + diff)^2 + 4 * p1 * p0 * f * (1 - 
  #                                                       f)) - (true + diff))/(2 * p0 * f)
  # est.Evalue = threshold(est.BF)
  
  ### Corrected point estimate
  # corrected RD for X=1 (women) stratum
  RDtW = ( pw_1 - pw_0 * .maxB ) * ( fw + ( 1 - fw ) / .maxB )
  # corrected RD for X=0 (men) stratum
  RDtM = ( pm_1 * .maxB - pm_0 ) * ( fm + ( 1 - fm ) / .maxB )
  
  RDt = RDtW - RDtM
  
  # sanity check
  # should recover uncorrected RDs when there is no bias
  if ( .maxB == 1 ) expect_equal(RDtW, pw_1 - pw_0)
  if ( .maxB == 1 ) expect_equal(RDtM, pm_1 - pm_0)
  
  ### Corrected confidence interval
  # calculate var for each stratum (W and M) separately
  # as in Ding & VanderWeele, eAppendix 5 (delta method)
  VarRDtW = RDt_var(f = fw,
                    p1 = pw_1,
                    p0 = pw_0,
                    n1 = nw_1,
                    n0 = nw_0,
                    .maxB = .maxB )
  
  VarRDtM = RDt_var(f = fm,
                    p1 = pm_1,
                    p0 = pm_0,
                    n1 = nm_1,
                    n0 = nm_0,
                    .maxB = .maxB )
  
  VarRDt = VarRDtW + VarRDtM
  
  
  ### Organize and return results
  res = data.frame( stratum = c("1", "0", "effectMod"),
                    
                    RD = c( RDtW, RDtM, RDt),
                    
                    se = c( sqrt(VarRDtW), sqrt(VarRDtM), sqrt(VarRDt) ) )
  
  crit = qnorm( 1 - alpha/2 )
  res$lo = res$RD - crit * res$se
  res$hi = res$RD - crit * res$se
  res$pval = 2 * ( 1 - pnorm( abs( res$RD / res$se ) ) )
  
  return(res)
}


# returns distance from RDt
# var: "RD" or "lo" for which to set to 0
# ...: args to be passed to RDt_bound
RD_distance = function(stratum,
                       varName,
                       true,
                       ...) {
  .RDs = RDt_bound(...)
  
  abs( .RDs[[varName]][ RDs$stratum == stratum ] - true )
}


# RD_distance( .maxB = 1,
# 
#              stratum = "effectMod",
#              varName = "RD",
# 
#              pw_1 = pw_1,
#              pw_0 = pw_0,
#              nw_1 = nw_1,
#              nw_0 = nw_0,
#              fw = fw,
# 
#              pm_1 = pm_1,
#              pm_0 = pm_0,
#              nm_1 = nm_1,
#              nm_0 = nm_0,
#              fm = fm,
# 
#              true = 0,
#              alpha = alpha )


RDc_evalue = function( stratum,
                       varName,
                       true = 0,
                       
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
  
  #@SAVE THIS IN R NOTES!
  # HOW TO PASS ALL (NAMED) ARGS TO ANOTHER FN:
  # https://stackoverflow.com/questions/29327423/use-match-call-to-pass-all-arguments-to-other-function
  # "-1" removes the name of the fn that was called ("RDc_evalue")
  .args = as.list(match.call()[-1])
  
  # # test match.call situation
  # .args$.maxB = 5
  # do.call(RD_distance, .args)
  
  #@ revisit upper bound of search space (500)
  opt = optimize( f = function(x){ .args$.maxB = x
                                    do.call( RD_distance, .args ) },
  interval = c(0, 500),
  maximum = FALSE )
  
  return( data.frame( evalue = g(opt$minimum),
                      bound = opt$objective ) )
  
}



( evalueEst = RDc_evalue( stratum = "1",
                       varName = "RD",
                       true = 0,
                       
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
                       
                       alpha = 0.05 ) )


( evalueCI = RDc_evalue( stratum = "1",
                          varName = "lo",
                          true = 0,
                          
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
                          
                          alpha = 0.05 ) )


# now try against package:
evalueOld = evalues.RD( n11 = dw$n[ dw$L == 1 & dw$Y == 1 ],
            n10 = dw$n[ dw$L == 1 & dw$Y == 0 ],
            n01 = dw$n[ dw$L == 0 & dw$Y == 1 ],
            n00 = dw$n[ dw$L == 0 & dw$Y == 0 ],
            true = 0,
            alpha = 0.05)


# they agree! :D
# woohoo!!!!!
evalueOld$est.Evalue; evalueEst$evalue
evalueOld$lower.Evalue; evalueCI$evalue

# # from package
# evalues.RD = function( n11, n10, n01, n00,  
#                        true = 0, alpha = 0.05, grid = 0.0001, ... ) {
#   
#   # sanity check
#   if ( any( c(n11, n10, n01, n00) < 0 ) ) stop("Negative cell counts are impossible.")
#   
#   # sample sizes
#   N = n10 + n11 + n01 + n00
#   N1 = n10 + n11  # total X=1
#   N0 = n00 + n01  # total X=0
#   
#   # compute f = P(X = 1)
#   f = N1 / N
#   
#   # P(D = 1 | X = 1)
#   p1  = n11 / N1
#   
#   # P(D = 1 | X = 0)
#   p0  = n01 / N0
#   
#   if( p1 < p0 ) stop("RD < 0; please relabel the exposure such that the risk difference > 0.")
#   
#   
#   # standard errors
#   se.p1 = sqrt( p1 * ( 1-p1 ) / N1 )
#   se.p0 = sqrt( p0 * ( 1-p0 ) / N0 )
#   
#   # back to Peng's code
#   s2.f   = f*( 1-f )/N
#   s2.p1  = se.p1^2
#   s2.p0  = se.p0^2
#   diff   = p0*( 1-f ) - p1*f
#   
#   # bias factor and E-value for point estimate
#   est.BF = ( sqrt( ( true + diff )^2 + 4 * p1 * p0 * f * ( 1-f )  ) - ( true + diff ) ) / ( 2 * p0 * f )
#   est.Evalue    = threshold(est.BF)   
#   if( p1 - p0 <= true ) stop("For risk difference, true value must be less than or equal to point estimate.")
#   
#   # compute lower CI limit
#   Zalpha        = qnorm( 1-alpha/2 )  # critical value
#   lowerCI       = p1 - p0 - Zalpha*sqrt( s2.p1 + s2.p0 )
#   
#   # check if CI contains null
#   if ( lowerCI <= true ) {
#     
#     # warning( "Lower CI limit of RD is smaller than or equal to true value." )
#     return( list( est.Evalue = est.Evalue, lower.Evalue = 1 ) )
#     
#   } else {
#     # find E-value for lower CI limit
#     # we know it's less than or equal to E-value for point estimate
#     BF.search = seq( 1, est.BF, grid )
#     
#     # population-standardized risk difference
#     RD.search = p1 - p0 * BF.search
#     f.search  = f + ( 1-f )/BF.search
#     
#     # using equation for RD^true on pg 376, compute the lower CI limit for these parameters
#     # RD.search * f.search is exactly the RHS of the inequality for RD^true ( population )
#     Low.search = RD.search * f.search -
#       Zalpha * sqrt( ( s2.p1 + s2.p0 * BF.search^2 ) * f.search^2 +
#                        RD.search^2 * ( 1 - 1 / BF.search )^2 * s2.f )
#     
#     # get the first value for BF_u such that the CI limit hits the true value
#     Low.ind    = ( Low.search <= true )
#     Low.no     = which( Low.ind==1 )[1]
#     lower.Evalue = threshold( BF.search[Low.no] )
#     
#     
#     return(list(est.Evalue   = est.Evalue,
#                 lower.Evalue = lower.Evalue))
#   }
#   
# }


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
