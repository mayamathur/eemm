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

# inverse of g (E-value -> bias factor)


# BIG STATS FUNCTIONS ----------------------


# quick variance of proportion
var_prop = function(p, n) ( p * (1 - p) ) / n

#@same for lower and upper?
# bias-corrected variance for one stratum
# Ding & VanderWeele, eAppendix 5 (delta method)
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

# # sanity check for symmetry
# v1 = RDt_var( f = .25,
#               p1 = 0.3,
#               p0 = 0.1,
#               n1 = 50,
#               n0 = 100,
#               .maxB = 2 )
# 
# # here I've changed argument order to reverse sign of RD
# v2 = RDt_var( f = 1-.25,
#               p0 = 0.3,
#               p1 = 0.1,
#               n0 = 50,
#               n1 = 100,
#               .maxB = 2 )
# 
# expect_equal(v1, v2)


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
                      
                      pm_1,
                      pm_0,
                      nm_1,
                      nm_0,
                      fm,
                      
                      alpha = 0.05,
                       
                      .maxB ) {
  
  #if ( pw_1 - pw_0 < 0 ) stop("Preventive RDw case not handled")
  #if ( pm_1 - pm_0 > 0 ) stop("Causative RDm case not handled")
  
  ### Corrected point estimate
  # corrected RD for X=1 (women) stratum (shift upward) - pg 376
  RDtW = ( pw_1 - pw_0 * .maxB ) * ( fw + ( 1 - fw ) / .maxB )
  # corrected RD for X=0 (men) stratum (shift downward) - pg 376
  RDtM = ( pm_1 * .maxB - pm_0 ) * ( fm + ( 1 - fm )/.maxB )  # without recoding f
  
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

# # sanity check for symmetry
# # RDw and RDm should match here by symmetry
# x1 = RDt_bound( pw_1 = 0.6,
#                 pw_0 = 0.4,
#                 nw_1 = 100,
#                 nw_0 = 10,
#                 fw = 0.25,
#                 
#                 pm_1 = 0.4,
#                 pm_0 = 0.6,
#                 nm_1 = 10,
#                 nm_0 = 100,
#                 fm = .25,
#                 
#                 alpha = 0.05,
#                 
#                 .maxB = 2 )
# 
# expect_equal( x1$RD[1], -x1$RD[2] )
# expect_equal( x1$se[1], x1$se[2] )



# NOT IN USE:
# p1 = 0.6
# p0 = 0.4
# f = 0.25
# n = 1000
# X1 = rbinom( n = n, size = 1, prob = f )
# Y = rep(NA, n)
# Y[ X1 == 1 ] = rbinom( n = sum(X1 == 1), size = 1, prob = p1 )
# Y[ X1 == 0 ] = rbinom( n = sum(X1 == 0), size = 1, prob = p0 )
# 
# # risk difference
# mean(Y[ X1 == 1 ]) - mean(Y[ X1 == 0 ]); p1 - p0
# 
# X2 = (X1 == 0)  # recoded exposure
# 
# # bm: Come back. Not understanding the recoding issue...
# RDt_bound( pw_1 = mean(Y[ X1 == 1 ]),
#            pw_0 = mean(Y[ X1 == 0 ]),
#            nw_1 = sum(X1==1),
#            nw_0 = sum(X1==0),
#            fw = mean(X1==1),
#            
#            pm_1 = mean(Y[ X2 == 1 ]),
#            pm_0 = mean(Y[ X2 == 0 ]),
#            nm_1 = sum(X2==1),
#            nm_0 = sum(X2==0),
#            fm = mean(X2==1),
#            
#            alpha = 0.05,
#            
#            .maxB = 2 )
# 










# returns distance of RDt for a given bias factor from desired true value ("true")
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
  
  # prepare to pass all arguments to another fn
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
                      bias = opt$minimum,
                      bound = opt$objective ) )
  
}


# # ~ Sanity checks to save -----------------
# ( evalueEst = RDc_evalue( stratum = "1",
#                        varName = "RD",
#                        true = 0,
#                        
#                        pw_1,
#                        pw_0,
#                        nw_1,
#                        nw_0,
#                        fw,
#                        
#                        pm_1,
#                        pm_0,
#                        nm_1,
#                        nm_0,
#                        fm,
#                        
#                        alpha = 0.05 ) )
# 
# 
# ( evalueCI = RDc_evalue( stratum = "1",
#                           varName = "lo",
#                           true = 0,
#                           
#                           pw_1,
#                           pw_0,
#                           nw_1,
#                           nw_0,
#                           fw,
#                           
#                           pm_1,
#                           pm_0,
#                           nm_1,
#                           nm_0,
#                           fm,
#                           
#                           alpha = 0.05 ) )
# 
# 
# # now try against package:
# evalueOld = evalues.RD( n11 = dw$n[ dw$L == 1 & dw$Y == 1 ],
#             n10 = dw$n[ dw$L == 1 & dw$Y == 0 ],
#             n01 = dw$n[ dw$L == 0 & dw$Y == 1 ],
#             n00 = dw$n[ dw$L == 0 & dw$Y == 0 ],
#             true = 0,
#             alpha = 0.05)
# 
# 
# # they agree! :D
# # woohoo!!!!!
# evalueOld$est.Evalue; evalueEst$evalue
# evalueOld$lower.Evalue; evalueCI$evalue
# 
# # end sanity checks



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
