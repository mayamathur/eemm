
# PRELIMINARIES -------------------------------------------------------

# This script uses renv to preserve the R environment specs (e.g., package versions.)
library(renv)
# run this if you want to reproduce results using the R environment we had:
# renv::restore()

library(EValue)
library(MetaUtility)
library(msm)
library(here)
library(testthat)
library(ggplot2)
library(dplyr)
library(tidyr)

# run this only if you want to update the R environment specs
# renv::snapshot()

# no sci notation
options(scipen=999)

# directory to save results
results.dir = here("Applied example/Results from R")

# for writing results to Overleaf paper
overleaf.dir = "~/Dropbox/Apps/Overleaf/EEMM: E-values for effect modification and interaction/R_objects"

# get helper fnss
setwd( here("Applied example/Code") )
source("helper.R")


# ENTER LETTENEUR DATA (CELL COUNTS) ----------------------

# only used for additive EMM

# L = low education (1=low vs. 0=high, excluding medium)
# G = sex (1=women, 0=men)
# Y = dementia
# C = age, age^2, study ID (because this was a pooled analysis)

# ~ Women ----------------------
# cell counts are from Table 1 (totals) and Table 4 (cases)
# women

# from Table 1:
# n with L=1 and with L=0 for women
nw_1 = 2988
nw_0 = 364

dw = data.frame(  Y = c(1, 1, 0, 0),
                  L = c(1, 0, 1, 0),
                  n = c( 158, 6, nw_1-158, nw_0-6 ) )

# P(Y = 1 | L=1) and P(Y = 1 | L=0) for women
( pw_1 = dw$n[ dw$L == 1 & dw$Y == 1 ] / sum(dw$n[ dw$L == 1 ]) )
( pw_0 = dw$n[ dw$L == 0 & dw$Y == 1 ] / sum(dw$n[ dw$L == 0 ]) )

# a bit different from 3.78 because mine is crude 
( crudeRRw = pw_1/pw_0 )

# HR using person-years from table: 3.01 (still different)
(158/6920) / (6/792)

# prevalence of low education among women
fw = nw_1 / (nw_1 + nw_0)

# total N for women
update_result_csv( name = "n women",
                   value = nw_1 + nw_0 )


# ~ Men ----------------------
# cell counts are from Table 1 (totals) and Table 4 (cases)
# men

# from Table 1:
# n with L=1 and with L=0 for women
nm_1 = 1790
nm_0 = 605

dm = data.frame(  Y = c(1, 1, 0, 0),
                  L = c(1, 0, 1, 0),
                  n = c( 64, 17, nm_1-64, nm_0-17 ) )

# P(Y = 1 | L=1) and P(Y = 1 | L=0) for mean
( pm_1 = dm$n[ dm$L == 1 & dm$Y == 1 ] / sum(dm$n[ dm$L == 1 ]) )
( pm_0 = dm$n[ dm$L == 0 & dm$Y == 1 ] / sum(dm$n[ dm$L == 0 ]) )

# a bit different from 1.09 because mine is crude 
( crudeRRm = pm_1/pm_0 )

# prevalence of low education among men
fm = nm_1 / (nm_1 + nm_0)

# total N for men
update_result_csv( name = "n men",
                   value = nm_1 + nm_0 )


# MULTIPLICATIVE EMM ----------------------

# checked this section 2021-3-8

# ~ Confounded point estimates and inference ----------------------

# *adjusted* RR from Table 1
RRw = 3.78
RRm = 1.09
( RRc = RRw/RRm )

# approximate CI limits for RRc since paper only gives CI limits by stratum
# these are on log-RR scale
VarLogRRw = scrape_meta(est = RRw, hi = 8.72 )$vyi
VarLogRRm = scrape_meta(est = RRm, hi = 1.94 )$vyi

VarLogRRc = VarLogRRw + VarLogRRm


# write to results
# sanity check: they report p=0.02 for interaction (pg 4)
resRRc = write_est_inf( est = log(RRc), 
                        var = VarLogRRc,
                        prefix = "RRc",
                        takeExp = TRUE )

# ~ E-values ----------------------

# ~~ Non-monotonic confounding ----------------------
# take square roots
( Emult = evalue( RR( sqrt(RRc) ),
                  lo = sqrt( resRRc$lo ) )["E-values", c("point", "lower") ] )


# ~~ Monotonic confounding ----------------------
# (regular E-value transformation)
( Emult.mono = evalue( RR(RRc),
                       lo = resRRc$lo )["E-values", c("point", "lower") ] )


update_result_csv( name = "RRc est evalue",
                   value = round( Emult[1], 2) )
update_result_csv( name = "RRc lo evalue",
                   value = round( Emult[2], 2) )

update_result_csv( name = "RRc est evalue mono",
                   value = round( Emult.mono[1], 2) )
update_result_csv( name = "RRc lo evalue mono",
                   value = round( Emult.mono[2], 2) )




# ADDITIVE EMM ----------------------


# ~ Confounded point estimates and inference ----------------------

# we don't have adjusted probabilities, so instead use 
#  crude ones

# get confounded estimates for each stratum and for EMM
# by passing bias factor of 1
RDs = RDt_bound( pw_1 = pw_1,
                 pw_0 = pw_0,
                 nw_1 = nw_1,
                 nw_0 = nw_0,
                 fw = fw,
                 biasDir_w = "positive",  # irrelevant because maxB = 0
                 maxB_w = 1,
                 
                 pm_1 = pm_1,
                 pm_0 = pm_0,
                 nm_1 = nm_1,
                 nm_0 = nm_0,
                 fm = fm,
                 biasDir_m = "positive",  # irrelevant because maxB = 0
                 maxB_m = 1 )

RDc = RDs$RD[ RDs$stratum == "effectMod" ]
RDw = RDs$RD[ RDs$stratum == "1" ]
RDm = RDs$RD[ RDs$stratum == "0" ]


# inference for risk differences
VarRDw = var_RD(p1 = pw_1,
                p0 = pw_0,
                n1 = nw_1,
                n0 = nw_0 )

VarRDm = var_RD(p1 = pm_1,
                p0 = pm_0,
                n1 = nm_1,
                n0 = nm_0 )

( VarRDc = VarRDw + VarRDm )


# write results for each stratum and for EMMs
resRDw = write_est_inf( est = RDw, 
                        var = VarRDw,
                        prefix = "RDw",
                        takeExp = FALSE )

resRDm = write_est_inf( est = RDm, 
                        var = VarRDm,
                        prefix = "RDm",
                        takeExp = FALSE )

resRDc = write_est_inf( est = RDc, 
                        var = VarRDc,
                        prefix = "RDc",
                        takeExp = FALSE )


# ~ E-values ----------------------

# ~~ Non-monotonic confounding ----------------------

### point estimate
( Eadd.est = IC_evalue( stratum = "effectMod",
                        varName = "RD",
                        true = 0,
                        monotonicBias = "no",
                        
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


### CI limit
( Eadd.CI = IC_evalue( stratum = "effectMod",
                       varName = "lo",
                       true = 0,
                       monotonicBias = "no",
                       
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

update_result_csv( name = "RDc est evalue",
                   value = round( Eadd.est$evalue, 2) )
update_result_csv( name = "RDc lo evalue",
                   value = round( Eadd.CI$evalue, 2) )


# ~~ Monotonic confounding ----------------------

### point estimate
( Eadd.est.mono = IC_evalue_outer( varName = "RD" ) )

# which bias direction is the winner (minimizes the E-value)?
Eadd.est.mono$evalueBiasDir

update_result_csv( name = "RDc est evalue mono",
                   value = round( Eadd.est.mono$evalue, 2) )

### CI limit
( Eadd.CI.mono = IC_evalue_outer( varName = "lo" ) )

# which bias direction is the winner (minimizes the E-value)?
Eadd.est.mono$evalueBiasDir

update_result_csv( name = "RDc lo evalue mono",
                   value = round( Eadd.CI.mono$evalue, 2) )



# MISC ----------------------


# ~ Plot RDt bound as a function of bias factor ----------------------

# sanity check: plot the bound as fn of max
library(ggplot2)

dp = data.frame( bias = seq(1,5,.001) )

dp = dp %>% rowwise() %>%
  mutate( boundW = RDt_bound( pw_1 = pw_1,
                              pw_0 = pw_0,
                              nw_1 = nw_1,
                              nw_0 = nw_0,
                              fw = fw,
                              
                              pm_1 = pm_1,
                              pm_0 = pm_0,
                              nm_1 = nm_1,
                              nm_0 = nm_0,
                              fm = fm,
                              
                              .max = bias )[1,"RD"],
          
          boundM = RDt_bound( pw_1 = pw_1,
                              pw_0 = pw_0,
                              nw_1 = nw_1,
                              nw_0 = nw_0,
                              fw = fw,
                              
                              pm_1 = pm_1,
                              pm_0 = pm_0,
                              nm_1 = nm_1,
                              nm_0 = nm_0,
                              fm = fm,
                              
                              .max = bias )[2,"RD"],
          
          boundEMM = RDt_bound( pw_1 = pw_1,
                                pw_0 = pw_0,
                                nw_1 = nw_1,
                                nw_0 = nw_0,
                                fw = fw,
                                
                                pm_1 = pm_1,
                                pm_0 = pm_0,
                                nm_1 = nm_1,
                                nm_0 = nm_0,
                                fm = fm,
                                
                                .max = bias )[3,"RD"] )

xmax = max(dp$bias)
#xmax = 1.5

#@fix this
# reshape for plotting joy
dp2 = dp %>% pivot_longer( !bias,
                           names_to = "Estimate",
                           values_to = "Value" )


# **here, boundEMM is when bias is non-monotonic confounding
#  and boundM and boundW are monotonic confounding (i.e., that stratum is the only one affected)
ggplot( data = dp2 ) +
  
  # line for the observed RDc
  geom_hline( yintercept = RDc, 
              lty = 2,
              color = "red") +
  geom_line(size = 1,  aes( x = bias,
                            y = Value,
                            color = Estimate)) +
  
  theme_classic() +
  
  scale_x_continuous( limits = c(1, xmax ),
                      breaks = seq( 1, xmax, .5) ) +
  
  # scale_y_continuous( limits = c(1,max(dp$bound)),
  #                     breaks = seq(1,max(dp$bound),.5)) +
  
  #xlab(bquote(RR[XY]^c)) +
  xlab("Bias factor") +
  ylab("Bound on RDt")



# 
# 
# 
# # plot it along with the normal E-value
# library(ggplot2)
# 
# 
# dp = data.frame( RR = seq(1,5,.001) )
# dp$evalue = dp$RR + sqrt( dp$RR * (dp$RR - 1) )
# dp$evalueInt = sqrt(dp$RR) + sqrt( sqrt(dp$RR) * ( sqrt(dp$RR) - 1) )
# 
# ggplot( data = dp ) +
#   
#   geom_abline(intercept = 0,
#               slope = 1,
#               lty = 2,
#               color = "gray") +
#   
#   geom_line(size = 1,  aes( x = RR,
#                             y = evalue)) +
#   geom_line(size = 1,  aes( x = RR,
#                             y = evalueInt)) +
#   theme_classic() +
#   
#   scale_x_continuous( limits = c(1, max(dp$RR)),
#                       breaks = seq(1,max(dp$RR),.5)) +
#   
#   scale_y_continuous( limits = c(1,5),
#                       breaks = seq(1,5,.5)) +
#   
#   xlab(bquote(RR[XY]^c)) + 
#   ylab("E-value")
# 
# 










