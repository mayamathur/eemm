

# CHECK SMALLER HELPER FNS ----------------------

# ~ RDt_var ----------------------
# sanity check for symmetry
v1 = RDt_var( f = .25,
              p1 = 0.3,
              p0 = 0.1,
              n1 = 50,
              n0 = 100,
              .maxB = 2 )

# here I've changed argument order to reverse sign of RD
v2 = RDt_var( f = .25,
              p0 = 0.3,
              p1 = 0.1,
              n0 = 50,
              n1 = 100,
              .maxB = 2 )

expect_equal(v1, v2)


# ~ RDt_bound ----------------------
# sanity check: symmetry when shifting strata in opposite directions
# RDw and RDm should match here by symmetry
x1 = RDt_bound( pw_1 = 0.6,
                pw_0 = 0.4,
                nw_1 = 100,
                nw_0 = 10,
                fw = 0.25,
                maxB_w = 2,
                biasDir_w = "positive",

                pm_1 = 0.4,
                pm_0 = 0.6,
                nm_1 = 10,
                nm_0 = 100,
                fm = .25,
                biasDir_m = "negative",

                alpha = 0.05 )

expect_equal( x1$RD[1], -x1$RD[2] )
expect_equal( x1$se[1], x1$se[2] )
expect_equal( x1$lo[1], -x1$hi[2] )
expect_equal( x1$lo[2], -x1$hi[1] )
expect_equal( x1$pval[1], x1$pval[2] )

# sanity check: equality when both strata are the same and are positively biased
x1 = RDt_bound( pw_1 = 0.6,
                pw_0 = 0.4,
                nw_1 = 100,
                nw_0 = 10,
                fw = 0.25,
                maxB_w = 1.6,
                biasDir_w = "positive",

                pm_1 = 0.6,
                pm_0 = 0.4,
                nm_1 = 100,
                nm_0 = 10,
                fm = .25,
                maxB_m = 1.6,
                biasDir_m = "positive",

                alpha = 0.05 )

expect_equal( x1$RD[1], x1$RD[2] )
expect_equal( x1$se[1], x1$se[2] )
expect_equal( x1$lo[1], x1$lo[2] )
expect_equal( x1$hi[1], x1$hi[2] )
expect_equal( x1$pval[1], x1$pval[2] )


# sanity check: equality when both strata are the same and are NEGATIVELY biased
x1 = RDt_bound( pw_1 = 0.6,
                pw_0 = 0.4,
                nw_1 = 100,
                nw_0 = 10,
                fw = 0.25,
                maxB_w = 1.6,
                biasDir_w = "negative",

                pm_1 = 0.6,
                pm_0 = 0.4,
                nm_1 = 100,
                nm_0 = 10,
                fm = .25,
                maxB_m = 1.6,
                biasDir_m = "negative",

                alpha = 0.05 )

expect_equal( x1$RD[1], x1$RD[2] )
expect_equal( x1$se[1], x1$se[2] )
expect_equal( x1$lo[1], x1$lo[2] )
expect_equal( x1$hi[1], x1$hi[2] )
expect_equal( x1$pval[1], x1$pval[2] )


# ADDITIVE EMM ----------------------

# ~ Confounded point estimates and inference ----------------------

# ~~ Point estimates from RDt_bound with B=1 should match naive RDs ----------------------

# check point estimates
expect_equal( RDw,
              pw_1 - pw_0,
              158/nw_1 - 6/nw_0 )

expect_equal( RDm,
              pm_1 - pm_0,
              158/nm_1 - 6/nm_0 )

expect_equal( RDc,
              (pw_1 - pw_0) - (pm_1 - pm_0) )


# check CI limit
expect_equal( RDw - 1.96 * sqrt(VarRDw),
              resRDw$lo,
              tol = 0.001 )

# ~ E-values ----------------------

# ~~ Non-monotonic confounding ----------------------

# ~~~ Bound from RDt_bound should agree with closed form in paper ----------------------

B = 2.3
term1 = (pw_1 - pw_0 * B) * ( fw + (1 - fw) / B )
term2 = (pm_1 * B - pm_0) * ( fm + (1 - fm) / B )

mine = term1 - term2

x = RDt_bound( pw_1 = pw_1,
               pw_0 = pw_0,
               nw_1 = nw_1,
               nw_0 = nw_0,
               fw = fw,
               biasDir_w = "positive",
               maxB_w = B,
               
               pm_1 = pm_1,
               pm_0 = pm_0,
               nm_1 = nm_1,
               nm_0 = nm_0,
               fm = fm,
               biasDir_m = "negative",
               maxB_m = B )
expect_equal( x$RD[ x$stratum == "effectMod" ], mine, tol = 0.0001 )



# ~~~ E-value for 1 stratum from IC_evalue should match R package -----------------
( evalueEst = IC_evalue( stratum = "1",
                         varName = "RD",
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
                         
                         alpha = 0.05 ) )


( evalueCI = IC_evalue( stratum = "1",
                        varName = "lo",
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
                        
                        alpha = 0.05 ) )


# now try against package:
evalueOld = EValue::evalues.RD( n11 = dw$n[ dw$L == 1 & dw$Y == 1 ],
                                n10 = dw$n[ dw$L == 1 & dw$Y == 0 ],
                                n01 = dw$n[ dw$L == 0 & dw$Y == 1 ],
                                n00 = dw$n[ dw$L == 0 & dw$Y == 0 ],
                                true = 0,
                                alpha = 0.05)


# they agree! :D
# woohoo!!!!!
expect_equal( evalueOld$est.Evalue, evalueEst$evalue, tol = 0.001 )
expect_equal( evalueOld$lower.Evalue, evalueCI$evalue, tol = 0.001 )

# ~~~ E-value from IC_evalue should be the solution to RDt_bound ----------------------
# pass the same bias factor to each stratum
x = RDt_bound( pw_1 = pw_1,
               pw_0 = pw_0,
               nw_1 = nw_1,
               nw_0 = nw_0,
               fw = fw,
               biasDir_w = "positive",
               maxB_w = Eadd.est$biasFactor,
               
               pm_1 = pm_1,
               pm_0 = pm_0,
               nm_1 = nm_1,
               nm_0 = nm_0,
               fm = fm,
               biasDir_m = "negative",
               maxB_m = Eadd.est$biasFactor )
expect_equal( x$RD[ x$stratum == "effectMod" ], 0, tol = 0.0001 )



##@TEMP ONLY
#bm: see if the RRs are biased by same factor when the E-value is attained
RRw = pw_1/pw_0
RRm = pm_1/pm_0

# corrected probab

( RRwt = RRw/Eadd.est$biasFactor )
( RRmt = RRm*Eadd.est$biasFactor )

# corrected probabilities
pw0.corr = pw_1 / RRwt
pw_1 - pw0.corr

pm1.corr = pm_0 * RRmt
pm1.corr - pm_0



pw1_corr = pw_0 * (RRw/Eadd.est$biasFactor)
pw1_corr - pw_0

pm0_corr = pm_1 / (RRm*Eadd.est$biasFactor)
pm_1 - pm0_corr
# bm: think about this


# try with different stratum probabilities


##@END OF TEMP

# ~~~ E-value from IC_evalue (grid search) should match closed form in paper ----------------------

gamma = pw_0 * (1 - fw) - pw_1 * fw + pm_1 * (1 - fm) - pm_0 * fm
term1 = 1 / ( 2 * (fw * pw_0 + fm * pm_1) )
term2 = 4 * (pw_0 * fw + pm_1 * fm) * ( pw_1 * (1 - fw) + pm_0 * (1-fm) )

# bias factor
( my.Badd.est = term1 * ( sqrt( gamma^2 + term2 ) - gamma ) )
expect_equal( Eadd.est$biasFactor, my.Badd.est, tol = 0.001 )

# E-value
expect_equal( Eadd.est$evalue, g(my.Badd.est), tol = 0.001 )


# # try the initial Wolfram form
# term1 = 1 / ( 2 * (fw * pw_0 + fm * pm_1) )
# term2 = fw * pw_0 + fw * pw_1 + fm * pm_0 + fm * pm_1 - pw_0 - pm_1
# term3 = -4 * ( fw * pw_0 + fm * pm_1 ) * ( fw * pw_1 + fm * pm_0 - pw_1 - pm_0 )
# 
# mine = term1 * ( sqrt( (-term2)^2 + term3 ) + term2 )


# ~~ Monotonic confounding ----------------------

# ~~~ Evalue candidate #1 (positive bias) should successfully move RDw down to RDm ----------------------
x = RDt_bound( pw_1 = pw_1,
           pw_0 = pw_0,
           nw_1 = nw_1,
           nw_0 = nw_0,
           fw = fw,
           biasDir_w = "positive",
           maxB_w = Eadd.est.mono$candidates$biasFactor[ Eadd.est.mono$candidates$biasDir == "positive" ],

           pm_1 = pm_1,
           pm_0 = pm_0,
           nm_1 = nm_1,
           nm_0 = nm_0,
           fm = fm,
           biasDir_m = "positive",
           maxB_m = 1 )
expect_equal( x$RD[ x$stratum == "1" ], RDm, tol = 0.0001 )

# and likewise for CI limit
( x = RDt_bound( pw_1 = pw_1,
                 pw_0 = pw_0,
                 nw_1 = nw_1,
                 nw_0 = nw_0,
                 fw = fw,
                 biasDir_w = "positive",
                 maxB_w = Eadd.CI.mono$candidates$biasFactor[ Eadd.CI.mono$candidates$biasDir == "positive" ],
                 
                 pm_1 = pm_1,
                 pm_0 = pm_0,
                 nm_1 = nm_1,
                 nm_0 = nm_0,
                 fm = fm,
                 biasDir_m = "positive",
                 maxB_m = 1 ) )
expect_equal( x$lo[ x$stratum == "effectMod" ], 0, tol = 0.0001 )

# ~~~ Evalue candidate #2 (negative bias) should successfully move RDm up to RDw ----------------------
( x = RDt_bound( pw_1 = pw_1,
               pw_0 = pw_0,
               nw_1 = nw_1,
               nw_0 = nw_0,
               fw = fw,
               biasDir_w = "negative",
               maxB_w = 1,

               pm_1 = pm_1,
               pm_0 = pm_0,
               nm_1 = nm_1,
               nm_0 = nm_0,
               fm = fm,
               biasDir_m = "negative",
               maxB_m = Eadd.est.mono$candidates$biasFactor[ Eadd.est.mono$candidates$biasDir == "negative" ] ) )
expect_equal( x$RD[ x$stratum == "0" ], RDw, tol = 0.0001 )

# and likewise for CI limit
( x = RDt_bound( pw_1 = pw_1,
                 pw_0 = pw_0,
                 nw_1 = nw_1,
                 nw_0 = nw_0,
                 fw = fw,
                 biasDir_w = "negative",
                 maxB_w = 1,
                 
                 pm_1 = pm_1,
                 pm_0 = pm_0,
                 nm_1 = nm_1,
                 nm_0 = nm_0,
                 fm = fm,
                 biasDir_m = "negative",
                 maxB_m = Eadd.CI.mono$candidates$biasFactor[ Eadd.CI.mono$candidates$biasDir == "negative" ] ) )
expect_equal( x$lo[ x$stratum == "effectMod" ], 0, tol = 0.0001 )



# ~~~ E-values from IC_evalue (grid search) should match closed form in paper ----------------------

### check ??
B = 4
true = ( pm_1 * B - pm_0 ) * ( fm + (1-fm) / B )

# suggestively name terms as in quadratic formula
termA = fm * pm_1
termB = pm_1 * ( 1 - fm ) - fm * pm_0 - true
termC = -pm_0 * (1 - fm)

# this is the polynomial that is to be solved for E-value
( mine = termA * B^2 + termB * B + termC )
expect_equal( mine, 0, tol = 0.0001 )

# check E-value for negatively biased stratum in paper
#  this is the one that arises from reversing roles and signs in the existing E-value on 
#  Ding Appendix, pg 18


### check E-value for positive bias (shift only stratum W)
# check against Ding Appendix, pg 18 (Prop A.11)
targetB = Eadd.est.mono$candidates$biasFactor[ Eadd.est.mono$candidates$biasDir == "positive" ]

lambda = pw_0 * (1 - fw) - pw_1 * fw
term1 = 1 / ( 2 * (fw * pw_0) )
term2 = 4 * pw_1 * pw_0 * fw * (1 - fw)
true = RDm
term3 = true + lambda

# bias factor
( myB = term1 * ( sqrt( term3^2 + term2 ) - term3 ) )
expect_equal( targetB, myB, tol = 0.001 )

### check E-value for negative bias (shift only stratum M)
targetB = Eadd.est.mono$candidates$biasFactor[ Eadd.est.mono$candidates$biasDir == "negative" ]

lambda = pm_1 * (1 - fm) - pm_0 * fm 
term1 = 1 / ( 2 * (fm * pm_1) )
term2 = 4 * pm_1 * pm_0 * fm * (1 - fm)
true = -RDw
term3 = true + lambda
  
# bias factor
( myB = term1 * ( sqrt( term3^2 + term2 ) - term3 ) )
expect_equal( targetB, myB, tol = 0.001 )












