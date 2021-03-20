

Eadd.est$evalue

gamma = pw_0 * (1 - fw) - pw_1 * fw + pm_1 * (1 - fm) - pm_0 * fm
term1 = 1 / ( 2 * (fw * pw_0 + fm * pm_1) )
term2 = 4 * (pw_0 * fw + pm_1 * fm) * ( pw_1 * (1 - fw) - pm_0 * (1-fm) )

# bias factor
( my.Badd.est = term1 * ( sqrt( gamma^2 + term2 ) - gamma ) )

Eadd.est$biasFactor

g(my.Badd.est)
