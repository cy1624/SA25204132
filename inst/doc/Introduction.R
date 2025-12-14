## -----------------------------------------------------------------------------
require(microbenchmark)
require(SA25204132)

## -----------------------------------------------------------------------------
data("CSI300_2")
microbenchmark(FactorNum(CSI300_2), FactorNum_cpp(CSI300_2), times = 20L)

microbenchmark(Loadings_Factor(CSI300_2, r=c(1,13)), 
               Loadings_Factor_cpp(CSI300_2, r=c(1,13)), 
               times = 20L)

