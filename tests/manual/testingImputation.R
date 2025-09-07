load("/Users/alessandrogiordano/Desktop/AnaConesa/NPLSDAPackage_Support/prepareOmics/Obj/TEDDY/arrayGCTOFX.RData")
library(tictoc)
dim(arrayGCTOFX)
tic("Exec Time BestFit GCTOFX ( 136 364   5")
bestGCTOFX = bestfittedmodel(X = arrayGCTOFX, centering = 0, top_fit_n = TRUE, verbose = FALSE)
bestGCTOFX$top_by_fit
toc()
tic("Exec Time Imputation GCTOFX ( 136 364   5")
fullarrayGCTOFX = Imputemethod(arrayGCTOFX,fac = c(4,2,4),max.iter = 1000,verbose = FALSE,seed =  123)

158.121/60
496.166/60
summary(fullarrayGCTOFX)
dim(fullarrayGCTOFX)
anyNA(fullarrayGCTOFX)
anyNA(arrayGCTOFX)
b

toc()
arrStategra = get(load("/Users/alessandrogiordano/Desktop/AnaConesa/NPLSDAPackage_Support/prepareOmics/Obj/Stategra/arrayStategraGE.RData"))
dim(arrStategra)
bestArrStategra = bestfittedmodel(X = arrStategra, centering = 0,min_comp = 2 ,max_comp =  4,top_fit_n = 5)
bestArrStategra$summary
bestArrStategracenter2 = bestfittedmodel(X = arrStategra, centering = 2,min_comp = 2 ,max_comp =  4,top_fit_n = 5)

bestArrStategracenter3 = bestfittedmodel(X = arrStategra, centering = 1,min_comp = 2 ,max_comp =  4,top_fit_n = 5)
