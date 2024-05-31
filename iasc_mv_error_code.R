library(Matrix)
library(matrixcalc)
library(truncnorm)

iasc_mv.error = function(matrix, matrixSE, rownum, colnum){
  
  svl.male = rtruncnorm(1, a = (matrix[1,1] - (2*matrixSE[1,1])), b = Inf, mean = matrix[1,1], sd = matrixSE[1,1])
  area.male = rtruncnorm(1, a = (matrix[2,2] - (1.5*matrixSE[2,2])), b = Inf, mean = matrix[2,2], sd = matrixSE[2,2])
  hue.male = rtruncnorm(1, a = (matrix[3,3] - (3.5*matrixSE[3,3])), b = Inf, mean = matrix[3,3], sd = matrixSE[3,3])
  bri.male = rtruncnorm(1, a = (matrix[4,4] - (2*matrixSE[4,4])), b = Inf, mean = matrix[4,4], sd = matrixSE[4,4])
  sat.male = rtruncnorm(1, a = (matrix[5,5] - (2.5*matrixSE[5,5])), b = Inf, mean = matrix[5,5], sd = matrixSE[5,5])
  svl.female = rtruncnorm(1, a = (matrix[6,6] - (2*matrixSE[6,6])), b = Inf, mean = matrix[6,7], sd = matrixSE[6,6])
  area.female = rtruncnorm(1, a = (matrix[7,7] - (2*matrixSE[7,7])), b = Inf, mean = matrix[7,7], sd = matrixSE[7,7])
  hue.female = rtruncnorm(1, a = (matrix[8,8] - (3.5*matrixSE[8,8])), b = Inf, mean = matrix[8,8], sd = matrixSE[8,8])
  bri.female = rtruncnorm(1, a = (matrix[9,9] - (1.5*matrixSE[9,9])), b = Inf, mean = matrix[9,9], sd = matrixSE[9,9])
  sat.female = rtruncnorm(1, a = (matrix[10,10] - (2.5*matrixSE[10,10])), b = Inf, mean = matrix[10,10], sd = matrixSE[10,10])
  
  rangeB = matrix[rownum,colnum]
  negrangeB = rangeB * -1
  
  if(rangeB > 0){
    blist = seq(negrangeB, rangeB, abs((2*rangeB)/101))
  } else {
    blist = seq(rangeB, negrangeB, abs((2*rangeB)/101))
  }
  
  
  for(i in 0:100){
    B = blist[i + 1]
    
    GBmat = matrix
    GBmat[1,1] = svl.male
    GBmat[2,2] = area.male
    GBmat[3,3] = hue.male
    GBmat[4,4] = bri.male
    GBmat[5,5] = sat.male
    GBmat[6,6] = svl.female
    GBmat[7,7] = area.female
    GBmat[8,8] = hue.female
    GBmat[9,9] = bri.female
    GBmat[10,10] = sat.female
    
    GBmat[rownum,colnum] = B
    GBmat[colnum,rownum] = B
    
    for(j in 1:length(s.area.m)){
      if(is.positive.definite(GBmat) == T){
        res = 0.5 * GBmat %*% dew.sel.t[,j]
        rm[j] = res[colnum,1]
        rf[j] = res[rownum,1]
      } else {
        GBmat = nearPD(GBmat)
        GBmat = matrix(GBmat$mat,nrow = 10)
        res = 0.5 * GBmat %*% dew.sel.t[,j]
        rm[j] = res[colnum,1]
        rf[j] = res[rownum,1]
      }
    }
    res.cor = cor.test(rm,rf)
    Br[i+1,1] = blist[i+1]
    Br[i+1,2] = res.cor$estimate
    Br[i+1,3] = area.male   
    Br[i+1,4] = area.female
    iasc.res <<- Br
  }
}

wittmat = matrix(NA,10,10)
wittmat[lower.tri(wittmat,diag = T)] = c(0.00048, 
                                         0.00150, -0.01410,
                                         -0.02956, 0.01389, 0.00004,
                                         0.00171, -0.00001, -0.01530, 0.03682,
                                         0.01308, -0.19268, -0.01104, 0.26991, 0.00063,
                                         0.00400, -0.13012, -0.26488, 0.15978, 7.4638, 2.3422,
                                         -11.2354, 0.00370, -0.08715, 7.8539, 8.5939, -3.8454, 16.1486,
                                         16.2746, -0.01022, -0.13771, -0.43307, 2.9513, -0.77655, 47.9309, 0.02370,
                                         0.18751, -13.3216, -12.1292, 15.9731, 0.00041, 0.00112, 0.02845, 00.00456, 0.07169,
                                         0.01806, -0.05958, -0.10137, 0.50522, 11.2918, 10.3346, -3.7493, 13.5126, -11.8242, 43.4240)
wittmat = Matrix::forceSymmetric(wittmat,uplo = "L")
wittmat = matrix(wittmat, nrow = 10)

witterrormat = matrix(NA,10,10)
witterrormat[lower.tri(witterrormat,diag = T)] = c(0.00022,
                                                   0.00108, 0.01505,
                                                   0.02812, 0.04285, 0.00014,
                                                   0.00091, 0.01821, 0.02410, 0.03877,
                                                   0.00808, 0.09450, 0.16092, 0.25569, 0.00079,
                                                   0.00504, 0.11271, 0.16366, 0.22873, 1.9603, 2.6368,
                                                   4.4424, 0.01213, 0.08358, 2.0522, 2.5040, 3.7739, 6.5648,
                                                   8.3152, 0.02458, 0.15757, 2.94525, 3.7458, 7.03649, 17.0944, 0.03972,
                                                   0.24898, 5.13778, 6.02776, 11.2576, 0.00017, 0.00082, 0.01575, 0.01843, 0.03623,
                                                   0.00738, 0.10551, 0.12630, 0.24476, 3.2128, 3.4264, 4.6739, 5.0548, 6.3608, 14.5926)
witterrormat = Matrix::forceSymmetric(witterrormat,uplo = "L")
witterrormat = matrix(witterrormat, nrow = 10)

s.svl.m = runif(10000, -1, 1)
s.area.m = runif(10000, -1, 1)
s.hue.m = runif(10000, -1, 1)
s.bri.m = runif(10000, -1, 1)
s.sat.m = runif(10000, -1, 1)
s.svl.f = runif(10000, -1, 1)
s.area.f = runif(10000, -1, 1) 
s.hue.f = runif(10000, -1, 1)
s.bri.f = runif(10000, -1, 1)
s.sat.f = runif(10000, -1, 1)

for(i in 1:length(s.sat.m)){ #change this to be the trait you want for antagonism
  if (s.sat.m[i] > 0) {
    s.sat.f[i] = runif(1, -1, 1)
    if(s.sat.f[i] > 0) {
      s.sat.f[i] = s.sat.f[i] * -1
    }
  } else {
    s.sat.f[i] = runif(1, -1, 1)
    if(s.sat.f[i] < 0) {
      s.sat.f[i] = s.sat.f[i] * -1
    }
  }
}
dew.sel = cbind(s.svl.m,s.area.m,s.hue.m,s.bri.m,s.sat.m,s.svl.f,s.area.f,s.hue.f,s.bri.f,s.sat.f)
dew.sel.t = t(dew.sel)


Br = matrix(ncol = 4, nrow = 101)
rm = vector()
rf = vector()

iasc_mv.error(wittmat, witterrormat, 10, 5)
plot(iasc.res[,1],iasc.res[,2])

fit = nls(iasc.res[,2]+1 ~ SSlogis(iasc.res[,1],Asym,xmid,scal))
fit2 = lm(iasc.res[,2] ~ poly(iasc.res[,1],2, raw = T))
fit3 = lm(iasc.res[,2] ~ iasc.res[,1])
aic.sig = AIC(fit)
aic.poly2 = AIC(fit2)
aic.lin = AIC(fit3)


if(aic.sig < aic.poly2 & aic.sig < aic.lin) {
  mean = coef(fit)[2]
} else {
  if(aic.poly2 < aic.sig & aic.poly2 < aic.lin){
    a = fit2$coefficients[3]
    b = fit2$coefficients[2]
    c = fit2$coefficients[1]
    
    mean = (-b - sqrt(b^2 - 4*a*c))/(2*a)
    print((-b + sqrt(b^2 - 4*a*c))/(2*a))
  } else {
    if(aic.lin < aic.sig & aic.lin < aic.poly2){
      a = fit3$coefficients[1]
      b = fit3$coefficients[2]
      
      mean =(-a)/b
    }
  }
}
out = cbind(mean, iasc.res[1,3], iasc.res[1,4])
write.csv(out, paste0("Simulation_", "replace_this", "_result.csv"))