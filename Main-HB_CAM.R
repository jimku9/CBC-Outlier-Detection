rm(list=ls())

#==================================================================#
# Copyright (c) 2015. Yu-Cheng Ku. All Rights Reserved.            #
# Bayesian Multinomial Probit Estimation for Choice-Base Conjoint  #
# Written by Yu-Cheng Ku, September 18, 2015                       #
#                                                                  #
# This code is the main function to handle input, initial settings #
# and output results                                               #
#==================================================================#



#===========#
#  Headers  # 
#===========#

require(MCMCpack)
require(bayesm)
require(mvtnorm)
require(truncnorm)


data.dir     = "C:/Users/Yu-Cheng/CBC/Data/"
all.dir      = "C:/Users/Yu-Cheng/CBC/Data/"
function.dir = "C:/Users/Yu-Cheng/CBC/Data/"
HB           = "HB_ALL.R"
IDConvert    = "ID Convert.R"
DummyKminus1 = "Dummy k-1_none.R"
ZGen         = "ZGen1.R"
Data         = paste("cam_input.csv")

source(paste(all.dir,      HB,           sep=""))
source(paste(function.dir, IDConvert,    sep=""))
source(paste(function.dir, DummyKminus1, sep=""))
source(paste(function.dir, ZGen,         sep=""))


#===================#
#  Data processing  #
#===================#

raw = read.csv(paste(data.dir, Data, sep=""), header=T)
id  = raw$Respondent
PIN = IDConvert(id)

N   = length(PIN)
H   = length(unique(PIN))
M   = length(unique(raw[,2]))
N.m = length(unique(raw[,3]))


raw[,1] = PIN
Brand   = raw[,4]
none    = rep(0, N)
none[which(Brand==0)] = 1

X.LEV = raw[,4:(ncol(raw)-1)]
X.DUM = DummyKminus1(as.matrix(X.LEV))

X.beta = cbind(X.DUM, none)
y      = raw[,ncol(raw)]

n.beta = ncol(X.beta)


#=================#
#  MCMC Settings  #
#=================#

n.iter = 10000
burn   = 5000
mcmc   = list(n.iter = n.iter, burn = burn)



#==================#
#  Prior Settings  #
#==================#

# The setting is for V^-1 ~ W(nu.0, V.0)
nu.0  = n.beta + 5
V.0   = (nu.0) * diag(n.beta)

# The seeting is for Sigma_mu_beta
Sig_mu.0 = 100 * diag(nrow(V.0))

hyper = list(nu = nu.0, V = V.0, Sig_mu = Sig_mu.0)




#==================#
#  Initial Values  #
#==================#

beta   = apply(X.beta, 2, mean)
V      = cov(X.beta)
beta.h = t(matrix(rep(beta, H), nrow = H, byrow = T)) # beta.h is n.beta by H
start  = list(beta = beta, V = V, beta.h = beta.h)


X = array( t(X.beta), dim = c(n.beta, N.m, M, H))	# X[,,i,j] is n.beta by N.m
Y = array( t(y), dim = c(N.m, M, H))


Xh = array(0, dim = c(n.beta, N.m*M, H))
for(h in 1:H)
{
	Xh[,,h] = matrix(as.vector(X[,,,h]), n.beta, N.m*M)
}


time.start = proc.time()
output = MCMC(Y, X, Xh, mcmc, hyper, start)
cat(proc.time() - time.start, "\n")


save.image(paste("C:/Users/Yu-Cheng/Dropbox/CBC/Submit to Marketing Letters/Revise/empirical/probit_result_everything.RData"))

beta = t(output$beta)

None = beta[,ncol(beta)]

estbetas = cbind(
	beta[,1:2],  0-apply(beta[,1:2],1,sum),
	beta[,3:4],  0-apply(beta[,3:4],1,sum),
	beta[,5:9],  0-apply(beta[,5:9],1,sum),
	beta[,10],   0-beta[,10],
	beta[,11:12],0-apply(beta[,11:12],1,sum),
	beta[,13],   0-beta[,13],
	beta[,14:16],0-apply(beta[,14:16],1,sum),
	beta[,17:19],0-apply(beta[,17:19],1,sum),
	None
)

colnames(estbetas) = c("Body_L","Body_M","Body_H","Chng_N","Chng_M","Chng_A",
				"Anno_0","Anno_1","Anno_2","Anno_3","Anno_4","Anno_5","Optn_N","Optn_Y","Zoom_0","Zoom_2","Zoom_4",
				"View_R","View_L","Set_N","Set_L","Set_V","Set_B","Prce_1","Prce_2","Prce_3","Prce_4","None")

round(as.matrix(apply(estbetas, 2, mean)), 3)



par(mar = c(4.5,4.5,2,2))
plot.ts(t(output$beta), plot.type = "single", col = c(1:19), ylab = "Part-worth Utilities", xlab = "Iteration after Burn-in", main = "")





