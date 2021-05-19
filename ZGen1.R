#==================================================================#
# Copyright (c) 2015. Yu-Cheng Ku. All Rights Reserved.            #
# Bayesian Multinomial Probit Estimation for Choice-Base Conjoint  #
# Written by Yu-Cheng Ku, September 18, 2015                       #
#                                                                  #
# This code is for the "Gibbs Thru" of the utility sampling        #
#==================================================================#


ZGen = function(Y.h, X.h, beta.h, Zh)
{
	# X.h[j,n,m]
	# x_ij is when in task m, for each att level, i is picked
	# Zh is N.m by M

	V.choi = rep(0, M)		 # representative util of choice
	V.base = rep(0, M)		 # representative util of second best
	Z.choi = rep(0, M)		 # Sampled utility of choice
	Z.base = rep(0, M)		 # Sampled value of second best
	eps    = rep(0, M)

	for(m in 1:M)
	{
		ym   = Y.h[,m]
		util = t(X.h[,,m]) %*% beta.h 
		zhm  = Zh[,m]

		ym.choi   = which(ym == 1)
		util.temp = util
		util.temp[ym.choi] = min(util.temp[-ym.choi])-1

		V.base[m] = max(util.temp)
		Zbound    = zhm[ym.choi]

		for(n in 1:N.m)
		{
			muu     = util[n]
			l.bound = max(zhm[-n])

			if(ym[n]==1)
			{
				Zh[n,m]   = rtrun(muu, 1, l.bound, 100)
				Z.choi[m] = Zh[n,m]
				V.choi[m] = muu
				eps[m]    = Zh[n,m] - muu

			} else
				{
					Zh[n,m] = rtrun(muu, 1, -100, min(l.bound, Zbound) )
				}			
		}

		Zh.temp = Zh[,m]
		Zh.temp[ym.choi] = min(Zh.temp[-ym.choi])-1

		Z.base[m] = max(Zh.temp)
	}
           	
	z.result = rbind(Zh, V.choi, V.base, Z.choi, Z.base, eps)

	return(z.result)
}


