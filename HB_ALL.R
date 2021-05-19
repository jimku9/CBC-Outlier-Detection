#==================================================================#
# Copyright (c) 2015. Yu-Cheng Ku. All Rights Reserved.            #
# Bayesian Multinomial Probit Estimation for Choice-Base Conjoint  #
# Written by Yu-Cheng Ku, September 18, 2015                       #
#                                                                  #
# This code is where the MCMC is performed and residuals are saved #
#==================================================================#



MCMC = function(y, X, Xh, mcmc, hyper, start)
{
	set.seed(201203)

	div = 5

	## MCMC Settings  
	n.iter = mcmc$n.iter
	burn   = mcmc$burn

	## Dimension information
	H   = dim(X)[4]
	M   = dim(X)[3]
	N.m = dim(X)[2]

	## Prior Settings
	nu.0     = hyper$nu
	V.0      = hyper$V
	Sig_mu.0 = hyper$Sig_mu

	## Initial Values
	beta   = start$beta
	V.inv  = solve(start$V)
	beta.h = start$beta.h
	n.beta = length(beta)
	z.h    = array(0.2, dim = c(N.m, M, H))	# initialize Z


	## Storage 
	L = n.iter - burn
	keep.beta   = matrix(0, n.beta, L/div)
	keep.betah  = array(0, dim = c(n.beta, H, L/div))
	keep.V_choi = array(0, dim = c(M, H, L/div))
	keep.V_base = array(0, dim = c(M, H, L/div))
	keep.Z_choi = array(0, dim = c(M, H, L/div))
	keep.Z_base = array(0, dim = c(M, H, L/div))
	keep.zh     = array(0, dim = c(M*N.m, H, L/div))
	keep.Sig_Bh = array(0, dim = c(n.beta*n.beta, H, L/div))
	eps2        = matrix(0, M, H)

	for(l in 1:n.iter)
	{
		## Generate h specific parameters
		for(h in 1:H)
		{
			X.h = X[,,,h]	# X.h is an array
			Y.h = Y[,,h]

			# Step 1a: Update z for h
			z.out    = ZGen(Y.h, X.h, beta.h[,h], z.h[,,h])

			# Step 1b: Update the entire z_h
			z.h[,,h] = z.out[1:N.m, ]
			zh       = as.vector(z.h[,,h])

			## Step 2: Generate beta_h
			#X.h = Xh[,,h]	# Xh is n.beta by N.m*M
			var.b      = solve(Xh[,,h] %*% t(Xh[,,h]) + V.inv)
			mean.b     = var.b %*% (Xh[,,h] %*% zh + V.inv %*% beta)
			beta.h[,h] = as.vector(rmvnorm(1, mean.b, var.b))


			if(l > burn & l %% div == 0)
			{
				eps2[,h]                       = eps2[,h] + z.out[(N.m+5),]^2	
				keep.V_choi[, h, (l-burn)/div] = z.out[(N.m+1),]
				keep.V_base[, h, (l-burn)/div] = z.out[(N.m+2),]
				keep.Z_choi[, h, (l-burn)/div] = z.out[(N.m+3),]
				keep.Z_base[, h, (l-burn)/div] = z.out[(N.m+4),]
				keep.beta[ , (l-burn)/div]     = beta
				keep.betah[, h, (l-burn)/div]  = beta.h[,h]
				keep.zh[ ,h, (l-burn)/div]     = zh
				keep.Sig_Bh[ , h, (l-burn)/div] = as.vector(var.b)
			} 			
		}

		## Step 3: Update beta_bar
		post.var.beta  = solve(H * V.inv + solve(Sig_mu.0))
		post.mean.beta = post.var.beta %*% V.inv %*% apply(beta.h, 1, sum)
		beta = as.vector(rmvnorm(1, post.mean.beta, post.var.beta))

		## Step 4: Update Sigma_beta		
		# problem could be the scale matrix
		V.inv = rwish( nu.0 + H, solve(V.0 + (beta.h - beta) %*% t(beta.h - beta)) )

		if(l %% 500 == 0)
		{
				plot.ts(t(keep.beta), plot.type = "single", col = c(1:9))
		}
	}

	eps2  = eps2/(L/div)	

	out = list(beta = keep.beta, betah =keep.betah, keep.V_choi = keep.V_choi, keep.V_base = keep.V_base, 
		     keep.Z_choi = keep.Z_choi, keep.Z_base = keep.Z_base, keep.zh = keep.zh, eps2 = eps2, Sig_Bh.mat = keep.Sig_Bh)

	return(out)
}




