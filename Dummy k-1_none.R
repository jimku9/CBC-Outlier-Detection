DummyKminus1 = function(Level)
{
	B     = Level
	N     = nrow(B)
	n_att = ncol(B)

	end    = 0
	Design = rep(0, N)
	for(i in 1:n_att)
	{
		att = B[,i]
		Design1 = matrix(0, N, length(unique(att))-1)
	
		end = end + length(unique(att))

		for(k in 1:max(unique(att)))
		{
			for(t in 1:N)
			{
				if(B[t,i] == k)
				{
					Design1[t,k] = 1
				}
			}
		}

		Design = cbind(Design, Design1[,-ncol(Design1)])
	}

	Design = Design[,-1]
	return(Design)
}



