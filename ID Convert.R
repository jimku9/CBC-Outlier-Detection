
####  This is to convert the respondent numbers to ordered numbers  ####

IDConvert = function(id)
{
	q = 1

	ID = NULL
	for(i in 1:(length(id)-1))
	{
		ID[i] = q

		if(id[i+1] != id[i])
		{
			q = q + 1
		}
	}

	ID = c(ID, ID[length(ID)])
	return(ID)
}