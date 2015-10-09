
manifoldLearning = function (method = "tsne", distances = D, 
		k = 2, 
		perplexity = neighborhoodSize, 
		outputDimension = 2,
		epoch = 25) 
{
	# sanity checks
	
	method = tolower (method)
	
	# do we have an x-construction?
	if (substr(method, 1, 2) == "x-") {
		# find neighbors first
		neighborList = knnx.index(data$X, data$X, k=neighborhoodSize, algorithm="cover_tree")
		dXpt = estimateNbdTopology (data, neighborList)
		D = computeXDistance (data)
		
		# remove x- from method
		method = substr (method, 3, nchar(method))
	} else {
		D = dist(data$X, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
	}

	# now apply the method
	if (method == "tsne") {
		Y = tsne (D, k = 2, perplexity = neighborhoodSize, epoch = 25)
	} else if (method == "lle") {
		Y = lleWrapper(data, D, outputDimension = 2, neighborhoodSize = neighborhoodSize)
	} else if (method == "isomap") {
		Y = isomap(D, neighborhoodSize, ndim = 2)
	} else if (method == "lem") {
		Y = spec.emb(as.matrix(D), outputDimension, norm = TRUE)
	} else {
		stop ("Unkown method! Please refer to the documentation.")
	}
	
	# finally...
	retObj = list(Y = Y)
	return (retObj)
}
	