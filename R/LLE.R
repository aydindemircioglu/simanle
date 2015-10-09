## Computes the Locally Linear Embedding by Roweis, Saul and Lawrence (2000).
## It is based on the Matlab implementation by Roweis.
##
## input: 
##	data: 		NxD matrix (N samples, D features)
##	dim:		dimension of the target space
##	neighbours:	number of neighbours (default 5 or number of samples)
##
## output: 
##	Nxdim matrix (N samples, dim features) with the reduced input data

LLE = function(data,dim=2,k){

  
   	## catch missing or bad input
  	if(missing(data))
          stop("data argument missing")
        else
          if(!is.matrix(data))
            stop("invalid argument: data argument is required to be a N x D matrix (N samples, D features)")
        if(!all(is.numeric(dim)) | dim < 1)
          stop("invalid argument: target dimension is required to be a positive integer value")
        if(missing(k))
                k = min(nrow(data),5)
        else{
          if(k >= nrow(data))
            stop("invalid argument: more neighbours than samples")
          if(!is.numeric(k) |k <= 1)
            stop("invalid argument: neighbour parameter is required to be an integer value >= 2")
        }
        k = round(k)
        dim = round(dim)

	num_samples = nrow(data)
	num_features = ncol(data)

	## compute pairwise distances
        message("Computing distance matrix ... ", appendLF = FALSE)
	d = t(rowSums(data^2))
	d = d[rep(1, each = num_samples),]
	d = d + t(d) - 2*data%*%t(data)
	diag(d) = 0
        d = sqrt(d)

	## determine the k nearest neighbours
	sort_idx = apply(d,2,order)
	neighbours = sort_idx[2:(k+1),]
        message("done")

	## build the weights
        message("Computing low dimensional embedding (using ", k, " nearest neighbours)... ", appendLF = FALSE)
	W = matrix(0,k,num_samples)
	for(i in 1:num_samples){
	
		## compute covariance matrix of the neighbours of the ith sample (center and square samples)
		N = t(data[neighbours[,i],]) - data[i,]
		Cov = t(N) %*% N

		## weights W are solution of CW=1
		W[,i] = solve(Cov,rep(1,k))
		## rescale weights that rows sum to one (-> invariance to translations)
		W[,i] = W[,i] / sum(W[,i])
	}

	## build the cost matrix M = (I-W)'(I-W)
	M = diag(1,num_samples)
	for(i in 1:num_samples){
		w = W[,i]
		n = neighbours[,i]
		M[i,n] = M[i,n] - t(w)
		M[n,i] = M[n,i] - w
		M[n,n] = M[n,n] + w %*% t(w)
	}
	
	## low dimensional embedding Y given by the eigenvectors belonging to the (dim+1) smallest eigenvalues (first eigenvalue is trivial)
	eig_M = eigen(M)
	sweep(eig_M$vectors, 2, sqrt(colSums(eig_M$vectors^2)), "/")    ## normalize eingenvectors
	Y = eig_M$vectors[,(num_samples-1):(num_samples-dim)] * sqrt(num_samples)	

        message("done")

	return(Y)
}
