
source("./tsne-internal.R")

computeQDistribution <- function (ydata, model)
{
    sum_ydata = apply(ydata^2, 1, sum)
    if (model == "ssne") {
        # compute exp( - || y_i - y_j ||^2) = nominator of q
        num = sum_ydata + sweep(-2 * ydata %*% t(ydata), 2, -t(sum_ydata))
        
        # numerical stability!!
#        num = num/sum(num)
        
        num = exp(-num)
    }
    
    if (model == "tsne") {
        # compute 1/(1+|| y_i - y_j ||^2) = nominator of q
        num =  1/(1 + sum_ydata + sweep(-2 * ydata %*% t(ydata), 2, -t(sum_ydata))) 
    }

    # make sure q_ii = 0
        diag(num)=0
        
    # simple sanity check
    if (any(is.nan(num))) 
    {
        message ('NaN in grad. descent')
        print (num)
    }
    # normalize the distribution
        Q = num / sum(num) 
    
    if (model == "ssne") {
        num = matrix(1, dim(num)[1], dim(num)[2])
    }
    
    # make sure no value is below eps (machine precision)
        eps = 2^(-52) 
        Q[Q < eps] = eps

    return (list("Q" = Q, "num" = num))
}    


# symmetric sne cost
CKL <- function (P, Q)
{
    # never divide by zero, please.
    eps = 2^(-52)
    
    if ( (abs(sum(P) - 1) > 10^-15 ) || (abs(sum(Q) - 1) > 10^-15 )) {
        P = P/sum(P)
        Q = Q/sum(Q)
    }

    # make sure no $Q$ or $P$ is == 0.
    cost = sum(apply(P * log((P+eps)/(Q+eps)),1,sum))
    return (cost)
}


sne <- function(X, 
                            k=2, 
                            initial_dims=30, 
                            perplexity=30, 
                            max_iter = 1000, 
                            min_cost=0*10^{-18}, 
                            epoch_callback=NULL,
                            whiten=TRUE,
                            epsilon = 100,
                            epoch = 25,
                 model = "tsne")
{
    # number of data points 
    n = attr(X,'Size')
 
    # parameters
	momentum = .5
	final_momentum = .8
	mom_switch_iter = 250

	min_gain = .01
	initial_P_gain = 4

	# typical machine precision
	eps = 2^(-52) 

	# initialize the solution randomly or the given config
    ydata = matrix(rnorm(k * n),n)
	
	
    # determine the probablity by finding the variance
    # that corresponds to the perplexity we specified
	r = .x2p (X, perplexity, 1e-5)
    P = r$P

    # p is not symmetric, it is the 'true' P
    P = (P + t(P))/(2)
    
    # make sure no P is smaller than eps
    # should not happen anyway.
	P[P < eps]<-eps

    # P will now sum to 1, it is now a probability distribution on $I \times I$
    # but we stil make sure it does what it should
	P = P/sum(P)
    

    # prepare stochastic gradient descent
	P = P * initial_P_gain
	grads =  matrix(0,nrow(ydata),ncol(ydata))
	incs =  matrix(0,nrow(ydata),ncol(ydata))
	gains = matrix(1,nrow(ydata),ncol(ydata))
    
	# big gd loop
	for (iter in 1:max_iter)
    {
        # report every epoch
        if ((iter %% epoch == 0) && (iter > 1)) 
        { 
            # compute current cost, using KL-divergence
            cost = CKL (P,Q)
            message("Epoch: Iteration #",iter," error is: ", cost)
		
            # check if we reached the specified min cost, in that case we end 
            if (cost < min_cost) 
                break
            
            # epoch callback
			if (!is.null(epoch_callback)) 
                epoch_callback(ydata)
		}

        # compute Q
		r = computeQDistribution(ydata, model)
        Q = r$Q
        num = r$num
        
        # compute gradient
        # first stiffnessmatrix, i.e. the 'fixed' part
        stiffnesses = 4 * (P-Q) 
        
		for (i in 1:n){
		    # compute y_i - y_. = vector
		    ydiff = sweep(-ydata, 2, -ydata[i,])

            # muliply both
			grads[i,] = apply(ydiff * stiffnesses[,i],2,sum)
		}

        # if the update goes into the same direction as the previous,
        # we add more speed to it (factor 1.2). if not, we probably crossed
        # the zero-point, and should decrease the speed (factor 0.8)
		    gains = (gains + .2) * abs(sign(grads) != sign(incs)) 
			    	+ gains * .8 * abs(sign(grads) == sign(incs))

        # make sure we have everywhere some gain
	        gains[gains < min_gain] = min_gain

        # do gradient step and add momentum  
    #messagef("%f g -- %f gr", sum(abs(gains)), sum(abs(grads)))
	        incs = momentum * incs - epsilon * (gains * grads)
    	    ydata = ydata + incs
		
        # substract mean of each point from data, to center it
            ymean = apply(ydata, 2, mean)
            ydata = sweep(ydata, 2, ymean)
		
        # after specified iteration, we only use the given momentum.
        if (iter == mom_switch_iter) 
            momentum = final_momentum
		
        # decrease the speed after 100 iterations, i.e. undo 'early exaggeration'
		if (iter == 100) 
            P = P/initial_P_gain
	}
	ydata
}




.checkgrad = function() {
	    # -- check gradient.
        # finite differences will give us the components of the gradient.
        approxgrads = grads
		costG = CKL (P, Q)
		for (i in 1:n)
        {
            eps = 1/100000000000000
            for (t in 1:2)
            {
                gdata = ydata
                gdata[i,t] = gdata[i,t] + eps
                Qe = computeQDistribution(gdata)
                costFD = CKL (P, Qe)
                approxgrads[i,t] = (costFD-costG)/eps
#            messagef( "%f vs %f", costG, costFD)
            }
        }
		# --- 
}        
