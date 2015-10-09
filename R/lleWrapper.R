#
# lle wrapper.
#   (c) 2014, Aydin Demircioglu
#
# stolen from the RDRToolbox.
#
#
#
# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This software is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this software. If not, see <http://www.gnu.org/licenses/>.
#


lleWrapper <- function (data, dist, outputDimension = 2, neighborhoodSize = 11, ...)
{
    # assume k > 2!
    k = neighborhoodSize
    dim = outputDimension

    # get a full symmetric matrix out of dist
    d = as.matrix(dist)

    num_samples = nrow(data)
    num_features = ncol(data)
    
    ## determine the k nearest neighbours
    sort_idx = apply(d,2,order)
    if (k > dim(sort_idx)[1] - 1)
    {
        k = dim(sort_idx)[1]-1
        print (k)
        }
    neighbours = sort_idx[2:(k+1),]
    message("done")
    
    ## build the weights
    message("  LLE: Computing low dimensional embedding (using ", k, " nearest neighbours)... ", appendLF = FALSE)
    W = matrix(0,k,num_samples)
    for(i in 1:num_samples){
        
        ## compute covariance matrix of the neighbours of the ith sample (center and square samples)
        N = t(data[neighbours[,i],]) - data[i,]
        G = t(N) %*% N

        #-- from lle package
        #regularisation
        delta <- 0.1
        #calculate eigenvalues of G
        e <- eigen(G, symmetric=TRUE, only.values=TRUE)$values
        #skip if all EV are null
        if( all( e==0 ) ) next
        
        #choose regularisation method  by tr(gram)
        r <- delta^2/k*sum(diag(G))
        
        #use regularisation if more neighbourse than dimensions!
        if( k>num_features ) alpha <- r else alpha <- 0
        
        #regularisation
        G <- G + alpha*diag(1,k)
        
        #calculate weights
        #using pseudoinverse ginv(A): works better for bad conditioned systems
        W[, i] = t(ginv(G)%*%rep(1,k)) 
        
        #--
        
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

