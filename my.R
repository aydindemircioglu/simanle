#
# Symmetric SNE example
#   (c) 2014, Aydin Demircioglu
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


library (methods)
library (SparseM)
library(BBmisc)
library(e1071)
library(tsne)
library(FNN)



source ("./divergence.R")
source ("./helper.R")
source ("./plotHelper.R")

source ("./estimateNbdTopology.R")



numberToColor <- function(p)
{
    cols <-  colorRamp(c("#000099", "#00FEFF", "#009E00", "#FCFF00", "#FF9400", "#FF3100"))
    np = (p-min(p))/(max(p)-min(p))
    return(cols(np)/256)
}

poseToColor <- function(p, minC, maxC)
{
    cols <-  colorRamp(c(minC, maxC))
    np = (p-min(p))/(max(p)-min(p))
    t = cols(np)
    t[,1] = (t[,1]-min(t[,1]))/(max(t[,1])-min(t[,1]))
    t[,2] = (t[,2]-min(t[,2]))/(max(t[,2])-min(t[,2]))
    t[,3] = (t[,3]-min(t[,3]))/(max(t[,3])-min(t[,3]))
    t[is.nan(t)] = 0
    return(t)
}


loadData <- function () {
  data1 = list()
	dataset = read.matrix.csr ("../data_3k_train.sparse")
	data1$X = as.matrix(dataset$x)
	data1$Y = as.matrix(as.numeric(as.character(dataset$y)))

	data2 = list()
	dataset = read.matrix.csr ("../data_3k_test.sparse")
	data2$X = as.matrix(dataset$x)
	data2$Y = as.matrix(as.numeric(as.character(dataset$y)))

	data = list()
	data$X = rbind(data1$X, data2$X)
	data$Y = rbind(data1$Y, data2$Y)
	data$col = numberToColor(data$Y)

  return (data)
}


cb <- function (ydata) {
    plotDR(ydata, data$X, col = data$col, name = model, text = heading)
}


# global variables
    model = "tsne"
    
# number of neighbors
    neighborhoodSize = 12
    outputDimension = 2

if (file.exists("XDist.Rdata") == FALSE) {
		
	# load data
		data = loadData ()
		nPoints = dim(data$X)[1]  

		
	# find neighbors first
		neighborList = knnx.index(data$X, data$X, k=neighborhoodSize, algorithm="cover_tree")

	# estimate topology
		messagef ("  Estimating topology.")
		dXpt = estimateNbdTopology (data, neighborList)


	# compute the new distance matrix    
		messagef ("  Computing distance matrix in X-space.")
		data$X = data$X[1:nPoints,]
		D = matrix(-1, nPoints, nPoints)
		for (i in 1:nPoints) {
			for (j in 1:i) {
				# get both covariance matrices and points
				S1 = dXpt$S[[i]]
				S2 = dXpt$S[[j]]
				mu1 = data$X[i, ]
				mu2 = data$X[j, ]
		
				# compute their distance
				D[i, j] = bhattacharyyaRiemannDivergence (mu1, S1, mu2, S2)
			}
		}

	# finish other half of matrix
		for (i in 1:nPoints) {
			for (j in i:nPoints) {
				D[i,j] = D[j,i]
			}
		}

	# convert it to a dist object
		XDist = as.dist(D, diag = FALSE, upper = FALSE)

	# cache results
		save(XDist, file = "XDist.Rdata")
} else {    
	# load from cache
		load ("XDist.Rdata")
    
	# apply tSNE to original data first
		messagef ("  Applying %s to data", model)
		XDist = dist(data$X, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
		Y = tsne (XDist, k = 2, perplexity = neighborhoodSize, epoch = 25)

		pdf("./notrandomized.pdf")
		plotDR(Y[1:1050,], data$X[1:1050,], col = data$col[1:1050,], name = model)
		plotDR(Y[1051:2100,], data$X[1051:2100,], col = data$col[1051:2100,], name = model)
		dev.off()

		pdf("./randomized.pdf")
		idx = sample(2100)
		plotDR(Y[idx[1:1050],], data$X[idx[1:1050],], col = data$col[idx[1:1050],], name = model)
		plotDR(Y[idx[1051:2100],], data$X[idx[1051:2100],], col = data$col[idx[1051:2100],], name = model)
		dev.off()
}






##########3

	# data is just data
    neighborhoodSize = 12
    outputDimension = 2
	
		
	Y = manifoldLearning (method = "tsne", distances = D, 
		k = 2, 
		perplexity = neighborhoodSize, 
		outputDimension = 2,
		epoch = 25)

	Y = manifoldLearning (method = "x-tsne", 
		divergence = "bhattacharyyaRiemann",
		k = 2, 
		perplexity = neighborhoodSize, 
		outputDimension = 2,
		epoch = 25)

	# plot things
	


computeXDistance = function (data = NULL, divergence = "bhattacharyyaRiemann", verbose = FALSE) {

	if (verbose == TRUE)
		messagef ("Computing distance matrix in X-space.")

	# sanity checks
	
	# compute the new distance matrix    
	nPoints = nrow(data)
	D = matrix(-1, nPoints, nPoints)
	for (i in 1:nPoints) {
		for (j in 1:i) {
			# get both covariance matrices and points
			S1 = dXpt$S[[i]]
			S2 = dXpt$S[[j]]
			mu1 = data$X[i, ]
			mu2 = data$X[j, ]
	
			# compute their distance
			D[i, j] = bhattacharyyaRiemannDivergence (mu1, S1, mu2, S2)
		}
	}

	# finish other half of matrix
	for (i in 1:nPoints) {
		for (j in i:nPoints) {
			D[i,j] = D[j,i]
		}
	}

	# convert it to a dist object
	XDist = as.dist(D, diag = FALSE, upper = FALSE)
	return (XDist)
}


