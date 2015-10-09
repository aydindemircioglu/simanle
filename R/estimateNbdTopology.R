#
# Neighborhood Topology Estimation
#   (c) 2014, Aydin Demircioglu
#
#
#
# This is an implemention of 
#   Learning a Metric Space for Neighbourhood Topology Estimation: 
#   Application to Manifold Learning
#   by Karim T. Abou-Moustafa, Dale Schuurmans and Frank Ferrie
#   JMLR 29 (2013), pp 341-356
#
# The implementation is a bit rough at the edges, as it does not 
# follow the paper for now. Will hopefully change in future.
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

library(robust)
library(covRobust)
library(corpcor)


estimateNbdTopology <- function (data, neighborList)
{
    # estimate gaussian in neighbors of every point
    Xpt <- data.frame(i=integer(), S=I(list()))
    for (currentIndex in 1:dim(data$X)[1])
    {
        if (currentIndex %% 100 == 1)
            messagef("Estimation at point %d", currentIndex - 1)
        
        point = data$X[currentIndex,]
        neighborIndices = neighborList[currentIndex,]
        neighbors = data$X[neighborIndices,]
        centeredNeighbors = sweep(neighbors, 2, t(point), "-")
        
        # finding covariance matrix
        #C = 1/length(neighborhoodIndices) * neighborhood' * neighborhood
        lambda = 10^{-6}
        covariance = cov.shrink(centeredNeighbors, lambda, verbose = FALSE)
        covariance = matrix(as.double(covariance), dim(covariance))
#        covariance = cov(centeredNeighbors)
#        covariance = covRob(centeredNeighbors, estim = "mcd")$cov
#        covariance = cov.nnve(as.matrix(centeredNeighbors), k = 5, extension = FALSE)$cov
        #print(str(covariance))
        
        # now add the variables to the frame
        Xpt[[currentIndex, 'S']] = covariance
        Xpt[[currentIndex, 'i']] = currentIndex
    }
    
    return (Xpt)
}
