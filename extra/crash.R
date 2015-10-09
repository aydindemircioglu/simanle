
library(simanle)

set.seed (5)


# take 0815 iris set, make it binary
d = iris[sample(nrow(iris)),]
x = as.matrix(d[,1:4])
y = as.matrix(as.numeric(d[,5]))
y[y==3] = 1
y[y==2] = -1

s = sample(nrow(iris))
p = round(runif(1)*(nrow(iris)-10))+5
trIdx = s[1:p]
testIdx = s[p+1:nrow(iris)]

degree = 60
coef0 = 1
cost = 30000
kernel = 1
selection = 1
	
model = lasvmTrain (x[trIdx,], y[trIdx,], degree = degree, coef0 = coef0, cost = cost, kernel = kernel , selection = selection, verbose = FALSE)
predictions = lasvmPredict (x[testIdx,], model, verbose = FALSE)



