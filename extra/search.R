
library(simanle)


# take 0815 iris set
set.seed(32)
d = iris[sample(nrow(iris)),]
x = as.matrix(d[,1:4])
y = as.matrix(as.numeric(d[,5]))
y[y==3] = 1
y[y==2] = -1

xt = as.matrix(x[101:150,])
yt = as.matrix(as.numeric(y[101:150,]))

model = lasvmTrain (x[1:100,], y[1:100], gamma = 0.1, cost = 1, epochs = 5, optimizer = 0, kernel = 2, selection = 1, verbose = FALSE)
predictions = lasvmPredict (xt, model, verbose = FALSE)



# test the same for polynomial kernel
error = 1
	set.seed(21)
while (error != 0) {
	degree = runif(1)*12
	coef0 = runif(1)*300
	cost = runif(1)*1000
	if (runif(1) < 0.5)
		coef0 = 1/coef0
	if (runif(1) < 0.5)
		cost = 1/cost
	if (runif(1) < 0.5)
		degree = 1/degree
	
	print (degree)
	print (coef0)
	print (cost)
	
	model = lasvmTrain (x[1:100,], y[1:100,], degree = degree, coef0 = coef0, cost = cost, kernel = 1, selection = 2, verbose = FALSE)
	predictions = lasvmPredict (xt, model, verbose = FALSE)
	error = sum(abs(predictions$predictions - yt))
}
