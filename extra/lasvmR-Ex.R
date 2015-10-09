pkgname <- "simanle"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "simanle-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('simanle')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("lasvmPredict")
### * lasvmPredict

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: lasvmPredict
### Title: lasvmPredict
### Aliases: lasvmPredict

### ** Examples

model = simanle::lasvmTrain (x = as.matrix(iris[seq(1,150,2),1:4]),
	y = (as.numeric(iris[seq(1,150,2),5]) %% 2)*2-1,
	gamma = 1, 
	cost = 1, 
	kernel = 2)
ytrue = (as.numeric(iris[seq(2,150,2),5]) %% 2)*2-1
result = lasvmPredict (x = as.matrix(iris[seq(2,150,2),1:4]), model)
ypred = result$predictions
error = sum(abs(ypred - ytrue))/length(ytrue)
cat ("Error rate =", error*100)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("lasvmPredict", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("lasvmTrain")
### * lasvmTrain

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: lasvmTrain
### Title: lasvmTrain
### Aliases: lasvmTrain

### ** Examples

model = simanle::lasvmTrain (x = as.matrix(iris[seq(1,150,2),1:4]),
			y = as.matrix(iris[seq(1,150,2),5]),
				gamma = 0.1, cost = 1.0, kernel = 2)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("lasvmTrain", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
