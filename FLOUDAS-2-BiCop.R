
########
######## Import the libraries
########

library(PearsonDS)
library(e1071)
# Copula package
library(copula)
library(VineCopula)
# Fancy 3D plain scatterplots
library(scatterplot3d)
# ggplot2
library(ggplot2)
# Useful package to set ggplot plots one next to the other
library(grid)





########
######## Generate the distributions of the uncertain parameters using statistical moments
########

set.seed(2)

n <- 1000

moments7 <- c (mean=16, variance=2.5, skewness= -0.8, kurtosis=3.8)
moments8 <- c (mean=11, variance=2.5, skewness= 0.8, kurtosis=3.3)

v7 <- rpearson (n, moments=moments7)
v8 <- rpearson (n, moments=moments8)

#data.min <- min(v7)
#data.max <- max(v7)
#data.mean <- mean(v7)
#data.var <- var(v7)
#data.dskew <- skewness(v7)*data.var^1.5 # Denormalized skewness
#data.dkurt <- (kurtosis(v7) + 3)*data.var^2 # Denormalized kurtosis

data <- data.frame(v7,v8)

write.table(data,file='exportdata.csv',sep=',')





########
######## Manipulation of data and histograms
########

par(mfrow = c(1, 2))

# Estimate x  gamma distribution parameters and visually compare simulated vs observed data
x_mean <- mean(data$v7)
x_var <- var(data$v7)
x_rate <- x_mean / x_var
x_shape <- ( (x_mean)^2 ) / x_var

hist(data$v7, breaks = 20, col = "green", density = 20)
hist(rgamma( nrow(data), rate = x_rate, shape = x_shape), breaks = 20,col = "blue", add = T, density = 20, angle = -45)

# Estimate y gamma distribution parameters and visually compare simulated vs observed data
y_mean <- mean(data$v8)
y_var <- var(data$v8)
y_rate <- y_mean / y_var
y_shape <- ( (y_mean)^2 ) / y_var

hist(data$v8, breaks = 20, col = "green", density = 20)
hist(rgamma(nrow(data), rate = y_rate, shape = y_shape), breaks = 20, col = "blue", add = T, density = 20, angle = -45)





########
######## Measure association using Kendall's Tau
########
cor(data, method = "kendall")




########
######## Selecting the BIVARIATE copula using VineCopula package
########

####### function pobs compute the peudo-observations of a given data-matrix
var_a <- pobs(data)[,1]
var_b <- pobs(data)[,2]

####### BiCopSelect  accepts pseudo observations as parameters
selectedCopula <- BiCopSelect(var_a, var_b, familyset = NA)
selectedCopula
fam <- selectedCopula$family
par <- selectedCopula$par
par2 <- selectedCopula$par2

####### BiCop creates the object given the selected copula family and parameters
cop_model <- BiCop(fam, par, par2, tau = NULL, check.pars = TRUE)




########
######## Generation of data using the selected BIVARIATE copula
########
df <- BiCopSim(1000,cop_model)
df7 <- df[,1]
df8 <- df[,2]



########
######## Calculation of the joint CDF of the bivariate distribution for the generated data
########

copulacdf <- BiCopCDF(df7,df8,cop_model)



########
######## Linear interpolation to find the values of the generated data given their CDF values
########

f7 <- ecdf(v7)
y7 <- f7(df7)
inv_ecdf <- function(f7){
  x <- environment(f7)$x
  y <- environment(f7)$y
  approxfun(y, x)
}
g7 <- inv_ecdf(f7)

newv7 <- g7(df7)

newv7[is.na(newv7)] <- min(newv7, na.rm = TRUE)


f8 <- ecdf(v8)
y8 <- f8(df8)
inv_ecdf <- function(f8){
  x <- environment(f8)$x
  y <- environment(f8)$y
  approxfun(y, x)
}
g8 <- inv_ecdf(f8)

newv8 <- g8(df8)

newv8[is.na(newv8)] <- min(newv8, na.rm = TRUE)

newdata <- data.frame(newv7,newv8)




########
######## Graphic plots (2-D) of the original and the generated Data
########

par(mfrow = c(1, 2))
# Plot the data for a visual comparison
plot(data$v7, data$v8, main = 'Original vs MVD simulated', col = "blue")
points(newdata[,1], newdata[,2], col = 'red')
legend('bottomleft', c('Observed', 'Simulated'), col = c('blue', 'red'), pch=21)




########
######## Export all data
########

alldata <- data.frame(data, df, newdata,copulacdf)

write.table(alldata,file='alldata.csv',sep=',')




scatterplot3d(v[,1],v[,2], cdf_mvd, color="red", main="CDF", xlab = "u1", ylab="u2", zlab="pMvdc",pch=".")

