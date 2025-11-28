df <- with(mtcars, data.frame(y=mpg, x1=disp, x2=hp, x3=wt))

nll_lm <- function(data, par){
  n <- nrow(data)
  y <- data$y
  x <- cbind(1, data$x1, data$x2, data$x3)
  beta <- par[1:4]
  sigma <- par[5]

  resid <- y - x %*% beta
  nll <- (n/2) * log(2 * pi * sigma^2) + (1/(2 * sigma^2)) * sum(resid^2)

  return(nll)
}

initial <- c(mean(df$y), 0, 0, 0, sd(df$y))

est <- optim(par = initial,
             fn = nll_lm,
             data = df,
             method = "L-BFGS-B", # ensures that sigma is positive
             lower = c(-Inf, -Inf, -Inf, -Inf, 0.001))

y <- df$y
x <- cbind(1, df$x1, df$x2, df$x3) # design matrix
beta_matrix <- as.vector(solve(crossprod(x), crossprod(x,y)))

n <- nrow(df)
p <- 4 # number of beta parameters

residual_matrix <- y-x %*% beta_matrix
RSS <- sum(residual_matrix^2) # residual sum of squares
sqrt(RSS/(n-p))



est <- optim(par = initial,
             fn = nll_lm,
             data = df,
             method = "L-BFGS-B",
             lower = c(-Inf, -Inf, -Inf, -Inf,0.001),
             hessian = TRUE)
varcov <- solve(est$hessian)
sqrt(diag(varcov))[1:4]
