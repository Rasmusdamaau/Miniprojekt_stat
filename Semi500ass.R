n <- 500
x_i <- 1:n
waldtest <- c(seq(1,n,1))
conff1 <- matrix(nrow=n, ncol=2)
limits <- c(seq(0,n,1))*0
errors <- c(seq(0,n,1))*0
for (ii in 1:500) {
  yi <- rnorm(n, mean=0, sd=1)
  beta <- function(b,c1,c2 = 0) {
    b ^ (c1 * x_i + c2)
  }
  
  i <- function(b) {
    a <- sum(yi)/sum(beta(b,1))
    ma <- sum(a ^ 2 * 2 * (x_i-1) * beta(b,2,-2) * x_i)-sum(a ^ 2 * x_i * (x_i-1) * beta(b,2,-2))
    ma^(-1)
  }
  
  intim <- function(b, funk1, funk2) {
    score <- funk1 
    fish <- funk2 
    bk <- b
    itt <-0
    itt1 <- 0
    itt2 <- 0 
    while(abs(score(bk) * fish(bk))>0.001) {
      itt<- itt + 1
      bk = bk + fish(bk) * score(bk)
      if (itt %% 500 == 0)
        print(bk)
      if (itt>5000) {
        cat("break", "\n")
        itt1 <- 1
        break
      }
      if (!is.finite(score(bk))==TRUE | !is.finite((fish(bk)))==TRUE) {
        cat("ERROR!", "\n")
        itt2 <- 1
        break
      }
    }
    c(bk,itt, itt1, itt2)
  }
  betahat <-  intim(2, s, i)[1]
  ahat <- a <- sum(yi)/sum(beta(betahat,1,0))
  j22y <- 1 / (sum(ahat * (2 * x_i - 1) * beta(betahat,2,-1) * x_i) -sum(ahat * x_i *(x_i - 1)* beta(betahat, 1, 0 ) * yi))
  waldtest[ii] <- (betahat-1) ^ 2 / j22y
  conff1[ii,1] <- betahat - 1.96 * sqrt(i(betahat))
  conff1[ii,2] <- betahat + 1.96 * sqrt(i(betahat))
  limits[ii] <- intim(2, s, i)[3]
  errors[ii] <- intim(2, s, i)[4]
  waldtest
  conff1
  cat("Iteration=", ii, "Limits=", sum(limits),  "Errors=",sum(errors) , "\n") 
}
