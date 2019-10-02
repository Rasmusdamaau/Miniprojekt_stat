library("tidyverse")
x_icount <- 10
n <- 50
x_i <- 1:x_icount
waldtestb <- rep(0,n)
conff1b <- matrix(nrow=n, ncol=2)
waldtesta <- rep(0,n)
conff1a <- matrix(nrow=n, ncol=2)
limits <- rep(0,n)
errors <- rep(0,n)
betahat <- rep(0,n)
ahat <- rep(0,n)
for (ii in 1:n) {
  yi <- rnorm(x_icount, mean=0, sd=1)
  beta <- function(b,c1,c2 = 0) {
    b ^ (c1 * x_i + c2)
  }
  
  s <- function(b) {
    a <- sum(yi)/sum(beta(b,1,0))
    - sum( a^2 * beta(b,2,-1) * x_i)  +  sum(a * x_i * beta(b,1,-1) * yi)
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
  betahat[ii] <-  intim(2, s, i)[1]
  ahat[ii] <- sum(yi)/sum(beta(betahat[ii],1,0))
  j22y <- 1 / (sum(ahat[ii] * (2 * x_i - 1) * beta(betahat[ii],2,-1) * x_i) -sum(ahat[ii] * x_i *(x_i - 1)* beta(betahat[ii], 1, 0 ) * yi))
  j11y <- 1 / (sum(beta(betahat[ii],2,0)))
  waldtesta[ii]  <- ((ahat[ii]-0)^2)/ j11y
  waldtestb[ii] <- (betahat[ii]-1) ^ 2 / j22y
  conff1b[ii,1] <- betahat[ii] - 1.96 * sqrt(abs(j22y))
  conff1b[ii,2] <- betahat[ii] + 1.96 * sqrt(abs(j22y))
  conff1a[ii,1] <- ahat[ii] - 1.96 * sqrt(abs(j11y))
  conff1a[ii,2] <- ahat[ii] + 1.96 * sqrt(abs(j11y))
  limits[ii] <- intim(2, s, i)[3]
  errors[ii] <- intim(2, s, i)[4]
  cat("Iteration=", ii, "Limits=", sum(limits),  "Errors=",sum(errors) , "\n") 
}
cat("Hvilke ahat waldtest ligger inde for confidence intervallet", "\n")
confresula <- which(conff1a[,1] < waldtesta & conff1a[,2] > waldtesta); confresula
length(confresula)
cat("Hvilke betahat waldtest ligger inde for confidence intervallet", "\n")
confresulb <- which(conff1b[,1] < waldtestb & conff1b[,2] > waldtestb); confresulb
length(confresulb)
ahat[which(ahat<=1e-16)]=0
ggplot(as.data.frame(ahat), aes(y=as.data.frame(ahat))) +
  geom_area(stat="bin") +
  scale_y_continuous()

view(ahat)
