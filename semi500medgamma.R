library("tidyverse")
x_icount <- 50
n <- 100
x_i <- 1:x_icount
waldtestb <- rep(0,n)
conff1b <- matrix(nrow=n, ncol=2)
waldtesta <- rep(0,n)
conff1a <- matrix(nrow=n, ncol=2)
limits <- rep(0,n)
errors <- rep(0,n)
betahat <- rep(0,n)
ahat <- rep(0,n)
waldtestg <- rep(0,n)
gammahat <- rep(0,n)
conff1g <- matrix(nrow=n, ncol=2)
resulthypotese <- rep(0,n)
crit = pchisq(0.05,df=2,lower.tail = FALSE)
beta <- function(b,c1 = 1,c2 = 0) {
  b ^ (c1 * x_i + c2)
}

for (ii in 1:n) {
  yi <- rnorm(x_icount, mean=0, sd=1)
  likelihood <- function(a,b) {
    1 / sqrt(2*pi)^x_icount * prod(exp (( a * b ^x_i -yi)^2)/2)
  }
  
  s <- function(b) {
    a_s <- sum(beta(b))
    b_s <- (sum(beta(b))) / n * sum(beta(b)) 
    c_s <- (sum(beta(b))) / n * sum(yi)
    d_s <- sum(yi*beta(b))
    sahat <- (c_s-d_s)/(b_s-a_s)
    gamma <- 1 / n * sum(sahat*beta(b)-yi) 
    - sum((sahat * beta(b) + gamma - yi ) * beta(b,1,-1) * x_i* sahat)
  }
  
  i <- function(b) {
    a_s <- sum(beta(b))
    b_s <- (sum(beta(b))) / n * sum(beta(b)) 
    c_s <- (sum(beta(b))) / n * sum(yi)
    d_s <- sum(yi*beta(b))
    sahat <- (c_s-d_s)/(b_s-a_s)
    gamma <- 1 / n * sum(sahat*beta(b))- sum(yi)
    imi <- sum(x_i * sahat * beta(b,1,-2) * ((2*x_i-1) * sahat *beta(b) + 
             gamma * (x_i-1) - ((sahat * beta(b) + gamma)) *(x_i -1)  ))
    imi^-1
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
      if (!is.finite(score(bk))==TRUE | !is.finite((fish(bk)))==TRUE) {
        cat("ERROR!", "\n")
        itt2 <- 1
        break
      }
      if (itt>5000) {
        cat("break", "\n")
        itt1 <- 1
        break(intim)
      }
    }
    c(bk,itt, itt1, itt2)
  }
  betahat[ii] <-  intim(1.2, s, i)[1] 
  ahat_s <- sum(beta(betahat[ii]))
  bhat_s <- (sum(beta(betahat[ii]))) / n * sum(beta(betahat[ii])) 
  chat_s <- (sum(beta(betahat[ii]))) / n * sum(yi)
  dhat_s <- sum(yi*beta(betahat[ii]))
  ahat[ii] <- (chat_s-dhat_s)/(bhat_s-ahat_s)
  gammahat[ii] <- 1 / n * sum(ahat[ii]*beta(betahat[ii]))- sum(yi)
  lmle <- likelihood(ahat[ii],betahat[ii])
  lognorm <- function(a,b) {
    -2*log(lmle/likelihood(0,1))
  }
  
  resulthypotese[ii] <- lognorm(ahat[ii], betahat[ii])<=crit
  j22y <- 1 / (sum(ahat[ii] * (2 * x_i - 1) * beta(betahat[ii],2,-1) * x_i) -sum(ahat[ii] * x_i *(x_i - 1)* beta(betahat[ii], 1, 0 ) * yi))
  j11y <- 1 / (sum(beta(betahat[ii],2,0)))
  j33y <- 1 / x_icount
  waldtesta[ii]  <- ((ahat[ii]-0)^2)/ j11y
  waldtestb[ii] <- (betahat[ii]-1) ^ 2 / j22y
  waldtestg[ii] <- (gammahat[ii]-0)^2/ j33y
  conff1b[ii,1] <- betahat[ii] - 1.96 * sqrt(abs(j22y))
  conff1b[ii,2] <- betahat[ii] + 1.96 * sqrt(abs(j22y))
  conff1a[ii,1] <- ahat[ii] - 1.96 * sqrt(abs(j11y))
  conff1a[ii,2] <- ahat[ii] + 1.96 * sqrt(abs(j11y))
  conff1g[ii,1] <- gammahat[ii] - 1.96 * sqrt(abs(j33y))
  conff1g[ii,2] <- gammahat[ii] + 1.96 * sqrt(abs(j33y))
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
cat("Hvilke gammahat waldtest ligger inde for confidence intervallet", "\n")
confresulg <- which(conff1g[,1] < waldtestg & conff1g[,2] > waldtestg); confresulg
length(confresulg)
cat("Number of accepted H0 hypothesis", "\n")
sum(resulthypotese)
ahat_plot <- data.frame(x=1:n, y= ahat)
betahat_plot <- data.frame(x=1:n, y= betahat)
ahat_betahat <- data.frame(x= ahat, y=betahat)
ggplot(ahat_plot, aes(x,y)) +
  geom_point() + 
  scale_y_continuous(limits=c(-1,1))
ggplot(betahat_plot, aes(x,y)) +
  geom_point() +
  scale_y_continuous(limits=c(0,5))
ggplot(ahat_betahat, aes(x,y)) +
  geom_point(size = 0.7) + 
  scale_y_continuous(limits=c(0.6,2.5)) +
  scale_x_continuous(limits=c(-0.8,0.8))

