library("tidyverse")
x_icount <- 20
n <- 1
x_i <- 1:x_icount
a_nul <- 0
beta_nul <- 1
gamma_nul <- 0
p_waldtest_a <- rep(0,n)
p_waldtest_b <- rep(0,n)
p_waldtest_g <- rep(0,n)
conff1b <- matrix(nrow=n, ncol=2)
conff1a <- matrix(nrow=n, ncol=2)
conff1g <- matrix(nrow=n, ncol=2)
limits <- rep(0,n)
errors <- rep(0,n)
betahat <- rep(0,n)
ahat <- rep(0,n)
gammahat <- rep(0,n)
resulthypotese <- rep(0,n)
bstar <- 1.1
crit = pchisq(0.05,df=3,lower.tail = FALSE)
beta <- function(b,c1 = 1,c2 = 0) {
  b ^ (c1 * x_i + c2)
}

for (ii in 1:n) {
  yi <- c(-0.8, 1.15, 1.82, 2.89, 1.43, 0.59, 1.09, 3.00, 2.81, 3.02, 2.39, 3.50, 3.06, 2.41, 2.81, 3.95, 4.19, 5.94, 7.06, 4.55)
  likelihood <- function(a,b,g) {
    1 / sqrt(2*pi)^x_icount * prod(exp (( a * b ^x_i -yi + g)^2)/2)
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
  resultatintim <- intim(bstar, s, i)
  betahat[ii] <-  resultatintim[1] 
  ahat_s <- sum(beta(betahat[ii]))
  bhat_s <- (sum(beta(betahat[ii]))) / n * sum(beta(betahat[ii])) 
  chat_s <- (sum(beta(betahat[ii]))) / n * sum(yi)
  dhat_s <- sum(yi*beta(betahat[ii]))
  ahat[ii] <- (chat_s-dhat_s)/(bhat_s-ahat_s)
  gammahat[ii] <- 1 / n * sum(ahat[ii]*beta(betahat[ii]))- sum(yi)
  lmle <- likelihood(ahat[ii],betahat[ii], gammahat[ii])
  lognorm <- function(a,b,g) {
    -2*log(likelihood(a,b,g)/lmle)
  }
  
  resulthypotese[ii] <- lognorm(ahat[ii], betahat[ii], gamma_nul)<=crit #Beta og alpha afhænger af 
  #gamma, derfor er dette ikke korrekt.
  j22y <- 1 / (sum(ahat[ii] * (2 * x_i - 1) * beta(betahat[ii],2,-1) * x_i) -sum(ahat[ii] * x_i *(x_i - 1)* beta(betahat[ii], 1, 0 ) * yi))
  j11y <- 1 / (sum(beta(betahat[ii],2,0)))
  j33y <- 1 / x_icount
  waldtesta <- ((ahat[ii]-a_nul)^2)/ j11y
  p_waldtest_a[ii] <- 2*pnorm(abs(waldtesta), lower.tail = FALSE)
  waldtestb <- (betahat[ii]-beta_nul) ^ 2 / j22y
  p_waldtest_b[ii] <- 2*pnorm(abs(waldtestb), lower.tail = FALSE)
  waldtestg <- (gammahat[ii]-gamma_nul)^2/ j33y
  p_waldtest_g[ii] <- 2*pnorm(abs(waldtestg), lower.tail = FALSE)
  conff1b[ii,1] <- betahat[ii] - 1.96 * sqrt(abs(j22y))
  conff1b[ii,2] <- betahat[ii] + 1.96 * sqrt(abs(j22y))
  conff1a[ii,1] <- ahat[ii] - 1.96 * sqrt(abs(j11y))
  conff1a[ii,2] <- ahat[ii] + 1.96 * sqrt(abs(j11y))
  conff1g[ii,1] <- gammahat[ii] - 1.96 * sqrt(abs(j33y))
  conff1g[ii,2] <- gammahat[ii] + 1.96 * sqrt(abs(j33y))
  limits[ii] <- resultatintim[3]
  errors[ii] <- resultatintim[4]
  cat("Iteration=", ii, "Limits=", sum(limits),  "Errors=",sum(errors) , "\n") 
}
cat("Hvilke ahat waldtest ligger inde for confidence intervallet", "\n")
confresula <- which(conff1a[,1] < a_nul & conff1a[,2] > a_nul); confresula
length(confresula)
cat("Hvilke betahat waldtest ligger inde for confidence intervallet", "\n")
confresulb <- which(conff1b[,1] < beta_nul & conff1b[,2] > beta_nul); confresulb
length(confresulb)
cat("Hvilke gammahat waldtest ligger inde for confidence intervallet", "\n")
confresulg <- which(conff1g[,1] < gamma_nul & conff1g[,2] > gamma_nul); confresulg
length(confresulg)
cat("Number of accepted H0 hypothesis", "\n")
sum(resulthypotese)
ahat_plot <- data.frame(x=1:n, y= ahat)
betahat_plot <- data.frame(x=1:n, y= betahat)
gammahat_plot <- data.frame(x=1:n, y=gammahat)
ahat_betahat <- data.frame(x= ahat, y=betahat)
ggplot(ahat_plot, aes(x,y)) +
  geom_point() + 
  labs(x = "index", y= "alpha")
ggplot(betahat_plot, aes(x,y)) +
  geom_point() +
  labs(x = "index", y= "beta")
ggplot(gammahat_plot, aes(x,y)) +
  geom_point() +
  labs(x = "index", y= "gamma")
ggplot(ahat_betahat, aes(x,y)) +
  geom_point(size = 0.7) + 
  labs(x = "alpha", y= "beta")