n <- 10
x_i <- seq(1,n,1)
yi <- c()
for (i in 1:n) {
  yi[i] <- rnorm(1, mean=runif(1, min=-1, max=1), sd=1)
}
sum(yi)
beta <- function(b,c1,c2) {
  b^(c1*x_i+c2)
}

s <- function(b) {
  a <- sum(yi)/sum(beta(b,1,0))
  -a^2*sum(beta(b,2,-1))*sum(x_i)+a*sum(x_i)*sum(beta(b,1,-1))*sum(yi)

}
s(2)
ma <- matrix(nrow=2, ncol=2)
i <- function(b) {
  a <- sum(yi)/sum(beta(b,1,0))
  ma[1,1] <- sum(beta(b,2,0))
  ma[1,2] <- a*2*sum(x_i)*sum(beta(b,2,-1))-sum(x_i)*sum(beta(b,2,-1))*a
  ma[2,1] <- a*2*sum(x_i)*sum(beta(b,2,-1))-a*sum(x_i)*sum(beta(b,2,-1))
  ma[2,2] <- a^2*sum(2*x_i-1)*sum(beta(b,2,-2))*sum(x_i)-a^2*sum(x_i)*sum(x_i-1)*sum(beta(b,2,-2))
  mainv <- ma^(-1)
  mainv[2,2]
}
i(10)
intim <- function(b) {
  bk <- b
  itt <-0
  while(abs(bk-(bk+i(bk)*s(bk)))>0.01) {
    itt<- itt+1
    bk=bk+i(bk)*s(bk)
    if (itt%%20==0)
    print(bk)
  }
  cat(bk,itt)
}
intim(5)
