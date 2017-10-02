#problem 4
A <- matrix(c(5, 5, 8, 5, 10, 9, 8, 9, 13), 3, 3, byrow=TRUE)
A2<- matrix(c(2.5,2.5,4,2.5,5,4.5,4,4.5,6.5),3,3, byrow=TRUE)

ev<- eigen(A2)
values <- ev$values
vectors <- ev$vectors
R(A)
A2.svd <- svd(A2)

#problem 9

A3 <- matrix(c(.25,.25,-.75,.25,-.75, 1.25, -1.75, 1.25, -1.25,-.25,1.75,-.25,1,0,1,-2,-1,1,-1,1,0,-1,0,1), 6, 4, byrow=TRUE)
A3.cov <- cov(A3)
eigenvalues <- eigen(A3.cov)
eigenvalues$values
eigenvalues$vectors

#problem 10
u1 <- matrix(c(.1641,.6278,-.2604,-.5389,.4637,.0752), 6, 1, byrow=TRUE)
u2 <- matrix(c(.2443,.107,-.8017,.4277,-.1373,-.2904),6,1,byrow=TRUE)
Y1 <- matrix(c(2,3,1,0,3,2),6,1)
Y2 <- matrix(c(-4,-5,0,3,1,-2),6,1)
Y3 <- matrix(c(2,3,0,1,3,2),6,1)
Y4 <- matrix(c(3,2,1,0,3,2),6,1)
delta <- matrix(c(-1.1069,1.2794,-2.68,2.5076,1.5480,0.5484,-1.2085,-.8879),2,4,byrow=TRUE)
mu <- matrix(c(7/4,7/4,5/4,2,2,1),6,1)
Y1bar <- Y1 - mu
W1 <- matrix(c(c(Y1bar) %*% c(u1),c(Y1bar) %*%c(u2)),2,1)
e1Y1 <- dist(cbind(W1,delta[,1]))
e2Y1 <- dist(cbind(W1,delta[,2]))
e3Y1 <- dist(cbind(W1,delta[,3]))
e4Y1 <- dist(cbind(W1,delta[,4]))
eY1 <- cbind(e1Y1,e2Y1,e3Y1,e4Y1)
mineY1 <- min(eY1)
Y2bar <- Y2 - mu
W2 <- matrix(c(c(Y2bar) %*% c(u1),c(Y2bar) %*%c(u2)),2,1)
e1Y2 <- dist(cbind(W2,delta[,1]))
e2Y2 <- dist(cbind(W2,delta[,2]))
e3Y2 <- dist(cbind(W2,delta[,3]))
e4Y2 <- dist(cbind(W2,delta[,4]))
eY2 <- cbind(e1Y2,e2Y2,e3Y2,e4Y2)
mineY2 <- min(eY2)
Y3bar <- Y3 - mu
W3 <- matrix(c(c(Y3bar) %*% c(u1),c(Y3bar) %*%c(u2)),2,1)
e1Y3 <- dist(cbind(W3,delta[,1]))
e2Y3 <- dist(cbind(W3,delta[,2]))
e3Y3 <- dist(cbind(W3,delta[,3]))
e4Y3 <- dist(cbind(W3,delta[,4]))
eY3 <- cbind(e1Y3,e2Y3,e3Y3,e4Y3)
mineY3 <- min(eY3)
Y4bar <- Y4 - mu
W4 <- matrix(c(c(Y4bar) %*% c(u1),c(Y4bar) %*%c(u2)),2,1)
e1Y4 <- dist(cbind(W4,delta[,1]))
e2Y4 <- dist(cbind(W4,delta[,2]))
e3Y4 <- dist(cbind(W4,delta[,3]))
e4Y4 <- dist(cbind(W4,delta[,4]))
eY4 <- cbind(e1Y4,e2Y4,e3Y4,e4Y4)
mineY4 <-min(eY4)

#problem 11
X1 <- matrix(c(2,-1,0,1,1,-3,5,2),8,1)
X2 <- matrix(c(-2,3,2,3,0,2,-1,1),8,1)
X3 <- matrix(c(-1,3,3,1,-1,4,5,2),8,1)
X4 <- matrix(c(3,-1,0,3,2,-1,3,0),8,1)
B <- cbind(X1,X2,X3,X4)
#generate the row means
mu2 <- rowMeans(B)
A <- B-mu2
#make covariance matrix from the A matrix
C <- cov(A)
eigenv <- eigen(C)
eigenv$vectors
eigenv$values
#calculate singular value decomposition of A, which gives us U, S, and V
svdA <- svd(A)
#multiply the transpose of U by A to obtain the scoring matrix. We only keep the top 3 columns because we are only interested in the 3 most significant eigenvectors of C
delta2 <- t(svdA$u) %*% A
delta2<-delta[-4,]

# problem 11 part b
y1 <- matrix(c(1,5,1,5,5,1,1,3),8,1)
y2 <- matrix(c(-2,3,2,3,0,2,-1,1),8,1)
y3 <- matrix(c(2,-3,2,3,0,0,2,-1),8,1)
y4 <- matrix(c(2,-2,2,2,-1,1,2,2),8,1)

y1bar <- y1 - mu2
y2bar <- y2 - mu2
y3bar <- y3 - mu2
y4bar <- y4 - mu2

w1 <- t(svdA$u) %*% y1bar
w2 <- t(svdA$u) %*% y2bar
w3 <- t(svdA$u) %*% y3bar
w4 <- t(svdA$u) %*% y4bar

#problem 14
U<-matrix(c(0.1641,0.2443,-.0710,.6278,.107,.2934,-.2604,-.8017,.3952,-.5389,.4277,.3439,.4637,-.1373,.3644,.0752,-.2904,-.7083),6,3,byrow=TRUE)
S<-matrix(c(4.0414,0,0,0,2.2239,0,0,0,1.7237), 3, 3, byrow=TRUE)
V<-matrix(c(-.2739,.6961,-.4364,.3166,.2466,.7674,-.6631,-.5434,.1224,.6205,-.3993,-.4534),4,3,byrow=TRUE)

transposeV <- t(V)

A4 <- U %*% S %*% transposeV

#problem 17

u1 <- matrix(c(.1641,.6278,-.2604,-.5389,.4637,.0752), 6, 1, byrow=TRUE)
u2 <- matrix(c(.2443,.107,-.8017,.4277,-.1373,-.2904),6,1,byrow=TRUE)
u3 <- matrix(c(-.071,.2934,.3952,.3439,.3644,-.7083),6,1,byrow=TRUE)
lambda1 <- 4.0833
lambda2 <- 1.2364
lambda3 <- 0.7428

loadingu1 <- sqrt(lambda1) * u1
loadingu2 <- sqrt(lambda2) * u2
loadingu3 <- sqrt(lambda3) * u3




