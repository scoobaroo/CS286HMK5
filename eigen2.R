
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

svdA$u <- svdA$u[,-4]
delta2 <- t(svdA$u) %*% A
delta2<-delta2[-4,]

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

e1y1 <- dist(cbind(w1,delta2[,1]))
e2y1 <- dist(cbind(w1,delta2[,2]))
e3y1 <- dist(cbind(w1,delta2[,3]))
e4y1 <- dist(cbind(w1,delta2[,4]))
ey1 <- cbind(e1y1,e2y1,e3y1,e4y1)
miney1 <-min(ey1)

e1y2 <- dist(cbind(w2,delta2[,1]))
e2y2 <- dist(cbind(w2,delta2[,2]))
e3y2 <- dist(cbind(w2,delta2[,3]))
e4y2 <- dist(cbind(w2,delta2[,4]))
ey2 <- cbind(e1y2,e2y2,e3y2,e4y2)
miney2 <-min(ey2)

e1y3 <- dist(cbind(w3,delta2[,1]))
e2y3 <- dist(cbind(w3,delta2[,2]))
e3y3 <- dist(cbind(w3,delta2[,3]))
e4y3 <- dist(cbind(w3,delta2[,4]))
ey3 <- cbind(e1y3,e2y3,e3y3,e4y3)
miney3 <-min(ey3)

e1y4 <- dist(cbind(w4,delta2[,1]))
e2y4 <- dist(cbind(w4,delta2[,2]))
e3y4 <- dist(cbind(w4,delta2[,3]))
e4y4 <- dist(cbind(w4,delta2[,4]))
ey4 <- cbind(e1y4,e2y4,e3y4,e4y4)
miney4 <-min(ey4)