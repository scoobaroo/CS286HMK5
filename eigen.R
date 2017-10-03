#problem 4
A <- matrix(c(5,5,8,5,10,9,8,9,13),3,3, byrow=TRUE)
A2<- matrix(c(2.5,2.5,4,2.5,5,4.5,4,4.5,6.5),3,3,byrow=TRUE)

ev<- eigen(A2)
values <- ev$values
vectors <- ev$vectors
R(A)
A2.svd <- svd(A2)

#problem 9

A3 <- matrix(c(.25,.25,-.75,.25,-.75, 1.25, -1.75, 1.25, -1.25,-.25,1.75,-.25,1,0,1,-2,-1,1,-1,1,0,-1,0,1),6,4,byrow=TRUE)
A3.cov <- cov(A3)
eigenvalues <- eigen(A3.cov)
eigenvalues$values
eigenvalues$vectors

#problem 10
#create eigenvectors
u1 <- matrix(c(.1641,.6278,-.2604,-.5389,.4637,.0752),6,1,byrow=TRUE)
u2 <- matrix(c(.2443,.107,-.8017,.4277,-.1373,-.2904),6,1,byrow=TRUE)
#create Y vectors
Y1 <- matrix(c(2,3,1,0,3,2),6,1)
Y2 <- matrix(c(-4,-5,0,3,1,-2),6,1)
Y3 <- matrix(c(2,3,0,1,3,2),6,1)
Y4 <- matrix(c(3,2,1,0,3,2),6,1)
#create scoring matrix
delta <- matrix(c(-1.1069,1.2794,-2.68,2.5076,1.5480,0.5484,-1.2085,-.8879),2,4,byrow=TRUE)
mu <- matrix(c(7/4,7/4,5/4,2,2,1),6,1)
Y1bar <- Y1 - mu
#create weight vector
W1 <- matrix(c(c(Y1bar) %*% c(u1),c(Y1bar) %*%c(u2)),2,1)
#calculate epsilons for each Y by calculating distance from W1 to each column of delta
e1Y1 <- dist(cbind(W1,delta[,1]))
e2Y1 <- dist(cbind(W1,delta[,2]))
e3Y1 <- dist(cbind(W1,delta[,3]))
e4Y1 <- dist(cbind(W1,delta[,4]))
eY1 <- cbind(e1Y1,e2Y1,e3Y1,e4Y1)
#find the minimum epsilon
mineY1 <- min(eY1)
#same steps repeates for Y2bar,Y3bar,and Y4bar
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

#problem 11 (second version)

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

#problem 13
#create malware samples
M1 <- matrix(c(1,-1,1,-1,-1,1),6,1)
M2 <- matrix(c(-2,2,2,-1,-2,2),6,1)
M3 <- matrix(c(1,3,0,1,3,1),6,1)
M4 <- matrix(c(2,3,1,1,-2,0),6,1)
M1bar <- M1-mu
M2bar <- M2-mu
M3bar <- M3-mu
M4bar <- M4-mu
#project samples onto eigen vector space
w1new <- t(svdnewA$u) %*% M1bar
w2new <- t(svdnewA$u) %*% M2bar
w3new <- t(svdnewA$u) %*% M3bar

e1M1 <- dist(cbind(w1new,delta[,1]))
e2M1 <- dist(cbind(w1new,delta[,2]))
e3M1 <- dist(cbind(w1new,delta[,3]))
e4M1 <- dist(cbind(w1new,delta[,4]))
eM1new <- cbind(e1M1,e2M1,e3M1,e4M1)
minepsilonM1 <-min(eM1new)

e1M2 <- dist(cbind(w2new,delta[,1]))
e2M2 <- dist(cbind(w2new,delta[,2]))
e3M2 <- dist(cbind(w2new,delta[,3]))
e4M2 <- dist(cbind(w2new,delta[,4]))
eM2new <- cbind(e1M2,e2M2,e3M2,e4M2)
minepsilonM2 <-min(eM2new)

e1M3 <- dist(cbind(w3new,delta[,1]))
e2M3 <- dist(cbind(w3new,delta[,2]))
e3M3 <- dist(cbind(w3new,delta[,3]))
e4M3 <- dist(cbind(w3new,delta[,4]))
eM3new <- cbind(e1M3,e2M3,e3M3,e4M3)
minepsilonM3 <-min(eM3new)

e1M4 <- dist(cbind(w4new,delta[,1]))
e2M4 <- dist(cbind(w4new,delta[,2]))
e3M4 <- dist(cbind(w4new,delta[,3]))
e4M4 <- dist(cbind(w4new,delta[,4]))
eM4new <- cbind(e1M4,e2M4,e3M4,e4M4)
minepsilonM4 <-min(eM4new)
#project malware samples onto eigenvectors

#create benign samples
B1 <- matrix(c(-1,2,1,2,-1,0),6,1)
B2 <- matrix(c(-2,1,2,3,2,1),6,1)
B3 <- matrix(c(-1,3,0,1,3,-1),6,1)
B4 <- matrix(c(0,2,3,1,1,-2),6,1)
b <- cbind(B1,B2,B3,B4)
Bbar <- rowMeans(b)
B1bar <- B1 - Bbar
B2bar <- B2 - Bbar
B3bar <- B3 - Bbar
B4bar <- B4 - Bbar
Bnew <- cbind(B1bar,B2bar,B3bar,B4bar)
#perform SVD on Bnew to obtain eigenvectors and eigenvalues
svdBnew <- svd(Bnew)
#project benign samples onto eigenvectors to obtain weights
w1Bnew <- t(svdBnew$u) %*% B1bar
w2Bnew <- t(svdBnew$u) %*% B2bar
w3Bnew <- t(svdBnew$u) %*% B3bar
w4Bnew <- t(svdBnew$u) %*% B4bar
#obtain scoring matrix for benign samples
deltaBenign <- t(svdBnew$u) %*% Bnew

#create Y vectors for part c
Y1new <- matrix(c(1,5,1,5,5,1),6,1)
Y2new <- matrix(c(-2,3,2,3,0,2),6,1)
Y3new <- matrix(c(2,-3,2,3,0,0),6,1)
Y4new <- matrix(c(2,-2,2,2,-1,1),6,1)
#subtract malware means and benign means from each Y element
Y1newbar <- Y1new - mu
Y1newbarB <- Y1new - Bbar
Y2newbar <- Y2new - mu
Y2newbarB <- Y2new - Bbar
Y3newbar <- Y3new - mu
Y3newbarB <- Y3new - Bbar
Y4newbar <- Y4new - mu
Y4newbarB <- Y4new - Bbar
#project test samples onto eigenvectors for both benign and malware and score
w1possiblemalware <- t(svdnewA$u) %*% Y1newbar
w2possiblemalware <- t(svdnewA$u) %*% Y2newbar
w3possiblemalware <- t(svdnewA$u) %*% Y3newbar
w4possiblemalware <- t(svdnewA$u) %*% Y4newbar
w1possiblebenign <- t(svdBnew$u) %*% Y1newbarB
w2possiblebenign <- t(svdBnew$u) %*% Y2newbarB
w3possiblebenign <- t(svdBnew$u) %*% Y3newbarB
w4possiblebenign <- t(svdBnew$u) %*% Y4newbarB
#delete the last 2 rows from deltaBenign to match malware scoring matrix
deltaBenign <- deltaBenign[-4,]
deltaBenign <- deltaBenign[-3,]
#we are calculating the distances to delta and deltaBenign and taking the minimum to classifty
e1m1 <- dist(cbind(w1possiblemalware[1:2],delta[,1]))
e1B1 <- dist(cbind(w1possiblebenign[1:2],deltaBenign[,1]))
#min is e1m1 so it is classified as malware

e2m2 <- dist(cbind(w2possiblemalware[1:2],delta[,2]))
e2B2 <- dist(cbind(w2possiblebenign[1:2],deltaBenign[,2]))
#min is e2B2 so it is classfied as benign

e3m3 <- dist(cbind(w3possiblemalware[1:2],delta[,3]))
e3B3 <- dist(cbind(w3possiblebenign[1:2],deltaBenign[,3]))
#min is e3m3 so it is classified as malware

e4m4 <- dist(cbind(w4possiblemalware[1:2],delta[,4]))
e4B4 <- dist(cbind(w4possiblebenign[1:2],deltaBenign[,4]))
#min is e4m4 so it is classified as malware

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




