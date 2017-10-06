norm_vec <- function(x) sqrt(sum(x^2))
#problem 4 (for verification only)
#create a matrix
a <- matrix(c(1,2,3,1,2,3),3,2,byrow=TRUE)
#calculate covariance matrix
c <- 0.5 * a %*% t(a)
#this is hand-calculated covariance matrix, should match c
C2<- matrix(c(2.5,2.5,4,2.5,5,4.5,4,4.5,6.5),3,3,byrow=TRUE)
#values from svd and eigen should match
ev<- eigen(c)
values <- ev$values
vectors <- ev$vectors
a.svd <- svd(a)
#verify that each eigenvector is of length 1
norm_vec(a.svd$u[,2])
norm_vec(a.svd$u[,1])

#problem 9
A3 <- matrix(c(.25,.25,-.75,.25,-.75, 1.25, -1.75, 1.25, -1.25,-.25,1.75,-.25,1,0,1,-2,-1,1,-1,1,0,-1,0,1),6,4,byrow=TRUE)
p9A3.C <- 1/4 * A3 %*% t(A3)
sumC <- sum(diag(p9A3.C))
eigenvalues <- eigen(p9A3.C)
sumProjection <- sum(diag(eigenvalues$values))
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
#create average which is taken from text
mu <- matrix(c(7/4,7/4,5/4,2,2,1),6,1)
#subtract average from each Y to obtain Y bars
Y1bar <- Y1 - mu
Y2bar <- Y2 - mu
Y3bar <- Y3 - mu
Y4bar <- Y4 - mu
#create weight vector (only 2 rows since we are only interested in the 2 most significant principal components)
W1 <- matrix(c(c(Y1bar) %*% c(u1),c(Y1bar) %*% c(u2)),2,1)
#calculate epsilons for each Y by calculating distance from W1 to each column of delta
e1Y1 <- dist(cbind(W1,delta[,1]))
e2Y1 <- dist(cbind(W1,delta[,2]))
e3Y1 <- dist(cbind(W1,delta[,3]))
e4Y1 <- dist(cbind(W1,delta[,4]))
eY1 <- cbind(e1Y1,e2Y1,e3Y1,e4Y1)
#find the minimum epsilon
mineY1 <- min(eY1)
#repeat same steps for Y2bar,Y3bar,and Y4bar

W2 <- matrix(c(c(Y2bar) %*% c(u1),c(Y2bar) %*%c(u2)),2,1)
e1Y2 <- dist(cbind(W2,delta[,1]))
e2Y2 <- dist(cbind(W2,delta[,2]))
e3Y2 <- dist(cbind(W2,delta[,3]))
e4Y2 <- dist(cbind(W2,delta[,4]))
eY2 <- cbind(e1Y2,e2Y2,e3Y2,e4Y2)
mineY2 <- min(eY2)

W3 <- matrix(c(c(Y3bar) %*% c(u1),c(Y3bar) %*%c(u2)),2,1)
e1Y3 <- dist(cbind(W3,delta[,1]))
e2Y3 <- dist(cbind(W3,delta[,2]))
e3Y3 <- dist(cbind(W3,delta[,3]))
e4Y3 <- dist(cbind(W3,delta[,4]))
eY3 <- cbind(e1Y3,e2Y3,e3Y3,e4Y3)
mineY3 <- min(eY3)

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
#make A matrix by subtracting mu2 from B
A <- B-mu2
#make covariance matrix from the A matrix
p11C <- 1/4 * A %*% t(A)
EigenV <- eigen(p11C)
EigenV$vectors

#calculate singular value decomposition of A, which gives us U, S, and V
svdA <- svd(A)
#multiply the transpose of U by A to obtain the scoring matrix. We only keep the top 3 columns because we are only interested in the 3 most significant principal components
delta2 <- t(svdA$u) %*% A
delta2 <- delta2[-4,]

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

#removing last vector from weights since we are only interested in top 3 principal components
w1 <- w1[-4,]
w2 <- w2[-4,]
w3 <- w3[-4,]
w4 <- w4[-4,]

e1y1 <- dist(cbind(w1,delta2[,1]))
e2y1 <- dist(cbind(w1,delta2[,2]))
e3y1 <- dist(cbind(w1,delta2[,3]))
ey1 <- cbind(e1y1,e2y1,e3y1)
miney1 <-min(ey1)

e1y2 <- dist(cbind(w2,delta2[,1]))
e2y2 <- dist(cbind(w2,delta2[,2]))
e3y2 <- dist(cbind(w2,delta2[,3]))
ey2 <- cbind(e1y2,e2y2,e3y2)
miney2 <-min(ey2)

e1y3 <- dist(cbind(w3,delta2[,1]))
e2y3 <- dist(cbind(w3,delta2[,2]))
e3y3 <- dist(cbind(w3,delta2[,3]))
ey3 <- cbind(e1y3,e2y3,e3y3)
miney3 <-min(ey3)

e1y4 <- dist(cbind(w4,delta2[,1]))
e2y4 <- dist(cbind(w4,delta2[,2]))
e3y4 <- dist(cbind(w4,delta2[,3]))
ey4 <- cbind(e1y4,e2y4,e3y4)
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
#project samples onto eigen vector space using given eigen vectors in book( we are using these eigenvectors
# because we are using the same scoring matrix in 4.12)
u1 <- matrix(c(.1641,.6278,-.2604,-.5389,.4637,.0752), 6, 1, byrow=TRUE)
u2 <- matrix(c(.2443,.107,-.8017,.4277,-.1373,-.2904),6,1,byrow=TRUE)
u3 <- matrix(c(-.071,.2934,.3952,.3439,.3644,-.7083),6,1,byrow=TRUE)

#we are projecting the malware samples onto the eigenvectors
w1new <- matrix(c(t(u1) %*% M1bar , t(u2) %*% M1bar),2,1,byrow=TRUE)
w2new <- matrix(c(t(u1) %*% M2bar , t(u2) %*% M2bar),2,1,byrow=TRUE)
w3new <- matrix(c(t(u1) %*% M3bar , t(u2) %*% M3bar),2,1,byrow=TRUE)
w4new <- matrix(c(t(u1) %*% M4bar , t(u2) %*% M4bar),2,1,byrow=TRUE)

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
#perform SVD on Bnew to obtain eigenvectors and eigenvalues (Bnew is essentially A matrix for benign samples)
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
#project test samples onto respective eigenvectors for both benign and malware and score
w1possiblemalware <- matrix(c(t(u1) %*% Y1newbar, t(u2) %*% Y1newbar),2,1,byrow=TRUE)
w2possiblemalware <- matrix(c(t(u1) %*% Y2newbar, t(u2) %*% Y2newbar),2,1,byrow=TRUE)
w3possiblemalware <- matrix(c(t(u1) %*% Y3newbar, t(u2) %*% Y3newbar),2,1,byrow=TRUE)
w4possiblemalware <- matrix(c(t(u1) %*% Y4newbar, t(u2) %*% Y4newbar),2,1,byrow=TRUE)
# we are only looking at top 2 principal components
svdBnew$u <- svdBnew$u[,-3:-4]
w1possiblebenign <- t(svdBnew$u) %*% Y1newbarB
w2possiblebenign <- t(svdBnew$u) %*% Y2newbarB
w3possiblebenign <- t(svdBnew$u) %*% Y3newbarB
w4possiblebenign <- t(svdBnew$u) %*% Y4newbarB
#delete the last 2 rows from deltaBenign to match malware scoring matrix
deltaBenign <- deltaBenign[-4:-3,]

#we are calculating the distances to delta and deltaBenign and taking the minimum to classify
###FIRST SAMPLE
p13epsilon1m1a <- dist(cbind(w1possiblemalware[1:2],delta[,1]))
p13epsilon1m1b <- dist(cbind(w1possiblemalware[1:2],delta[,2]))
p13epsilon1m1c <- dist(cbind(w1possiblemalware[1:2],delta[,3]))
p13epsilon1m1d <- dist(cbind(w1possiblemalware[1:2],delta[,4]))
p13epsilonFirstSampleMalwareScore <- min(cbind(p13epsilon1m1a,p13epsilon1m1b,p13epsilon1m1c,p13epsilon1m1d))

p13epsilon1m1aB <- dist(cbind(w1possiblebenign[1:2],deltaBenign[,1]))
p13epsilon1m1bB <- dist(cbind(w1possiblebenign[1:2],deltaBenign[,2]))
p13epsilon1m1cB <- dist(cbind(w1possiblebenign[1:2],deltaBenign[,3]))
p13epsilon1m1dB <- dist(cbind(w1possiblebenign[1:2],deltaBenign[,4]))
p13epsilonFirstSampleBenignScore <- min(cbind(p13epsilon1m1aB,p13epsilon1m1bB,p13epsilon1m1cB,p13epsilon1m1dB))

#min is p13epsilonFirstSampleMalwareScore so it is classified as malware
###SECOND SAMPLE
p13epsilon1m1a2 <- dist(cbind(w2possiblemalware[1:2],delta[,1]))
p13epsilon1m1b2 <- dist(cbind(w2possiblemalware[1:2],delta[,2]))
p13epsilon1m1c2 <- dist(cbind(w2possiblemalware[1:2],delta[,3]))
p13epsilon1m1d2 <- dist(cbind(w2possiblemalware[1:2],delta[,4]))
p13epsilonSecondSampleMalwareScore <- min(cbind(p13epsilon1m1a2,p13epsilon1m1b2,p13epsilon1m1c2,p13epsilon1m1d2))

p13epsilon1m1aB2 <- dist(cbind(w2possiblebenign[1:2],deltaBenign[,1]))
p13epsilon1m1bB2 <- dist(cbind(w2possiblebenign[1:2],deltaBenign[,2]))
p13epsilon1m1cB2 <- dist(cbind(w2possiblebenign[1:2],deltaBenign[,3]))
p13epsilon1m1dB2 <- dist(cbind(w2possiblebenign[1:2],deltaBenign[,4]))
p13epsilonSecondSampleBenignScore <- min(cbind(p13epsilon1m1aB2,p13epsilon1m1bB2,p13epsilon1m1cB2,p13epsilon1m1dB2))

#min is p13epsilonSecondSampleMalwareScore so it is classfied as malware

###THIRD SAMPLE
p13epsilon1m1a3 <- dist(cbind(w3possiblemalware,delta[,1]))
p13epsilon1m1b3 <- dist(cbind(w3possiblemalware,delta[,2]))
p13epsilon1m1c3 <- dist(cbind(w3possiblemalware,delta[,3]))
p13epsilon1m1d3 <- dist(cbind(w3possiblemalware,delta[,4]))
p13epsilonThirdSampleMalwareScore <- min(cbind(p13epsilon1m1a3,p13epsilon1m1b3,p13epsilon1m1c3,p13epsilon1m1d3))

p13epsilon1m1aB3 <- dist(cbind(w3possiblebenign,deltaBenign[,1]))
p13epsilon1m1bB3 <- dist(cbind(w3possiblebenign,deltaBenign[,2]))
p13epsilon1m1cB3 <- dist(cbind(w3possiblebenign,deltaBenign[,3]))
p13epsilon1m1dB3 <- dist(cbind(w3possiblebenign,deltaBenign[,4]))
p13epsilonThirdSampleBenignScore <- min(cbind(p13epsilon1m1aB3,p13epsilon1m1bB3,p13epsilon1m1cB3,p13epsilon1m1dB3))

#min is p13epsilonThirdSampleBenignScore so it is classified as benign

###FOURTH SAMPLE
p13epsilon1m1a4 <- dist(cbind(w4possiblemalware,delta[,1]))
p13epsilon1m1b4 <- dist(cbind(w4possiblemalware,delta[,2]))
p13epsilon1m1c4 <- dist(cbind(w4possiblemalware,delta[,3]))
p13epsilon1m1d4 <- dist(cbind(w4possiblemalware,delta[,4]))
p13epsilonFourthSampleMalwareScore <- min(cbind(p13epsilon1m1a4,p13epsilon1m1b4,p13epsilon1m1c4,p13epsilon1m1d4))

p13epsilon1m1aB4 <- dist(cbind(w4possiblebenign,deltaBenign[,1]))
p13epsilon1m1bB4 <- dist(cbind(w4possiblebenign,deltaBenign[,2]))
p13epsilon1m1cB4 <- dist(cbind(w4possiblebenign,deltaBenign[,3]))
p13epsilon1m1dB4 <- dist(cbind(w4possiblebenign,deltaBenign[,4]))
p13epsilonFourthSampleBenignScore <- min(cbind(p13epsilon1m1aB4,p13epsilon1m1bB4,p13epsilon1m1cB4,p13epsilon1m1dB4))

#min is p13epsilonFourthSampleMalwareScore so it is classified as malware


#problem 14
U<-matrix(c(0.1641,0.2443,-.0710,.6278,.107,.2934,-.2604,-.8017,.3952,-.5389,.4277,.3439,.4637,-.1373,.3644,.0752,-.2904,-.7083),6,3,byrow=TRUE)
S<-matrix(c(4.0414,0,0,0,2.2239,0,0,0,1.7237), 3, 3, byrow=TRUE)
V<-matrix(c(-.2739,.6961,-.4364,.3166,.2466,.7674,-.6631,-.5434,.1224,.6205,-.3993,-.4534),4,3,byrow=TRUE)

transposeV <- t(V)

A4 <- U %*% S %*% transposeV

#problem 17

u1 <- matrix(c(.1641,.6278,-.2604,-.5389,.4637,.0752),6,1,byrow=TRUE)
u2 <- matrix(c(.2443,.107,-.8017,.4277,-.1373,-.2904),6,1,byrow=TRUE)
u3 <- matrix(c(-.071,.2934,.3952,.3439,.3644,-.7083),6,1,byrow=TRUE)
lambda1 <- 4.0833
lambda2 <- 1.2364
lambda3 <- 0.7428

loadingu1 <- sqrt(lambda1) * u1
loadingu2 <- sqrt(lambda2) * u2
loadingu3 <- sqrt(lambda3) * u3
norm_vec(loadingu1)
norm_vec(loadingu2)
norm_vec(loadingu3)


