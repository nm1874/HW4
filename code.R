#1.a 
library("pracma")
a1<-c(5,1,1)
a2<-c(-1,3,-1)
a3<-c(0,0,4)
A<-rbind(a1,a2,a3)
A
w <- c(1,0,0); w; 
A%*%w; A%*%A%*%w

#These three vectors in R^2 have to be linearly dependent.
#Confirm this by forming a matrix.
T <- cbind( w, A%*%w, A%*%A%*%w); T
#Row reduce to find the exact form of the dependence.
rref(T)

#So A^2w = -16Iw + 8Aw, or A^2w -8Aw + 16Iw = (A^2-8A+16I)w = 0.
#We have discovered the polynomial
w<-c(1,0,0)
p <- function(t) t^2 - 8*t + 16
#which in this case is the same as the characteristic polynomial.
#Since p(t)=(t-4)(t-4)
I <- diag(c(1,1,1))
(A - 4*I)%*%(A-4*I)%*%w    #has to be zero
#Now we can find two eigenvectors
v1 <- (A-4*I)%*%w     #eigenvector for eigenvalue 1
A%*%v1; 4*v1

rref(A-4*I)
#top row (x1, 1, 0) sets x1 = -1
#top row (x1, 0, 1) sets x1 = -1

#so we have two eigenvectors (two dim in the kernal of A-4*I)
v1<-c(-1,1,0)
v2<-c(1,-1,0)
eigen(A)$vectors

#1.b 
#A = D+N
d1<-c(4,0,0)
d2<-c(0,4,0)
d3<-c(0,0,4)
D<-rbind(d1,d2,d3); D
N<-A-D; N
N%*%N



#2
a_1<-c(4,-1,1)
a_2<-c(-1,3,2)
a_3<-c(1,2,-3)
A_<- rbind(a_1,a_2,a_3)
A_

identical(t(A_),A_)    #it is symmetric
w_ <- c(1,0,0)
T <- cbind(w_, A_%*%w_, A_%*%A_%*%w_, A_%*%A_%*%A_%*%w_); T
pCoef <- rref(T)[,4]; pCoef  #the last column
p <- function(t) t^3 - pCoef[3]*t^2 - pCoef[2]*t - pCoef[1] 
curve(p(x), from = -8, to =8); abline(h=0, col="red")

lam1 <- uniroot(p, c(4,5))$root;lam1; p(lam1)
lam2 <- uniroot(p, c(3,4))$root;lam2; p(lam2)
lam3 <- uniroot(p, c(-5,-1))$root;lam3; p(lam3)

I_ <- diag(c(1,1,1))
v_ <- (A_-lam2*I_ )%*%(A_-lam3*I_ )%*%w_; v1_<-v_/sqrt(sum(v_^2))
A_%*%v1_; lam1*v1_    #it checks!

v_  <- (A_-lam1*I_ )%*%(A_-lam3*I_ )%*%w_; v2_<-v_/sqrt(sum(v_^2))
v_  <- (A_-lam1*I_ )%*%(A_-lam2*I_ )%*%w_; v3_<-v_/sqrt(sum(v_^2))

P_ <- cbind(v1_,v2_,v3_); PInv_ <-solve(P_); D_ <- diag(c(lam1,lam2,lam3))
round(P_%*%D_%*%PInv_, digits = 4);A_  #we diagonalized A!!

round(t(P_) %*%P_, digits = 4)   #the columns are orthonormal
eigen(A_)$vectors; P_            #and we agree with the built-in function



#3

b1<-c(3,4,-4)
b2<-c(1,3,-1)
b3<-c(3,6,-4)
B<-rbind(b1,b2,b3); B

library("pracma")   #for rref()

#Any vector will do for w. An obvious choice is the first standard basis vector
w <- c(1,0,0)
T <- cbind(w, B%*%w, B%*%B%*%w, B%*%B%*%B%*%w); T
rref(T)
#So A^3 -2A^2-A+2 = 0
p <- function(t) t^3 - 2*t^2 -t +2
curve(p(x), from = -1, to = 4); abline(h=0, col = "red")

#Since P has three distinct real roots -1, 1, 2, we will get three eigenvectors.
I <- diag(c(1,1,1))
v0 <- (B - I)%*%(B - 2*I)%*%w; v0
B%*%v0; -1*v0      #eigenvector for eigenvalue -1
v1 <- (B + I)%*%(B - 2*I)%*%w; v1
B%*%v1; 1*v1      #eigenvector for eigenvalue 1
v2 <- (B + I)%*%(B - I)%*%w; v2
B%*%v2; 2*v2      #eigenvector for eigenvalue 2

c1<-c(1,0,0)
c2<-c(1,3,-1)
c3<-c(1,2,0)
C<-rbind(c1,c2,c3); C

w <- c(1,0,0)
T <- cbind(w, C%*%w, C%*%C%*%w); T
rref(T)
#So A^2 -3A + 2 = 0
p <- function(t) t^2 - 3*t +2
curve(p(x), from = -1, to = 4); abline(h=0, col = "red")

I <- diag(c(1,1,1))
w<-c(1,0,0)
v0 <- (C - 2*I)%*%w; v0
C%*%v0; 1*v0      #eigenvector for eigenvalue 1
v1 <- (C - I)%*%w; v1
C%*%v1; 2*v1      #eigenvector for eigenvalue 2
w<-c(0,1,0)
v2 <- (C - 2*I)%*%w; v2
C%*%v2; 1*v2 #second eigenvector for eigenvalue 1


