nsim = 10e3

a = c(5, 2, 3)

xd = t(replicate(nsim, rDirichlet(a)))
dim(xd)

(K = length(a))
gamdel = shape_Dir2genDir(a)
## both should produce Dirichlet draws
xsbm = t(replicate(nsim,
                   rSBM(K=K, p_small=0.0, p_large=0.0,
                        eta=1.0e3, gam=gamdel[,1], del=gamdel[,2])))
xsbm2 = t(replicate(nsim,
                    rPostSBM(x=rep(0,K), p_small=0.0, p_large=0.0,
                        eta=1.0e3, gam=gamdel[,1], del=gamdel[,2])$w))

a / sum(a)
colMeans(xd)
colMeans(xsbm)
colMeans(xsbm2)



### Explore prior attributes
K = 7
a = rep(1, K)
(gamdel = shape_Dir2genDir(a))
(gamdel = matrix(c(1.5, 1.5), ncol=2))

p_small = 0.5
p_large = 0.0

(delfact = delCorrection_SBM(p_small, p_large, K))
(delfact = K*p_small/(K-1)) # if p_large = 0

x1 = t(replicate(
  nsim,
  rSBM(
    K = K,
    p_small = p_small,
    p_large = p_large,
    eta = 1.0e3,
    gam = gamdel[, 1],
    del = gamdel[, 2]
  )
))

x2 = t(replicate(
  nsim,
  rSBM(
    K = K,
    p_small = p_small,
    p_large = p_large,
    eta = 1.0e3,
    gam = gamdel[, 1],
    del = delfact*gamdel[, 2]
  )
))

round(head(x1), 3)
round(head(x2), 3)

colMeans(x1)
colMeans(x2)

msg = "Dirichlet"
msg = "Fixed delta, gamma\np_small=0.5, p_large=0.2"
msg = "Tapered, adjusted delta\np_small=0.5, p_large=0"

boxplot(x1, ylim=c(0.0, 1.0), main=msg)
boxplot(x2, ylim=c(0.0, 1.0), main=msg)

sum(diag(cov(x1)))
sum(diag(cov(x2)))



### explore posterior attributes
rm(list=ls())

nsim = 10e3

K = 9
(a = rep(1.0/K, K))
(gamdel = shape_Dir2genDir(a))

eta = 5e2
p_small = 0.7
p_large = 0.1

(delfact = delCorrection_SBM(p_small, p_large, K))
# delfact = 1.0

x = c(1, 3, 0, 5, 1, 0, 0, 4, 0)

thD = t(replicate(nsim, rDirichlet(a+x)))

th = t(replicate(
  nsim,
  rPostSBM(
    x = x,
    p_small = p_small,
    p_large = p_large,
    eta = eta,
    gam = gamdel[, 1],
    del = gamdel[, 2] * delfact
  )$w
))
str(th)
head(th)

library("ggtern")
th1id = 2
th2id = 4
th3id = 8
th_df = data.frame(theta1=th[,th1id], theta2=th[,th2id], theta3=th[,th3id])
thD_df = data.frame(theta1=thD[,th1id], theta2=thD[,th2id], theta3=thD[,th3id])

ggtern(data=th_df, aes(theta1, theta2, theta3)) +
  stat_density_tern(geom="polygon", aes(fill=..level..), bins=50)

ggtern(data=thD_df, aes(theta1, theta2, theta3)) +
  stat_density_tern(geom="polygon", aes(fill=..level..), bins=50)

k = 9
hist(th[,k], xlim=c(0,1))
hist(thD[,k], xlim=c(0,1))

dev.off()

