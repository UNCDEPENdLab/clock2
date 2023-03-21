setwd("/Users/hallquist/Data_Analysis/clock2/code")
source("von_mises_basis.R")
##
#sig_spread=prop_spread*range(tvec); %determine SD of spread function

# should be same as basis elements
#sig_spread=sig;


##

# test basis set
vset <- vm_circle_set(n_basis=12, weights=1:12, width_sd=0.2)
plot(vset$get_pvec(), type="l")
plot(vset$get_vfunc(), type="l")
vset$get_weights()
vset$get_centers()

library(ggplot2)
library(dplyr)
df <- data.frame(pos=vset$get_pvec(), y=vset$get_vfunc())
df <- data.frame(pos=vset$get_pvec(), y=vset$get_vfunc())
ggplot(df, aes(x=pos, y=y)) + geom_line() + coord_polar()

x <- vset$get_basis()
vdf <- reshape2::melt(x, varnames=c("point", "basis"))
#bdf <- merge(vdf, xdf, by="point", all = TRUE)
ggplot(vdf, aes(x=point, y= value, color=factor(basis))) + geom_line()

dd <- vset$get_basis_df()
ggplot(dd, aes(x=loc, y=value, color=factor(basis))) + geom_line()

ee <- vset$eligibility
ee$center <- 0.0
plot(ee$get_vfunc())
sum(ee$get_vfunc())

plot(ee$basis_df$basis_norm)
vset$set_eligibility_center(0)
vset$eligibility$width_sd <- 0.2
plot(vset$get_eligibilities())
vset$get_eligibilities()

#em <- pracma::repmat(ee$get_vfunc(), 12, 1)
#em <- replicate(12, ee$get_vfunc())
em <- replicate(12, ee$basis_df$basis_norm)
b <- vset$get_basis()
ep <- em * b

# check overlap
plot(colSums(ep))

ep <- reshape2::melt(ep)
ggplot(ep, aes(x=Var1, y=weight, color=factor(Var2))) + geom_line()

plot(ee$get_vfunc())
plot(b[,7])

# verification that these are identical
cor(ee$get_vfunc(), b[,7])
summary(ee$get_vfunc() - b[,7])

# note that when we multiply two Gaussians together, the effective SD shrinks and even when the overlap is perfect,
# the maximum reaches a weight of .68.

# In the original SCEPTIC, we take the AUC of an untruncated Gaussian as a scaling factor:
refspread <- sum(gaussmf(1:100, sigma=10, mu=50))

# We then multiply this against any AUC=1.0 basis
auc1_example <- gaussmf(1:100, sigma=10, mu=30)
summary(auc1_example)
auc1_example <- auc1_example/sum(auc1_example)
summary(auc1_example)
auc1_example <- auc1_example*refspread # simple undoing of AUC scaling -- goes back to 1.0 max
summary(auc1_example)

# In the simple case (no truncation), this means that elig in h_sceptic_fixed_decay is a 1.0 max gaussmf centered on the RT
elig <- gaussmf(1:100, sigma=10, mu=50)

# This is the multiplied against a basis set where each element has AUC = 1.0
basis <- cbind(
  gaussmf(1:100, sigma=10, mu=40),
  gaussmf(1:100, sigma=10, mu=50),
  gaussmf(1:100, sigma=10, mu=60)
)

basis <- basis/colSums(basis)
apply(basis, 2, sum)
plot(basis[,2])

elig_mat <- replicate(3, elig)

elig_prod <- elig_mat * basis

df <- data.frame(
  t=1:100,
  p1=elig_prod[,1],
  p2=elig_prod[,2],
  p3=elig_prod[,3],
  e=elig/sum(elig)
)

# note that even with the perfect overlap of p2 and elig, the sum is ~.7.
# this is because the product of two Gaussians has a smaller SD even when the means are the same:
# https://rpsychologist.com/calculating-the-overlap-of-two-normal-distributions-using-monte-carlo-integration
colSums(df[,-1])

plot(elig_prod[,1], type="l")
lines(elig_prod[,2], type="l", col="blue")
lines(elig_prod[,3], type="l", col="orange")
lines(elig/sum(elig), type="l", col="green") #original elig -- rescaled to AUC1

# Critically, note that p1 -- the product of the eligbility function and the basis aligned to the same mean is more narrow
# than the eligibility function itself. This is a critical flaw in our products approach
dfm <- df %>% tidyr::gather(key="key", value = "value", -t)
ggplot(dfm, aes(x=t, y=value, color=key)) + geom_line()


sum((ee$basis_df$basis_norm*1.5) *b[,7])


sum((ee$basis_df$basis_norm*1.5) *b[,7])

sum(ee$basis_df$basis_norm *b[,7])

plot(ee$basis_df$basis_norm *b[,7], type="l")
lines(b[,7], type="l", col="blue")

plot(convolve(b[,7], ee$basis_df$basis_norm), type="l")

plot(b[,7], type="l")
plot(ee$basis_df$basis_norm, type="l")

# The correct approach is to normalize all functions to have AUC 1.0 and then to take the 
# pmin -- parallel minima -- to understand their overlap.
sum(b[,7])
sum(ee$basis_df$basis)


dd <- density(ee$basis_df$basis_norm, n = 360)
cbind(ee$basis_df$basis, dd$y)

sum(pmin(b[,7], ee$basis_df$basis))
sum(pmin(b[,8], ee$basis_df$basis))
sum(pmin(b[,6], ee$basis_df$basis))
sum(pmin(b[,5], ee$basis_df$basis))
sum(pmin(b[,1], ee$basis_df$basis))


cbind(pmin(b[,1], ee$basis_df$basis))
sum(apply(cbind(pmin(b[,7], ee$basis_df$basis)), 1, min))


sum(b[,7])
sum(ee$basis_df$basis_norm)


ee$basis_df$basis_norm
sum(ee$basis_df$basis_norm)



# curve 1
dvec <- dnorm(x=100:200, mean=150, sd=10)
dvec <- dvec/max(dvec)

plot(dvec, type="l")

d2 <- dnorm(x=100:200, mean=180, sd=10)
sum(pmin(dvec, d2))



lines(dvec*d2, type="l", col="blue")

int_f <- function(x, mu1, mu2, sd1, sd2) {
  f1 <- dnorm(x, mean=mu1, sd=sd1)
  f2 <- dnorm(x, mean=mu2, sd=sd2)
  pmin(f1, f2)
}
integrate(int_f, -Inf, Inf, mu1=0, mu2=0.0, sd1=1, sd2=1)

#e = sum(repmat(elig,nbasis,1).*inF.gaussmat_trunc, 2);
