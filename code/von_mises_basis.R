# von mises basis function class
vm_bf <- R6::R6Class(
  "vm_bf",
  private = list(
    pvt_center = 0, # center of basis
    pvt_width_sd = 1, # SD of basis
    pvt_weight = 5, # height of basis
    
    # von mises membership function when normalize=TRUE or PDF when normalize=FALSE
    vm_mf = function(x=NULL, mu, kappa, normalize=TRUE) {
      y <- exp(kappa * cos(x - mu))/(2*pi*besselI(kappa, nu = 0))
      if (isTRUE(normalize)) {
        y <- y/max(y)
      }
      return(y)
    },
    setup_bvals = function() {
      x <- seq(-pi, pi, length.out=self$n_points)
      mu <- self$center
      sd <- self$width_sd
      
      # convert to radians for basis calculation
      if (self$units == "degrees") {
        mu <- (pi/180) * mu
        sd <- (pi/180) * sd
        pvec <- x/(pi/180)
      } else {
        pvec <- x
      }
      
      vm <- private$vm_mf(x=x, mu=mu, kappa=1/sd, normalize=FALSE)
      nvm <- vm/max(vm) # normalized basis (membership function)
      avm <- vm/sum(vm) # AUC 1.0 basis (compute against eligibility)
      #weights <- nvm*self$weight
      weights <- avm*self$weight #use AUC version
      df <- data.frame(pvec = pvec, pdf=vm, basis_norm=nvm, basis=avm, weights=weights)
      self$basis_df <- df
    }
  ),
  active = list(
    center = function(value) {
      if (missing(value)) {
        return(private$pvt_center)
      } else {
        checkmate::assert_number(value)
        private$pvt_center <- value
        private$setup_bvals() # re-define basis based on updated center
      }
    },
    width_sd = function(value) {
      if (missing(value)) {
        return(private$pvt_width_sd)
      } else {
        checkmate::assert_number(value)
        private$pvt_width_sd <- value
        private$setup_bvals() # re-define basis based on updated center
      }
    },
    weight = function(v) {
      if (missing(v)) {
        return(private$pvt_weight)
      } else {
        checkmate::assert_number(v)
        private$pvt_weight <- v
        private$setup_bvals() # re-define basis based on updated center
      }
    }
  ),
  public = list(
    units="radians",
    #probability = 0.5,
    weight_sd = 1,
    min_weight = 0,
    n_points=360, # resolution of basis
    basis_df = NULL, # data.frame of positions (x) and basis weights (y)
    initialize = function(units = NULL, center = NULL, width_sd = NULL, weight = NULL, weight_sd = NULL, min_weight = NULL) {
      if (!is.null(units)) self$units <- units
      if (!is.null(center)) private$pvt_center <- center
      if (!is.null(width_sd)) private$pvt_width_sd <- width_sd
      if (!is.null(weight)) private$pvt_weight <- weight
      if (!is.null(weight_sd)) self$weight_sd <- weight_sd
      if (!is.null(min_weight)) self$min_weight <- min_weight
      private$setup_bvals() # setup basis based on inputs
    },
    get_pdf = function() {
      self$basis_df$pdf
    },
    # get the time
    get_pvec = function() {
      self$basis_df$pvec
    },
    get_basis = function() {
      self$basis_df$basis
    },
    get_vfunc = function() {
      self$basis_df$weights
    }
  )
)

# prototype Gaussian radial basis
rbf <- R6::R6Class(
  "rbf",
  public = list(
    units="degrees",
    center = 0,
    width_sd = 10,
    #probability = 0.5,
    weight = 5,
    weight_sd = 1,
    min_weight = 0,
    initialize = function(units = NULL, center = NULL, width_sd = NULL, weight = NULL, weight_sd = NULL, min_weight = NULL) {
      if (!is.null(units)) self$units <- units
      if (!is.null(center)) self$center <- center
      if (!is.null(width_sd)) self$width_sd <- width_sd
      if (!is.null(weight)) self$weight <- weight
      if (!is.null(weight_sd)) self$weight_sd <- weight_sd
      if (!is.null(min_weight)) self$min_weight <- min_weight
    },
    get_vfunc = function() {
      dvec <- dnorm(x=1:360, mean=self$center, sd=self$width_sd)
      dvec <- dvec/max(dvec)*self$weight #renormalize to max=1
      return(dvec)
    }
  )
)


# testing
# a <- vm_bf$new(center=0, width_sd=0.1, weight=5)
# dd <- a$basis_df
# ggplot(dd, aes(x=pvec, y=weights)) + geom_line()
# ggplot(dd, aes(x=pvec, y=pdf)) + geom_line()
# ggplot(dd, aes(x=pvec, y=basis)) + geom_line()


vm_circle_set <- function(n_basis=12, weights=0, width_sd=0.3) {
  # repeat scalar weight for each basis function
  if (length(weights) == 1L) weights <- rep(weights, n_basis)
  if (length(width_sd) == 1L) width_sd <- rep(width_sd, n_basis)
  
  stopifnot(length(weights) == n_basis)
  stopifnot(length(width_sd) == n_basis)
  
  d_theta <- (2*pi)/n_basis
  
  # since circle wraps, we want to avoid adding a redundant basis function at -pi vs. at pi
  loc <- seq(-pi, pi - d_theta, by=d_theta)
  
  vset <- lapply(seq_along(loc), function(ii) {
    vm_bf$new(center=loc[ii], width_sd=width_sd[ii], weight=weights[ii])
  })
  
  vset <- rbf_set$new(elements=vset)
  
  # vdf <- reshape2::melt(vset, varnames=c("point", "basis"))
  # bdf <- merge(vdf, xdf, by="point", all = TRUE)
  # 
  # ggplot(bdf, aes(x=loc_rad, y=weight, color=factor(basis))) +
  #   geom_line() + coord_polar(theta = "x")
  # 
  # ggplot(bdf, aes(x=loc_rad, y=weight, color=factor(basis))) +
  #   geom_line()
}

# class that wraps a set of basis functions
rbf_set <- R6::R6Class(
  "rbf_set",
  public = list(
    elements = list(),
    eligibility = NULL,
    initialize = function(elements = NULL, eligibility=NULL) {
      if (!is.null(elements)) {
        checkmate::assert_list(elements)
        sapply(elements, checkmate::assert_multi_class, classes=c("rbf", "vm_bf"))
        # should probably validate that the pvec is identical for each element in the list...
        self$elements <- elements
      }
      
      if (!is.null(eligibility)) {
        checkmate::assert_multi_class(eligibility, classes=c("rbf", "vm_bf"))
        self$eligibility <- eligibility
      } else {
        # use the first element of the basis set as the eligibility function. Note that the center of the eligibility
        # is updated based on the input, so the location doesn't matter.
        self$eligibility <- elements[[1]]$clone(deep=TRUE)
      }
    },
    get_eligibility = function() {
      self$eligibility$get_basis()
    },
    # get position vector
    get_pvec = function() {
      self$elements[[1]]$get_pvec() # assume that the first position vector is the same for all elements
    },
    get_centers = function() {
      sapply(self$elements, function(x) x$center)
    },
    get_weights = function() {
      sapply(self$elements, function(x) x$weight)
    },
    get_vfunc = function() {
      rowSums(sapply(self$elements, function(x) x$get_vfunc()))
    },
    get_basis = function() {
      sapply(self$elements, function(x) x$get_basis())
    },
    get_basis_df = function() {
      a <- sapply(self$elements, function(x) x$get_basis())
      a_df <- reshape2::melt(a, varnames=c("point", "basis"))
      pvec <- self$elements[[1]]$get_pvec()
      p_df <- data.frame(point=seq_along(pvec), loc=pvec)
      b_df <- merge(a_df, p_df, by="point", all=TRUE) %>%
        dplyr::arrange(basis, point)
    },
    apply_drift = function(step_size) {
      v <- rnorm(1, 0, sd=step_size)
      c <- self$get_centers()
      cnew <- c + v
      cnew <- sapply(cnew, function(x) {
        if (x > 360) x = x - 360
        else if (x < 0) x = 360 + x
        else x
      })
      # update centers
      sapply(seq_along(cnew), function(ii) self$elements[[ii]]$center <- cnew[ii])
      return(self)
    }
  )
)


##
#sig_spread=prop_spread*range(tvec); %determine SD of spread function

# should be same as basis elements
#sig_spread=sig;


##

# test basis set
vset <- vm_circle_set(n_basis=12, weights=1:12, width_sd=0.3)
plot(vset$get_pvec(), type="l")
plot(vset$get_vfunc(), type="l")
vset$get_weights()
vset$get_centers()

df <- data.frame(pos=vset$get_pvec(), y=vset$get_vfunc())
df <- data.frame(pos=vset$get_pvec(), y=vset$get_vfunc())
ggplot(df, aes(x=pos, y=y)) + geom_line() + coord_polar()

x <- vset$get_basis()
vdf <- reshape2::melt(x, varnames=c("point", "basis"))
bdf <- merge(vdf, xdf, by="point", all = TRUE)
ggplot(vdf, aes(x=point, y= weight, color=factor(basis))) + geom_line()

dd <- vset$get_basis_df()
ggplot(dd, aes(x=loc, y=weight, color=factor(basis))) + geom_line()

ee <- vset$eligibility
ee$center <- 0.0
plot(ee$get_vfunc())
sum(ee$get_vfunc())

plot(ee$basis_df$basis_norm)


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
