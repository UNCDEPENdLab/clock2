# von mises basis function class
vm_bf <- R6::R6Class(
  "vm_bf",
  private = list(
    pvt_center = 0, # center of basis
    pvt_width_sd = 1, # SD of basis
    pvt_weight = 5, # height of basis
    
    # von Mises: membership function when normalize=TRUE or PDF when normalize=FALSE
    vm_mf = function(x=NULL, mu, kappa, normalize=TRUE) {
      y <- exp(kappa * cos(x - mu))/(2*pi*besselI(kappa, nu = 0))
      if (isTRUE(normalize)) {
        y <- y/max(y)
      }
      return(y)
    },
    
    # setup data.frame containing basis function values across the circular space
    setup_bvals = function() {
      x <- seq(0, 2*pi, length.out=self$n_points)
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
      
      vm <- private$vm_mf(x=x, mu=mu, kappa=1/sd, normalize=FALSE) # probability density
      nvm <- vm/max(vm) # normalized max=1.0 basis (membership function)
      avm <- vm/sum(vm) # AUC 1.0 basis (used for update eligibility)
      weights <- nvm*self$weight
      #weights <- avm*self$weight #use AUC version
      df <- data.frame(pvec = pvec, pdf=vm, basis_norm=nvm, basis=avm, weights=weights)
      self$basis_df <- df
    }
  ),
  active = list(
    #' @field center The center point of the basis function
    center = function(value) {
      if (missing(value)) {
        return(private$pvt_center)
      } else {
        checkmate::assert_number(value)
        private$pvt_center <- value
        private$setup_bvals() # re-define basis based on updated center
      }
    },
    #' @field width_sd The width of the basis function specified as its standard deviation
    width_sd = function(value) {
      if (missing(value)) {
        return(private$pvt_width_sd)
      } else {
        checkmate::assert_number(value)
        private$pvt_width_sd <- value
        private$setup_bvals() # re-define basis based on updated center
      }
    },
    #' @field weight The weight (aka height) of the basis function
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
    initialize = function(n_points = NULL, units = NULL, center = NULL, width_sd = NULL, weight = NULL, weight_sd = NULL, min_weight = NULL) {
      if (!is.null(n_points)) self$n_points <- n_points
      if (!is.null(units)) self$units <- units
      if (!is.null(center)) private$pvt_center <- center
      if (!is.null(width_sd)) private$pvt_width_sd <- width_sd
      if (!is.null(weight)) private$pvt_weight <- weight
      if (!is.null(weight_sd)) self$weight_sd <- weight_sd
      if (!is.null(min_weight)) self$min_weight <- min_weight
      private$setup_bvals() # setup basis based on inputs
    },
    #' @description get probability density function for basis (sum=1.0)
    get_pdf = function() {
      self$basis_df$pdf
    },
    #' @description get the position vector for the basis function
    #' @return a vector of the positions tracked by the basis
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



# testing
# a <- vm_bf$new(center=0, width_sd=0.1, weight=5)
# dd <- a$basis_df
# ggplot(dd, aes(x=pvec, y=weights)) + geom_line()
# ggplot(dd, aes(x=pvec, y=pdf)) + geom_line()
# ggplot(dd, aes(x=pvec, y=basis)) + geom_line()

vm_circle_contingency <- function(centers, weights, widths, units="degrees") {
  gg <- lapply(seq_along(centers), function(ii) {
    r <- vm_bf$new(weight=weights[ii], center=centers[ii], width_sd = widths[ii], units=units)
  })
  contingency <- rbf_set$new(elements = gg)
  return(contingency)
}

vm_circle_set <- function(n_basis=12, weights=0, width_sd=0.3) {
  # repeat scalar weight for each basis function
  if (length(weights) == 1L) weights <- rep(weights, n_basis)
  if (length(width_sd) == 1L) width_sd <- rep(width_sd, n_basis)
  
  stopifnot(length(weights) == n_basis)
  stopifnot(length(width_sd) == n_basis)
  
  d_theta <- (2*pi)/n_basis
  
  # since circle wraps, we want to avoid adding a redundant basis function at 0 vs. at 2*pi
  loc <- seq(0, 2*pi - d_theta, by=d_theta)
  
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
  private = list(
    units="radians",
    
    # function that computes the proportion overlap between 
    compute_overlap = function(b1, b2) {
      e1 <- b1$get_basis()
      e2 <- b2$get_basis()
      stopifnot(abs(sum(e1) - 1) < 1e-5) # verify AUC = 1
      stopifnot(abs(sum(e2) - 1) < 1e-5)
      sum(pmin(e1, e2))
    }
  ),
  public = list(
    elements = list(),
    eligibility = NULL,
    initialize = function(elements = NULL, eligibility=NULL) {
      if (!is.null(elements)) {
        checkmate::assert_list(elements)
        sapply(elements, checkmate::assert_multi_class, classes=c("rbf", "vm_bf"))
        unit_vec <- sapply(elements, "[[", "units")
        if (length(unique(unit_vec)) > 1L) {
          print(unit_vec)
          stop("Units for elements do not match")
        } else {
          private$units <- elements[[1]]$units
        }
        
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
    set_eligibility_center = function(value) {
      checkmate::assert_number(value)
      self$eligibility$center <- value
    },
    get_eligibilities = function(efunc=NULL) {
      # allow for input of eligibility from outside
      if (is.null(efunc)) efunc <- self$eligibility
      checkmate::assert_multi_class(efunc, classes=c("rbf", "vm_bf"))
      # compute eligibility of each basis function against the eligibility
      sapply(self$elements, function(e) {
        private$compute_overlap(e, efunc)
      })
    },
    # get position vector
    get_pvec = function() {
      self$elements[[1]]$get_pvec() # assume that the first position vector is the same for all elements
    },
    get_centers = function() {
      sapply(self$elements, function(x) x$center)
    },
    set_weights = function(v) {
      checkmate::assert_numeric(v, len=length(self$elements))
      sapply(seq_along(self$elements), function(ii) self$elements[[ii]]$weight <- v[ii])
      return(self)
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
        if (private$units == "degrees") {
          if (x > 360) x <- x - 360
          else if (x < 0) x <- 360 + x
        } else {
          x <- x %% (2*pi) # wrap any radians onto [0, 2*pi] interval
          if (x > 2*pi) x <- x - 2*pi
          else if (x < 0) x <- 2*pi + x
        }
        return(x)
      })
      
      # update centers
      sapply(seq_along(cnew), function(ii) self$elements[[ii]]$center <- cnew[ii])
      return(self)
    }
  )
)





### 
# prototype Gaussian radial basis (not currently used)
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


scepticc <- R6::R6Class(
  "scepticc",
  private = list(
    pvt_bf_set = NULL,
    pvt_eligibility = NULL,
    pvt_alpha = 0.1,
    pvt_gamma = 0.3,
    pvt_beta = 2,
    pvt_n_points = 50,
    pvt_history = NULL,
    pvt_contingency = NULL
  ),
  active = list(
    alpha = function(v) {
      if (missing(v)) return(private$pvt_alpha)
      else {
        checkmate::assert_number(v, lower=0.001, upper=0.999)
        private$pvt_alpha <- v
      }
    },
    gamma = function(v) {
      if (missing(v)) return(private$pvt_gamma)
      else {
        checkmate::assert_number(v, lower=0.001, upper=0.999)
        private$pvt_gamma <- v
      }
    },
    #' @field beta temperature controlling softmax (higher values -> more stochastic)
    beta = function(v) {
      if (missing(v)) return(private$pvt_beta)
      else {
        checkmate::assert_number(v, lower=0.001, upper=1e4)
        private$pvt_beta <- v
      }
    }
  ),
  public = list(
    initialize = function(n_basis = 12, n_points = NULL, basis_sd = 0.3, weights_0 = 0, elig_sd = 0.3, 
                          alpha=NULL, gamma=NULL, beta=NULL, contingency=NULL) {
      checkmate::assert_number(n_basis, lower=2)
      if (checkmate::test_number(weights_0)) weights_0 <- rep(weights_0, n_basis)
      if (checkmate::test_number(basis_sd)) basis_sd <- rep(basis_sd, n_basis)
      stopifnot(length(weights_0) == n_basis)
      stopifnot(length(basis_sd) == n_basis)
      
      if (!is.null(n_points)) {
        checkmate::assert_number(n_points, lower=2)
        private$pvt_n_points <- n_points
      }
      
      # setup positions
      d_theta <- (2*pi)/n_basis
      
      # since circle wraps, we want to avoid adding a redundant basis function at 0 vs. at 2*pi
      loc <- seq(0, 2*pi - d_theta, by=d_theta)
      vset <- lapply(seq_along(loc), function(ii) {
        vm_bf$new(n_points = private$pvt_n_points, center=loc[ii], width_sd=basis_sd[ii], weight=weights_0[ii])
      })
      
      private$pvt_bf_set <- rbf_set$new(elements=vset)
      private$pvt_eligibility <- vm_bf$new(n_points = private$pvt_n_points, center=0, width_sd=elig_sd, weight=1)
      
      if (!is.null(alpha)) self$alpha <- alpha
      if (!is.null(beta)) self$beta <- beta
      if (!is.null(gamma)) self$gamma <- gamma
      
      if (!is.null(contingency)) {
        checkmate::assert_class(contingency, classes = "troll_world")
        private$pvt_contingency <- contingency # assign by reference -- will update original object
      }
    },
    get_weights = function() {
      private$pvt_bf_set$get_weights()
    },
    update_weights = function(tau, sd=NULL, outcome, model="decay") {
      checkmate::assert_number(tau)
      private$pvt_eligibility$center <- tau
      e <- private$pvt_bf_set$get_eligibilities(private$pvt_eligibility)
      w <- private$pvt_bf_set$get_weights()
      pe <-  outcome - w
      if (model == "decay") {
        decay <- -self$gamma * (1-e) * w  
      } else {
        decay <- 0
      }
      
      w_new <- w + self$alpha*e*pe + decay
      private$pvt_bf_set$set_weights(w_new)
      return(self)
    },
    get_choice_probs = function() {
      #v=x_t(1:nbasis)*ones(1,ntimesteps) .* gaussmat; %use vector outer product to replicate weight vector
      #v_func = sum(v); %subjective value by timestep as a sum of all basis functions
      
      v_func <- private$pvt_bf_set$get_vfunc()
      
      p_choice <- (exp((v_func-max(v_func))/private$pvt_beta)) / (sum(exp((v_func-max(v_func))/private$pvt_beta))) # Divide by temperature
      return(p_choice)
    },
    emit_choice = function() {
      # sample choice based on softmax probabilities
      sample(1:private$pvt_n_points, 1, prob = self$get_choice_probs())
    },
    
  )
)