scepticc <- R6::R6Class(
  "scepticc",
  private = list(
    pvt_bf_set = NULL,
    pvt_eligibility = NULL,
    pvt_alpha = 0.1,
    pvt_gamma = 0.3,
    pvt_beta = 2,
    pvt_n_points = 50,
    pvt_model = "decay", # what variant of sceptic to run
    pvt_history = NULL,
    pvt_contingency = NULL,
    pvt_seed = 1001, # random seed set for running contingency
    pvt_sample = 1
  ),
  active = list(
    alpha = function(v) {
      if (missing(v)) {
        return(private$pvt_alpha)
      } else {
        checkmate::assert_number(v, lower=0.001, upper=0.999)
        private$pvt_alpha <- v
      }
    },
    gamma = function(v) {
      if (missing(v)) {
        return(private$pvt_gamma)
      } else {
        checkmate::assert_number(v, lower=0.001, upper=0.999)
        private$pvt_gamma <- v
      }
    },
    #' @field beta temperature controlling softmax (higher values -> more stochastic)
    beta = function(v) {
      if (missing(v)) {
        return(private$pvt_beta)
      } else {
        checkmate::assert_number(v, lower=0.001, upper=1e4)
        private$pvt_beta <- v
      }
    }
  ),
  public = list(
    initialize = function(n_basis = 12, n_points = NULL, basis_sd = 0.3, weights_0 = 0, elig_sd = 0.3, 
                          alpha=NULL, gamma=NULL, beta=NULL, contingency=NULL, data = NULL, seed=NULL) {
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
      
      checkmate::assert_false(!is.null(contingency) & !is.null(data))
      if (!is.null(contingency)) {
        checkmate::assert_class(contingency, classes = "troll_world")
        private$pvt_contingency <- contingency # assign by reference -- will update original object
      }
      
      if (!is.null(data)) {
        checkmate::assert_class(data, classes = "data.frame")
        private$pvt_data <- data # assign by reference -- will update original object
      }
      

      if (!is.null(seed)) {
        checkmate::assert_integerish(seed, lower = 1, len = 1L)
        private$pvt_seed <- seed
      }
    },
    get_weights = function() {
      private$pvt_bf_set$get_weights()
    },
    update_weights = function(tau, sd=NULL, outcome, model="decay") {
      #cat("alpha: ", self$alpha, ", gamma: ", self$gamma, ", beta: ", self$beta, "\n")
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
      private$pvt_history[[private$pvt_sample]] <- w_new
      private$pvt_sample <- private$pvt_sample + 1
      return(self)
    },
    get_weight_history = function() {
      do.call(rbind, private$pvt_history)
    },
    get_entropy_history = function() {
      w <- do.call(rbind, private$pvt_history)
      wnorm <- w/rowSums(w)
      h <- NA[seq_len(nrow(wnorm))]
      for (t in seq_len(nrow(wnorm))) {
        h[t] <- -sum(wnorm[t,]*log2(wnorm[t,]))
      }
      return(h)
    },
    
    get_func_history = function() {
      tmp_bf <- private$pvt_bf_set$clone(deep=TRUE)
      h <- self$get_weight_history()
      vhist <- t(sapply(seq_len(nrow(h)), function(i) {
        tmp_bf$set_weights(h[i,])
        tmp_bf$get_wfunc()
      }))
      return(vhist)
    },
    reset_history = function() {
      private$pvt_history <- list()
      private$pvt_sample <- 1
    },
    get_choice_probs = function(model="sceptic") {
      #v=x_t(1:nbasis)*ones(1,ntimesteps) .* gaussmat; %use vector outer product to replicate weight vector
      #v_func = sum(v); %subjective value by timestep as a sum of all basis functions

      if (model == "sceptic") {
        v_func <- private$pvt_bf_set$get_wfunc()
        p_choice <- (exp((v_func - max(v_func)) / private$pvt_beta)) / (sum(exp((v_func - max(v_func)) / private$pvt_beta))) # Divide by temperature
      } else if (model == "random") {
        n_t <- length(private$pvt_bf_set$get_pvec())
        p_choice <- rep(1 /n_t, n_t)
      }

      return(p_choice)
    },
    emit_choice = function() {
      # sample choice based on softmax probabilities
      # print(private$pvt_bf_set$get_wfunc())
      # lookup choice against positions in radians
      s <- sample(seq_len(private$pvt_n_points), 1, prob = self$get_choice_probs())
      return(private$pvt_eligibility$get_pvec()[s])
    },
    run_contingency = function(pvec=NULL, optimize=FALSE) {
      if (is.null(private$pvt_contingency)) {
        stop("Cannot run contingency if it does not exist")
      }
      
      if (!is.null(pvec)) {
        private$pvt_alpha <- pvec["alpha"]
        private$pvt_beta <- pvec["beta"]
        private$pvt_gamma <- pvec["gamma"]
      }

      n_trials <- private$pvt_contingency$get_n_trials()
      
      private$pvt_contingency$reset_counter()
      private$pvt_bf_set$reset_weights() # set back to 0
      self$reset_history()
      set.seed(private$pvt_seed)
      
      if (!optimize) {
        history <- data.frame(trial=1:n_trials, choice=rep(NA, n_trials), outcome=rep(NA, n_trials))
        for (i in seq_len(n_trials)) {
          history$choice[i] <- choice <- self$emit_choice()
          history$outcome[i] <- private$pvt_contingency$harvest(choice)
          self$update_weights(choice, outcome=history$outcome[i], model="decay")
        }
        return(history)
      } else {
        ovec <- rep(NA, n_trials)
        for (i in seq_len(n_trials)) {
          choice <- self$emit_choice()
          ovec[i] <- private$pvt_contingency$harvest(choice)
          self$update_weights(choice, outcome=ovec[i], model="decay")
          #cat("c: ", choice, "\n")
        }
        #print(pvec)
        #cat("o: ", sum(ovec), "\n")
        return(-sum(ovec))
      }
    },
    optimize_returns = function() {
      pvec <- c(alpha = self$alpha, gamma = self$gamma, beta = self$beta)
      lb <- c(0.001, 0, 0.05) %>% setNames(c("alpha", "gamma", "beta"))
      ub <- c(0.999, 0.95, 20) %>% setNames(c("alpha", "gamma", "beta"))
      parscale <- c(.1, .1, 1) %>% setNames(c("alpha", "gamma", "beta"))#rough step sizes for parameters
      browser()
      # optResult <- nlminb(start=pvec, objective=self$run_contingency, optimize=TRUE,
      #                     lower=lb, upper=ub, control=list(eval.max=500, iter.max=500, step.min=0.1))
      
      
      avals <- seq(0.01, 0.95, by=0.001)
      rets <- sapply(avals, function(a) {
        self$alpha <- a
        self$run_contingency(optimize=TRUE)
      })
      
      optResult <- optim(par=pvec, fn=self$run_contingency, optimize=TRUE, method="L-BFGS-B", lower=lb, upper=ub,
                         control=list(parscale=c(.1, .1, 1), maxit=1e4))
      
      
      #scale=1/parscale, 
    },
    # run_data
    run_data = function() {
      if (is.null(private$pvt_data)) {
        stop("Need data to run")
      }
      n_trials <- length(unique(private$pvt_data$trial))
      history <- private$pvt_data %>% mutate(
        choice <- pos_shifted*pi/180,
        u_location <- u_location*pi/180,
      )
      # private$pvt_contingency$reset_counter()
      self$reset_history()
      # set.seed(private$pvt_seed)
      for (i in seq_len(n_trials)) {
        history$choice[i] <- choice <- private$pvt_data$pos_shifted[i]*pi/180
        history$outcome[i] <- private$pvt_data$Earnings[i]
        self$update_weights(choice, outcome=history$outcome[i], model="decay")
      }
      return(history)
    },
    get_contingency = function() {
      private$pvt_contingency
    }
    

  )
)