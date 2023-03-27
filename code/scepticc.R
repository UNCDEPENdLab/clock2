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
                          alpha=NULL, gamma=NULL, beta=NULL, contingency=NULL, seed=NULL) {
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

      if (!is.null(seed)) {
        checkmate::assert_integerish(seed, lower = 1, len = 1L)
        private$pvt_seed <- seed
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
      private$pvt_history[[private$pvt_sample]] <- w_new
      private$pvt_sample <- private$pvt_sample + 1
      return(self)
    },
    get_weight_history = function() {
      do.call(rbind, private$pvt_history)
    },
    get_func_history = function() {
      tmp_bf <- private$pvt_bf_set$clone(deep=TRUE)
      h <- self$get_weight_history()
      vhist <- t(sapply(seq_len(nrow(h)), function(i) {
        tmp_bf$set_weights(h[i,])
        tmp_bf$get_vfunc()
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
        v_func <- private$pvt_bf_set$get_vfunc()
        p_choice <- (exp((v_func - max(v_func)) / private$pvt_beta)) / (sum(exp((v_func - max(v_func)) / private$pvt_beta))) # Divide by temperature
      } else if (model == "random") {
        n_t <- length(private$pvt_bf_set$get_pvec())
        p_choice <- rep(1 /n_t, n_t)
      }

      return(p_choice)
    },
    emit_choice = function() {
      # sample choice based on softmax probabilities
      sample(1:private$pvt_n_points, 1, prob = self$get_choice_probs())
    },
    run_contingency = function() {
      if (is.null(private$pvt_contingency)) {
        stop("Cannot run contingency if it does not exist")
      }

      n_trials <- private$pvt_contingency$get_n_trials()
      history <- data.frame(trial=1:n_trials, choice=rep(NA, n_trials), outcome=rep(NA, n_trials))
      private$pvt_contingency$reset_counter()
      self$reset_history()
      set.seed(private$pvt_seed)
      for (i in seq_len(n_trials)) {
        history$choice[i] <- choice <- self$emit_choice()
        history$outcome[i] <- private$pvt_contingency$harvest(choice)
        self$update_weights(choice, outcome=history$outcome[i], model="decay")
      }
      return(history)
    }

  )
)