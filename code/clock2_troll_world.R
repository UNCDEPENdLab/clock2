troll_world <- R6::R6Class(
  "troll_world",
  private = list(
    pvt_pvec = NULL, # position of values around circle in units of measurement (radians)
    shift_vec = function(x, delta) {
      ## converts x[n] into x[n-1] by multiplying the Fourier
      ## transform on x[n] by z^-1 i.e. by e^-iw
      ## lifted from emuR package
      N <- length(x)
      delta <- delta %% N # wrap around if shift is larger than length
      if (delta < 0) delta <- N + delta
      h <- c(rep(0, delta), 1)
      stats::filter(x, h, sides = 1, circular = TRUE)[1:N]
    },
    # internal function to calculate drifting and flexed values
    calculate_values = function() {
      if (isTRUE(private$pvt_clean)) return(invisible(NULL))
      
      has_drift <- abs(private$pvt_drift_sd) > 1e-5
      has_flex <- !is.null(self$spread)
      has_jump <- !is.null(private$pvt_jump_vec)

      if (has_drift) {
        # use sweep for efficient addition of a vector to each row
        dv <- private$pvt_drift_vec
        private$pvt_drift_tvals <- t(sapply(seq_len(private$pvt_n_trials), function(i) {
          private$shift_vec(self$tvals[i, ], dv[i])
        }))
      } else {
        # drifting values are identical to regular
        private$pvt_drift_tvals <- self$tvals
        dv <- rep(0, private$pvt_n_trials)
      }

      if (has_flex) {
        private$pvt_flex_tvals <- t(sapply(seq_len(private$pvt_n_trials), function(i) {
          tval <- self$v_rescale(self$tvals[i, ], spread = self$spread[i])
        }))
      } else {
        # no flexing -- original values apply
        private$pvt_flex_tvals <- self$tvals
      }

      if (has_jump || has_drift) {
        # add slow drift to jump
        sv <- private$pvt_jump_vec + dv
        private$pvt_objective_tvals <- t(sapply(seq_len(private$pvt_n_trials), function(i) {
          private$shift_vec(private$pvt_flex_tvals[i, ], sv[i])
        }))
      } else {
        private$pvt_objective_tvals <- private$pvt_flex_tvals
      }

      # cached value matrices are now current
      private$pvt_clean <- TRUE
    },
    pvt_n_trials = 100,
    pvt_n_timesteps = 200,
    pvt_drift_sd = 0,
    pvt_drift_vec = NULL,
    pvt_drift_tvals = NULL, # drift only
    pvt_flex_tvals = NULL, # flex only
    pvt_objective_tvals = NULL, # objective values with all manipulations
    pvt_jump_vec = NULL, # shift vector for jump of location during high entropy (flex)
    pvt_clean = FALSE
  ),
  active = list(
    drift_sd = function(v) {
      if (missing(v)) {
        return(private$pvt_drift_sd)
      } else {
        checkmate::assert_integerish(v, lower = 0, len = 1L)
        private$pvt_drift_sd <- v

        # setup drift vector
        private$pvt_drift_vec <- cumsum(round(rnorm(private$pvt_n_trials, mean = 0, sd = v)))
      }
    }
  ),
  public = list(
    erased_segments = list(),
    cur_trial = 1,
    units="radians",
    tvals=matrix(0, nrow=100, ncol=200),
    spread = NULL, # tracks the point spread between min and max (entropy flexing)
    epoch = NULL, # tracks what phase we are in of the entropy manipulation
    p_reward = 0.7, # current fixed probability for reward

    #' @param values A vector of values for every timestep, or a trials x timesteps matrix
    initialize = function(n_trials = NULL, values = NULL, drift_sd = NULL) {
      if (!is.null(n_trials)) {
        checkmate::assert_integerish(n_trials, lower=1, len=1L)
        private$pvt_n_trials <- n_trials
      }

      if (!is.null(values)) {
        if (is.vector(values)) {
          self$tvals <- pracma::repmat(values, n = private$pvt_n_trials, m = 1)
        } else if (is.matrix(values)) {
          # allow trials x timesteps matrix of values as input
          stopifnot(nrow(values) == private$pvt_n_trials)
          self$tvals <- values
        }
      }

      private$pvt_n_timesteps <- ncol(self$tvals)
      
      if (self$units == "radians") {
        max_val <- 2*pi
      } else if (self$units == "degrees") {
        max_val <- 360
      } else { stop("wrong units") }
      
      # setup positions
      d_theta <- max_val/private$pvt_n_timesteps
      
      # since circle wraps, we want to avoid adding a redundant basis function at 0 vs. at 2*pi
      private$pvt_pvec <- seq(0, max_val - d_theta, by=d_theta)

      if (!is.null(drift_sd)) {
        checkmate::assert_integerish(drift_sd, lower=0, len=1L)
        self$drift_sd <- drift_sd
      }
    },

    erase_segment = function(size_deg = 30, trial = 1) {
      pos_start <- sample(1:360, size=1)
      pos_end <- pos_start + size_deg - 1
      pos_vec <- pos_start:pos_end

      pos_vec <- sapply(pos_vec, function(x) {
        if (x > 360) x = x - 360
        else x
      })

      browser()

      # tmp fill
      # to fill rows of a matrix with the same vector, we need to recall that matrices are filled column-wise, not row-wise.
      # Thus, we can either double transpose the matrix for the operation, or use rep() to replicate the vector for columnwise assignment
      # https://stackoverflow.com/questions/11424047/assigning-values-to-rows-from-a-vector
      self$tvals[trial:private$pvt_n_trials, pos_vec] <- runif(size_deg) #pracma::repmat(runif(size_deg), n=length(trial:private$pvt_n_trials), m=1)
      return(self)
    },
    get_n_trials = function() {
      private$pvt_n_trials
    },
    get_n_timesteps = function() {
      private$pvt_n_timesteps
    },
    get_starting_values = function() {
      self$tvals[1, ] #starting point
    },
    get_values_matrix = function(type="objective", quiet=TRUE) {
      private$calculate_values() # always ensure fresh values

      if (type=="original") {
        if (!quiet) message("Returning original values matrix")
        self$tvals
      } else if (type == "drift") {
        if (!quiet) message("Returning values matrix with GRW drift only")
        private$pvt_drift_tvals
      } else if (type == "flex") {
        if (!quiet) message("Returning values matrix with range flex only")
        private$pvt_flex_tvals
      } else {
        if (!quiet) message("Returning objective values matrix with all value manipulations")
        private$pvt_objective_tvals
      }
    },
    reset_counter = function() {
      self$cur_trial <- 1
      return(self)
    },
    apply_flex = function(low_avg=20, low_spread=5, decrease_avg=10, decrease_spread=5, high_avg=5, high_spread=2, 
                          increase_avg=10, increase_spread=5, spread_max=80, spread_min=10, jump_high=TRUE, start_low=TRUE) {

      t_counter <- 1 # how many trials have we generated
      spread <- rep(NA, private$pvt_n_trials)
      epoch <- rep(NA_character_, private$pvt_n_trials)
      while (t_counter < private$pvt_n_trials) {
        # for now, always start with low entropy
        # low entropy period
        low_duration <- low_avg + sample(-low_spread:low_spread, 1)
        spread[t_counter:(t_counter + low_duration - 1)] <- spread_max
        epoch[t_counter:(t_counter + low_duration - 1)] <- "low"
        t_counter <- t_counter + low_duration

        # increasing entropy
        increase_duration <- increase_avg + sample(-increase_spread:increase_spread, 1)
        # need to trim first and last elements in approximation so that we only keep the changing values
        ivals <- approx(x=c(1, increase_duration + 2), y=c(spread_max, spread_min), xout=2:(increase_duration+1))$y
        spread[t_counter:(t_counter + increase_duration - 1)] <- ivals
        epoch[t_counter:(t_counter + increase_duration - 1)] <- "increasing"
        t_counter <- t_counter + increase_duration

        # high entropy
        high_duration <- high_avg + sample(-high_spread:high_spread, 1)
        spread[t_counter:(t_counter + high_duration - 1)] <- spread_min
        epoch[t_counter:(t_counter + high_duration - 1)] <- "high"
        t_counter <- t_counter + high_duration

        # decreasing entropy
        decrease_duration <- decrease_avg + sample(-decrease_spread:decrease_spread, 1)
        dvals <- approx(x=c(1, decrease_duration + 2), y=c(spread_min, spread_max), xout=2:(decrease_duration+1))$y
        spread[t_counter:(t_counter + decrease_duration - 1)] <- dvals
        epoch[t_counter:(t_counter + decrease_duration - 1)] <- "decreasing"
        t_counter <- t_counter + decrease_duration
      }

      # chop the spread back to the number of trials
      self$spread <- spread[1:private$pvt_n_trials]
      self$epoch <- epoch[1:private$pvt_n_trials]

      # handle the jump of location during high entropy periods
      if (isTRUE(jump_high)) {
        high_pos <- which(dplyr::lag(self$epoch) != self$epoch & self$epoch == "high")
        n_jumps <- length(high_pos)
        jv <- rep(0, private$pvt_n_trials)
        jv[high_pos] <- (-1)^sample(1:2, n_jumps, replace=TRUE) * sample(seq(90, 180, by = 10), size = n_jumps, replace=TRUE)
        private$pvt_jump_vec <- cumsum(jv) # use cumsum to allow vectorized shifts of value matrix
      }

      # need to regenerate values when requested
      private$pvt_clean <- FALSE
    },
    #' @description rescale a vector to have a given 2*spread range around the mean and a designated mean
    v_rescale = function(x, spread = 20, mean_val=50) {
      #mx <- mean(x)
      y <- scales::rescale(x, to = c(mean_val - spread, mean_val + spread))
      adjust <- mean_val - mean(y)
      y <- y + adjust
      return(y)
    },
    # get_flex_values = function(jump_high=TRUE) {
    #   tmat <- self$tvals
    #   if (isTRUE(jump_high)) {
    #     high_pos <- which(dplyr::lag(self$epoch) != self$epoch & self$epoch == "high")
    #     for (h in high_pos) {
    #       # shift by 60-90 degrees CCW or CW
    #       new_vec <- private$shift_vec(tmat[h,], (-1)^sample(1:2, 1)*sample(seq(90, 180, by=10), size = 1))
    #       tmat[h:nrow(tmat), ] <- pracma::repmat(new_vec, n=nrow(tmat) - h + 1, m=1) # shift subsequent rows
    #     }
    #   }
    # 
    #   # apply the flex
    #   vv <- sapply(seq_len(private$pvt_n_trials), function(i) {
    #     tval <- self$v_rescale(tmat[i,], spread=self$spread[i])
    #     # if (isTRUE(jump_high) && i > 1 && self$epoch[i] == "high" && self$epoch[i-1] != "high") {
    #     #   tval <- private$shift_vec()
    #     # }
    #   })
    #   return(vv)
    # },
    get_next_values = function() {
      if (self$cur_trial > private$pvt_n_trials) {
        warning("You have already iterated through trials. Use $reset_counter to start over")
      }

      private$calculate_values()
      r <- private$pvt_objective_tvals[self$cur_trial, ]
      self$cur_trial <- self$cur_trial + 1

      return(r)
    },
    harvest = function(choice) {
      v <- self$get_next_values()
      p <- runif(1)
      # computationally expensive location of the closest matching position
      pos <- which.min(abs(choice - private$pvt_pvec))
      ifelse(p < self$p_reward, v[pos], 0)
    }
  )
)