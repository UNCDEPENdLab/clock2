if (!exists("rcpp_cache_dir")) rcpp_cache_dir <- "~/Rcpp_cache"
if (!exists("rcpp_source_dir")) rcpp_source_dir <- "~/Data_Analysis/clock2/code/cpp"

troll_world <- R6::R6Class(
  "troll_world",
  private = list(
    # private fields
    pvt_n_trials = 100,            # number of trials
    pvt_n_timesteps = 200,         # number of timesteps around circle
    pvt_drift_sd = 0,              # sd of slow drift of contingency across trials (Gaussian random walk)
    pvt_drift_vec = NULL,          # pre-calculated drift vector
    pvt_drift_tvals = NULL,        # contingency values with drift only
    pvt_flex_tvals = NULL,         # contingency values with flex only
    pvt_objective_tvals = NULL,    # objective values with all manipulations
    pvt_jump_vec = NULL,           # shift vector for jump of location during high entropy (flex)
    pvt_clean = FALSE,             # boolean indicating whether contingency values across trials are up to date
    pvt_choices = NULL,            # data.frame containing history of choices made by sampling agent
    pvt_pvec = NULL,               # position of values around circle in units of measurement (radians)
    pvt_has_erasures = FALSE,      # whether this world has erasure dynamics

    # private methods
    # determine end of current erasure condition/block
    get_block_end = function(trial=NULL) {
      if (is.null(trial)) trial <- self$cur_trial

      # show segment -- populate design file
      switch_pos <- which(self$erase_condition[(trial+1):private$pvt_n_trials] != self$erase_condition[trial])[1] - 1
      if (is.na(switch_pos)) {
        block_end <- private$pvt_n_trials # no change here to end
      } else {
        block_end <- trial + switch_pos  
      }
      return(block_end)
    },
    
    # helper function to determine whether a choice in radians fits within an arc whose boundaries are in radians
    seg_test = function(choice, min_rad = NULL, max_rad = NULL) {
      if (is.null(min_rad) || is.na(min_rad)) return(FALSE) # return FALSE if no valid segment
      
      # what is the CW rotation from minimum to maximum?
      seg_width_cw <- max_rad - min_rad
      
      # if difference is negative, we are spanning 0 radians and need to add 2*pi to max
      if (seg_width_cw < -1e-6) seg_width_cw <- max_rad + 2*pi - min_rad
      
      # rotation ccw from point to seg min
      to_min_ccw <- choice - min_rad
      if (to_min_ccw < -1e-6) to_min_ccw <- to_min_ccw + 2*pi # add 2pi to keep units clear
      
      # rotation cw from point to seg min
      to_max_cw <- max_rad - choice
      if (to_max_cw < -1e-6) to_max_cw <- to_max_cw + 2*pi
      
      # use floating point math check since there can be small imprecision problems at the boundary
      in_range <- ifelse(to_min_ccw + to_max_cw - seg_width_cw < 1e-6, TRUE, FALSE)
      return(in_range)
    },

    # helper to shift a vector x circularly by a number of positions delta
    # positive values shift the vector left (early values wrap to the end)
    # negative values shift the vector right (late values wrap to the beginning)
    shift_vec = Rcpp::cppFunction(code = readLines(file.path(rcpp_source_dir, "shift_vec.cpp")), rebuild=FALSE, cacheDir=rcpp_cache_dir),
    
    # internal function to calculate drifting and flexed values
    calculate_values = function() {
      if (isTRUE(private$pvt_clean)) return(invisible(NULL))

      has_drift <- abs(private$pvt_drift_sd) > 1e-5
      has_flex <- !is.null(self$spread)
      has_jump <- sum(abs(private$pvt_jump_vec) > 0)

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
    }
    
  ),
  
  # active bindings
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
    # erasure properties
    erasure_segments = list(),
    disappear_clicks = NULL, # how many clicks until segment disappears
    timeout_trials = NULL, # how many trials without a segment before one can occur again
    erasure_width = pi/6, # fixed width of erasures
    erasure_elig_v_qs = c(.25, .75), # values that can be chosen for the erased value
    erasure_elig_pos_qs = c(.25, .75), # quantiles of positions that can be erased
    
    
    cur_trial = 1,
    units="radians",
    tvals = NULL,
    spread = NULL, # tracks the point spread between min and max (entropy flexing)
    epoch = NULL, # tracks what phase we are in of the entropy manipulation
    erase_condition = NULL, # tracks what erasure phase is active
    trial_type = NULL, # tracks what kind of trial is present on the screen (e.g., erasure)
    p_reward = 0.7, # current fixed probability for reward

    #' @param values A vector of values for every timestep, or a trials x timesteps matrix
    initialize = function(n_trials = NULL, n_timesteps = NULL, values = NULL, drift_sd = NULL) {
      if (!is.null(n_trials)) {
        checkmate::assert_integerish(n_trials, lower=1, len=1L)
        private$pvt_n_trials <- n_trials
      }
      
      if (!is.null(n_timesteps)) {
        checkmate::assert_integerish(n_timesteps, lower=10, len=1L)
        private$pvt_n_timesteps <- n_timesteps
      }

      if (!is.null(values)) {
        if (is.vector(values)) {
          self$tvals <- pracma::repmat(values, n = private$pvt_n_trials, m = 1)
        } else if (is.matrix(values)) {
          # allow trials x timesteps matrix of values as input
          if(nrow(values) != private$pvt_n_trials) {
            warning(glue("Setting n_trials to {nrow(values)} based on values matrix input"))
            private$pvt_n_trials <- nrow(values)
          }
          
          if(ncol(values) != private$pvt_n_timesteps) {
            warning(glue("Setting n_timesteps to {ncol(values)} based on values matrix input"))
            private$pvt_n_timesteps <- ncol(values)
          }
          
          self$tvals <- values
        }
      } else {
        warning("Initializing troll world with all zeros for value")
        self$tvals <- matrix(0, nrow=private$pvt_n_trials, ncol = private$pvt_n_timesteps)
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
      
      # populate 0 jump vector if there are no jumps or apply_flex has not been run
      private$pvt_jump_vec <- rep(0, private$pvt_n_trials)

      if (!is.null(drift_sd)) {
        checkmate::assert_integerish(drift_sd, lower=0, len=1L)
        self$drift_sd <- drift_sd
      }
      
      private$pvt_choices <- data.frame(
        trial = 1:private$pvt_n_trials,
        trial_type = NA_character_,
        choice_rad = NA_real_,
        outcome = NA_real_,
        in_segment = NA,
        segment_shown = NA,
        segment_min = NA_real_,
        segment_max = NA_real_
      )
      
      self$erasure_segments <- data.frame(
        trial = 1:private$pvt_n_trials,
        trial_type = NA_character_,
        segment_shown = NA,
        segment_min = NA_real_,
        segment_max = NA_real_,
        clicks_remain = NA_integer_,
        timeouts_remain = NA_integer_
      )
    },
    
    # method to erase a segment on a given trial
    erase_segment = function(erase = TRUE, trial=NULL) {
      # default to current trial if not specified
      if (is.null(trial)) trial <- self$cur_trial
      
      # get untainted value vector for this trial -- other manipulations (e.g., flex) will be added later
      v <- self$get_values_matrix(type = "original")[trial,]
      loc <- private$pvt_pvec # positions of each timestep around circle
      
      stopifnot(length(v) == length(loc))
      v_new <- v # value function to be updated by erasure
      
      # Examine which positions are eligible for erasure based on their quantiles in the value vector
      qs <- quantile(v, self$erasure_elig_pos_qs)
      pp <- which(v >= qs[1] & v <= qs[2])
      
      # Choose the mid-point of erased segment at random from eligible positions
      mid_pos <- sample(pp, 1) # vector index of midpoint
      mid_rad <- loc[mid_pos]  # position of midpoint in radians
      
      # Find the extent of the segment by going CCW and CW by width/2
      low_rad <- (mid_rad - self$erasure_width/2) %% (2*pi) # location of start point in radians based on location
      high_rad <- (mid_rad + self$erasure_width/2) %% (2*pi) # wrap back onto circle by modulus
      low_pos <- which.min(abs(loc - low_rad)) # vector index of starting point
      high_pos <- which.min(abs(loc - high_rad)) # vector index of ending point
      
      # handle scenario where we wrap around 0/2*pi boundary
      if (low_pos > high_pos) {
        wrap <- TRUE
        erase_indices <- c(1:high_pos, low_pos:length(loc))
      } else {
        wrap <- FALSE
        erase_indices <- low_pos:high_pos
      }

      # erase segment with NAs
      v_new[erase_indices] <- NA

      # find the value for the midpoint of the erasure based on the eligible quantiles
      v_qs <- quantile(v, self$erasure_elig_v_qs)
      v_elig <- v[v > v_qs[1L] & v < v_qs[2L] ]
      v_point <- sample(v_elig, 1) # sample only from eligible quantiles
      # place this value sample/point at the mid-point location for interpolation
      v_new[mid_pos] <- v_point

      # use cubic spline interpolation to connect the dots
      if (wrap) {
        # for spline interpolation, we need to put erasure away from boundaries since na.spline doesn't think circularly
        v_new <- private$shift_vec(v_new, length(v_new)/2)
        v_new <- zoo::na.spline(v_new)
        v_new <- private$shift_vec(v_new, -length(v_new)/2)
      } else {
        v_new <- zoo::na.spline(v_new)
      }
      
      attr(v_new, "ecent") <- mid_pos # center of erasure as vector index
      attr(v_new, "delta") <- v[mid_pos] - v_point # change in value at midpoint
      attr(v_new, "pct_orig") <- (v_point/v[mid_pos]) * 100
      
      # apply erasure to original value vector from here to the end.
      if (isTRUE(erase)) {
        replace_trials <- trial:private$pvt_n_trials
        self$tvals[replace_trials,] <- matrix(v_new, nrow=length(replace_trials), ncol = private$pvt_n_timesteps, byrow=TRUE)
      }
      
      # need to regenerate values when requested since erasure altered tvals
      private$pvt_clean <- FALSE
      
      # show segment -- populate design file
      block_end <- private$get_block_end(trial)
      
      self$erasure_segments$segment_shown[trial:block_end] <- TRUE
      self$erasure_segments$segment_min[trial:block_end] <- low_rad
      self$erasure_segments$segment_max[trial:block_end] <- high_rad
      self$erasure_segments$clicks_remain[trial:block_end] <- self$disappear_clicks
      self$erasure_segments$timeouts_remain[trial:block_end] <- self$timeout_trials
      
      return(v_new)
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
    setup_erasure_blocks = function(erase_prop = 1/3, attention_prop = 1/3, block_length = 15, disappear_clicks = 2, timeout_trials = 0) {
      private$pvt_has_erasures <- TRUE
      self$disappear_clicks <- disappear_clicks
      self$timeout_trials <- timeout_trials

      checkmate::assert_number(erase_prop, lower=0, upper=1)
      checkmate::assert_number(attention_prop, lower = 0, upper = 1)
      stopifnot(erase_prop + attention_prop <= 1)
      no_erasure_prop <- 1 - erase_prop - attention_prop

      div <- ifelse(erase_prop > 0 && attention_prop > 0, 3, 2)

      # allow for blocks to be shortened by up to 5 trials to reduce get proportions right
      # this logic is a bit flawed since if the attention and erasure props diverge a lot, we
      # need something more like the the GCF to determine number of blocks
      possible_lengths <- block_length - 0:5

      rem <- (private$pvt_n_trials / possible_lengths) %% div

      best <- which.min(rem)
      # if a length other than what specified is best, notify user
      if (best != 1) {
        block_length <- block_length - best + 1 # +1 to handle 0 in the first position
        message(sprintf("Adjusting block_length to %d to better balance phases", block_length))
      }

      n_blocks <- floor(private$pvt_n_trials / block_length)
      erase_blocks <- round(erase_prop * n_blocks)
      attention_blocks <- round(attention_prop * n_blocks)
      no_blocks <- round(no_erasure_prop * n_blocks)
      if (erase_blocks + attention_blocks + no_blocks != n_blocks) {
        stop("My rounding is not working!")
      }

      # easy programming for balanced condition
      if (abs(erase_prop - 1/3) < 1e-5 && abs(attention_prop - 1/3) < 1e-5) {
        conditions <- c("no erasure", "erasure", "attention")
        cvec <- as.vector(replicate(n_blocks/3, sample(conditions, 3, replace = FALSE)))
      } else {
        # conditions <- c("no erasure")
        # probs <- no_erasure_prop
        # if (erase_prop > 0) {
        #   conditions <- c(conditions, "erasure")
        #   probs <- c(probs, erase_prop)
        # }
        # if (attention_prop > 0) {
        #   conditions <- c(conditions, "attention")
        #   probs <- c(probs, attention_prop)
        # }

        cvec <- c(rep("erasure", erase_blocks), rep("attention", attention_blocks), rep("no erasure", no_blocks))

        # this is flawed since we have no assurance that the empirical frequencies are close to the target
        #e <- sample(conditions, n_blocks, replace = TRUE, prob = probs) 

        #ec <- rep(e, each=block_length)

        cvec <- sample(cvec, length(cvec))

      }

      self$erase_condition <- rep(cvec, each=block_length)
      aa <- rep(cvec, each=block_length)
      
      if (length(aa) < private$pvt_n_trials) {
        # lengthen last phase if it is not a multiple of the block length
        aa <- c(aa, rep(aa[length(aa)], private$pvt_n_trials - length(aa)))
      } else if (length(aa) > private$pvt_n_trials) {
        # shorten last phase if it goes beyond number of trials
        aa <- aa[1:private$pvt_n_trials]
      }

      self$erase_condition <- aa
      self$erasure_segments$trial_type <- aa

      # seed erased segments that start each attention and erasure phase
      pshift <- which(aa != dplyr::lag(aa, default = "FIRST") & aa != "no erasure")
      for (pp in pshift) {
        e <- ifelse(aa[pp] == "erasure", TRUE, FALSE)
        self$erase_segment(erase = e, trial=pp)
      }
    },

    apply_flex = function(low_avg=20, low_spread=5, decrease_avg=10, decrease_spread=5, high_avg=5, high_spread=2, 
                          increase_avg=10, increase_spread=5, spread_max=80, spread_min=10, jump_high=TRUE, start_low=TRUE) {

      t_counter <- 1 # how many trials have we generated
      spread <- rep(NA, private$pvt_n_trials)
      epoch <- rep(NA_character_, private$pvt_n_trials)
      while (t_counter <= private$pvt_n_trials) {
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

      # handle the jump of location during high entropy periods -- overrides default 0 jump vector
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
    #' @param spread the range of rescaled values
    #' @param mean_val the average of rescaled values
    #' @param force_min if > 0, the minimum rescaled value (can undermine mean)
    v_rescale = Rcpp::cppFunction(code = readLines(file.path(rcpp_source_dir, "v_rescale.cpp")), rebuild=FALSE, cacheDir=rcpp_cache_dir),

    get_cur_values = function() {
      if (self$cur_trial > private$pvt_n_trials) {
        warning("You have already iterated through trials. Use $reset_counter to start over")
      }

      private$calculate_values()
      r <- private$pvt_objective_tvals[self$cur_trial, ]

      return(r)
    },
    get_cur_segment = function() {
      as.list(self$erasure_segments[self$cur_trial,])
    },
    harvest = function(choice, choose_erasure = TRUE) {
      v <- self$get_cur_values()
      s <- self$get_cur_segment()
      
      # computationally expensive location of the closest matching position to choice
      pos <- which.min(abs(choice - private$pvt_pvec))
      
      # force lots of sampling to sort out bugs
      # if (runif(1) < .9) choice <- s$segment_min + 0.05
      
      # test whether choice is in the current erasure/attention segment. Will be FALSE if no segment
      in_seg <- private$seg_test(choice, s$segment_min, s$segment_max)
      
      # harvest probabilistic outcome
      p <- runif(1)
      outcome <- ifelse(p < self$p_reward, v[pos], 0)
      
      # populate info for this trial
      private$pvt_choices$trial_type[self$cur_trial] <- s$trial_type
      private$pvt_choices$choice_rad[self$cur_trial] <- choice
      private$pvt_choices$outcome[self$cur_trial] <- outcome
      private$pvt_choices$in_segment[self$cur_trial] <- in_seg
      private$pvt_choices$segment_shown[self$cur_trial] <- s$segment_shown
      private$pvt_choices$segment_min[self$cur_trial] <- s$segment_min
      private$pvt_choices$segment_max[self$cur_trial] <- s$segment_max
      
      # if the next trial is in the same erase condition, handle dynamics of segments disappearing and appearing
      if (private$pvt_has_erasures && self$cur_trial < private$pvt_n_trials && 
          self$erasure_segments$trial_type[self$cur_trial] == self$erasure_segments$trial_type[self$cur_trial + 1]) {
        
        generate_new_erasure <- FALSE
        clicks_remain <- self$erasure_segments$clicks_remain[self$cur_trial] # how many clicks remain as of this trial
        timeouts_remain <- self$erasure_segments$timeouts_remain[self$cur_trial] # how many timeouts remain as of this trial
        block_end <- private$get_block_end()
        
        # small bug in logic here. The last trial in each block will not be decremented as expected because the opening
        # if block has the cur trial type = trial type + 1. Leaving the bug for now since it doesn't affect dynamics
        # decrement timeouts_remain if we have already hit clicks_remain == 0
        if (!is.na(timeouts_remain) && clicks_remain == 0L) {
          self$erasure_segments$timeouts_remain[(self$cur_trial):block_end] <- timeouts_remain - 1
          
          # if the timeouts have elapsed, sample new erasure
          if (timeouts_remain - 1 == 0) generate_new_erasure <- TRUE
        }
        
        # if person clicks in segment, we should decrement the clicks that remain and potentially resample an erasure if we hit 0
        if (in_seg) {
          # this needs to fire after timeouts collapse, even if they click outside segment
          if (clicks_remain <= 1 && timeouts_remain == 0) {
            # if this was the last click before the segment disappears and there is no timeout, erase a new segment
            generate_new_erasure <- TRUE
          } else if (clicks_remain > 0) {
            # decrement clicks remain
            self$erasure_segments$clicks_remain[(self$cur_trial + 1):block_end] <- clicks_remain - 1
            
            # if we've made the last click, turn off segment display
            if (clicks_remain - 1 == 0) self$erasure_segments$segment_shown[(self$cur_trial + 1):block_end] <- FALSE
          }
        }
        
        if (generate_new_erasure) {
          e <- ifelse(self$erasure_segments$trial_type[self$cur_trial] == "erasure", TRUE, FALSE)
          self$erase_segment(erase = e, trial = self$cur_trial + 1)
        }
      }
      
      # increment trial counter
      self$cur_trial <- self$cur_trial + 1
      
      return(outcome)
    },
    get_pvec = function() {
      private$pvt_pvec
    },
    get_choices = function() {
      private$pvt_choices
    }
  )
)


##
#  # tmp fill
# to fill rows of a matrix with the same vector, we need to recall that matrices are filled column-wise, not row-wise.
# Thus, we can either double transpose the matrix for the operation, or use rep() to replicate the vector for columnwise assignment
# https://stackoverflow.com/questions/11424047/assigning-values-to-rows-from-a-vector
# self$tvals[trial:private$pvt_n_trials, pos_vec] <- runif(size_deg) #pracma::repmat(runif(size_deg), n=length(trial:private$pvt_n_trials), m=1)
# return(self)



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


# handle erasure condition
# if (!is.null(self$erase_condition)) {
#   cur_condition <- self$erase_condition[self$cur_trial]
#   
#   # for now an NA for segment shown indicates that we need to do something
#   if (is.na(self$erasure_segments$segment_shown[self$cur_trial])) {
#     if (cur_condition == "no erasure") {
#       self$erasure_segments$segment_shown[self$cur_trial] <- FALSE # never show
#     } else if (cur_condition == "erasure") {
#       if (self$erasure_segments$segment)
#       self$erase_segment(erase = TRUE)
#         
#     }
#   }
#   ##clicks_remain = NA_integer_,
#   ##timeouts_remain = NA_integer_
#   
# }