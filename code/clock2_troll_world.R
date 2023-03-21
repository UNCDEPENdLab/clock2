
floor_val <- 10


ncenters <- 6
mean_val <- 10
sd_val <- 2
centers <- sample(seq(0, 360, by=10), ncenters, replace=FALSE)
values <- sample(truncnorm::rtruncnorm(ncenters, a=0, mean=mean_val, sd=sd_val))
width_sd <- 0.2 # fixed

bump_prominence <- 10
bump_value <- mean_val*bump_prominence
bump_center <- sample(seq(0, 360, by=10), 1, replace=FALSE)

# gg <- lapply(seq_len(ncenters), function(ii) {
#   r <- rbf$new(value=values[ii], value_sd=0, center=centers[ii], width_sd = width_sd)
# })
# 
# bump_rbf <- rbf$new(value=bump_value, value_sd=0, center=bump_center, width_sd = width_sd)
# 
# gg <- c(gg, bump_rbf)


# VM version
contingency <- vm_circle_contingency(centers = c(centers, bump_center), weights=c(values, bump_value), widths = rep(width_sd, ncenters+1))



# a <- rbf$new(value=10, value_sd=1, center=100, width_sd = 20)
# vv <- a$get_tvec()

#contingency <- rbf_set$new(elements = gg)
plot(contingency$get_vfunc())

contingency$get_centers()
contingency$get_weights()


#contingency$apply_drift(2)
#contingency$get_centers()

for (i in 1:10) {
  contingency$apply_drift(2)
  print(contingency$get_centers())
  plot(contingency$get_vfunc())
  Sys.sleep(.3)
}


troll_world <- R6::R6Class(
  "troll_world",
  private = list(
    shift_vec = function(x, delta) {
      ## converts x[n] into x[n-1] by multiplying the Fourier
      ## transform on x[n] by z^-1 i.e. by e^-iw
      ## lifted from emuR package
      N <- length(x)
      if (delta < 0) delta <- N + delta
      h <- c(rep(0, delta), 1)
      stats::filter(x, h, sides = 1, circular = TRUE)[1:N]
    }
  ),
  public = list(
    erased_segments = list(),
    n_trials = 100,
    cur_trial = 1,
    units="degrees",
    tvals=matrix(0, nrow=100, ncol=360),
    drift_sd = 0,
    spread = NULL, # tracks the point spread between min and max (entropy flexing)
    epoch = NULL, # tracks what phase we are in of the entropy manipulation
    
    initialize = function(n_trials = NULL, tvec = NULL, drift_sd = NULL) {
      if (!is.null(n_trials)) {
        checkmate::assert_integerish(n_trials, lower=1, len=1L)
        self$n_trials <- n_trials
      }
      
      if (!is.null(tvec)) {
        self$tvals <- pracma::repmat(tvec, n=self$n_trials, m=1)
      }
      
      if (!is.null(drift_sd)) {
        checkmate::assert_integerish(drift_sd, lower=1, len=1L)
        self$drift_sd <- drift_sd
      }
      
      # controls spread of values for entropy manipulation
      self$spread <- rep(NA, self$n_trials)
    },
    
    erase_segment = function(size_deg = 30, trial = 1) {
      pos_start <- sample(1:360, size=1)
      pos_end <- pos_start + size_deg - 1
      pos_vec <- pos_start:pos_end
      
      pos_vec <- sapply(pos_vec, function(x) {
        if (x > 360) x = x - 360
        else x
      })
      
      # tmp fill
      self$tvals[trial:self$n_trials, pos_vec] <- pracma::repmat(runif(size_deg), n=length(trial:self$n_trials), m=1)
      return(self)
    },
    get_original_values = function() {
      self$tvals[1,] #starting point
    },
    get_values_matrix = function() {
      self$tvals
    },
    reset_counter = function() {
      self$cur_trial <- 1
      return(self)
    },
    apply_flex = function(low_avg=20, low_spread=5, decrease_avg=10, decrease_spread=5, high_avg=5, high_spread=2, 
                          increase_avg=10, increase_spread=5, spread_max=80, spread_min=10, flip_high=TRUE, start_low=TRUE) {

      t_counter <- 1 # how many trials have we generated
      spread <- rep(NA, self$n_trials)
      epoch <- rep(NA_character_, self$n_trials)
      while (t_counter < self$n_trials) {
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
      self$spread <- spread[1:self$n_trials]
      self$epoch <- epoch[1:self$n_trials]
    },
    v_rescale = function(x, spread = 20, mean_val=50) {
      #mx <- mean(x)
      y <- scales::rescale(x, to = c(mean_val - spread, mean_val + spread))
      adjust <- mean_val - mean(y)
      y <- y + adjust
      return(y)
    },
    get_flex_values = function(flip_high=TRUE) {
      tmat <- self$tvals
      if (isTRUE(flip_high)) {
        high_pos <- which(dplyr::lag(self$epoch) != self$epoch & self$epoch == "high")
        for (h in high_pos) {
          # shift by 60-90 degrees CCW or CW
          new_vec <- private$shift_vec(tmat[h,], (-1)^sample(1:2, 1)*sample(seq(90, 180, by=10), size = 1))
          tmat[h:nrow(tmat),] <- pracma::repmat(new_vec, n=nrow(tmat) - h + 1, m=1) # shift subsequent rows
        }
      }
        
      # apply the flex
      vv <- sapply(seq_len(self$n_trials), function(i) {
        tval <- self$v_rescale(tmat[i,], spread=self$spread[i])
        # if (isTRUE(flip_high) && i > 1 && self$epoch[i] == "high" && self$epoch[i-1] != "high") {
        #   tval <- private$shift_vec()
        # }
      })
      return(vv)
    },
    get_next_values = function() {
      if (self$cur_trial > self$n_trials) {
        warning("You have already iterated through trials. Use $reset_counter to start over")
      }
      r <- self$tvals[self$cur_trial,]
      
      # apply drift
      if (self$drift_sd > 0 && self$cur_trial < self$n_trials) {
        v <- round(rnorm(1, 0, self$drift_sd))
        self$tvals[self$cur_trial + 1,] <- private$shift_vec(r, v)
      }
      
      self$cur_trial <- self$cur_trial + 1
      
      return(r)
    }
  )
  
)

tt <- troll_world$new(n_trials=100, tvec=contingency$get_tvec(), drift_sd=10)
tt$apply_flex(high_avg = 1, high_spread = 0)
plot(tt$spread)


mm <- tt$get_flex_values()
for (ii in 1:ncol(mm)) {
  plot(mm[,ii], type="l", main=paste("Trial", ii, "epoch", tt$epoch[ii]), ylim = range(mm))
  Sys.sleep(.1)
}


tt$erase_segment(30, trial=10)

for (i in 1:10) {
  plot(tt$get_next_values())
  Sys.sleep(.3)
}

# placement of erasures
# - to maintain entropy manipulation, preserve top and bottom 5% of the unerased segments
# - erasure indicators fade instantly (one trial) because the values then enter the overall value distribution
# and become a mix of softmax stochastic and uncertainty-driven. Only the first trial is clear on process.



xx <- tt$get_values_matrix()

plot(xx[15,])


exp_renorm <- function(x, beta=1, min_val, max_val) {
  trans <- exp((x-max(x))/beta) #avoid floating point overflow
  
  #https://www.ibm.com/support/pages/transforming-different-likert-scales-common-scale
  trans <- (max_val - min_val)*(trans - min(trans))/(max(trans) - min(trans)) + min_val
  
  # minshift <- min(trans) - min(x)
  # trans <- trans - minshift
  # rpeak <- (max(trans) - min(trans))/(max(x) - min(x))
  return(trans)
}

#plot(exp_renorm(EV, beta=.1), type="l", main="SCEPTIC value compression", xlab = "Time", ylab = "Value")
plot(EV, type="l", main="SCEPTIC value compression", xlab = "Time", ylab = "Value")
lines(exp_renorm(EV, beta=1, min_val=10, max_val=30), type="l", col="red")
#lines(exp_renorm(EV, beta=3), type="l", col="green")
lines(exp_renorm(EV, beta=10, min_val=15, max_val=20), type="l", col="gray")
# lines(exp_renorm(EV, beta=20), type="l", col="blue")
# lines(exp_renorm(EV, beta=30), type="l", col="orange")
# lines(exp_renorm(EV, beta=1000), type="l", col="orange")

x1 <- exp_renorm(EV, beta=10, min_val=15, max_val=20) + 33.40002
x2 <- exp_renorm(EV, beta=1, min_val=10, max_val=30) + 38.68209

v_rescale <- function(x, spread = 20, mean_val=50) {
  #mx <- mean(x)
  y <- scales::rescale(x, to = c(mean_val - spread, mean_val + spread))
  adjust <- mean_val - mean(y)
  y <- y + adjust
  return(y)
}

aa <- v_rescale(EV, spread=20)
bb <- v_rescale(EV, spread=40)
a2 <- v_rescale(bb, spread=20)

summary(aa)
summary(bb)

plot(bb, type="l", col="blue")
lines(aa, type="l")

v_rescale <- function(x, dynamic_range=1, mean_val=50) {
  #mx <- mean(x)
  min_val <- min(x)
  max_val <- min_val*dynamic_range
  y <- scales::rescale(x, to = c(min_val, max_val))
  adjust <- mean_val - mean(y)
  y <- y + adjust
  return(y)
}

aa <- v_rescale(EV, dynamic_range = 10)
bb <- v_rescale(EV, dynamic_range = 6)
cc <- v_rescale(EV, dynamic_range = 3)
dd <- v_rescale(EV, dynamic_range = 1.5)

# aa <- v_rescale(EV, dynamic_range = 1.5)
# bb <- v_rescale(aa, dynamic_range = 10)
# cc <- v_rescale(EV, dynamic_range = 10)

plot(aa, type="l")
lines(bb, type="l", col="blue")
lines(cc, type="l", col="red")
lines(dd, type="l", col="green")

summary(aa)
summary(bb)



plot(aa, type="l")
plot(bb, type="l")


hist(aa, breaks = 10)
hist(bb, breaks = 10)



# m <- mean(df$value)
# gaussmat <- sapply(seq_along(centers), function(x) {
#   dvec <- dnorm(x=1:360, mean=centers[x], sd=15)
#   dvec <- dvec/max(dvec)*floor_weights[x] #renormalize to max=1
# })
