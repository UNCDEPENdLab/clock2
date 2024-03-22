# scepticc inversion in stan
ncenters <- 9 # how many gaussians there are
mean_val <- 10 # mean reward rate
sd_val <- 2 # standard deviation of reward / range of rewards
centers <- sample(seq(0, 360, by = 10), ncenters, replace = FALSE) # line up gaussians here
values <- sample(truncnorm::rtruncnorm(ncenters, a = 0, mean = mean_val, sd = sd_val))
width_sd <- 20 # fixed, how wide are the underlying Gaussians

bump_prominence <- 8 # bump will always be higher, but it will change
bump_value <- mean_val * bump_prominence

centers <- (pi / 180) * centers
width_sd <- (pi / 180) * width_sd
bump_center <- (pi / 180) * bump_center
# 
contingency <- vm_circle_contingency(centers = c(centers, bump_center), weights = c(values, bump_value), widths = rep(width_sd, ncenters + 1), units = "radians")
ttd <- troll_world$new(n_trials=300, values=contingency$get_wfunc(), drift_sd=5)
ttd$apply_flex(high_avg = 1, high_spread = 0, low_avg = 10, spread_max = 150, jump_high = T)
ttd$setup_erasure_blocks(disappear_clicks = 2, timeout_trials = 1)
plot(ttd$get_starting_values())
plot(ttd$spread)



generate_scepticc_data <- function(n=25, contingency = NULL) {
  require(truncnorm)
  if (!inherits(contingency, "troll_world")) stop("contingency must be troll_world object")
  
  seeds <- sample(1:1e4, size=n, replace=FALSE)
  
  alpha_m <- 0.1
  alpha_sd <- 0.1
  alphas <- rtruncnorm(n, a=0, b=1, mean=alpha_m, sd=alpha_sd)
  
  beta_m <- 0.1
  beta_sd <- 2
  betas <- rtruncnorm(n, a=0, mean=beta_m, sd=beta_sd)
  
  gamma_m <- 0.3
  gamma_sd <- 0.1
  gammas <- rtruncnorm(n, a=0, b=1, mean=gamma_m, sd=gamma_sd)
  
  epsilon_u_m <- 1/12 # 8%
  epsilon_u_sd <- 1/8
  epsilons <- rtruncnorm(n, a=0, b=1, mean=epsilon_u_m, sd=epsilon_u_sd)
  
  epsilon_u_m <- 1/12 # 8%
  epsilon_u_sd <- 1/8
  epsilon_us <- rtruncnorm(n, a=0, b=1, mean=epsilon_u_m, sd=epsilon_u_sd)
    
  epsilon_a_m <- 1/12 # 8%
  epsilon_a_sd <- 1/8
  epsilon_as <- rtruncnorm(n, a=0, b=1, mean=epsilon_a_m, sd=epsilon_a_sd)
  
  dlist <- lapply(seq_len(n), function(i) {
    ag <- scepticc$new(n_basis=12, n_points=360, contingency=contingency,
                       alpha = alphas[i], beta = betas[i], gamma = gammas[i],
                       epsilon_u = epsilon_us[i], epsilon_a = epsilon_as[i], seed=seeds[i])
    
    learning_history <- ag$run_contingency(optimize = FALSE)
    world <- ag$get_contingency()
    df <- world$get_choices()
    return(df)
    
  })
  
  ddf <- tibble(id=1:n, alphas, gammas, betas, epsilon_us, epsilon_as, data=dlist) %>%
    unnest(data)
  
}

ff <- generate_scepticc_data(n=25, contingency = ttd)



stan_str <- "
functions {
  ## parallel minimum function: DONE
  vector pmin(vector a, vector b) {
    if (num_elements(b) != num_elements(a))
      reject('Number of elements in a and b vectors to pmin are unequal');
      
    vector[num_elements(b)] res;
    for(i in 1:num_elements(b)) {
      res[i] = min({a[i], b[i]});
    }

    return res;
  }
  
  # compute proportion overlap between two vectors: DONE
  real compute_overlap(vector b1, vector b2) {
    real s1 = sum(b1);
    if (abs(s1) - 1.0 > 1e-5)
      reject('b1 does not sum to 1');
    
    real s2 = sum(b2);
    if (abs(s2) - 1.0 > 1e-5) 
      reject('b2 does not sum to 1');
    
    real res = sum(pmin(b1, b2));
    return(res);
  }
  
  # von mises basis in radians
  vector vm_pdf(int n, real mu, real sd) {
  
    // generate positions between 0 and 2*pi
    vector x = linspaced_vector(n, low, 2*pi);

    real kappa = 1/sd; // concentration parameter
    
    // von Mises probability density function (pdf)
    // note that modified_bessel_first_kind takes order (0), then argument (kappa)
    pdf = exp(kappa * cos(x - mu)) / (2 * pi * modified_bessel_first_kind(0, kappa))
    
    // need sum to normalize pdf to AUC 1.0
    real s = sum(pdf);

    pdf = pdf / s
    
    return(pdf)
    }
    
}
          private$compute_overlap(e$get_basis(), efunc$get_basis())

  #### update vm basis weights: tau is chosen point on circle
  vector update_weights(real tau, real elig_sd, real outcome, vector w, ) {
    int n_w = num_elements(w);
    
    // create eligibility function centering on tau
    
    
    vector e[n_w];
    for (i in 1:n_w) {
      e[i] = compute_overlap()
    }
    
            private$compute_overlap(e$get_basis(), efunc$get_basis())
  
      checkmate::assert_number(tau)
      private$pvt_eligibility$center <- tau
      vector e
      e <- private$pvt_bf_set$get_eligibilities(private$pvt_eligibility)
      w <- private$pvt_bf_set$get_weights()
      pe <-  outcome - w
      
      int model = 1;
      if (model == 1) {
        decay <- -self$gamma * (1-e) * w  
      } else {
        decay <- 0
      }
      
      w_new <- w + self$alpha*e*pe + decay
      private$pvt_bf_set$set_weights(w_new)
      private$pvt_history[[private$pvt_sample]] <- w_new
      private$pvt_sample <- private$pvt_sample + 1
      return(self)
    }

}
"



