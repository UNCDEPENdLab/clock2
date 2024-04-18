// Basis function matrix could be passed in as data: 
// https://discourse.mc-stan.org/t/how-do-i-make-a-basis-function-matrix/27293/2
// On further examination, however, it works well to compute it once in transformed data, which makes the code more standalone

functions {
  // adaptation of .bincode from base R for binning continuous choices
  int[] bincode(vector x, vector breaks, int right, int include_border) {
    int n = num_elements(x);
    int nb = num_elements(breaks);
    int code[n]; // n element integer vector return

    int lo, hi, nb1 = nb - 1, new;
    int lft = !right;

    // This relies on breaks being sorted, so wise to check that
    for (i in 2:nb) 
	    if (breaks[i-1] > breaks[i]) reject("'breaks' is not sorted");

    for (i in 1:n) {
	    // code[i] = nan(); // not used in stan -- default is NaN anyhow
      if (!is_nan(x[i])) {
	      lo = 1;
	      hi = nb;

	      if (x[i] <  breaks[lo] || breaks[hi] < x[i] ||
	        (x[i] == breaks[lft ? hi : lo] && ! include_border));
	      else {
		      while (hi - lo >= 2) {
		        new = (hi + lo)/2;
		        if (x[i] > breaks[new] || (lft && x[i] == breaks[new]))
			        lo = new;
		        else
			        hi = new;
          }
        }
        code[i] = lo;
      }
    }

    return(code);
  }

  // parallel minimum function used by compute_overlap
  vector pmin(vector a, vector b) {
    if (num_elements(b) != num_elements(a))
      reject("Number of elements in a and b vectors to pmin are unequal. a=", num_elements(a), ", b=", num_elements(b));
      
    vector[num_elements(b)] res;
    for (i in 1:num_elements(b)) {
      res[i] = min({a[i], b[i]});
    }

    return res;
  }
  
  // compute proportion overlap between two vectors
  real compute_overlap(vector b1, vector b2) {
    real s1 = sum(b1);
    if (abs(s1) - 1.0 > 1e-5)
      reject("b1 does not sum to 1. Sum=", s1);
    
    real s2 = sum(b2);
    if (abs(s2) - 1.0 > 1e-5) 
      reject("b2 does not sum to 1. Sum=", s2);
    
    real res = sum(pmin(b1, b2));
    return(res);
  }
  
  // von mises basis in radians
  vector vm_pdf(int n, real mu, real sd) {
  
    // generate positions between 0 and 2*pi
    vector[n] x = linspaced_vector(n, 0.0, 2.0*pi());

    real kappa = 1.0/sd; // concentration parameter
    
    // von Mises probability density function (pdf)
    // note that modified_bessel_first_kind takes order (0), then argument (kappa)
    vector[n] pdf = exp(kappa * cos(x - mu)) / (2.0 * pi() * modified_bessel_first_kind(0, kappa));
    
    // need sum to normalize pdf to AUC 1.0
    real s = sum(pdf);
    pdf = pdf / s;
    
    return(pdf);
  }

}

// basis is B functions x P points
data {
  int<lower=2> B; // number of basis functions
  int<lower=2> P; // number of evaluation points around the circle
 
  int<lower=1>  N; // number of subjects
  int<lower=1>  T; // number of trials (max/complete)
  int<lower=1, upper=T> Tsubj[N]; // number of trials for each subject (if incomplete)

  real<lower=0.01> sb; // standard deviation of basis functions
  real<lower=0.01> sg; // standard deviation of generalization function
  
  real outcomes[N,T];     // continuous outcomes by trial
  real choices_rad[N,T];  // choices matrix in radians
  int trial_types[N,T];   // trial types: 1=no erasure, 2=erasure, 3=attention
  int<lower=0,upper=1> segment_shown[N,T]; // whether the segment is shown: 0=no, 1=yes
  real segment_min[N,T];  // segment minimum in radians
  real segment_max[N,T];  // segment maximum in radians
}

/* 
for epsilon_u and epsilon_int, we need
a) a 0/1 indicator function S of the same size as V (P points) for each trial, where 1 indicates the erased or attentional segment
b) parameter epsilon for general segment preference
c) parameter omega for erasure preference relative to attention
d) a mixing choice rule choices[i,t] ~ categorical_logit(V / beta[i] + (epsilon + omega) * S);

To calculate S, we need the segment_min and segment_max locations and to bin them identically to choices
*/

transformed data {
  // initialize basis weights at 0
  vector[B] initW;
  initW  = rep_vector(0.0, B);
 
  // setup basis
  matrix[B, P] Phi;   // von mises basis matrix
  vector[B] centers;  // of each basis function in radians
  vector[P] pvec;     // vector of positions (in radians) at which value function is evaluated
  vector[P+1] breaks; // breaks for cutting observed choices (in radians) to P bins
  real d_theta;       // spacing in radians

  // define positions of P evaluation points around the circle
  d_theta = 2.0*pi()/P;

  // since circle wraps, we want to avoid having redundant evaluation points at 0 and 2*pi
  // these positions will form the centers of binned data, and to avoid 0/2pi wrapping in the bin-cutting,
  // start first points at 1/2 d_theta above 0 radians and end last point 1/2 d_theta below 2*pi
  pvec = linspaced_vector(P, 0.5*d_theta, 2.0*pi() - 0.5*d_theta); // N.B. As in C++, don't use 1/2 or it converts to int!
  
  breaks[1:P] = pvec  - 0.5 * d_theta; // shift centers ccw by 1/2 d_theta to start at 0
  breaks[P+1] = 2.0*pi(); // last cut sits at 2*pi to complete the circle

  // discretize choices into bins for multinomial fitting using evaluation points
  int choices[N,T]; // choices matrix as integer bins

  for (i in 1:N) {
    //choices[i,] = to_row_vector(bincode(to_vector(choices_rad[i,]), breaks, 1, 1));
    choices[i,] = bincode(to_vector(choices_rad[i,]), breaks, 1, 1);
  }

  // setup positions of B basis functions around the circle
  d_theta = 2.0*pi()/B;

  // since circle wraps, we want to avoid adding a redundant basis function at 0 vs. at 2*pi
  centers = linspaced_vector(B, 0.0, 2.0*pi() - d_theta);

  // populate basis matrix
  for (i in 1:B) {
    Phi[i,] = to_row_vector(vm_pdf(P, centers[i], sb));
  }

  // precalculate all basis eligibilities at each possible choice bin around circle
  matrix[P, B] EB; // here, P is the chosen point around circle
  vector[P] tmp_g; // generalization function, used only for populating EB

  for (p in 1:P) {
      // create eligibility function centered on current position p with SD sg
      tmp_g = vm_pdf(P, pvec[p], sg);

      // populate eligibility matrix for this point and basis
      for (b in 1:B) {
        EB[p,b] = compute_overlap(tmp_g, to_vector(Phi[b,]));
      }
  }

  // precalculate segment indicator matrix
  // simplest approach to avoid lots of type conversion complaints is to create a 2D array with P-length vectors, rather than a 3-D array
  array[N,T] vector[P] S; // subjects x trials array where each element is a P-length vector
  
  int min_pos[N,T]; // segment minimum in bins
  int max_pos[N,T]; // segment maximum in bins

  // convert segment positions into bins, mirroring choice bins
  for (i in 1:N) {
    //if (segment_min[i,] < 0.0) segment_min[i,] = nan();
    //if (segment_max[i,] < 0.0) segment_max[i,] = nan();
    min_pos[i,] = bincode(to_vector(segment_min[i,]), breaks, 1, 1);
    max_pos[i,] = bincode(to_vector(segment_max[i,]), breaks, 1, 1);
  }

  // now populate 0/1 vectors for each trial based on bins
  vector[P] svec; // P-length indicator vector
  int l;
  int m;
  for (i in 1:N) {
    for (t in 1:Tsubj[i]) {
      svec = rep_vector(0.0, P);
      
      // negative values in segment min denote non-segment (no erasure) trials
      if (segment_min[i,t] > 0.0) {
        l = min_pos[i,t];
        m = max_pos[i,t];
        if (m < l) {
          // handle wrapping at 0
          svec[l:num_elements(svec)] = rep_vector(1.0, num_elements(svec) - l + 1);
          svec[1:m] = rep_vector(1.0, m);
        } else {
          svec[l:m] = rep_vector(1.0, m - l + 1);
        }
      }
      
      S[i,t] = svec;
    }
  }
}

parameters {
  // group hyperparameters
  // 1 = alpha
  // 2 = beta
  // 3 = gamma
  // 4 = epsilon_u
  // 5 = epislon_a
  vector[5] mu_pr; 
  vector<lower=0>[5] sigma;

  vector[N] alpha_pr;     // learning rate
  vector[N] beta_pr;      // temperature
  vector[N] gamma_pr;     // decay
  vector[N] epsilon_u_pr; // uncertainty selection preference
  vector[N] epsilon_int_pr; // overall (intercept) segment selection preference
}

transformed parameters {
  vector<lower=0, upper=1>[N] alpha;
  vector<lower=0>[N] beta;
  vector<lower=0, upper=1>[N] gamma;
  // vector<lower=0, upper=1>[N] epsilon_u;
  // vector<lower=0, upper=1>[N] epsilon_int;
  vector[N] epsilon_u;
  vector[N] epsilon_int;

  // sigmoid transform for gamma and alpha
  // info on transformation here: https://mc-stan.org/docs/2_18/stan-users-guide/reparameterization-section.html
  for (i in 1:N) {
    alpha[i] = Phi_approx(mu_pr[1] + sigma[1] * alpha_pr[i]);
    gamma[i] = Phi_approx(mu_pr[3] + sigma[3] * gamma_pr[i]);
    // epsilon_u[i] = Phi_approx(mu_pr[4] + sigma[4] * epsilon_u_pr[i]);
    // epsilon_int[i] = Phi_approx(mu_pr[5] + sigma[5] * epsilon_int_pr[i]);

    // switch to continuous normal case (allowing negatives)
    epsilon_u[i] = mu_pr[4] + sigma[4] * epsilon_u_pr[i];
    epsilon_int[i] = mu_pr[5] + sigma[5] * epsilon_int_pr[i];
  }

  // exponential transform on temperature
  beta = exp(mu_pr[2] + sigma[2] * beta_pr);
}

model {
  // hyperparameters
  mu_pr[1:5] ~ std_normal(); // currently all ~N(0,1) on hyperparameters
  sigma[1:5] ~ normal(0, 0.2); // individual variances in hyperparameters

  // untransformed individual parameters
  alpha_pr ~ std_normal();
  beta_pr ~ std_normal();
  gamma_pr ~ std_normal();
  epsilon_u_pr ~ std_normal();
  epsilon_int_pr ~ std_normal();

  // loop over subjects and trials
  for (i in 1:N) {
    vector[P] V;  // value vector
    vector[P] st;  // segment vector on this trial
    vector[B] e;  // eligibility vector (for each basis)
    vector[B] w;  // weights vector
    vector[B] pe; // basis-wise PEs
    vector[B] decay; // basis-wise decay
    vector[P] g; // generalization function
    real mV;    // mean value
    real epsilon_t; // current trial epsilon (based on trial type and whether segment is highlighted)

    w = initW;
    
    //print("beta: ", beta[i]);
    //print("V: ", V);

    for (t in 1:Tsubj[i]) { // trial types: 1=no erasure, 2=erasure, 3=attention
      if (trial_types[i,t] == 1) {
        epsilon_t = 0.0;
      } else if (trial_types[i,t] == 2) {
        epsilon_t = epsilon_int[i] + epsilon_u[i];
      } else if (trial_types[i,t] == 3) {
        epsilon_t = epsilon_int[i];
      }

      // compute evaluated value function
      V = to_vector(to_row_vector(w) * Phi); // w is B and Phi is B x P -> 1 x P
      mV = mean(V);
      // st = S[i,t] .* mV; // original parameterization, scaling by mean value
      st = S[i,t]; // no scaling by mean value
      V = V - max(V); // subtract max to make categorical_logit easier
      
      // V = rep_vector(0.0, P); // debug with a forced value vector of 0
      // print("V: ", V);

      // predicted softmax choice -- multinomial
      // choices[i,t] ~ categorical_logit(V / beta[i] + epsilon_t .* st); // original parameterization
      choices[i,t] ~ categorical_logit((V + epsilon_t .* st) / beta[i]);

      // look up eligibility for each basis using the position of the choice
      e = to_vector(EB[choices[i,t],]);

      // prediction error for each basis
      pe = outcomes[i,t] - w;

      // decay by basis
      decay = gamma[i] * (1 - e) .* w;
      // decay = alpha[i] * gamma[i] * (1 - e) .* w;
      // update weights based on PE and decay
      w = w + alpha[i] * e .* pe - decay;
    }
  }

}

generated quantities {

  // placing code within a local section using curly braces removes these objects from the output
  
  // debug basis setup
  vector[P] p_pvec; // positions vector
  p_pvec = pvec;
  matrix[P, B] p_EB; // here, P is the chosen point around circle
  p_EB = EB;
  vector[P] p_tmp_g; // generalization function, used only for populating EB
  p_tmp_g = tmp_g;
  vector[B] p_centers; // centers of basis functions
  p_centers = centers;
  matrix[B, P] p_Phi; // von mises basis matrix
  p_Phi = Phi;

  int p_choices[N,T]; // binned choices
  p_choices = choices;
  vector[P+1] p_breaks; // breaks for cutting observed choices (in radians) to P bins
  p_breaks = breaks;

  // array[N,T] vector[P] p_S; // subjects x trials x points
  /*
  array[5,T] vector[P] p_S; // subjects x trials x points
  for (i in 1:5) {
    for (j in 1:T) {
      //print("S[i,j]: ", S[i,j]);
      p_S[i,j] = S[i,j];
    }
  }
  */

  // p_S = S; // subjects x trials x points
  
  // run sceptic model loop to obtain predictions for auditing (in progress)
  /*
  {
    for (i in 1:N) {
      vector[P] V;  // value vector
      vector[B] e;  // eligibility vector (for each basis)
      vector[B] w;  // weights vector
      vector[B] pe; // basis-wise PEs
      vector[B] decay; // basis-wise decay
      vector[P] g; // generalization function

      w = initW;
      
      //print("beta: ", beta[i]);
      //print("V: ", V);

      for (t in 1:Tsubj[i]) {
        // compute evaluated value function
        V = to_vector(to_row_vector(w) * Phi); // w is B and Phi is B x P -> 1 x P
        // V = rep_vector(0.0, P);
        // print("V: ", V);

        // predicted softmax choice -- multinomial
        choices[i,t] ~ categorical_logit(beta[i] * V);

        // look up eligibility for each basis using the position of the choice
        e = to_vector(EB[choices[i,t],]);

        // prediction error for each basis
        pe = outcomes[i,t] - w;

        // decay by basis
        // decay = gamma[i] * (1-e) * w;
        decay = gamma[i] * (1 - e) .* w;

        // update weights based on PE and decay
        w = w + alpha[i] * e .* pe - decay;
      }
    }

  }
  */

}
