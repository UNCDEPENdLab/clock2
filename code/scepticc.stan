// it looks like the basis should be passed into stan so that it is static
// https://discourse.mc-stan.org/t/how-do-i-make-a-basis-function-matrix/27293/2
// on further examination, it works well to compute it once in transformed data

functions {
  // parallel minimum function
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
    vector[n] x = linspaced_vector(n, 0, 2*pi());

    real kappa = 1/sd; // concentration parameter
    
    // von Mises probability density function (pdf)
    // note that modified_bessel_first_kind takes order (0), then argument (kappa)
    vector[n] pdf = exp(kappa * cos(x - mu)) / (2 * pi() * modified_bessel_first_kind(0, kappa));
    
    // need sum to normalize pdf to AUC 1.0
    real s = sum(pdf);

    pdf = pdf / s;
    
    return(pdf);
  }


  // update vm basis weights: tau is chosen point on circle
  // basis is functions x positions
  // vector update_weights(real alpha, real gamma, real tau, real elig_sd, real outcome, matrix Phi, vector w) {
  //   int p = ncol(Phi); // number of points around circle
  //  
  //   // create eligibility function centered on tau
  //   vector elig = vm_pdf(p, tau, elig_sd);
  //  
  //   vector e[p];
  //   for (i in 1:n_w) {
  //     e[i] = compute_overlap(elig, Phi[i,]);
  //   }
  //
  //   vector pe = outcome - w;
  //   real decay = 0;
  //      
  //   int model = 1; // decay = 1, full = 0 for now
  //   if (model == 1) {
  //     decay = -gamma * (1-e) * w;
  //   }
  //  
  //   w_new = w + alpha * e * pe + decay;
  //   return(w_new);
  // }
 
}

// basis is B functions x P points, computed and passed in
data {
  int<lower=2> B; // number of basis functions
  int<lower=2> P; // number of evaluation points around the circle
 
  int<lower=1>  N; // number of subjects
  int<lower=1>  T; // number of trials (max/complete)
  int<lower=1, upper=T> Tsubj[N]; // number of trials for each subject (if incomplete)

  real<lower=0.01> sb; // standard deviation of basis functions
  real<lower=0.01> sg; // standard deviation of generalization function
  
  //matrix w[T, B]; // basis weights by trial

  real outcomes[N,T]; // continuous outcomes by trial
  real choices[N,T]; // choices matrix with integer value represent point around circle
  int trial_types[N,T]; // trial types: 1=no erasure, 2=erasure, 3=attention
  int<lower=0,upper=1> segment_shown[N,T]; // whether the segment is shown: 0=no, 1=yes
  
}
transformed data {
  // initialize weights at 0
  vector[B] initW;
  initW  = rep_vector(0.0, B);

  // define positions of P evaluation points around the circle
  real d_theta;
  vector[P] pvec; // positions vector
  d_theta = 2*pi()/P;
      
  // since circle wraps, we want to avoid adding a redundant basis function at 0 vs. at 2*pi
  pvec = linspaced_vector(P, 0, 2*pi() - d_theta);

  // setup basis
  matrix[B, P] Phi; // von mises basis matrix

  // setup positions of each basis funciton around the circle
  d_theta = (2*pi())/B;
  
  // since circle wraps, we want to avoid adding a redundant basis function at 0 vs. at 2*pi
  vector[B] centers;
  centers = linspaced_vector(B, 0, 2*pi() - d_theta);

  // populate basis matrix
  // row_vector[P] tt; // temporary row vector allow matrix assignment
  for (i in 1:B) {
    Phi[B,] = to_row_vector(vm_pdf(P, centers[i], sb));
  }

  // precalculate all eligibilities at each point around circle
  matrix[P, B] EB; // here, P is the chosen point around circle
  vector[P] tmp_g; // generalization function, used only for populating EB

  for (p in 1:P) {
      // create eligibility function centered on current position with SD sg
      tmp_g = vm_pdf(P, pvec[p], sg);

      // populate eligibility matrix for this point and basis
      for (b in 1:B) {
        EB[p,b] = compute_overlap(tmp_g, to_vector(Phi[b,]));
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
  vector[N] epsilon_a_pr; // attention selection preference
}
transformed parameters {
  vector<lower=0, upper=1>[N] alpha;
  vector<lower=0>[N] beta;
  vector<lower=0, upper=1>[N] gamma;
  vector<lower=0, upper=1>[N] epsilon_u;
  vector<lower=0, upper=1>[N] epsilon_a;

  // sigmoid transform for gamma and alpha
  // info on transformation here: https://mc-stan.org/docs/2_18/stan-users-guide/reparameterization-section.html
  for (i in 1:N) {
    alpha[i] = Phi_approx(mu_pr[1] + sigma[1] * alpha_pr[i]);
    gamma[i] = Phi_approx(mu_pr[3] + sigma[3] * gamma_pr[i]);
    epsilon_u[i] = Phi_approx(mu_pr[4] + sigma[4] * epsilon_u_pr[i]);
    epsilon_a[i] = Phi_approx(mu_pr[5] + sigma[5] * epsilon_a_pr[i]);
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
  epsilon_a_pr ~ std_normal();

  // loop over subjects and trials

  for (i in 1:N) {
    vector[P] V;  // value vector
    vector[B] e;  // eligibility vector (for each basis)
    vector[B] w;  // weights vector
    vector[B] pe; // basis-wise PEs
    vector[B] decay; // basis-wise decay
    vector[P] g; // generalization function

    w = initW;
    V = to_vector(to_row_vector(w) * Phi); // w is B and Phi is B x P -> 1 x P

    for (t in 1:Tsubj[i]) {
      // predicted softmax choice -- multinomial
      choices[i,t] ~ categorical_logit(beta[i] * V);

      // deprecated: use precalculated eligibilities

      // create eligibility function centered on choice with sg SD
      // g = vm_pdf(P, choices[i,t], sg);

      // for (b in 1:B) {
      //   e[b] = compute_overlap(g, Phi[b,]);
      // }

      e = to_vector(EB[choices[i,t],]);

      // prediction error for each basis
      pe = outcomes[i,t] - w;

      // decay by basis
      // decay = gamma[i] * (1-e) * w;
      decay = gamma[i] * (1 - e) .* w;

      // update weights based on PE and decay
      w = w + e * alpha[i] .* pe - decay;
    }
  }

}

