## it looks like the basis shoudl be passed into stan so that it is static
## https://discourse.mc-stan.org/t/how-do-i-make-a-basis-function-matrix/27293/2

# basis is B functions x P points, computed and passed in
data {
  int<lower=1> i; // number of trials
  int<lower=2> b; // number of basis functions
  int<lower=2> p; // number of evaluation points around the circle

  int<lower=1>  N; //number of subjects
  int<lower=1, upper=i> isubj[N]; // number of trials for each subject (if incomplete)

  matrix Phi[b, p]; // von mises basis (passed in)
  matrix w[i,b]; // basis weights by trial

  vector outcomes[i];
  vector choices[i];
  
}
parameters {
  vector[N] alpha;     // learning rate
  vector[N] beta;      // temperature
  vector[N] gamma;     // decay
  vector[N] epsilon_u; // uncertainty selection preference
  vector[N] epsilon_a; // attention selection preference
}

model {
  for ()


}

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
  
  #### von mises basis in radians
  vector vm_pdf(int n, real mu, real sd) {
  
    // generate positions between 0 and 2*pi
    vector x = linspaced_vector(n, 0, 2*pi);

    real kappa = 1/sd; // concentration parameter
    
    // von Mises probability density function (pdf)
    // note that modified_bessel_first_kind takes order (0), then argument (kappa)
    vector pdf = exp(kappa * cos(x - mu)) / (2 * pi * modified_bessel_first_kind(0, kappa));
    
    // need sum to normalize pdf to AUC 1.0
    real s = sum(pdf);

    pdf = pdf / s;
    
    return(pdf);
  }


  #### update vm basis weights: tau is chosen point on circle
  // basis is functions x positions
  vector update_weights(real alpha, real gamma, real tau, real elig_sd, real outcome, matrix Phi, vector w) {
    int p = ncol(Phi); // number of points around circle
    
    // create eligibility function centered on tau
    vector elig = vm_pdf(p, tau, elig_sd);
    
    vector e[p];
    for (i in 1:n_w) {
      e[i] = compute_overlap(elig, Phi[i,]);
    }

    vector pe = outcome - w;
    real decay = 0;
        
    int model = 1; // decay = 1, full = 0 for now
    if (model == 1) {
      decay = -gamma * (1-e) * w;
    }
    
    w_new = w + alpha * e * pe + decay;
    return(w_new);
  }
 
}
 