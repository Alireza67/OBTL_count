data {
    int<lower=0> n_t;
    int<lower=1> d;
    int<lower=2> nu_mu;
    int<lower=2> nu_r;
    real<lower=0> M_t[d];
    real<lower=0> S_t[d];
    int x_t[n_t,d];
}
parameters {
    real<lower=0> chi_mu_t[d];
    real<lower=0> chi_r_t[d];
}
transformed parameters {
    real mu_t[d];
    real r_t[d];
    for (i in 1:d){
        mu_t[i] = M_t[i] * chi_mu_t[i];
        r_t[i] = S_t[i] * chi_r_t[i];
    }
}
model {
    for (i in 1:d){
        chi_mu_t[i] ~ chi_square(nu_mu);
        chi_r_t[i] ~ chi_square(nu_r);
    }
    for (i in 1:n_t){
        for (j in 1:d){
            x_t[i,j] ~ neg_binomial_2(mu_t[j], r_t[j]);
        }
    }
}