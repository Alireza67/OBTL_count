data {
    int<lower=0> n_t;
    int<lower=0> n_s;
    int<lower=1> d;
    int<lower=2> nu_mu;
    int<lower=2> nu_r;
    real<lower=0> M_t[d];
    real<lower=0> M_s[d];
    real<lower=0> S_t[d];
    real<lower=0> S_s[d];
    real<lower=0> rho_mu;
    real<lower=0> rho_r;
    int x[n_t+n_s,d];
}
transformed data{
    matrix[2, 2] L_mu[d];
    matrix[2, 2] L_r[d];
    matrix[2, 2] M[d];
    matrix[2, 2] S[d];
    for (i in 1:d) {
        M[i][1,1] = M_t[i];
        M[i][2,2] = M_s[i];
        M[i][1,2] = sqrt(rho_mu * M_t[i] * M_s[i]);
        M[i][2,1] = sqrt(rho_mu * M_t[i] * M_s[i]);
        S[i][1,1] = S_t[i];
        S[i][2,2] = S_s[i];
        S[i][1,2] = sqrt(rho_r * S_t[i] * S_s[i]);
        S[i][2,1] = sqrt(rho_r * S_t[i] * S_s[i]);
        L_mu[i] = cholesky_decompose(M[i]);
        L_r[i] = cholesky_decompose(S[i]);
    }
}
parameters {
    real<lower=0> c1_mu[d];
    real<lower=0> c2_mu[d];
    real z_mu[d];
    real<lower=0> c1_r[d];
    real<lower=0> c2_r[d];
    real z_r[d];
}
transformed parameters{
     matrix[2,2] A_mu[d];
     matrix[2,2] A_r[d];
     matrix[2,2] MU[d];
     matrix[2,2] R[d];
     real mu_t[d];
     real mu_s[d];
     real r_t[d];
     real r_s[d];
     for (i in 1:d){
         A_mu[i][1, 2] = 0.0;
         A_r[i][1, 2] = 0.0;
         A_mu[i][2, 1] = z_mu[i];
         A_r[i][2, 1] = z_r[i];
         A_mu[i][1, 1] = sqrt(c1_mu[i]);
         A_r[i][1, 1] = sqrt(c1_r[i]);
         A_mu[i][2, 2] = sqrt(c2_mu[i]);
         A_r[i][2, 2] = sqrt(c2_r[i]);
         MU[i] = L_mu[i] * A_mu[i] * A_mu[i]' * L_mu[i]';
         R[i] = L_r[i] * A_r[i] * A_r[i]' * L_r[i]';
         mu_t[i] = MU[i][1,1];
         mu_s[i] = MU[i][2,2];
         r_t[i] = R[i][1,1];
         r_s[i] = R[i][2,2];
     }
}
model {
     for (i in 1:d){
          z_mu[i] ~ normal(0, 1);
          z_r[i] ~ normal(0, 1);
          c1_mu[i] ~ chi_square(nu_mu);
          c1_r[i] ~ chi_square(nu_r);
          c2_mu[i] ~ chi_square(nu_mu-1);
          c2_r[i] ~ chi_square(nu_r-1);
     }
    for (i in 1:n_t){
        for (j in 1:d){
            x[i,j] ~ neg_binomial_2(mu_t[j], r_t[j]);
        }
    }
    for (i in n_t+1:n_t+n_s){
        for (j in 1:d){
            x[i,j] ~ neg_binomial_2(mu_s[j], r_s[j]);
        }
    }
}