function log_pp = log_posterior_predictive_gamma_NB(x,r_t,mu_t)

M = size(r_t,1);
X = repmat(x,M,1);

log_p = gammaln(X + r_t) - gammaln(r_t) - gammaln(X + 1) + X .* log(mu_t./(mu_t + r_t)) + r_t .* log(r_t./(r_t + mu_t));
C = sum(log_p,2);
C_max = max(C);
log_pp = C_max + log((1/M) * sum(exp(C - C_max)));