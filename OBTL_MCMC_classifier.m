function [error, accuracy] = OBTL_MCMC_classifier(fit_OBTL, X, Label, param, n_t, n_s)


L = param.L;
d = param.d;

r_t_post = cell(1,L);
mu_t_post = cell(1,L);

for j = 1:L
    
    OBTL_data = struct('n_t',n_t{j},'n_s',n_s{j},'d',d,'nu_mu',param.nu_mu,'nu_r',param.nu_r,'M_t',param.M_t{j},'S_t',param.S_t{j}, 'M_s',...
        param.M_s{j},'S_s', param.S_s{j},'rho_mu',param.rho_mu,'rho_r',param.rho_r,'x',[X.t{j};X.s{j}]);
    fit = stan('fit',fit_OBTL,'data',OBTL_data, 'chains',4,'iter',500,'thin',1);
    fit.block();
    
    r_t_post{j} = fit.extract('permuted',true).r_t;
    mu_t_post{j} = fit.extract('permuted',true).mu_t;
    
end

X_test = X.test_all;
Label_test = Label.test_all;
n_test = length(Label_test);

Label_pred = zeros(1,n_test);
for i = 1:n_test
    x = X_test(i,:);
    log_pp = zeros(1,L);
    for j = 1:L
        log_pp(j) = log_posterior_predictive_gamma_NB(x,r_t_post{j},mu_t_post{j});
    end
    [~,I] = max(log_pp);
    Label_pred(i) = I;
end

accuracy = sum(Label_pred == Label_test)/n_test;
error = 1 - accuracy;


