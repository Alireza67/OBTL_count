clear 
clc

% Add the path of MatlabStan in your computer:

addpath('/scratch/user/karb1367/MatlabStan-2.15.1.0')

% MatlabStan code for OBC:

OBC_code = {
'data {'
'    int<lower=0> n_t;'
'    int<lower=1> d;'
'    int<lower=2> nu_mu;'
'    int<lower=2> nu_r;'
'    real<lower=0> M_t[d];'
'    real<lower=0> S_t[d];'
'    int x_t[n_t,d];'
'}'
'parameters {'
'    real<lower=0> chi_mu_t[d];'
'    real<lower=0> chi_r_t[d];'
'}'
'transformed parameters {'
'    real mu_t[d];'
'    real r_t[d];'
'    for (i in 1:d){'
'        mu_t[i] = M_t[i] * chi_mu_t[i];'
'        r_t[i] = S_t[i] * chi_r_t[i];'
'    }'
'}'
'model {'
'    for (i in 1:d){'
'        chi_mu_t[i] ~ chi_square(nu_mu);'
'        chi_r_t[i] ~ chi_square(nu_r);'
'    }'
'    for (i in 1:n_t){'
'        for (j in 1:d){'
'            x_t[i,j] ~ neg_binomial_2(mu_t[j], r_t[j]);'
'        }'
'    }'
'}'
 };


% MatlabStan code for OBTL:

OBTL_code = {
'data {'
'    int<lower=0> n_t;'
'    int<lower=0> n_s;'
'    int<lower=1> d;'
'    int<lower=2> nu_mu;'
'    int<lower=2> nu_r;'
'    real<lower=0> M_t[d];'
'    real<lower=0> M_s[d];'
'    real<lower=0> S_t[d];'
'    real<lower=0> S_s[d];'
'    real<lower=0> rho_mu;'
'    real<lower=0> rho_r;'
'    int x[n_t+n_s,d];'
'}'
'transformed data{'
'    matrix[2, 2] L_mu[d];'
'    matrix[2, 2] L_r[d];'
'    matrix[2, 2] M[d];'
'    matrix[2, 2] S[d];'
'    for (i in 1:d) {'
'        M[i][1,1] = M_t[i];'
'        M[i][2,2] = M_s[i];'
'        M[i][1,2] = sqrt(rho_mu * M_t[i] * M_s[i]);'
'        M[i][2,1] = sqrt(rho_mu * M_t[i] * M_s[i]);'
'        S[i][1,1] = S_t[i];'
'        S[i][2,2] = S_s[i];'
'        S[i][1,2] = sqrt(rho_r * S_t[i] * S_s[i]);'
'        S[i][2,1] = sqrt(rho_r * S_t[i] * S_s[i]);'
'        L_mu[i] = cholesky_decompose(M[i]);'
'        L_r[i] = cholesky_decompose(S[i]);'
'    }' 
'}'
'parameters {'
'    real<lower=0> c1_mu[d];'
'    real<lower=0> c2_mu[d];'
'    real z_mu[d];'
'    real<lower=0> c1_r[d];'
'    real<lower=0> c2_r[d];'
'    real z_r[d];'
'}'
'transformed parameters{'
'     matrix[2,2] A_mu[d];'
'     matrix[2,2] A_r[d];'
'     matrix[2,2] MU[d];'
'     matrix[2,2] R[d];'
'     real mu_t[d];'
'     real mu_s[d];'
'     real r_t[d];'
'     real r_s[d];'
'     for (i in 1:d){'
'         A_mu[i][1, 2] = 0.0;'
'         A_r[i][1, 2] = 0.0;'
'         A_mu[i][2, 1] = z_mu[i];'
'         A_r[i][2, 1] = z_r[i];'
'         A_mu[i][1, 1] = sqrt(c1_mu[i]);'
'         A_r[i][1, 1] = sqrt(c1_r[i]);'
'         A_mu[i][2, 2] = sqrt(c2_mu[i]);'
'         A_r[i][2, 2] = sqrt(c2_r[i]);'
'         MU[i] = L_mu[i] * A_mu[i] * A_mu[i]'' * L_mu[i]'';'
'         R[i] = L_r[i] * A_r[i] * A_r[i]'' * L_r[i]'';'
'         mu_t[i] = MU[i][1,1];'
'         mu_s[i] = MU[i][2,2];'
'         r_t[i] = R[i][1,1];'
'         r_s[i] = R[i][2,2];'
'     }'
'}'
'model {'
'     for (i in 1:d){'
'          z_mu[i] ~ normal(0, 1);'
'          z_r[i] ~ normal(0, 1);'
'          c1_mu[i] ~ chi_square(nu_mu);'
'          c1_r[i] ~ chi_square(nu_r);'
'          c2_mu[i] ~ chi_square(nu_mu-1);'
'          c2_r[i] ~ chi_square(nu_r-1);'
'     }'
'    for (i in 1:n_t){'
'        for (j in 1:d){'
'            x[i,j] ~ neg_binomial_2(mu_t[j], r_t[j]);'
'        }'
'    }'
'    for (i in n_t+1:n_t+n_s){'
'        for (j in 1:d){'
'            x[i,j] ~ neg_binomial_2(mu_s[j], r_s[j]);'
'        }'
'    }'
'}'
 };

% Building the STAN model for OBC and OBTL:

OBC_data = struct('n_t',2,'d',2,'nu_mu',5,'nu_r',5,'M_t',5 * ones(1,2),'S_t', ones(1,2), 'x_t',[30 50;40 60]);
fit_OBC_initial = stan('model_code',OBC_code,'data',OBC_data, 'model_name','OBC_model','chains',1,'iter',100,'file_overwrite',true);
fit_OBC_initial.block();

OBTL_data = struct('n_t',2,'n_s',2,'d',2,'nu_mu',5,'nu_r',5,'M_t',10 * ones(1,2),'S_t', ones(1,2), 'M_s',10 * ones(1,2),'S_s', ones(1,2),'rho_mu',.9,'rho_r',.9,'x',[30 50;40 60;130 150;140 160]);
fit_OBTL_initial = stan('model_code',OBTL_code,'data',OBTL_data, 'model_name','OBTL_model','chains',1,'iter',100,'file_overwrite',true);
fit_OBTL_initial.block();

load('data_target_source.mat')   % TCGA data for two types of lung cancers LUAD and LUSC for target and source domains. Names of 2857 genes are in gene_names.
load('partition_of_50_LUAD_LUSC_d_10_n_t_5_features_11_to_20.mat')  % First feature set
%load('partition_of_50_LUAD_LUSC_d_10_n_t_5_features_51_to_60.mat')  % Second feature set

d = 10;           % number of features (genes)
L = 2;            % number of classes
n_s{1} = 100;     % number of training source samples in class 1
n_s{2} = 100;     % number of training source samples in class 2
n_t{1} = 5;       % number of training target samples in class 1
n_t{2} = 5;       % number of training target samples in class 2
n_test{1} = 100;  % number of test target samples in class 1
n_test{2} = 100;  % number of test target samples in class 2
partition = 50;   % number of partitions


rho_mu = 0.9;
rho_r = 0.9;
error_OBTL = zeros(1,partition);
error_OBC = zeros(1,partition);
paramchoose = choose_parameters(d, L, rho_mu, rho_r);  % Setting the values of hyperparameters

parfor i = 1:partition
    i
    
    
    % This part is for paralellization 
    
    if ~exist([tempdir 'worker_' num2str(i)], 'dir')
        mkdir([tempdir 'worker_' num2str(i)]);
    end
    
    OBC_data(i) = struct('n_t',2,'d',2,'nu_mu',5,'nu_r',5,'M_t',10 * ones(1,2),'S_t', ones(1,2), 'x_t',[30 50;40 60]);
    fit_OBC(i) = stan('fit',fit_OBC_initial,'data',OBC_data(i),'chains',1,'iter',100,'working_dir',[tempdir 'worker_' num2str(i)]);
    fit_OBC(i).block();
    
    OBTL_data(i) = struct('n_t',2,'n_s',2,'d',2,'nu_mu',5,'nu_r',5,'M_t',10 * ones(1,2),'S_t', ones(1,2), 'M_s',10 * ones(1,2),'S_s', ones(1,2),'rho_mu',.9,'rho_r',.9,'x',[30 50;40 60;130 150;140 160]);
    fit_OBTL(i) = stan('fit',fit_OBTL_initial,'data',OBTL_data(i),'chains',1,'iter',100,'working_dir',[tempdir 'worker_' num2str(i)]);
    fit_OBTL(i).block();
    
    
    % Classifiers (OBTL and OBC):
    
    [error_OBTL(i), ~] = OBTL_MCMC_classifier(fit_OBTL(i), X{i}, Label{i}, paramchoose, n_t, n_s);
    [error_OBC(i), ~] = OBC_MCMC_classifier(fit_OBC(i), X{i}, Label{i}, paramchoose, n_t);
    
end

% Average classification error (over 50 partitions) of OBTL and OBC:
error_OBTL_mean = mean(error_OBTL);
error_OBC_mean = mean(error_OBC);

fprintf(' The average error of OBTL is: %5.3f\n', error_OBTL_mean)
fprintf(' The average error of OBC is: %5.3f\n', error_OBC_mean)

