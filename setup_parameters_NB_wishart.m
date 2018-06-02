function param = setup_parameters_NB_wishart(d, L)

param.d = d;
param.L = L;
nu_mu = 6;
param.nu_mu = nu_mu;
nu_r = 6;
param.nu_r = nu_r;

for i = 1:d
    param.rho_mu{1}(i) = .9;   % correlation in [0 1)
    param.M_t{1}(i) = 1000/nu_mu;  % expected = 1000
    param.M_s{1}(i) = 5000/nu_mu;
    param.M_ts{1}(i) = sqrt(param.rho_mu{1}(i) * param.M_t{1}(i) * param.M_s{1}(i));
    param.M{1}{i} = [param.M_t{1}(i) param.M_ts{1}(i); param.M_ts{1}(i) param.M_s{1}(i)];
    
    param.rho_mu{2}(i) = .9;   % correlation in [0 1)
    param.M_t{2}(i) = 1500/nu_mu;
    param.M_s{2}(i) = 6000/nu_mu;
    param.M_ts{2}(i) = sqrt(param.rho_mu{2}(i) * param.M_t{2}(i) * param.M_s{2}(i));
    param.M{2}{i} = [param.M_t{2}(i) param.M_ts{2}(i); param.M_ts{2}(i) param.M_s{2}(i)];
    
    
    param.rho_r{1}(i) = .9;   % correlation in [0 1)
    param.S_t{1}(i) = 1/nu_r;  % expected = .1
    param.S_s{1}(i) = .5/nu_r;
    param.S_ts{1}(i) = sqrt(param.rho_r{1}(i) * param.S_t{1}(i) * param.S_s{1}(i));
    param.S{1}{i} = [param.S_t{1}(i) param.S_ts{1}(i); param.S_ts{1}(i) param.S_s{1}(i)];
    
    param.rho_r{2}(i) = .9;   % correlation in [0 1)
    param.S_t{2}(i) = 1/nu_r;
    param.S_s{2}(i) = .5/nu_r;
    param.S_ts{2}(i) = sqrt(param.rho_r{2}(i) * param.S_t{2}(i) * param.S_s{2}(i));
    param.S{2}{i} = [param.S_t{2}(i) param.S_ts{2}(i); param.S_ts{2}(i) param.S_s{2}(i)];
    
end





