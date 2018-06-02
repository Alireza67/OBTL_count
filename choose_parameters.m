function paramchoose = choose_parameters(d, L, rho_mu, rho_r)


paramchoose.d = d;
paramchoose.L = L;
paramchoose.nu_mu = 2;
paramchoose.nu_r = 2;
paramchoose.rho_mu = rho_mu;
paramchoose.rho_r = rho_r;

for i = 1:L
    paramchoose.M_t{i} = 2000 * ones(1,d)/paramchoose.nu_mu;
    paramchoose.M_s{i} = 500 * ones(1,d)/paramchoose.nu_mu;
    paramchoose.S_t{i} = 4 * ones(1,d)/paramchoose.nu_r;
    paramchoose.S_t_OBC{i} = 4 * ones(1,d)/paramchoose.nu_r;
    paramchoose.S_s{i} = 4 * ones(1,d)/paramchoose.nu_r;
end

