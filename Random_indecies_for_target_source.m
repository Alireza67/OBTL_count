function [X, Label] = Random_indecies_for_target_source(X_tr_0, X_tr_1, X_sr_0, X_sr_1, d, n_t, n_s, n_test, partition, feature_set)

N_t_0 = size(X_tr_0,2);
N_t_1 = size(X_tr_1,2);
N_s_0 = size(X_sr_0,2);
N_s_1 = size(X_sr_1,2);

idx_t_0_set = 1:N_t_0;
idx_t_1_set = 1:N_t_1;
idx_s_0_set = 1:N_s_0;
idx_s_1_set = 1:N_s_1;

n_t_0 = n_t{1};
n_t_1 = n_t{2};
n_s_0 = n_s{1};
n_s_1 = n_s{2};
n_test_0 = n_test{1};
n_test_1 = n_test{2};

X = cell(1,partition);
Label = cell(1,partition);
features = feature_set(1:d);

for i = 1:partition
        %features = datasample(1:m, d, 'Replace', false);
        
        idx_t_0 = datasample(idx_t_0_set, n_t_0, 'Replace', false);
        X{i}.t{1} = X_tr_0(features, idx_t_0)';
        temp_0 = idx_t_0_set;
        temp_0(idx_t_0) = [];
        idx_test_0 = datasample(temp_0, n_test_0, 'Replace', false);
        X{i}.test{1} = X_tr_0(features, idx_test_0)';
        Label{i}.t{1} = ones(1, n_t_0);
        Label{i}.test{1} = ones(1, n_test_0);
        
        
        idx_t_1 = datasample(idx_t_1_set, n_t_1, 'Replace', false);
        X{i}.t{2} = X_tr_1(features, idx_t_1)';
        temp_1 = idx_t_1_set;
        temp_1(idx_t_1) = [];
        idx_test_1 = datasample(temp_1, n_test_1, 'Replace', false);
        X{i}.test{2} = X_tr_1(features, idx_test_1)';
        Label{i}.t{2} = 2 * ones(1, n_t_1);
        Label{i}.test{2} = 2 * ones(1, n_test_1);
        
        Label{i}.test_all = [Label{i}.test{1} Label{i}.test{2}];
        X{i}.test_all = [X{i}.test{1}; X{i}.test{2}];
        
        
        idx_s_0 = datasample(idx_s_0_set, n_s_0, 'Replace', false);
        X{i}.s{1} = X_sr_0(features, idx_s_0)';
        Label{i}.s{1} = ones(1, n_s_0);
        
        
        idx_s_1 = datasample(idx_s_1_set, n_s_1, 'Replace', false);
        X{i}.s{2} = X_sr_1(features, idx_s_1)';
        Label{i}.s{2} = 2 * ones(1, n_s_1);
end
        
        
        
        
        