function [model_pars, n_dat, cov_matrix] = fit_model_reg(vel, varargin)
% This function fits the given model to the existing data

n_dat = numel(vel); % number of velocity data
model_mat = [varargin{:}];  % model matrix
n_pars = size(model_mat, 2); % number of model parameters
% handle rank deficient model matrices
model_pars = nan(n_pars, 1);
cov_matrix = nan(n_pars, n_pars);

if rank(model_mat) < n_pars
    return
end
% fit model tot data
cn = cond(model_mat);
p=0;
if cn>100
    disp(strcat('n_data = ', num2str(n_dat), 'while cond(A) = ', num2str(cn)))
    p = 2;
end
G = p*eye(n_pars);
model_pars = (model_mat'*model_mat + G'*G)\model_mat'*vel;
end