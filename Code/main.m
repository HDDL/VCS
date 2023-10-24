clc,clear;
close all;
%% import true_flow_matrix
true_file_name = './true_flow_matrix.txt';
true_flow_matrix = importdata(true_file_name);
%% import error_flow_matrix
error_file_name = './error_flow_matrix.txt';
error_flow_matrix = importdata(error_file_name);
%% import error p
p_file_name = './true_p.txt';
P = importdata(p_file_name);
%% import exo_var
x_file_name = './exo_var.txt';
X = importdata(x_file_name);
%% import network 
network_folder = './';
A = get_network_adjacent_matrix(network_folder);
A(14:19,:)=[];
%% test the flow-conservation
conservation_result = A*true_flow_matrix; % due to the decimal, some results may be 1, -1;
%% gurobi model
[node_number, link_number] = size(A);
[~, time_interval_number] = size(true_flow_matrix);
ae_origin = abs(error_flow_matrix-true_flow_matrix);
mae_origin = mean(ae_origin(:));
ape_origin = ae_origin./(true_flow_matrix+1); % avoid zero
mape_origin = mean(ape_origin(:));
error_link_index_origin = round(double(ae_origin~=0));
%% estimate p from data
[mape_p, estimated_p, estimated_b] = ep_from_data(A, error_flow_matrix, X, P, conservation_result);
link_conservation_flag = get_link_conservation_flag(A, error_flow_matrix, conservation_result);


%% decision variables
[M, N] = size(A);
ele_number = link_number*time_interval_number;
q_hat_number = ele_number; % recovered q
z_number = ele_number; % free variables z
gamma_number = ele_number; % occurrences of errors \gamma
var_number = q_hat_number + z_number + gamma_number; % all vars
%% coefficent of primary objective function
obj_coeff = zeros(1, var_number); 
for ii = 1:ele_number
    %obj_coeff(ii) = -1;
    link_index = floor((ii-1)/time_interval_number) + 1;
    if estimated_p(link_index) <= 0
        obj_coeff(q_hat_number + z_number + ii) = 10;
    else
        obj_coeff(q_hat_number + z_number + ii) = (log(1-estimated_p(link_index))-log(estimated_p(link_index))); % 网络中只使用了流量数据log(1-estimated_p(link_index))-log(estimated_p(link_index)) 
    end
end
%% con1: \hat_q = (1-gamma)*q + z
A1 = sparse(ele_number, var_number);
B1 = zeros(ele_number, 1); % A1对应的B1
for jj = 1:link_number
    for tt = 1:time_interval_number
        ii = (jj-1)*time_interval_number + tt;
        column_index = ii;
        column_index_z = q_hat_number + column_index;
        column_index_gamma = q_hat_number + z_number + column_index;
        A1(ii, column_index) = 1;
        A1(ii, column_index_z) = -1;
        A1(ii, column_index_gamma) = error_flow_matrix(jj, tt);
        B1(ii) = error_flow_matrix(jj, tt);
    end
end
%% con2: z \leq gamma*M
A2 = sparse(ele_number, var_number);
B2 = zeros(ele_number, 1);
BIG_NUMBER = 5000;
for ii = 1:ele_number
    column_index = ii;
    column_index_z = q_hat_number + column_index;
    column_index_gamma = q_hat_number + gamma_number + column_index;
    A2(ii, column_index_z) = 1;
    A2(ii, column_index_gamma) = -BIG_NUMBER;
end
%% con3: system constraints |\sum q+ - \sum q-| \leq \sigma
A3 = sparse(node_number*time_interval_number, var_number);
temp_result = conservation_result';
B3 = temp_result(:);
for ii = 1:node_number
    for tt = 1:time_interval_number
        row_index = (ii-1)*time_interval_number + tt;
        for jj = 1:link_number
            column_index = (jj-1)*time_interval_number + tt;
            A3(row_index, column_index) = A(ii, jj);
        end
    end
end

%% con4: \gamma=0 when data conforms flow balance law
observed_conservation_result = A*error_flow_matrix;
node_conservation_flag = observed_conservation_result == conservation_result;
A4 = sparse(link_number*time_interval_number, var_number);
B4 = ones(link_number*time_interval_number, 1);
for jj = 1:link_number
    for tt = 1:time_interval_number
        row_index = (jj-1)*time_interval_number + tt;
        column_index_gamma = q_hat_number + z_number + row_index;
        A4(row_index, column_index_gamma) = 1;
        B4(row_index, 1) = (1-link_conservation_flag(jj, tt));
    end
end


%% combination of constaints
Aeq = [A1;A3];
beq = [B1;B3];
Aineq = [A2;A4];%
bineq = [B2;B4];%

% clear -regexp ^B\d+
%% build Gurobi model
params.outputflag = 1; 
params.LogFile = 'log.lp';
%params.MIPGap = 5; % minimum errors
model.obj = obj_coeff; % primary obj
model.A = [Aeq; Aineq;]; % sparse matrix
n = size(model.A, 2); % number of vars
vtype = [repmat('C', q_hat_number, 1); repmat('C', z_number, 1); repmat('B', gamma_number, 1)];
model.vtype = vtype;
sense = [repmat('=',size(Aeq, 1), 1); repmat('<',size(Aineq, 1), 1);];
model.sense = sense;
model.rhs = [beq; bineq];
model.Params.IntegralityFocus=1;
gurobi_write(model, 'primary_model.lp');
result = gurobi(model, params);
%% build cvx model 
potential_mae = zeros(1, length(result.pool));
potential_mape = zeros(1, length(result.pool));
potential_detected_link_error_index = zeros(link_number, time_interval_number, length(result.pool));
potential_q = zeros(link_number, time_interval_number, length(result.pool));
% second obj: take the first one as an example
var_value = (result.pool(1).xn)';
estimated_q_value = var_value(1:ele_number);
z_value = var_value(ele_number+1:2*ele_number);
gamma_value = var_value(2*ele_number+1:3*ele_number);
gamma = reshape(gamma_value, time_interval_number, link_number);
gamma = gamma';
recovered_q = reshape(estimated_q_value, time_interval_number, link_number);
recovered_q = recovered_q';
Omega = 1 - gamma;
temp_error_flow_matrix = (error_flow_matrix).*Omega;
cvx_begin sdp
variable Q(link_number,time_interval_number)
subject to
    Omega.*Q == temp_error_flow_matrix;
minimize(norm_nuc(Q))
cvx_end
% our error 
ae_result = abs(Q-true_flow_matrix);
mae_result = mean(ae_result(:));
ape_result = ae_result./(true_flow_matrix+1);
mape_result = mean(ape_result(:));
% results
disp(['MAPE of estimated error probabilities:', num2str(mape_p)])
disp(['Origin MAPE:', num2str(mape_origin)]);
disp(['Residual MAPE:', num2str(mape_result)]);





