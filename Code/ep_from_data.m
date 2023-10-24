function [mape_p_link, estimated_p_link, estimated_b, estimated_p_node] = ep_from_data(A, error_data, X, true_p_link, conservation_result)
%  estimate error probabilities from error flow matrix and exogenous variables
[node_number, link_number] = size(A);
[~, time_interval_number] = size(error_data);
error_node_flow = A*error_data;
flow_diff = error_node_flow ~= conservation_result;
flow_diff = sum(flow_diff, 2);
estimated_p_node = flow_diff/time_interval_number;
% ape_node = abs(true_p_node - estimated_p_node);
% mape_node = mean(ape_node(:));
% disp(mape_node)
node_y = -log(1-estimated_p_node);
node_X = abs(A)*X;
%% gls
mdl = fitglm(node_X,node_y,'intercept', false);
regress(node_y, node_X)
b = mdl.Coefficients;
b  = table2array(b(:,1));
estimated_p_link = 1 - exp(-X*b);
mape_p_link = mean(abs(estimated_p_link-true_p_link)./true_p_link);
disp('mape_p_link_from_data')
disp(mape_p_link)
estimated_b = b;
end

