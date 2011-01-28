function conditional_x_a_given_x_b = condition(g,indexes, vals)
%GAUSSIAN Gaussian conditional distribution returns the Gaussian
%distribution resulting from conditioning on a subset of the variables the
%subset of variables whose value is given is indicated by the first
%argument which is a logical array of the same dimension as the mean of the
%original gaussian with 1's for the variables not specified and 0's for the
%variables that are specified.  The second argument is an array of
%conditioning values with number of elements equal to sum(indexes==0)

joint_precision = inv(g.c);
joint_mean = g.m;

x_a_indexes = logical(indexes);
x_b_indexes = ~x_a_indexes;

x_b = vals;

mean_x_a_given_x_b = joint_mean(x_a_indexes) - inv(joint_precision(x_a_indexes,x_a_indexes))*joint_precision(x_a_indexes,x_b_indexes)*(x_b-joint_mean(x_b_indexes));
covariance_x_a_given_x_b = inv(joint_precision(x_a_indexes,x_a_indexes));

conditional_x_a_given_x_b = gaussian(mean_x_a_given_x_b,covariance_x_a_given_x_b);

% subplot(2,1,2)
% plot(conditional_x_a_given_x_b);

