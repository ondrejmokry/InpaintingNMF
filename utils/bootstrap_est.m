function [ means, lower, upper ] = bootstrap_est(data, alpha, draws)
% BOOTSTRAP_EST computes the interval estimation for the mean of given data
% computed as the bootstrap-t interval [p. 160, 1].
%
% It operates along the first dimension of the data.
%
% [1] Efron, B. & Tibshirani, R. J.
%     An Introduction to the Bootstrap 
%     Chapmann & Hall, 1994
%
% input arguments
%   data      the data
%   alpha     significance level
%   draws     number of draws for the bootstrap estimate
%
% output arguments
%   means     mean values of the data along the first dimension
%   lower     lower (1-alpha) interval estimate for the mean of the data
%   upper     upper (1-alpha) interval estimate for the mean of the data
%
% Date: 03/09/2021
% By Ondrej Mokry
% Brno University of Technology
% Contact: ondrej.mokry@mensa.cz

if nargin < 3
  draws = 1000;
end

if nargin < 2
  alpha = 0.05;
end

% reshape the data to a 2D array
dims   = size(data);
data2D = reshape(data, dims(1), []);   
[m, n] = size(data2D);

% initialize
t      = NaN(draws, 1); % the stuff computed for each draw
lower  = NaN(n, 1);
upper  = NaN(n, 1);

% simple mean and variance
means = mean(data2D, 1, 'omitnan');
vars  = var(data2D, [], 1, 'omitnan');

% loop over the 2nd dimension of data2D
for i = 1:n

    % bootstrap upper and lower bound
    s = sqrt(vars(i));
    for j = 1:draws
        draw = randsample(data2D(:,i),m,true);
        t(j) = (mean(draw,'omitnan') - means(i))/sqrt(var(draw,'omitnan'))*sqrt(m-1);
    end
    t_lower = quantile(t, 1-alpha/2);
    t_upper = quantile(t, alpha/2);
    lower(i) = means(i) - t_lower*s/sqrt(m-1);
    upper(i) = means(i) - t_upper*s/sqrt(m-1);

end

% reshape to correspond to the original shape (if necessary)
if length(dims) > 2
    means = reshape(means, dims(2:end));
    lower = reshape(lower, dims(2:end));
    upper = reshape(upper, dims(2:end));
end

end

