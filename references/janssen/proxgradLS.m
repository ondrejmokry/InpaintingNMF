function [ x_hat, objective, times ] = proxgradLS(prox, grad, startpoint, iterations, g, h, saveall)
% PROXGRAD Proximal gradient algorithm with backtracking linesearch. It can
% be used as the projected gradient algorithm when the input proximal
% operator is a projection.
%
% The problem to solve is
%
%      x_hat = arg min g(x) + h(x)
%
% where we know the proximal operator of h and the gradient of g.
% The algorithm is taken from
% Ryan Tibshirani, Proximal Gradient Descent(and Acceleration), p. 10
%
% Input parameters
%       prox           proximal operator of c*h(x) used as prox(x, c) where
%                      c > 0 is a constant
%       grad           gradient of g(x), assumed to be linear
%       startpoint     starting point of the algorithm
%       iterations     number of iterations of the algorithm
%       g, h           the two parts of the objective function
%       saveall        switch for saving the solution during iterations
%
% Output parameters
%       x_hat          solution
%       objective      values of the objective function during iterations
%       times          cumulative computation time during iterations
%
% Date: 21/05/2021
% By Ondrej Mokry
% Brno University of Technology
% Contact: ondrej.mokry@mensa.cz

% initialization
x      = startpoint;
t_init = 1;
beta   = 0.5;
if nargout > 1
    objective = NaN(iterations,1);
end
if nargout > 2
    times = NaN(iterations,1);
    tic
end
if saveall
    x_hat = NaN(length(startpoint),iterations);
end

% define the function Gt
Gt = @(x, t) 1/t * (x - prox(x-t*grad(x), t));

% iterate
for i = 1:iterations

    % backtracking line search
    t     = t_init;
    Gtx   = Gt(x, t);
    gradx = grad(x);
    while g(x - t*Gtx) > g(x) - t*gradx'*Gtx + t/2*norm(Gtx)^2
        t = beta*t;
        Gtx = Gt(x, t);
    end

    % perform proximal gradient update
    x = x - t*Gtx;

    % saving the objective
    if nargout > 1
        objective(i) = g(x) + h(x);
    end

    % saving the time
    if nargout > 2
        times(i) = toc;
    end

    % update the solution
    if saveall
        x_hat(:,i) = x;
    else
        x_hat = x;
    end
end

end