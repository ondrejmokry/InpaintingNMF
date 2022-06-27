function [ x_hat, obj_val, rel_norm ] = DouglasRachford(param,paramsolver)
% Douglas-Rachford algorithm for solving
%
%      x_hat = arg min f(x) + g(x)
%
% algorithm taken from
%     Combettes, Patrick & Pesquet, Jean-Christophe. (2009). Proximal
%     Splitting Methods in Signal Processing. Fixed-Point Algorithms for
%     Inverse Problems in Science and Engineering. 49.
%     10.1007/978-1-4419-9569-8_10.
%
% inputs:
%      param
%         .f ........ function f(x)
%         .g ........ function g(x)
%         .prox_f ... proximal operator of f(x), used as
%                     prox_f(arg, parameter)
%         .prox_g ... proximal operator of g(x), used as
%                     prox_g(arg, parameter)
%         .dim ...... length of x
%      paramsolver
%         .lambda ... step size
%         .gamma .... parameter of the proximal operators
%         .y0 ....... starting point
%         .maxit .... maximal number of iterations
%         .tol ...... tolerance of the relative norm of x between two
%                     iterations
%
% outputs:
%     x_hat ......... solution
%     obj_val ....... value of f(x) + g(x) during the iterations
%     rel_norm ...... relative norm of x wrt previous iteration
%
% Date: 06/01/2021
% By Ondrej Mokry
% Brno University of Technology
% Contact: ondrej.mokry@mensa.cz

%% default settings
if ~isfield(paramsolver,'lambda')
    paramsolver.lambda = 1;
end
if ~isfield(paramsolver,'gamma')
    paramsolver.gamma = 1;
end
if ~isfield(paramsolver,'y0')
    paramsolver.y0 = zeros(param.dim,1);
end
if ~isfield(paramsolver,'maxit')
    paramsolver.maxit = 500;
end
if ~isfield(paramsolver,'tol')
    paramsolver.tol = 5e-4;
end

%% initialization
obj_val  = Inf(paramsolver.maxit,1);
rel_norm = Inf(paramsolver.maxit,1);
x_old    = paramsolver.y0;
y        = paramsolver.y0;

%% algorithm
i = 1;
while i <= paramsolver.maxit && rel_norm(i) > paramsolver.tol
    x_new = param.prox_g(y, paramsolver.gamma);
    y     = y + paramsolver.lambda*(param.prox_f(2*x_new - y, paramsolver.gamma) - x_new);
    i     = i + 1;
    if nargout > 1
        obj_val(i)  = param.f(x_new) + param.g(x_new);
        rel_norm(i) = norm(x_new-x_old)/norm(x_old);
    end
    x_old = x_new;
end

%% output
x_hat = x_new;
if nargout > 1
    obj_val  = obj_val(1:i-1);
    rel_norm = rel_norm(1:i-1);
end