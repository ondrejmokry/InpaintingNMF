function [ x_hat, obj_val, rel_norm ] = ChambollePock(param,paramsolver)
% Chambolle-Pock algorithm for solving
%
%      x_hat = arg min f(Kx) + g(x)
%
% algorithm taken from
%     Malitsky, Yura & Pock, Thomas. (2016). A First-Order Primal-Dual
%     Algorithm with Linesearch. SIAM Journal on Optimization. 28.
%     10.1137/16M1092015. 
%
% inputs:
%      param
%         .f ........ function f(x)
%         .g ........ function g(x)
%         .K ........ linear function K(x)
%         .K_adj .... linear function K*(x)
%         .prox_f ... proximal operator of f(x), used as
%                     prox_f(arg, parameter)
%         .prox_g ... proximal operator of g(x), used as
%                     prox_g(arg, parameter)
%         .dim ...... length of x
%      paramsolver
%         .sigma .... parameter of prox_f*
%         .tau ...... parameter of proc_g
%         .theta .... step size
%         .x0 ....... starting point of primal variable
%         .y0 ....... starting point of dual variable
%         .maxit .... maximal number of iterations
%         .tol ...... tolerance of the relative norm of x between two
%                     iterations
%
% outputs:
%     x_hat ......... solution
%     obj_val ....... value of f(Kx) + g(x) during the iterations
%     rel_norm ...... relative norm of x wrt previous iteration
%
% Date: 08/07/2020
% By Ondrej Mokry
% Brno University of Technology
% Contact: ondrej.mokry@mensa.cz

%% default settings
if ~isfield(paramsolver,'sigma')
    paramsolver.sigma = 5;
end
if ~isfield(paramsolver,'tau')
    paramsolver.tau = 0.2;
end
if ~isfield(paramsolver,'theta')
    paramsolver.theta = 1;
end
if ~isfield(paramsolver,'x0')
    paramsolver.x0 = zeros(param.dim,1);
end
if ~isfield(paramsolver,'y0')
    paramsolver.y0 = zeros(size(param.K(zeros(param.dim,1))));
end
if ~isfield(paramsolver,'maxit')
    paramsolver.maxit = 500;
end
if ~isfield(paramsolver,'tol')
    paramsolver.tol = 5e-4;
end

%% inicialization
obj_val = Inf(paramsolver.maxit,1);
rel_norm = Inf(paramsolver.maxit,1);
x_old = paramsolver.x0;
x_line = paramsolver.x0;
y = paramsolver.y0;

%% algorithm
i = 1;
while i <= paramsolver.maxit && (rel_norm(i) > paramsolver.tol || i < 10)
    arg = y + paramsolver.sigma*param.K(x_line);
    y = arg - paramsolver.sigma*param.prox_f(arg/paramsolver.sigma, 1/paramsolver.sigma);
    x_new = param.prox_g(x_old - paramsolver.tau*param.K_adj(y), paramsolver.tau);
    x_line = x_new + paramsolver.theta*(x_new-x_old);    
    i = i + 1;
    obj_val(i) = param.f(param.K(x_new)) + param.g(x_new);
    rel_norm(i) = norm(x_new-x_old)/norm(x_old);
    x_old = x_new;
end

%% output
x_hat = x_new;
obj_val = obj_val(1:i-1);
rel_norm = rel_norm(1:i-1);