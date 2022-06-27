function [ outsig ] = reweighted( insig, mask, param, CPparamsolver, DRparamsolver, RWparamsolver )
% general function for l1-relaxed audio inpainting with reweighting
%
% insig ........... input (degraded) signal
% mask ............ logical vector indicating the missing samples
% param
%   .type ......... switch between 'analysis' and 'synthesis'
%   .F ............ Parseval-tight Gabor frame
%   .offset ....... 'full' / 'half' / 'none'
%   .weighting .... 'none' / 'supp' / 'abs' / 'norm' / 'energy'
%   .reweighting .. do we want to use the iterative reweighting?
% CPparamsolver
%   .sigma .... parameter of prox_f*
%   .tau ...... parameter of prox_g
%   .theta .... step size
%   .x0 ....... starting point of primal variable
%   .y0 ....... starting point of dual variable
%   .maxit .... maximal number of iterations
%   .tol ...... tolerance of the relative norm of x between two iterations
% DRparamsolver
%   .lambda ... step size
%   .gamma .... parameter of the proximal operators
%   .y0 ....... starting point
%   .maxit .... maximal number of iterations
%   .tol ...... tolerance of the relative norm of x between two iterations
% RWparamsolver
%   .maxit .... maximum number of iterations of the outer cycle
%   .epsilon .. parameter of the reweighting step
%   .delta .... stopping criterion for the outer cycle
%
% Date: 08/07/2020
% By Ondrej Mokry
% Brno University of Technology
% Contact: ondrej.mokry@mensa.cz

%% input signal shortening
a = param.F.a;
M = param.F.M;
g = param.F.g; g = g{2}; w = g{2};
N = length(insig);
s = find(~mask,1,'first');
f = find(~mask,1,'last');
[q,~,~,~,~,~,~,~,U,V,L] = min_sig_supp_2(w,a,M,s,f,N,1,offset(s,f,a,param.offset));
if L < framelength(param.F,L)
    L = framelength(param.F,L);
end
if q+L-1 <= N && q >= 1
    signal = insig(q:q+L-1);
    mask = mask(q:q+L-1);
else
    q = max(q,1);
    padding = zeros(L-(N-q+1),1);
    signal = [ insig(q:end); padding ];
    mask = [ mask(q:end); true(L-(N-q+1),1) ];
end
param.F = frameaccel(param.F,L); % acceleration of snythesis and analysis

%% weighting
if strcmp(param.weighting,'none')
    wts = 1;
else
    wts = weights( param.F, L, mask, U, V, syn_atoms(param.F), param.weighting);
end

%% reconstruction
soft = @(x,gamma) sign(x) .* max(abs(x)-gamma, 0);
if ~param.reweighting
    if strcmp(param.type,'analysis')
        CPparam.f = @(x) norm(wts.*x,1);
        CPparam.g = @(x) 0; % it is actually the indicator function...
        CPparam.prox_f = @(x,gamma) soft(x,wts*gamma);
        CPparam.prox_g = @(x,gamma) proj_time(x,mask,signal.*mask);
        CPparam.dim = L;
        CPparam.K = @(x) frana(param.F,x);
        CPparam.K_adj = @(x) frsyn(param.F,x);

        % algorithm
        [ restored, ~, ~ ] = ChambollePock(CPparam,CPparamsolver);
        restored = real(restored(1:L));
    else
        DRparam.f = @(x) norm(wts.*x,1);
        DRparam.g = @(x) 0; % it is actually the indicator function...
        DRparam.prox_f = @(x,gamma) soft(x,wts*gamma);
        c = frana(param.F,signal.*mask);
        DRparam.prox_g = @(x,gamma) x - frana(param.F,mask.*frsyn(param.F,x)) + c;
        DRparam.dim = frameclength(param.F,L);

        % algorithm
        [ x_hat, ~, ~ ] = DouglasRachford(DRparam,DRparamsolver);
        x_hat = DRparam.prox_g(x_hat);
        restored = frsyn(param.F,x_hat);
        restored = real(restored(1:L));
    end
else
    if strcmp(param.type,'analysis')
        rwts = 1;
        CPparam.g = @(x) 0; % it is actually the indicator function...
        CPparam.prox_g = @(x,gamma) proj_time(x,mask,signal.*mask);
        CPparam.dim = L;
        CPparam.K = @(x) frana(param.F,x);
        CPparam.K_adj = @(x) frsyn(param.F,x);
        y_hat = zeros(CPparam.dim,1); % initial solution
        
        % the outer cycle
        for iteration = 1:RWparamsolver.maxit
            y_old = y_hat;
            CPparam.f = @(x) norm(rwts.*wts.*x,1);
            CPparam.prox_f = @(x,gamma) soft(x,rwts.*wts*gamma);
            [ y_hat, ~, ~ ] = ChambollePock(CPparam,CPparamsolver);
            if norm(y_old - y_hat) < RWparamsolver.delta
                break;
            end
            rwts = 1./(abs(frana(param.F,y_hat)) + RWparamsolver.epsilon);  
        end
        
        % solution
        restored = real(y_hat(1:L));
    else
        rwts = 1;
        DRparam.g = @(x) 0; % it is actually the indicator function...      
        c = frana(param.F,signal.*mask);
        DRparam.prox_g = @(x,gamma) x - frana(param.F,mask.*frsyn(param.F,x)) + c;
        DRparam.dim = frameclength(param.F,L);        
        x_hat = zeros(DRparam.dim,1); % initial solution
        
        % the outer cycle
        for iteration = 1:RWparamsolver.maxit
            x_old = x_hat;
            DRparam.f = @(x) norm(rwts.*wts.*x,1);
            DRparam.prox_f = @(x,gamma) soft(x,rwts.*wts*gamma);
            [ x_hat, ~, ~ ] = DouglasRachford(DRparam,DRparamsolver);
            if norm(x_old - x_hat) < RWparamsolver.delta
                break;
            end
            rwts = 1./(abs(x_hat) + RWparamsolver.epsilon);
        end
        
        % solution
        restored = frsyn(param.F,x_hat);
        restored = real(restored(1:L));
    end
end

%% putting the restored signal together
outsig = insig;
outsig(q:q+L-1) = restored;
outsig = outsig(1:N);