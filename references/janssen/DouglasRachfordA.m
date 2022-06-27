function [ x, objval ] = DouglasRachfordA(y, h, P, proxP, Nx, MAX_ITER, varargin)
% solves
% min_x 0.5 || y - H x ||_2^2 + P(x)
% where H is a Toeplitz matrix derived from the causal filter h of size
% length(y) by Nx, using the Douglas-Rachford algorithm as described in 
% 
% I. Bayram, "Proximal Mappings Involving Almost Structured Matrices,"
% in IEEE Signal Processing Letters, vol. 22, no. 12, pp. 2264-2268,
% Dec. 2015, doi: 10.1109/LSP.2015.2476381.
%
% This is an implementation of Algorithm 2 from the manuscript.
%
% Code written by Ilker Bayram, 
% ibayram@itu.edu.tr
% Istanbul Technical University, 2015
%
% Modified by Ondrej Mokry
% ondrej.mokry@mensa.cz
% Brno University of Technology, 2021
%
%   - some notation and comments
%   - generalized the proximal step
%   - added computation of the objective function

%% parse the inputs
% create the parser
p = inputParser;

% add optional name-value pairs
addParameter(p,'alpha',0.1)
addParameter(p,'lambda',0.1)

% parse
parse(p, varargin{:})

% save the parsed results to nice variables
alpha  = p.Results.alpha;
lambda = p.Results.lambda;

%% initialization
m  = length(y);
mp = length(h) + Nx - 1;
Q  = 128;
mp = Q * ceil(mp / Q);
if nargout > 1
    objval = NaN(MAX_ITER, 1);
end
u      = zeros(mp,1);
t      = zeros(mp,1);
t(1:m) = y;
H      = fft(h,mp);
beta   = alpha/(1+alpha); % beta in the manuscript
invH   = 1./(1 + beta * abs(H).^2); % inverse in the FFT domain

%% the algorithm
% wb = waitbar(0,'Algorithm 2 running');
for iter = 1 : MAX_ITER
    % waitbar(iter/MAX_ITER,wb);
    
    % update utilde
    T = fft(t);
    T = invH .* (fft(u) + beta * conj(H) .* T);
    ut = ifft(T);
    
    % update ttilde
    tt = (1-alpha) * t / (1+alpha) + 2 * beta * ifft( H .* T);
        
    % update x
    u(1:Nx) = u(1:Nx) + 2 * (1-lambda) * (proxP(2*ut(1:Nx) - u(1:Nx), alpha) - ut(1:Nx));
    
    % update z
    u(Nx+1:end) = u(Nx+1:end) - 2 * (1-lambda) * ut(Nx+1:end) ;
    
    % update t
    t(1:m) = lambda * t(1:m) + (1-lambda) * (2*y - tt(1:m));
    t(m+1:end) = lambda * t(m+1:end) + (1-lambda) * tt(m+1:end);    
    
    % compute the objective
    if nargout > 1
        z = conv(ut(1:Nx),h);
        objval(iter) = 0.5*norm(y - z(1:length(y)))^2 + P(ut(1:Nx));
    end
end
% close(wb);

% the final prox
T = fft(t);
T = invH .* (fft(u) + beta * conj(H) .* T);
ut = ifft(T);
x = ut(1:Nx);