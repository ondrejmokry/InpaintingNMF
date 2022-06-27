function [ outsig ] = tdc( insig, mask, param, CPparamsolver, DRparamsolver, TDCparamsolver )
% function for direct time-domain compensation of energy decrease inside
% the filled gap
%
% inputs: 
% insig ........... input (degraded) signal
% mask ............ logical vector indicating the missing samples, note
%                   that only single gap needs to be indicated
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
% TDCparamsolver
%   .segs ..... number of segments from which the energy is computed
%   .lens ..... length of the segments from which the energy is computed
%   .gaps ..... number of additional gaps
%   .shift .... shift between adjacent artificial gaps
%
% Date: 08/07/2020
% By Ondrej Mokry
% Brno University of Technology
% Contact: ondrej.mokry@mensa.cz

s = find(~mask,1,'first'); % first sample in the gap
f = find(~mask,1,'last'); % last sample in the gap
h = f - s + 1; % number of missing samples (in case of only one gap)
L = length(insig);
g = param.F.g; g = g{2}; w = g{2};

%% computing default parameters
if ~isfield(TDCparamsolver,'segs')
    TDCparamsolver.segs = 10;
end
if ~isfield(TDCparamsolver,'gaps')
    TDCparamsolver.gaps = 4;
end
if ~isfield(TDCparamsolver,'shift')
    TDCparamsolver.shift = w/2;
end

%% setting the first samples of additional gaps
starts = zeros(1,TDCparamsolver.gaps);
for i = 1:ceil(TDCparamsolver.gaps/2)
    starts(i) = s - w + 1 - h - TDCparamsolver.shift*(i-1);
end
for i = ceil(TDCparamsolver.gaps/2)+1:TDCparamsolver.gaps
    starts(i) = f + w + TDCparamsolver.shift*(i-ceil(TDCparamsolver.gaps/2)-1);
end
ends = starts + h - 1;

%% initialization of the arrays for computing energy curves
E_orig = zeros(TDCparamsolver.segs,TDCparamsolver.gaps); % energy of segments in original signal
E_reco = E_orig; % energy of segments in inpainted signal
    
%% segment parameters
seg_len = max(TDCparamsolver.lens, h/TDCparamsolver.segs);
seg_shift = (h - seg_len)/(TDCparamsolver.segs - 1);
    
for i = 1:TDCparamsolver.gaps
    %% reconstruction of an additional hole
    iter_mask = true(L,1); iter_mask(starts(i):ends(i)) = false;
    iter_gapped = insig.*iter_mask;
    iter_restored = reweighted( iter_gapped, iter_mask, param, CPparamsolver, DRparamsolver, [] );
        
    %% computation of energy
    for j = 1:TDCparamsolver.segs
        E_orig(j,i) = norm(...
            insig        (starts(i)+round((j-1)*seg_shift):starts(i)-1+round((j-1)*seg_shift+seg_len))...
            )^2;
        E_reco(j,i) = norm(...
            iter_restored(starts(i)+round((j-1)*seg_shift):starts(i)-1+round((j-1)*seg_shift+seg_len))...
            )^2;
    end
end

%% least squares computation of the nodes of weighting curve
m = ones(TDCparamsolver.segs,1);
for i = 1:TDCparamsolver.segs
    m(i) = E_orig(i,:)*E_reco(i,:)'/norm(E_reco(i,:))^2;
end
n = sqrt(m);

%% computation of the weighting curve
x = linspace(seg_len/2,h-seg_len/2,TDCparamsolver.segs);
y = transpose(n);

% ensuring that the first and last index of the rates correspond to the
% first and last index of the gap and the rates are equal 1 there

x = [ 1 x h ];
y = [ 1 y 1 ];

I = logical(x);
for i = 2:length(x)
    if x(i) == x(i-1)
        I(i) = 0;
    end
end
x = x(I); y = y(I);

% ensuring the supposed shape
y = (y + flip(y))/2; % symmetrization
for i = 2:ceil(length(y)/2) % y increasing in the first half...
    y(i) = max(y(i-1),y(i));
end
y(ceil(length(y)/2) + 1:end) = y(floor(length(y)/2):-1:1); % ...and decreasing in the second half

% spline computation, derivative in 1 and h is zero
cs = spline(x,[0 y 0]);

% weighting vector for the time domain
q = ppval(cs,1:h);

%% updating the reconstructed signal
outsig = insig;
outsig(s:f) = outsig(s:f).*q';

