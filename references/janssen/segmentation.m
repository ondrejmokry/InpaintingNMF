function restored = segmentation(method,signal,masks,lambda,p,maxit,varargin)
% segmentation performs segment-wise processing of audio using the Janssen
% algorithm or its variations:
% (1) inpainting using the Janssen algorithm [1] with sparse regularization
%     of the AR coefficients
% (2) declipping using the Janssen algorithm [1] with sparse regularization
%     of the AR coefficients and consistency constraints on the signal
% (3) declipping using the generalized linear prediction [2] with sparse
%     regularization of the AR coefficients
%
% implementation of the signal estimation step in (1) and (3) is taken from
% the Audio Inpainting Toolbox as described in [3]
%
% the subproblems of estimating the AR coefficients in (1)-(3) and
% estimating the signal in (2) are solved using the Douglas-Rachford
% algorithm [4], which is possibly accelerated as described in [5]
%
% [1] A. Janssen, R. Veldhuis and L. Vries, "Adaptive interpolation of
%     discrete-time signals that can be modeled as autoregressive
%     processes," in IEEE Transactions on Acoustics, Speech, and Signal
%     Processing, vol. 34, no. 2, pp. 317-330, 1986, doi:
%     10.1109/TASSP.1986.1164824.
% [2] L. Atlas and C. P. Clark, "Clipped-waveform repair in acoustic
%     signals using generalized linear prediction," US patent US8126578B2,
%     2007.
% [3] A. Adler, V. Emiya, M. G. Jafari, M. Elad, R. Gribonval and M. D.
%     Plumbley, "Audio Inpainting," in IEEE Transactions on Audio, Speech,
%     and Language Processing, vol. 20, no. 3, pp. 922-932, 2012, doi:
%     10.1109/TASL.2011.2168211.
% [4] P. Combettes and J.-C. Pesquet, "Proximal Splitting Methods in Signal
%     Processing," in Fixed-Point Algorithms for Inverse Problems in
%     Science and Engineering, vol. 49, pp. 185-212, 2009, doi:
%     10.1007/978-1-4419-9569-8_10.
% [5] I. Bayram, "Proximal Mappings Involving Almost Structured Matrices,"
%     in IEEE Signal Processing Letters, vol. 22, no. 12, pp. 2264-2268,
%     2015, doi: 10.1109/LSP.2015.2476381.
%
% input arguments
%   method        switch between the algorithms (1)-(3)
%                 (1) 'inpainting'
%                 (2) 'declipping'
%                 (3) 'glp'
%   signal        the input (degraded) signal
%   masks.R       mask of the reliable samples
%        .U       mask of the samples clipped on the upper clipping level,
%                 needed only for (2), (3)
%        .L       mask of the samples clipped on the lower clipping level,
%                 needed only for (2), (3)
%   lambda        regularization parameters:
%                 lambda(1) ... for the AR model estimation, if
%                               lambda(1) = 0, no regularization is used
%                 lambda(2) ... for the signal estimation, if lambda(2) = 0
%                               or length(lambda) = 1, the indicator
%                               function is used instead of the distance
%                               function
%   p             order of the AR model
%   maxit         number of iterations of the whole Janssen algorithm
%   varargin      name-value pairs
%                 'wtype' ('sine')     window shape for overlap-add
%                 'w' (4096)           window length for overlap-add
%                 'a' (2048)           window shift for overlap-add
%                 'coefaccel' (false)  accelerate the coef. estimation
%                 'coefextra' (false)  accelerate the signal estimation
%                 'sigaccel' (false)   extrapolate the signal
%                 'sigextra' (false)   extrapolate the coefficients
%                 'linesearch' (false) search for the coefficients and the
%                                      signal using the linesearch strategy
%                                      and the current and previous
%                                      solutions
%                 'saveall' (false)    save the solution during iterations
%                 'DRmaxit' (1000)     number of iterations for the solver
%                                      of the subproblems; if the parameter
%                                      is set to non-numeric value,
%                                      progressive strategy is applied with
%                                      100 iterations of DR in the first
%                                      iteration of Janssen and 1000
%                                      iterations of DR in the last
%                                      iteration of Janssen
%                  'mat' ('toeplitz')  how to build the matrices needed for
%                                      the subproblems, accepted values are
%                                      'toeplitz', 'xcorr', 'conv', 'fft'
%                  'decompose' (true)  use the Cholesky decomposition to
%                                      substitute for the multiplication
%                                      with matrix inversion in the
%                                      non-accelerated DR algorithm
%                  'gammaC' (0.1)      parameter of the DR algorithm for
%                                      coef. estimation
%                  'gammaS' (10)       parameter of the DR algorithm for
%                                      signal estimation
%                  'verbose' (true)    print number of the segment being
%                                      processed
%
% output arguments
%   restored      the solution; if saveall is true, restored is of size
%                 length(signal) x maxit, or length(signal) x 1 otherwise
%
% Date: 21/07/2021
% By Ondrej Mokry
% Brno University of Technology
% Contact: ondrej.mokry@mensa.cz

%% parse the inputs for overlap-add
% create the parser
pars = inputParser;
pars.KeepUnmatched = true;

% add optional name-value pairs
addParameter(pars,'wtype','sine')
addParameter(pars,'w',4096)
addParameter(pars,'a',2048)
addParameter(pars,'saveall',false)
addParameter(pars,'verbose',true)

% parse
parse(pars, varargin{:})

% save the parsed results to nice variables
wtype      = pars.Results.wtype;
w          = pars.Results.w;
a          = pars.Results.a;
saveall    = pars.Results.saveall;
verbose    = pars.Results.verbose;

%% padding the signal
% L is divisible by a and minimum amount of zeros equals gl (window length)
% zeros will be appended to avoid periodization of nonzero samples
L       = ceil(length(signal)/a)*a + (ceil(w/a)-1)*a;
S       = L/a; % number of signal segments
data    = [signal;  zeros(L-length(signal),1)];
masks.R = [masks.R; true(L-length(signal), 1)];
masks.U = [masks.U; false(L-length(signal),1)];
masks.L = [masks.L; false(L-length(signal),1)];

%% initializing the solution array
if saveall
    mrestored = zeros(w,maxit,S);
else
    mrestored = zeros(w,S);
end

%% construction of analysis and synthesis windows
g    = gabwin(wtype, a, w, L);
gana = normalize(g,'peak'); % peak-normalization of the analysis window
gana = fftshift(gana);
gsyn = gabdual(gana, a, w)*w; % computing the synthesis window

%% segment settings
mdata = NaN(w,S);
mR = false(w,S);
mL = false(w,S);
mU = false(w,S);
for s = 1:S
    % defining the indices of the current block
    indices = 1 + (s-1)*a - floor(w/2) : (s-1)*a + ceil(w/2);
    indices = 1 + mod(indices-1, L);
    
    % defining the segment data and masks
    mdata(:,s) = data(indices) .* gana;
    mR(:,s) = masks.R(indices);
    mL(:,s) = masks.L(indices);
    mU(:,s) = masks.U(indices);
end

%% segment processing via parfor
if saveall
    parfor s = 1:S
        if verbose
            fprintf('Processing segment %d of %d...\n',s,S)
        end
        smasks = struct('R',mR(:,s),'L',mL(:,s),'U',mU(:,s));   
        mrestored(:,:,s) = janssen(method,mdata(:,s),smasks,lambda,p,maxit,...
            varargin{:}); %#ok<*PFBNS>
    end
else
    parfor s = 1:S
        if verbose
            fprintf('Processing segment %d of %d...\n',s,S)
        end
        smasks = struct('R',mR(:,s),'L',mL(:,s),'U',mU(:,s));   
        mrestored(:,s) = janssen(method,mdata(:,s),smasks,lambda,p,maxit,...
            varargin{:});
    end
end

%% overlap-add
if saveall
    restored = zeros(L,maxit);
else
    restored = zeros(L,1);
end
for s = 1:S
    indices = 1 + (s-1)*a - floor(w/2) : (s-1)*a + ceil(w/2);
    indices = 1 + mod(indices-1, L);
    if saveall
        restored(indices,:) = restored(indices,:) + mrestored(:,:,s).*repmat(gsyn,1,maxit);
    else
        restored(indices) = restored(indices) + mrestored(:,s).*gsyn;
    end
end

%% cropping the solution to the original length
restored = restored(1:length(signal),:);

end