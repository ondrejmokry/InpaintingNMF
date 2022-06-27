function [restored, W, H, relnorms, objectives] = ainmf(method,signal,mask,K,maxit,varargin)
% AINMF wraps three closely related reconstruction algorithms for audio
% inpainting based on nonnegative matrix factorization (NMF):
% (EM-tf) the approach from [1], i.e. the generalized EM algorithm where
%         the complete data is the spectrogram of the clean signal, and the
%         individual time-frequency coefficients are Gaussian, independent
%         and with variances corresponding to a matrix suitable for NMF
% (EM-t)  same as (EM-tf), except for the complete data being defined as
%         the clean time-domain frame*
% (AM)    alternating minimization approach to solve the maximum likelihood
%         problem
%
% input arguments
%   method        switch between the algorithms 'EM-tf', 'EM-t', 'AM', or
%                 set an integer < maxit to identify the switching mode
%                 from 'AM' to 'EM-tf'
%   signal        the input (degraded) signal
%   mask          mask of the reliable samples
%   K             number of components for NMF
%                 (20 in the declippind survey [2])
%   maxit         number of iterations of the whole algorithm
%                 (50 in declipping survey [2])
%   varargin      name-value pairs
%                 'nmfit'   (1)     number of inner iterations of NMF
%                 'saveall' (false) save the solution during iterations
%                 'M'       (2048)  length of the window
%                 'a'       (1024)  window shift
%                 'F'       (2048)  number of frequency channels
%                 'T'       ([])    the synthesis operator (applied to each
%                                   time frame), it can be either anonymous
%                                   function or matrix
%                 'U'       ([])    the analysis operator (applied to each
%                                   time frame), it can be either anonymous
%                                   function or matrix
%                 'verbose' (false) controls the outputs to the command
%                                   window
%                 'drawing' (false) controls plotting additional data
%                                   during iterations (Itakura-Saito
%                                   divergence progression during the
%                                   multiplicative updates, current
%                                   solution etc.)
%                 'likelihood' ([]) allows to change the objective
%                                   computation, choices are 'full' or
%                                   'observation'
%                 'epsilon' (1e-4)  additive constant for NMF to avoid
%                                   zeros
%                 'tolsol'  (0)     stopping value for the relative
%                                   norm of the solution update
%                 'tolobj'  (0)     stopping value for the relative
%                                   objective change
%
% output arguments
%   restored      the solution; if saveall is true, restored is of size
%                 length(signal) x maxit, or length(signal) x 1 otherwise
%   W, H          the resulting nonnegative matrix fatorization
%   relnorms      relative norms of the solution updates
%   objective     values of the objective function during iterations
%
%
% [1] C. Bilen, A. Ozerov, and P. Perez, "Audio declipping via nonnegative
%     matrix factorization," in 2015 IEEE Workshop on Applications of
%     Signal Processing to Audio and Acoustics (WASPAA), 2015.
% [2] P. Zaviska, P. Rajmic, A. Ozerov, and L. Rencker, "A survey and
%     an extensive evaluation of popular audio declipping methods," in IEEE
%     Journal of Selected Topics in Signal Processing, 2021.
%
%
% Date: 27/06/2022
% By Ondrej Mokry
% Brno University of Technology
% Contact: ondrej.mokry@mensa.cz
%
% _______________________
% * note that (EM-tf) and (EM-t) are equivalent if analysis after synthesis
%   in each signal frame equlas identity

%% parsing the inputs
% create the parser
pars = inputParser;
pars.KeepUnmatched = true;

% add optional name-value pairs
addParameter(pars,'nmfit',1)
addParameter(pars,'saveall',false)
addParameter(pars,'M',2048)
addParameter(pars,'a',1024)
addParameter(pars,'F',2048)
addParameter(pars,'T',[])
addParameter(pars,'U',[])
addParameter(pars,'verbose',false)
addParameter(pars,'drawing',false)
addParameter(pars,'likelihood',[])
addParameter(pars,'epsilon',1e-4)
addParameter(pars,'tolsol',0)
addParameter(pars,'tolobj',0)

% parse
parse(pars, varargin{:})

% save the parsed results to nice variables
nmfit      = pars.Results.nmfit;
saveall    = pars.Results.saveall;
M          = pars.Results.M;
a          = pars.Results.a;
F          = pars.Results.F;
T          = pars.Results.T;
U          = pars.Results.U;
verbose    = pars.Results.verbose;
drawing    = pars.Results.drawing;
likelihood = pars.Results.likelihood;
epsilon    = pars.Results.epsilon;
tolsol     = pars.Results.tolsol;
tolobj     = pars.Results.tolobj;

% handle the method parameter
if isnumeric(method)
    switchit = method;
else
    switchit = [];
end

if drawing
    % initialize the figure
    fig = figure;
    tiles = tiledlayout('flow');

    % initialize the tiles
    if nmfit > 1
        tilenmf = nexttile;
    end
    if nargout > 3 || tolsol > 0
        tilesol = nexttile;
    end
    if nargout > 4 || tolobj > 0
        tileobj = nexttile;
    end
    tileV = nexttile;
end

%% constructing the operators T and U
% if T and U are not set, we use the (potentionally oversampled) FFT
if isempty(T) && isempty(U)
    U = @(x) fft(x,F)/sqrt(F);
    T = @(x) crop(ifft(x),M)*sqrt(F);
end

% if T or U are anonymous functions, convert them to matrices
if isa(T, 'function_handle')
    T = T(eye(F));
end
if isa(U, 'function_handle')
    U = U(eye(M));
end

% precompute the inversion of T, if it exists
if M == F && rank(T) == M
    invT = inv(T);
    invertible = true;
else
    invT = T';
    invertible = false;
end

%% initialization
% signal padding
L    = ceil(length(signal)/a)*a; % + (ceil(M/a)-1)*a;
N    = L/a; % number of signal segments
data = [signal; zeros(L-length(signal),1)];
mask = [mask; true(L-length(signal),1)]; % logical mask of the reliable samples

% random initialization of the parameters W, H
rng(0)
W = 0.5*(1.5*abs(randn(F,K))+0.5); % taken from betaNTF
H = 0.5*(1.5*abs(randn(K,N))+0.5); % taken from betaNTF except for the shape

% initialize the solution matrix
if saveall
    mrestored_all = zeros(M,maxit,N);
else
    mrestored = zeros(M,N);
end

% construct the analysis and synthesis windows
g    = gabwin('sine', a, M, L);
gana = normalize(g,'peak'); % peak-normalization of the analysis window
gana = fftshift(gana);
gsyn = gabdual(gana, a, M)*M; % computing the synthesis window

% segment settings
mdata = NaN(M,N);
mmask = false(M,N);
for n = 1:N
    % defining the indices of the current block
    indices = 1 + (n-1)*a - floor(M/2) : (n-1)*a + ceil(M/2);
    indices = 1 + mod(indices-1, L);
    
    % defining the block data and masks
    mdata(:,n) = data(indices) .* gana;
    mmask(:,n) = mask(indices);
end

% saving some time by precomputing U*T for EM-t
UT = U*T;

% if desired, initialize the relative norm and objective values
if nargout > 3 || tolsol > 0
    relnorms    = NaN(maxit,1);
    presolution = NaN(L,1);
end
if nargout > 4 || tolobj > 0
    objectives = NaN(maxit,1);
end

% prepare the low-rank representation
V = W*H + epsilon;

%% main iterations
if verbose
    str = [];
end
for i = 1:maxit
    
    if ~isempty(switchit)
        if i <= switchit
            method = 'AM';
        else
            method = 'EM-tf';
        end
    end

    if verbose
        fprintf(repmat('\b',1,length(str)))
        str = sprintf('%-4s Iteration %d of %d.\n', [method,':'], i, maxit);
        fprintf('%s', str)
    end

    % prepare the matrix for the power spectrum
    P = NaN(F,N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                 e-step (EM-tf, EM-t) / signal update (AM)                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if strcmpi(method,'EM-tf')
        parfor n = 1:N
            if sum(mmask(:,n)) == M && invertible
                s         = invT*mdata(:,n);
                diagSigma = 0;
            else
                MT        = T(mmask(:,n),:); %#ok<*PFBNS>
                DTM       = MT'.*V(:,n);
                MTDTM     = MT*DTM;
                DTMMTDTM  = DTM/MTDTM;
                s         = DTMMTDTM*mdata(mmask(:,n),n);
                diagSigma = V(:,n) - sum(DTMMTDTM .* conj(DTM), 2);
            end
            
            % saving the current solution
            if saveall
                mrestored_all(:,i,n) = real(T*s); 
            else
                mrestored(:,n) = real(T*s);
            end
            
            % posterior power spectrum
            P(:,n) = abs(s).^2 + abs(diagSigma);
        end
    end
    if strcmpi(method,'EM-t')
        parfor n = 1:N
            if sum(mmask(:,n)) == M
                s         = U*mdata(:,n);
                diagSigma = 0;
            else
                MT        = T(mmask(:,n),:); %#ok<*PFBNS>
                DTM       = MT'.*V(:,n);
                MTDTM     = MT*DTM;
                DTMMTDTM  = DTM/MTDTM;
                s         = UT*(DTMMTDTM*mdata(mmask(:,n),n));
                diagSigma = sum((UT*(diag(V(:,n)) - DTMMTDTM*DTM')) .* conj(UT), 2);
            end
            
            % saving the current solution
            if saveall
                mrestored_all(:,i,n) = real(T*s); 
            else
                mrestored(:,n) = real(T*s);
            end
            
            % posterior power spectrum
            P(:,n) = abs(s).^2 + abs(diagSigma);
        end
    end
    if strcmpi(method,'am')
        parfor n = 1:N
            % signal update
            if sum(mmask(:,n)) == M
                xhat = mdata(:,n);
                s    = invT*xhat;
            else
                MT    = T(mmask(:,n),:); %#ok<*PFBNS>
                DTM   = MT'.*V(:,n);
                MTDTM = MT*DTM;
                MTDTMx = MTDTM\mdata(mmask(:,n),n);
                s     = DTM*MTDTMx; % useful later for power spectrum
                xhat  = T*s;
            end

            % saving the current solution
            if saveall
                mrestored_all(:,i,n) = real(xhat);
            else
                mrestored(:,n) = real(xhat);
            end

            % power spectrum (for the model update)
            P(:,n) = abs(s).^2;
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                  m-step (EM-tf, EM-t) / model update (AM)                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if drawing
        dis = zeros(nmfit,1);
    end

    % iterate
    for j = 1:nmfit
        
        % multiplicative updates
        W = W .* ((V.^(-2) .* P) * H') ./ (V.^(-1) * H');
        V = W*H + epsilon;
        H = H .* (W' * (V.^(-2) .* P)) ./ (W' * V.^(-1));
        V = W*H + epsilon;
        
        % normalization
        scale = sum(W,1);
        W = W ./ scale;
        H = H .* scale';

        % computing the Itakura-Saito divergence
        if drawing
            dis(j) = sum(P./V - log(P./V) - 1,'all');
        end
        
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                          relative norm update                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargout > 3 || tolsol > 0
        
        % read only the current matrix solution
        if saveall
            mrestored = squeeze(mrestored_all(:,i,:));
        end
        
        % overlap-add
        postsolution = zeros(L,1);
        for n = 1:N
            indices = 1 + (n-1)*a - floor(M/2) : (n-1)*a + ceil(M/2);
            indices = 1 + mod(indices-1, L);
            postsolution(indices) = postsolution(indices) + mrestored(:,n).*gsyn;
        end
        
        % compute the relative norm
        presolution = presolution(1:length(signal));
        postsolution = postsolution(1:length(signal));
        relnorms(i) = norm(postsolution - presolution)/norm(presolution);
        
        % update presolution for the next iteration
        presolution = postsolution;
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                            objective update                             %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if nargout > 4 || tolobj > 0
        if isempty(likelihood)
            if strcmpi(method,'am')
                likelihood = 'full';
            else
                likelihood = 'observation';
            end
        end
        if strcmpi(likelihood,'observation')
            objectives(i) = objectiveObservation(V,T,mmask,mdata);
        else
            if saveall
                mrestored = squeeze(mrestored_all(:,i,:));
            end
            objectives(i) = objectiveFull(V,T,mmask,mrestored);
        end
        
        % stop execution if objective computation fails
        if isnan(objectives(i))
            break
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                  plot                                   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if drawing
        if nmfit > 1
            semilogy(tilenmf,dis)
            xlabel(tilenmf,'iteration of NMF')
            ylabel(tilenmf,'DIS(P,V)')
            title(tilenmf,'Itakura-Saito divergence')
        end
        if nargout > 3 || tolsol > 0
            semilogy(tilesol,relnorms)
            xlabel(tilesol,['iteration of ', method])
            ylabel(tilesol,'relative solution change')
            title(tilesol,'relative solution change')
        end
        if nargout > 4 || tolobj > 0
            plot(tileobj,objectives)
            xlabel(tileobj,['iteration of ', method])
            ylabel(tileobj,'objective function')
            title(tileobj,'objective function')
        end
        imagesc(tileV,log10(abs(V)))
        set(gca,'YDir','normal')
        colorbar
        xlabel(tileV,'time frame')
        ylabel(tileV,'frequency bin')
        title(tileV,'magnitude of V = W*H + epsilon')
        title(tiles,sprintf('%s iteration %d',method,i))
        drawnow
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                       check the stopping criteria                       %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % relative norm of the solution update
    if tolsol > 0 && i > 1
        if relnorms(i) < tolsol
            break
        end        
    end

    % relative objective change
    if tolobj > 0 && i > 1
        % check for decreasing objective (due to numerical errors or
        % non-default objective computation setting)
        if objectives(i) < objectives(i-1)
            % then check the relative change
            if abs(objectives(i)-objectives(i-1))/abs(objectives(i-1)) < tolobj
                break
            end
        end
    end
end

%% overlap-add
if saveall
    restored = zeros(L,maxit);
    for n = 1:N
        indices = 1 + (n-1)*a - floor(M/2) : (n-1)*a + ceil(M/2);
        indices = 1 + mod(indices-1, L);
        restored(indices,:) = restored(indices,:) + mrestored_all(:,:,n).*repmat(gsyn,1,maxit);
    end
else
    restored = zeros(L,1);
    for n = 1:N
        indices = 1 + (n-1)*a - floor(M/2) : (n-1)*a + ceil(M/2);
        indices = 1 + mod(indices-1, L);
        restored(indices) = restored(indices) + mrestored(:,n).*gsyn;
    end
end

%% output
% crop the solution back to the original length
restored = crop(restored,length(signal));

% if the iteration stopped before maxit, crop also along the iteration
% dimension
if saveall
    restored = restored(:,1:i);
end

% crop the convergence profiles
if nargout > 3
    relnorms = relnorms(1:i);
end
if nargout > 4
    objectives = objectives(1:i);
end

% close the figure showing Itakura-Saito divergence
if drawing
    close(fig)
end

end

%% functions
% the objective function computed only from the observed samples
% (default in EM-tf or EM-t)
function val = objectiveObservation(V,T,mmask,mdata)
    N = size(V,2);
    val = NaN(N,1);
    parfor n = 1:N
        MT      = T(mmask(:,n),:);
        DTM     = MT'.*V(:,n);
        MTDTM   = MT*DTM;
        MTDTM   = (MTDTM + MTDTM')/2; % ensure real diagonal
        [R, er] = chol(MTDTM); % useful for stable determinant computation
        if ~er
            val(n) = sum(mmask(:,n)) * log(pi)...
                 + 2 * sum(log(diag(R)))...
                 + mdata(mmask(:,n),n)'*(R\(R'\mdata(mmask(:,n),n)));      
        end
    end
    val = real(sum(val));
end

% the objective function computed from the whole recovered signal
% (default in AM)
function val = objectiveFull(V,T,mmask,mrestored)
    N = size(V,2);
    val = NaN(N,1);
    parfor n = 1:N
        xhat    = mrestored(:,n);
        TDT     = T*(T'.*V(:,n));
        TDT     = (TDT + TDT')/2; % ensure real diagonal
        [R, er] = chol(TDT); % useful for stable determinant computation
        if ~er
            val(n) = length(mmask(:,n)) * log(pi)...
                 + 2 * sum(log(diag(R)))...
                 + xhat'*(R\(R'\xhat));
        end
    end
    val = real(sum(val));
end