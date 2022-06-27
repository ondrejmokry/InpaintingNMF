function [restored, objective, times] = janssen(method,signal,masks,lambda,p,maxit,varargin)
% janssen is a function that wraps three closely related reconstruction
% algorithms for audio declipping/inpainting:
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
%
% output arguments
%   restored      the solution; if saveall is true, restored is of size
%                 length(signal) x maxit, or length(signal) x 1 otherwise
%   objective     values of the objective function during iterations
%   times         cumulative computation time during iterations
%
% Date: 02/07/2021
% By Ondrej Mokry
% Brno University of Technology
% Contact: ondrej.mokry@mensa.cz

%% parse the inputs
% create the parser
pars = inputParser;
pars.KeepUnmatched = true;

% add optional name-value pairs
addParameter(pars,'coefaccel',false)
addParameter(pars,'coefextra',false)
addParameter(pars,'sigaccel',false)
addParameter(pars,'sigextra',false)
addParameter(pars,'linesearch',false)
addParameter(pars,'saveall',false)
addParameter(pars,'DRmaxit',1000)
addParameter(pars,'mat','toeplitz')
addParameter(pars,'decompose',true)
addParameter(pars,'gammaC',0.1)
addParameter(pars,'gammaS',10)
addParameter(pars,'plotLS',false)

% parse
parse(pars, varargin{:})

% save the parsed results to nice variables
coefaccel  = pars.Results.coefaccel;
coefextra  = pars.Results.coefextra;
sigaccel   = pars.Results.sigaccel;
sigextra   = pars.Results.sigextra;
linesearch = pars.Results.linesearch;
saveall    = pars.Results.saveall;
DRmaxit    = pars.Results.DRmaxit;
mat        = pars.Results.mat;
decompose  = pars.Results.decompose;
gammaC     = pars.Results.gammaC;
gammaS     = pars.Results.gammaS;
plotLS     = pars.Results.plotLS;

% update plotLS
plotLS = linesearch && plotLS;

%% initialization
solution = signal;
N = length(signal);
if saveall
    restored = NaN(N,maxit);
end
oldcoef = zeros(p+1,1); % auxiliary variable for the extrapolation
oldsolution = zeros(N,1); % auxiliary variable for the extrapolation

% define some useful functions
soft = @(x,t) sign(x).*max(abs(x)-t,0);
proj = @(x) signal .* masks.R + ...
            max(x, signal) .* masks.U + ...
            min(x, signal) .* masks.L;
if nargout > 1 || linesearch
    if length(lambda) > 1 && strcmpi(method,'declipping') && lambda(2) < Inf
        Q = @(x,c) 0.5*norm(fft(c,N+p).*fft(x,N+p))^2 / (N+p)...
            + lambda(1)*norm(c,1)...
            + lambda(2)*0.5*norm(x-proj(x))^2;
    else
        Q = @(x,c) 0.5*norm(fft(c,N+p).*fft(x,N+p))^2 / (N+p)...
            + lambda(1)*norm(c,1);
    end
    objective = NaN(maxit,1);
end

% if desired, compute the progression of the Douglas-Rachford iterations
if isfloat(DRmaxit)
    if length(DRmaxit) == 1
        maxits = DRmaxit*ones(maxit,1);
    else
        maxits = DRmaxit;
    end
else
    maxits = round(logspace(2,3,maxit));
end

% prepare some matrices for the inpainting / glp case
if strcmpi(method,'inpainting') || strcmpi(method,'glp')
    indmiss  = find(~masks.R);
    indobs   = find(masks.R);
    IAA      = abs(repmat(indmiss,1,N)-repmat(1:N,length(indmiss),1));
    IAA1     = IAA <= p;
end

% if desired, initialize the figures for line search plotting
if plotLS
    fig1 = figure; % objective function
    fig2 = figure; % second derivative
    [rows, cols] = sbplts(maxit-2);
end

% if desired, start the timer
if nargout > 2
    times = NaN(maxit,1);
    tic
end

%% main iteration
for i = 1:maxit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                           AR model estimation                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if lambda(1) == 0        
        coef = lpc(solution,p)'; 
    else
        % prepare the matrix XX
        if strcmpi(mat,'toeplitz')
            X  = toeplitz([solution; zeros(p,1)], [solution(1), zeros(1,p)]);
            XX = X'*X;
            clear X
        elseif strcmpi(mat,'xcorr')
            x  = xcorr(solution,p)';
            XX = spdiags(ones(p+1,1)*x,-p:p,p+1,p+1);
        elseif strcmpi(mat,'conv')
            x  = conv(solution,flip(solution))';
            XX = spdiags(ones(p+1,1)*x(N-p:N+p),-p:p,p+1,p+1);
        elseif strcmpi(mat,'fft')
            x  = ifft(fft([solution' zeros(1,p)]) .* fft([flip(solution)' zeros(1,p)]));
            XX = spdiags(ones(p+1,1)*x(N-p:N+p),-p:p,p+1,p+1);
        end
        
        % check the positive definiteness using chol
        try
            dXX = decomposition(eye(p+1) + gammaC*XX, 'chol');
        catch
            break
        end
        
        % check the conditioning
        if isIllConditioned(dXX)
            break
        end
        
        if coefaccel
            % set the parameters of the accelerated Douglas-Rachford
            % algorithm
            P     = @(a) lambda(1)*norm(a, 1);
            proxP = @(a,t) [1; soft(a(2:end),lambda(1)*t)]; % proximal operator of the penalty

            % solve the task
            coef = DouglasRachfordA(zeros(N+p,1),solution,P,proxP,p+1,maxits(i),'alpha',gammaC);
        else
            % set the parameters of the Douglas-Rachford algorithm
            DR.lambda = 1;
            DR.gamma  = gammaC;
            DR.y0     = [1; zeros(p,1)];
            DR.maxit  = maxits(i);
            DR.tol    = -Inf;
            
            % set the parameters of the model
            DR.f  = @(x) 0.5*norm(fft(coef,N+p).*fft(x,N+p))^2 / (N+p);
            if decompose
                DR.prox_f = @(a,t) dXX\a;
            else
                invmat = inv(eye(p+1) + gammaC*XX);
                DR.prox_f = @(a,t) invmat*a; %#ok<MINV>
            end
            DR.g      = @(a) 0.5*norm(fft(solution,N+p).*fft(a,N+p))^2 / (N+p);
            DR.prox_g = @(a,t) [1; soft(a(2:end),lambda(1)*t)]; 
            DR.dim = p+1;

            % solve the task
            coef = DouglasRachford(DR,DR);
        end
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                          AR model extrapolation                         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if coefextra && i > 1
        % this formula is heuristic
        coef = coef + 2*(maxit-i)/(maxit-1) * (coef - oldcoef);
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                            signal estimation                            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if strcmpi(method,'inpainting') || strcmpi(method,'glp')
        AA       = zeros(size(IAA));
        b        = coef'*hankel(coef',[coef(end),zeros(1,p)]);
        AA(IAA1) = b(IAA(IAA1)+1);
        [R, er]  = chol(AA(:,indmiss));
        if er
            break
        else
            solution(~masks.R) = -R\(R'\(AA(:,indobs)*signal(indobs)));
        end
    else
        % prepare the matrix AA
        if strcmpi(mat, 'toeplitz')
            A  = toeplitz([coef; zeros(N-1,1)], [coef(1), zeros(1,N-1)]);
            AA = A'*A;
            clear A
        elseif strcmpi(mat, 'xcorr')
            b  = xcorr(coef,p)';
            AA = spdiags(ones(N,1)*b,-p:p,N,N);
        elseif strcmpi(mat, 'conv')
            b  = conv(coef,flip(coef))';
            AA = spdiags(ones(N,1)*b,-p:p,N,N);
        elseif strcmpi(mat, 'fft')
            b  = ifft(fft([coef' zeros(1,p)]) .* fft([flip(coef') zeros(1,p)]));
            AA = spdiags(ones(N,1)*b,-p:p,N,N);
        end

        % check the positive definiteness using chol
        try
            dAA = decomposition(eye(N) + gammaS*AA, 'chol');
        catch
            break
        end
        
        % check the conditioning
        if isIllConditioned(dAA)
            break
        end

        if sigaccel
            % set the parameters of the accelerated Douglas-Rachford algorithm
            if length(lambda) > 1  && lambda(2) < Inf
                P     = @(x) lambda(2)*0.5*norm(x-proj(x))^2;
                proxP = @(x,t) lambda(2)*t/(lambda(2)*t+1)*proj(x) + 1/(lambda(2)*t+1)*x;
            else
                P     = @(x) 0; % this is actually the indicator function
                proxP = @(x,t) proj(x);
            end

            % solve the task
            solution = DouglasRachfordA(zeros(N+p,1),coef,P,proxP,N,maxits(i),'alpha',gammaS);
        else
            % set the parameters of the Douglas-Rachford algorithm
            DR.lambda = 1;
            DR.gamma  = gammaS;
            DR.y0     = solution;
            DR.maxit  = maxits(i);
            DR.tol    = -Inf;

            % set the parameters of the model
            DR.f = @(x) 0.5*norm(fft(coef,N+p).*fft(x,N+p))^2 / (N+p);
            if decompose
                DR.prox_f = @(x,t) dAA\x;
            else
                invmat = inv(eye(N) + gammaS*AA);
                DR.prox_f = @(x,t) invmat*x; %#ok<MINV>
            end
            if length(lambda) > 1  && lambda(2) < Inf
                DR.g      = @(x) lambda(2)*0.5*norm(x-proj(x))^2;
                DR.prox_g = @(x,t) lambda(2)*t/(lambda(2)*t+1)*proj(x) + 1/(lambda(2)*t+1)*x;
            else
                DR.g      = @(x) 0;
                DR.prox_g = @(x,t) proj(x);
            end
            DR.dim = N;

            % solve the task
            solution = DouglasRachford(DR,DR);
        end
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                    signal rectification (for GLP only)                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if strcmpi(method,'glp') && i < maxit        
        % signal segments to be checked
        solutionU = solution(masks.U);
        solutionL = solution(masks.L);
        signalU   = signal(masks.U);
        signalL   = signal(masks.L);
        
        % indicators of the samples to be rectified
        IU = solutionU < signalU;
        IL = solutionL > signalL;
                
        % rectification
        solutionU(IU)     = signalU(IU) + abs(solutionU(IU) - signalU(IU));
        solutionL(IL)     = signalL(IL) - abs(solutionL(IL) - signalL(IL));
        solution(masks.U) = solutionU;
        solution(masks.L) = solutionL;
    end
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                           signal extrapolation                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if sigextra && i > 1
        % this formula is heuristic
        solution = solution + (maxit-i)/(maxit-1) * (solution - oldsolution);
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                linesearch                               %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if linesearch && i > 1 && i < maxit
        
        % search for the solution and coefs on the line given the previous
        % and current solution and previous and current coefs
        steps = logspace(-4,2,100);
        QQ = zeros(length(steps), 1);
        for j = 1:length(steps)
            QQ(j) = Q(solution + steps(j)*(solution-oldsolution),...
                coef + steps(j)*(coef-oldcoef));
        end
        
        % plot the objective function and the second derivative along the
        % line
        if plotLS
            x1 = solution;
            x2 = x1 - oldsolution;
            A1 = toeplitz([coef; zeros(N-1,1)], [coef(1), zeros(1,N-1)]);
            A2 = A1 - toeplitz([oldcoef; zeros(N-1,1)], [oldcoef(1), zeros(1,N-1)]);
            C0 = x1'*(A1'*A1)*x1;
            C1 = 2*(x2'*A1' + x1'*A2')*A1*x1;
            C2 = x1'*(A2'*A2)*x1 + x2'*(A1'*A1)*x2 + 2*x2'*(A2'*A1 + A1'*A2)*x1;
            C3 = 2*(x1'*A2' + x2'*A1')*A2*x2;
            C4 = x2'*(A2'*A2)*x2;
            pol1 = @(t) 0.5*(C4*t.^4 + C3*t.^3 + C2*t.^2 + C1*t + C0);
            pol2 = @(t) 6*C4*t.^2 + 3*C3*t + C2;

            % plot the objective function
            figure(fig1)
            subplot(rows,cols,i-1)
            loglog(steps, pol1(steps))
            hold on
            loglog(steps, QQ)
            legend('quadratic part','objective')
            [m1, I] = min(QQ);
            x1 = steps(I);
            Q2min = @(step) Q(solution + step*(solution-oldsolution),...
                    coef + step*(coef-oldcoef)); 
            [x2, m2] = fminbnd(Q2min, steps(1), steps(end));
            title(sprintf(...
                'objective in iteration %d, no search: [0, %.3f]\nmin. point: [%.3f, %.3f]\nfminbnd: [%.3f, %.3f]',...
                i, Q(solution, coef), x1, m1, x2, m2))

            % plot the secont derivative
            figure(fig2)
            subplot(rows,cols,i-1)
            numerical = gradient(gradient(QQ));
            polynomial = pol2(steps);
            loglog(steps(numerical>0), polynomial(numerical>0))
            hold on
            loglog(steps(numerical>0),numerical(numerical>0))          
            legend('polynomial (only quadratic part)','numerical')
            title(sprintf('2nd derivative in iteration %d',i))
            xlim(steps([1,end]))
        end

        % update the solution
        [minimum, index] = min(QQ);
        if minimum < Q(solution, coef)
            solution = solution + steps(index)*(solution-oldsolution);
            coef = coef + steps(index)*(coef-oldcoef);
        end
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                       solution and objective update                     %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % update the solution
    if saveall
        restored(:,i) = solution;
    else
        restored = solution;
    end
    
	% compute the objective value
    if nargout > 1
        objective(i) = Q(solution,coef);
    end
    
    % save the current solution for the next iteration
    oldcoef = coef;
    oldsolution = solution;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                time update                              %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if nargout > 2
        times(i) = toc;
    end

end

% rename the figures
if plotLS
    figure(fig1)
    sgtitle({method;...
        sprintf('c. accel, s. accel = [%d %d]',...
        coefaccel,sigaccel)})
    set(gcf,'WindowState','maximized')

    figure(fig2)
    sgtitle(fig2,{method;...
        sprintf('c. accel, s. accel = [%d %d]',...
        coefaccel,sigaccel)})
    set(gcf,'WindowState','maximized')
end