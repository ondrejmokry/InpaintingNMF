close all
clear
clc
addpath('utils')
addpath(genpath('references'))

% future file name
filename = 'results/inpainting_comparison_01.mat';

%% settings
load('signals/EBU_SQAM.mat')
signals = { 'a08_violin',...
            'a16_clarinet',...
            'a18_bassoon',...
            'a25_harp',...
            'a35_glockenspiel',...
            'a41_celesta',...
            'a42_accordion',...
            'a58_guitar_sarasate',...
            'a60_piano_schubert',...
            'a66_wind_ensemble_stravinsky' };
methods = { 'EMtf',...
            'AM',...
            'Janssen',...
            'SPAIN',...
            'reweighted',...
            'SPAINmod',...
            'SPAINlearned',...
            'Janssenmod',...
            'LR',...
            'AMtoEMtf'};
methodnames = { 'EM-tf',...
            'AM'...
            'Janssen',...
            'SPAIN',...
            'reweighted',...
            'SPAINmod',...
            'SPAINlearned',...
            'Janssenmod',...
            'LR',...
            'AM-to-EM-tf',};

% signals used
signums = 1:10;

% gap lengths used (in ms)
glengths = 20:10:80;

% number of gaps per signal
gaps = 10;

%% parameters of all the methods
% NMF-based methods (EM-tf, EM-t, AM)
NMF.M = 4096;             % window length
NMF.a = 2048;             % window shift
NMF.F = 4096;             % number of frequency channels
NMF.K = 20;               % number of NMF components
NMF.maxit = 100;          % number of EM/AM iterations
NMF.nmfit = 10;           % number of inner loop iterations (multiplicative updates)
NMF.epsilon = 1e-6;       % additive constant for NMF to avoid zeros
NMF.tolsol = 1e-5;        % stopping value for the relative norm of the solution update

% Janssen
JAN.w = 4096;             % window length
JAN.a = 2048;             % window shift
JAN.p = 512;              % AR model order
JAN.wtype = 'sine';       % window shape
JAN.maxit = 100;          % number of iterations
JAN.lambda = 0;           % regularization parameter of the coef. sparsity

% SPAIN
SPA.algorithm = 'aspain'; % algorithm
SPA.w = 4096;             % window length
SPA.a = 2048;             % window shift
SPA.wtype = 'sine';       % window shape
SPA.F = frame('dft');     % frequency transform
SPA.F.redundancy = 2;     % frequency transform redundancy
SPA.F.frana = @(insig)dft([insig; zeros(length(insig)*(SPA.F.redundancy-1),1)]);
SPA.F.frsyn = @(insig)postpad(idft(insig),length(insig)/SPA.F.redundancy);
SPA.s = 1;                % increment of k
SPA.r = 1;                % every r-th iteration increment k by s   
SPA.epsilon = 0.001;      % stopping criterion of termination function
SPA.maxit = ceil(floor(SPA.w*SPA.F.redundancy/2+1)*SPA.r/SPA.s); % maximum number of iterations
SPA.store_snr = false;    % save SNR during iterations?
SPA.store_obj = false;    % save objective during iterations?

% (re)weighted l1 relaxation
REW.type = 'analysis';    % formulation
REW.F = frametight(frame('dgt',{'sine',4096},2048,4096,'timeinv')); % frame
REW.offset = 'half';      % offset
REW.weighting = 'energy'; % atom weighting
REW.reweighting = false;  % reweighting
REW.maxit = 2000;         % number of iterations
REW.tol = 1e-5;           % stopping value for the relative norm of the solution update
REW.tau = 1;              % parameter of the proximal step
REW.sigma = 1;            % parameter of the proximal step

% SPAIN modified
SPM.margin = 8192;        % margin around each gap
SPM.algorithm = 'aspain_mod'; % formulation
SPM.w = 4096;             % window length
SPM.a = 2048;             % window shift
SPM.wtype = 'sine';       % window shape
SPM.M = 8192;             % number of frequency channels
SPM.gwindow = gabwin(SPM.wtype,SPM.a,SPM.M);
SPM.gwindow = normalize(SPM.gwindow,'peak'); % peak-normalization of the analysis window
SPM.gdual = gabdual(SPM.gwindow,SPM.a,SPM.M); 
SPM.s = 1;                % increment of k
SPM.r = 1;                % every r-th iteration increment k by s   
SPM.epsilon = 0.001;      % stopping criterion of termination function
SPM.maxit = ceil(floor(SPM.w/2+1)*SPM.r/SPM.s)*2; % maximum number of iterations
SPM.store_snr = false;    % save SNR during iterations?
SPM.store_obj = false;    % save objective during iterations?

% SPAIN learned
SPL.margin = 8192;        % margin around each gap
SPL.algorithm = 'aspain_LEARNED'; % formulation
SPL.w = 4096;             % window length
SPL.a = 2048;             % window shift
SPL.wtype = 'sine';       % window shape
SPL.M = 8192;             % number of frequency channels
SPL.gwindow = gabwin(SPL.wtype, SPL.a, SPL.M);
SPL.gwindow = normalize(SPL.gwindow,'peak'); % peak-normalization of the analysis window
SPL.gdual = gabdual(SPL.gwindow,SPL.a,SPL.M); 
SPL.s = 1;                % increment of k
SPL.r = 1;                % every r-th iteration increment k by s   
SPL.epsilon = 0.001;      % stopping criterion of termination function
SPL.maxit = ceil(floor(SPL.w/2+1)*SPL.r/SPL.s)*2; % maximum number of iterations
SPL.store_snr = false;    % save SNR during iterations?
SPL.store_obj = false;    % save objective during iterations?

% Janssen modified
JAM.margin = 4096;        % margin around each gap
JAM.p      = 512;         % AR model order
JAM.maxit  = 100;         % number of iterations
JAM.lambda = 0;           % regularization parameter of the coef. sparsity

% LR
LRE.margin = 4096;        % margin around each gap
LRE.p      = 512;         % AR model order

%% tables for results
% initialization
types = {'cell','double','double','cell','double','cell'};
names = {'signal','fs','gap','relnorms','time','SNR'};
units = {'','Hz','ms','','s','dB'};
rows  = length(signums)*length(glengths);
for m = 1:length(methods)
    tables.(methods{m}) = table('Size',[rows, length(types)],...
        'VariableTypes',types,'VariableNames',names);
    tables.(methods{m}).Properties.VariableUnits = units;
end

% future folder for restored signals
% sigfold = ['results/signals/',filename(9:end-4)];
% [~,~] = mkdir(sigfold);

% if the experiment has been stopped, load the data
if isfile(filename)
    load(filename)
end
skipped = 0;

% save time for global timer
cinit = clock;

%% processing
row = 0;
for signum = signums
    for glength = glengths

        row = row + 1;
        
        % if the experiment has been stopped, skip what has been already done
        if ~isempty(tables.(methods{length(methods)}).SNR{row})
            skipped = skipped + 1;
            continue
        end
        
        % command window output
        str = sprintf('Signal: %s',signals{signum});
        fprintf(repmat('=',length(str),1))
        fprintf('\n')
        fprintf(str)
        fprintf('\nGap length: %d ms\n',glength)
        fprintf(repmat('=',length(str),1))
        fprintf('\n')

        % estimate the remaining time
        c = clock;
        d = etime(c, cinit); % elapsed time in seconds
        d = d/3600; % elapsed time in hours
        fprintf('Elapsed time: %d hours\n',round(d))
        fprintf('Estimated remaining time: %d hours\n',...
            round(d*((rows-skipped)/(row-skipped)-1)))
        
        %% loading and degrading the signal
        % load the signal
        signal = eval(signals{signum});
        L      = length(signal);

        % compute gap parameters
        h      = round(fs*glength/1000); % gap length in samples 
        margin = 0.25; % length of reliable part at the start / end of the signal (is seconds)

        % set the mask
        mask    = true(L,1);
        slength = round((length(signal)-2*margin*fs)/gaps);
        w       = 2048; % half of the longest window length (used as a margin)
        rng(0)
        for i = 1:gaps
            idxs = round(margin*fs)+((i-1)*slength+1:i*slength);
            s = (w+1) + rand()*(slength-2*w-h); % gap start in samples
            s = round(s);
            f = s + h - 1; % gap end in samples
            smask = true(slength,1); % indicator of the reliable samples
            smask(s:f) = false;
            mask(idxs) = smask; 
        end
        gapped = signal.*mask;

        %% inpainting
        for m = 1:length(methods)
            fprintf('%-15s ',[methodnames{m},'...'])
            tic
            if strcmpi(methods{m},'EMtf') || strcmpi(methods{m},'EMt') || strcmpi(methods{m},'am')
                restored = ainmf(...
                    methods{m},... % method
                    gapped,... % degraded signal
                    mask,... % mask of the reliable samples
                    NMF.K,... % number of NMF components
                    NMF.maxit,... % number of outer iterations
                    'nmfit',NMF.nmfit,... % number of NMF iterations
                    'M',NMF.M,'a',NMF.a,'F',NMF.F,... % transform parameters
                    'epsilon',NMF.epsilon,... % parameter of NMF
                    'tolsol',NMF.tolsol,... % stopping criterion
                    'saveall',true,'verbose',false); % outputs
            elseif strcmpi(methods{m},'AMtoEMtf')
                restored = ainmf(...
                    5,... % switching iteration from AM to EM-tf
                    gapped,... % degraded signal
                    mask,... % mask of the reliable samples
                    NMF.K,... % number of NMF components
                    NMF.maxit,... % number of outer iterations
                    'nmfit',NMF.nmfit,... % number of NMF iterations
                    'M',NMF.M,'a',NMF.a,'F',NMF.F,... % transform parameters
                    'epsilon',NMF.epsilon,... % parameter of NMF
                    'tolsol',NMF.tolsol,... % stopping criterion
                    'saveall',true,'verbose',false); % outputs
            elseif strcmpi(methods{m},'janssen')
                masks.R = mask;
                masks.U = false(L,1);
                masks.L = false(L,1);
                restored = segmentation(...
                    'inpainting',... % method
                    gapped,... % degraded signal
                    masks,... % mask of the reliable samples
                    JAN.lambda,... % regularization parameter of the coef. sparsity
                    JAN.p,... % AR model order
                    JAN.maxit,... % number of iterations
                    'w',JAN.w,'a',JAN.a,'wtype',JAN.wtype,... % window parameters
                    'saveall',true,'verbose',false); % outputs
            elseif strcmpi(methods{m},'spain')
                SPA.mask = mask;
                SPA.Ls = L;
                tic
                restored = spain_segmentation(...
                    gapped,... % degraded signal
                    SPA,... % model parameters
                    SPA,... % algorithm parameters
                    signal); % clean signal (to compute SNR)
            elseif strcmpi(methods{m},'reweighted')
                restored = reweighted(...
                    gapped,... % degraded signal
                    mask,... % mask of the reliable samples...
                    REW,... % transform parameters
                    REW,... % Chambolle-Pock parameters
                    [],... % Douglas-Rachford parameters (irrelevant)
                    []); % reweighting parameters (irrelevant)
            elseif strcmpi(methods{m},'spainmod')
                % find missing indices
                missing = find(~mask);
                
                % prepare solution
                restored = gapped;

                % iterate over all the gaps
                for i = 1:gaps
                    marge   = SPM.margin;
                    idxs    = missing(1+(i-1)*h) - marge : missing(i*h) + marge;
                    smask   = mask(idxs);
                    sdata   = signal(idxs);
                    sgapped = gapped(idxs);
                    
                    % inpaint
                    SPM.Ls = length(smask);
                    SPM.mask = smask;
                    if contains(SPM.algorithm,'aspain')
                        srestored = a_spain_learned(sgapped, SPM, SPM, sdata, eye(SPM.M/2+1));
                    else
                        srestored = s_spain_learned(sgapped, SPM, SPM, sdata, eye(SPM.M/2+1));
                    end
                    restored(idxs) = srestored(1:length(idxs));
                end
            elseif strcmpi(methods{m},'spainlearned')
                % find missing indices
                missing = find(~mask);
                
                % prepare solution
                restored = gapped;

                % iterate over all the gaps
                for i = 1:gaps
                    marge   = SPL.margin;
                    idxs    = missing(1+(i-1)*h) - marge : missing(i*h) + marge;
                    smask   = mask(idxs);
                    sdata   = signal(idxs);
                    sgapped = gapped(idxs);
                    
                    % basis optimization starts...
                    q  = idxs(1);
                    Q  = idxs(end);
                    q1 = q - SPL.w*5;
                    Q1 = Q + SPL.w*5;

                    % save some data
                    sdata_store   = sdata;
                    sgapped_store = sgapped;
                    if Q1 > length(sgapped_store) 
                        sgapped_store = [sgapped_store; zeros(Q1-length(sgapped_store) ,1)]; %#ok<AGROW>
                        sdata_store   = [sdata; zeros(Q1-length(sdata) ,1)];
                    end

                    % prepare learning neighborhood
                    coeff_TF   = dgtreal(sgapped_store,SPL.gwindow,SPL.a,SPL.M,'timeinv');
                    Mprime     = size(coeff_TF,1);
                    Nprime     = size(coeff_TF,2);
                    center     = round((Q1-q1)/(2*SPL.a) + min(q1,0)/SPL.a);
                    offs       = round((Q-q)/(2*SPL.a));
                    learn_offs = (2*SPL.w/SPL.a); % determines the size of learning neighborhood
                    X_tr       = [coeff_TF(2:(Mprime-1),(max((center-offs-learn_offs),1)):(max((center-offs),1))) coeff_TF(2:(Mprime-1),(center+offs):(center+offs+learn_offs))];
                    
                    % compute optimized basis
                    [Basis_small, sparsity_init, sparsity_final] = basis_opt_new(X_tr,1,2^(-20));
                    Basis = eye(Mprime);
                    Basis(2:(Mprime-1),2:(Mprime-1)) = Basis_small;

                    % inpaint
                    SPL.Ls = length(smask);
                    SPL.mask = smask;
                    if contains(SPL.algorithm,'aspain')
                        srestored = a_spain_learned(sgapped, SPL, SPL, sdata, Basis);
                    else
                        srestored = s_spain_learned(sgapped, SPL, SPL, sdata, Basis);
                    end
                    restored(idxs) = srestored(1:length(idxs));
                end
            elseif strcmpi(methods{m},'janssenmod')
                % find missing indices
                missing = find(~mask);
                
                % prepare solution
                restored = repmat(gapped,1,JAM.maxit);

                % iterate over all the gaps
                for i = 1:gaps
                    marge    = JAM.margin;
                    idxs     = missing(1+(i-1)*h) - marge : missing(i*h) + marge;
                    smask    = mask(idxs);
                    smasks.R = smask;
                    smasks.U = false(length(smask),1);
                    smasks.L = false(length(smask),1);
                    
                    % inpaint
                    restored(idxs,:) = janssen('inpainting',... % method
                        gapped(idxs),... % degraded signal
                        smasks,... % mask of the reliable samples
                        JAN.lambda,... % regularization parameter of the coef. sparsity
                        JAM.p,... % AR model order
                        JAM.maxit,... % number of iterations
                        'saveall',true,'verbose',false); % outputs
                end
            elseif strcmpi(methods{m},'lr')
                % find missing indices
                missing = find(~mask);
                
                % prepare solution
                restored = gapped;

                % iterate over all the gaps
                for i = 1:gaps
                    marge   = LRE.margin;
                    idxs    = missing(1+(i-1)*h) - marge : missing(i*h) + marge;
                    smask   = mask(idxs);
                    sgapped = gapped(idxs);
                    presig  = sgapped(1:marge);
                    postsig = sgapped(marge+h+1:end);
                    t       = linspace(0, pi/2, h)';
                    sqCos   = cos(t).^2;
                    
                    % forward prediction
                    af = lpc(presig, LRE.p);
                    zf = filtic(1, af, presig(end-(0:(LRE.p-1))));
                    prediction = filter(1,af,zeros(1,h),zf)';

                    % backward prediction
                    postsig = flip(postsig);
                    ab = lpc(postsig, LRE.p);
                    zb = filtic(1, ab, postsig(end-(0:(LRE.p-1))));
                    postdiction = flip(filter(1,ab,zeros(1,h),zb)');

                    % composition
                    restored(missing(1+(i-1)*h) : missing(i*h)) = sqCos.*prediction + flip(sqCos).*postdiction;
                end
            end
            t = toc;
                
            % relative norms
            relnorms = NaN(size(restored,2),1);
            for i = 2:size(restored,2)
                relnorms(i) = norm(restored(:,i)-restored(:,i-1))...
                    /norm(restored(:,i-1));
            end
                
            % SNR
            SNR = NaN(size(restored,2),1);
            for i = 1:size(restored,2)
                SNR(i) = snr(signal(~mask),signal(~mask)-restored(~mask,i));
            end
            
            % command window output
            fprintf('time: %4d s, SNR end: %6.2f dB, SNR max: %6.2f dB\n',round(t),SNR(end),max(SNR))

            % write the data
            tables.(methods{m}).signal{row} = signals{signum};
            tables.(methods{m}).fs(row) = fs;
            tables.(methods{m}).gap(row) = glength;
            tables.(methods{m}).relnorms{row} = relnorms;
            tables.(methods{m}).time(row) = t;
            tables.(methods{m}).SNR{row} = SNR;
            
            % save the signals
%             dt = datetime(clock);
%             save(sprintf('%s/%s_%02d_%s.mat',sigfold,signals{signum}(1:3),glength,methods{m}),...
%                 'fs','glength','signal','mask','restored','SNR')
         end
        
        %% save the data
        save(filename,...
            'glengths',...
            'gaps',...
            'NMF','JAN','SPA','REW','SPM','SPL','JAM','LRE',...
            'methodnames',...
            'methods',...
            'signals',...
            'signums',...
            'tables',...
            '-v7.3')
        
    end % gapnum

end % signum