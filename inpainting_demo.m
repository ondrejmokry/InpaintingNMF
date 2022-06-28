clear
clc
close all
rng(0)
addpath('utils')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                               the inputs                                %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% create the signal and mask
% signal
[y, fs] = audioread('signals/mamavatu.wav');
y  = resample(y,16e3,fs);
fs = 16e3;
L  = length(y);

% mask
mask = rand(L,1) > 0.6;

%% settings
% algorithm
maxit = 100; % (50 in declipping survey)
nmfit = 1;   % (1 in declipping survey)
epsilon = 1e-4;

% segmentation
M = 1024; % segment length (2048 in declipping survey)
a = 512;  % segment shift  (1024 in declipping survey)

% model
F = 2*M;  % number of frequency channels (M in declipping survey)
K = 10;   % number of components for NMF (20 in declipping survey)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                              the algorithm                              %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

methods = {'AM','EMtf',5,'EMtf','EMt'};
methodnames = {...
    'AM (F = M)',...
    'EM-tf (F = M)',...
    'AM-to-EM-tf (F = M)',...
    'EM-tf (F = 2M)',...
    'EM-t (F = 2M)'};

results.signals    = cell(length(methods),1);
results.Ws         = cell(length(methods),1);
results.Hs         = cell(length(methods),1);
results.SNRs       = cell(length(methods),1);
results.objectives = cell(length(methods),1);
results.relnorms   = cell(length(methods),1);

for m = 1:length(methods)
    
    tic
    if m <= 3 % the case with F = M: AM, EM-tf (equivalent with EM-t), switching
        [restored, W, H, ~, objectives] = ainmf(methods{m},y,mask,K,maxit,...
            'M',M,'a',a,'F',M,'saveall',true,'verbose',true,'nmfit',nmfit,...
            'drawing',true,...
            'epsilon',epsilon);
    else % the case with F = 2*M: EM-tf, EM-t
        [restored, W, H, ~, objectives] = ainmf(methods{m},y,mask,K,maxit,...
            'M',M,'a',a,'F',F,'saveall',true,'verbose',true,'nmfit',nmfit,...
            'drawing',true,...
            'epsilon',epsilon);
    end
    t = toc;
 
    % compute SNRs
    SNRs = NaN(size(restored,2),1);
    for i = 1:size(restored,2)
        SNRs(i) = snr(y(~mask),y(~mask)-restored(~mask,i));
    end
    
    % compute relative norms
    relnorms = Inf(size(restored,2),1);
    for i = 2:size(restored,2)
        relnorms(i) = norm(restored(:,i)-restored(:,i-1))...
            /norm(restored(:,i-1));
    end
    
    fprintf('     Elapsed time is %.1f seconds.\n',t)
    fprintf('     SNR is %.2f dB.\n',SNRs(end))
    
    % save
    results.signals{m}    = restored(:,end);
    results.Ws{m}         = W;
    results.Hs{m}         = H;
    results.SNRs{m}       = SNRs;
    results.objectives{m} = objectives;
    results.relnorms{m}   = relnorms;
    
end

%% save the data
dt = datetime(clock);
save(fname,'a','dt','epsilon','F','K','M','mask','maxit','methodnames','methods','nmfit','results','y')