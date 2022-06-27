function [data_rec_fin, snr_iter, obj_iter] = spain_segmentation(data_gapped, param, paramsolver, data_orig)
% This function performs input signal padding, computation of analysis and
% synthesis windows, windowing the signal and after processing each block by
% aspain.m or sspain.m, it folds the blocks back together.
%
% Input parameters
%       data_gapped    vector of gapped signal
%       param          structure of parameters containing:
%                            Ls         length of the original signal
%                            w          window length (in samples)
%                            a          window shift (in samples)
%                            wtype      string defining the type of the window, see http://ltfat.github.io/doc/sigproc/firwin.html
%                            F          frame (usually DFT frame)
%                            algorithm  string used to select A-SPAIN or S-SPAIN
%                            mask       logical vector indicating the reliable samples
%                            fs         sampling frequency of the signal
%       paramsolver    structure of parameters for SPAIN containing:
%                            s          relaxation stepsize
%                            r          relaxation steprate
%                            epsilon    stopping threshold of termination function
%                            maxit      maximal possible number of iterations with particular settings
%                            store_snr  switch to enable computing SNR in each iteration
%                            store_obj  switch to enable storing the value of termination function in each iteration
%                            f_update   switch between different approximations of f-update in case of S-SPAIN
%       data_orig      vector of original clean signal used to compute SNR during iterations
%
% Output parameters
%       data_rec_fin   vector of reconstructed signal after OLA
%       snr_iter       matrix of SNR values during iterations for all blocks of signal
%       obj_iter       matrix of termination function values during iterations for all blocks of signal
%
% Date: 08/07/2020
% By Ondrej Mokry, Pavel Zaviska
% Brno University of Technology
% Contact: ondrej.mokry@mensa.cz

% preparation for padding at the end of the signal
L = ceil(param.Ls/param.a)*param.a+(ceil(param.w/param.a)-1)*param.a; % L is divisible by a and minimum amount of zeros equals gl (window length). Zeros will be appended to avoid periodization of nonzero samples.
N = L/param.a; % number of signal blocks

% padding the signals and mask to length L
padding = zeros(L-param.Ls, 1);
data_gapped = [data_gapped; padding];
data_orig = [data_orig; padding];
param.mask = [param.mask; true(L-param.Ls,1)];

% construction of analysis and synthesis windows
g = gabwin(param.wtype, param.a, param.w, L);
gana = normalize(g,'peak'); % peak-normalization of the analysis window
gsyn = gabdual(gana, param.a, param.w)*param.w; % computing the synthesis window

% this is substituting fftshift (computing indexes to swap left and right half of the windows)
idxrange = [0:ceil(param.w/2)-1,-floor(param.w/2):-1];
idxrange2 = idxrange+abs(min(idxrange))+1;

% initialization of param_seg (parameters for one signal block)
param_seg = param;
param_seg.Ls = param.w;
param_seg.mask = true(param.w,1);

% initialization of signal blocks
data_block = zeros(param.w,1);
data_orig_block = zeros(param.w,1);
data_rec_fin = zeros(L,1);

% initialization of matrices for recording SNR and for termination function during iterations
snr_iter = NaN(paramsolver.maxit,N);
obj_iter = NaN(paramsolver.maxit,N);


for n=0:N-1
    % multiplying signal block with windows and choosing corresponding masks
    idx = mod(n*param.a + idxrange,L) + 1;
    data_block(idxrange2) = data_gapped(idx).*gana;
    data_orig_block(idxrange2) = data_orig(idx).*gana;
    param_seg.mask(idxrange2) = param.mask(idx);
    
    % if the block does not contain any samples to restore, skip it
    if sum(param_seg.mask) == param.w
        continue;
    end
    
    % perform SPAIN
    switch param.algorithm
        case 'sspain'
           [data_rec_block, snr_iter(:,n+1), obj_iter(:,n+1)] = sspain(data_block, param_seg, paramsolver, data_orig_block);
        case 'aspain'
           [data_rec_block, snr_iter(:,n+1), obj_iter(:,n+1)] = aspain(data_block, param_seg, paramsolver, data_orig_block);
        case 'pspain'
           [data_rec_block, snr_iter(:,n+1), obj_iter(:,n+1)] = pspain(data_block, param_seg, paramsolver, data_orig_block);
    end
    
    % folding blocks together using Overlap-add approach (OLA)
    data_rec_block = real(ifftshift(data_rec_block));
    data_rec_fin(idx) = data_rec_fin(idx) + data_rec_block.*gsyn;
end

% ensure equality of solution and data in reliable part
data_rec_fin(param.mask) = data_gapped(param.mask);

% crop the padding of reconstructed signal
data_rec_fin = data_rec_fin(1:param.Ls);