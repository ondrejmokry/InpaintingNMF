function [data_rec, snr_iter, obj_val] = aspain(data_gapped, param, paramsolver, data_orig)
% implementation of the A-SPAIN algorithm
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
%       data_rec       vector of restored signal
%       snr_iter       SNR values during iterations
%       obj_iter       termination function values during iterations
%
%   x_hat (reconstructed signal - waveform)
%   z_bar (signal coefficients after hard thresholding)
%   u (dual variable, vector of signal coefficients)
%   k (required sparsity)
%   zEst = A*x_hat
%
% Date: 08/07/2020
% By Ondrej Mokry, Pavel Zaviska
% Brno University of Technology
% Contact: ondrej.mokry@mensa.cz


% initialization of variables
x_hat = data_gapped;
zEst = frana(param.F, x_hat);  % zEst = A*x_hat
u = zeros(length(zEst), 1);
k = paramsolver.s;
cnt = 1;
bestObj = Inf;

% initialization of vectors for SNR and termination function during iterations
snr_iter = NaN(paramsolver.maxit, 1);
obj_val = NaN(paramsolver.maxit, 1);

while cnt <= paramsolver.maxit
    
    % set all but k largest coefficients to zero
    % (complex conjugated pairs are taken into consideration)
    z_bar = hard_thresholding(zEst + u, k);
    
    objVal = norm(zEst - z_bar); % update termination function
    
    % make a record if the objective function has decreased
    if objVal <= bestObj
        data_rec = x_hat;
        bestObj = objVal;
    end
    
    % termination step
    if objVal <= paramsolver.epsilon
        break
    end
    
    % projection onto the set of feasible solutions
    b = z_bar - u;
    syn = frsyn(param.F, b);
    x_hat = proj_time(syn, param.mask, data_gapped);
    
    % computation and storing of SNR and termination function if requested
    if paramsolver.store_snr
        snr_iter(cnt) = snr_n(data_orig,  x_hat);
    end
    if paramsolver.store_obj
        obj_val(cnt) = objVal;
    end
    
    % dual variable update
    zEst = frana(param.F, x_hat);
    u = u + zEst - z_bar;
    
    % iteration counter update
    cnt = cnt+1;
    
    % incrementation of variable k (require less sparse signal in next iteration)
    if mod(cnt,paramsolver.r) == 0
        k = k + paramsolver.s;
    end
    
end

end