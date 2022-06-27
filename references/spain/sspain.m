function [data_rec, snr_iter, obj_val] = sspain(data_gapped, param, paramsolver, data_orig)
% implementation of the S-SPAIN algorithm
%
% f-update approximation:
%    (H)    Hk(As) (hard thresholding) ( s = x_hat - u )
%    (OMP)  franamp(s,...), maxit = k,  conj_atoms = true
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
%       data_rec       vector of restored block of signal
%       snr_iter       SNR values during iterations
%       obj_iter       Termination function values during iterations
%
%   x_hat (reconstructed signal - waveform)
%   z_bar (signal coefficients after hard thresholding)
%   xEst (signal synthesized from z_bar)
%   u (dual variable)
%   k (required sparsity)
%
% Date: 08/07/2020
% By Ondrej Mokry, Pavel Zaviska
% Brno University of Technology
% Contact: ondrej.mokry@mensa.cz

% initialization of variables
x_hat = data_gapped;
u = zeros(length(x_hat), 1);
k = paramsolver.s;
cnt = 1;
bestObj = Inf;

% initialization of vectors for SNR and termination function during iterations
snr_iter = NaN(paramsolver.maxit, 1);
obj_val = NaN(paramsolver.maxit, 1);

while cnt <= paramsolver.maxit
    % sets all but k largest coefficients to zero
    %(complex conjugated pairs are taken into consideration)
    switch paramsolver.f_update
        case {'H','h'}
            z_bar = hard_thresholding(frana(param.F,x_hat-u), k);
        case {'OMP','omp'}
            z_bar = franamp(param.F,x_hat-u,'omp','qr','maxit',k);            
    end
    objVal = norm(frsyn(param.F,z_bar) - x_hat);
    
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
    xEst = real(frsyn(param.F,z_bar));
    x_hat = proj_time(xEst + u, param.mask, data_gapped);
    
    % computation and storing of SNR and termination function if requested
    if paramsolver.store_snr
        snr_iter(cnt) = snr_n(data_orig, x_hat);
    end
    if paramsolver.store_obj
        obj_val(cnt) = objVal;
    end
    
    % dual variable update
    u = u + xEst - x_hat;
    
    % iteration counter update
    cnt = cnt+1;
    
    % incrementation of variable k (require less sparse signal in next iteration)
    if mod(cnt,paramsolver.r) == 0
        k = k + paramsolver.s;
    end
    
end

end