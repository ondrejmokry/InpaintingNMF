function [data_rec, snr_iter, obj_val] = a_spain_learned(data_gapped, param, paramsolver,data_orig,Basis)
% 
% Input parameters
%       data_gapped    vector of gapped signal
%       param          structure of parameters containing:                 
%                            a          window shift (in samples)
%                            M          number of channels
%                            mask       logical vector indicating the reliable samples
%                            gwindow    Gabor window
%                            gdual	should be chosen as canonical dual window of gwindow 
%       paramsolver    structure of parameters for SPAIN containing:
%                            s          relaxation stepsize
%                            r          relaxation steprate
%                            epsilon    stopping threshold of termination function
%                            maxit      maximal possible number of iterations with particular settings
%                            store_snr  switch to enable computing SNR in each iteration
%                            store_obj  switch to enable storing the value of termination function in each iteration
%       data_orig      vector of original clean signal used to compute SNR during iterations
%       Basis	       sparsity enhancing unitary matrix
%
% Output parameters
%       data_rec       vector of restored signal
%       snr_iter       SNR values during iterations
%       obj_iter       termination function values during iterations
%

% initialization of variables
x_hat = data_gapped;
L = length(data_gapped);
BasisInv=Basis'; %unitary matrix
zEst = Basis*(dgtreal(x_hat,param.gwindow,param.a,param.M,dgtlength(L,param.a,param.M),'timeinv'));  % zEst = A*x_hat
dim_zEst=size(zEst);
u = zeros(dim_zEst);
k = paramsolver.s;
cnt = 1;
bestObj = Inf;

% initialization of vectors for SNR and termination function during iterations
snr_iter = NaN(paramsolver.maxit, 1);
obj_val = NaN(paramsolver.maxit, 1);

while cnt <= paramsolver.maxit
 
   z_bar = hard_thresholding_dgtreal(zEst + u, k); 

   objVal = norm((zEst - z_bar),'fro'); % update termination function
    
    % make a record if the objective function has decreased
    if objVal <= bestObj
        data_rec = x_hat;
        bestObj = objVal;
    end
    
    % termination step
    if objVal <= paramsolver.epsilon
        break
    end
    
    b = z_bar-u;
    syn = idgtreal(BasisInv*b,param.gdual,param.a,param.M,'timeinv');
    x_hat = syn;
    x_hat(param.mask) = data_gapped(param.mask);
     
    % computation and storing of SNR and termination function if requested
    if paramsolver.store_snr
        snr_iter(cnt) = snr_n(data_orig(~(param.mask)),  x_hat(~(param.mask)));
    end
    if paramsolver.store_obj
        obj_val(cnt) = objVal;
    end
    % dual variable update   
    zEst = Basis*(dgtreal(x_hat,param.gwindow,param.a,param.M,dgtlength(L,param.a,param.M),'timeinv'));
    u = u + zEst - z_bar;
    
    % iteration counter update
    cnt = cnt + 1;
    
    % incrementation of variable k (require less sparse signal in next iteration)
    if mod(cnt,paramsolver.r) == 0
        k = k + paramsolver.s;
    end
 
end
end
%end of function
