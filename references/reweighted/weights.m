function [ wt ] = weights( F, L, rel, U, V, atoms, opt)
% weighting atoms of frame of the type TIMEINV
%
% inputs:   F ...... frame, F = frame('type',{'window',w},a,M)
%           L ...... signal length
%           rel .... mask of the reliable samples
%           U ...... index of the first window overlapping with the gap
%                    (U = 1 for declipping)                     
%           V ...... index of the last window overlapping with the gap
%                    (V = L/a for declipping)
%           atoms .. matrix of frame atoms for zero shift
%           opt .... 'supp'  / 'abs' / 'norm' / 'energy'
% output:   wt ..... atom weights
%
% Date: 08/07/2020
% By Ondrej Mokry
% Brno University of Technology
% Contact: ondrej.mokry@mensa.cz

%% frame parameters
% m ... number of atoms per one shift
type = F.type;
a = F.a;
M = F.M;
if strcmp(type,'dgtreal')
        m = floor(M/2) + 1;
else
    if strcmp(type,'dgt')
        m = M;      
    else
        disp('Frame type not compatible. Use ''dgt'' or ''dgtreal''.');
    end
end
[w,~] = size(atoms);
n = frameclength(F,L)/m; % number of shifts
wt = ones(m*n,1); % initialization of the weights

%% weights of the "relevant" atoms
atoms = fftshift(atoms,1);
rel_shift = zeros(L+w,1);
rel_shift(1:w/2) = rel(L+1-w/2:L);
rel_shift(w/2+1:w/2+L) = rel;
rel_shift(w/2+L+1:end) = rel(1:w/2);
for i = U:V
    reliable = rel_shift((i-1)*a + 1:(i-1)*a + w);
    switch opt
        case 'supp' % based on atom support
            AE = logical(abs(atoms));
            Sum_E = ones(1,w)*AE;
            en = reliable'*AE;
            wt((i-1)*m+1:i*m) = en./Sum_E;
            
        case 'energy' % based on atom energy
            AE = abs(atoms).^2;
            Sum_E = ones(1,w)*AE;
            en = reliable'*AE;
            wt((i-1)*m+1:i*m) = en./Sum_E;

        case 'norm' % based on atom l2 norm
            AE = abs(atoms).^2;
            Sum_E = ones(1,w)*AE;
            en = reliable'*AE;
            wt((i-1)*m+1:i*m) = (en.^0.5)./(Sum_E.^0.5);

        case 'abs' % based on atom l1 norm
            AE = abs(atoms);
            Sum_E = ones(1,w)*AE;
            en = reliable'*AE;
            wt((i-1)*m+1:i*m) = en./Sum_E;
    end
end
end