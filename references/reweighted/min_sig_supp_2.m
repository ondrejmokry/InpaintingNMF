function [q,Q,p,P,S,F,u,v,U,V,L] = min_sig_supp_2(w,a,M,s,f,N,varargin)
% Find the minimum range of signal that still carries all info about the gap in DGT
%
%   Usage:  [q,Q] = min_sig_supp(w,a,M,s,f,N);
%              computes the first and last index of signal needed for reconstruction
%              and first and last atom of the dictionary needed for inpainting
%
%           [q,Q] = min_sig_supp(w,a,M,s,f,N,neigh,offset);
%
%           [q,Q,p,P,S,F,u,v,U,V,L] = min_sig_supp(w,a,M,s,f,N);
%
%
%   Input parameters:
%         w     : Window length.
%         a     : Time shift of window.
%         M     : Number of modulation channels.
%         s     : Left edge index of the gap.
%         f     : Right edge index of the gap.
%         N     : Length of the original signal.
%        offset : Position offset of the first overlapping window.
%                 Choose from [0,w-1]; optional, default 0
%         neig  : horizontal size of neighborhood in case of structured sparsity.
%                 Default 1.
%    
%   Output parameters:
%         q     : First index of shortened signal wrt. original signal.
%         Q     : Last index of shortened signal wrt. original signal.
%         p     : Middle index of the first "useful" window overlapping with gap.
%         P     : Middle index of the last "useful" window overlapping with gap.
%         S     : Position of first "useful" window in the series of DGT coefficients.
%         F     : Position of last "useful" window in the series of DGT coefficients.
%         u     : Left edge of the gap within shortened signal.
%         v     : Right edge of the gap within shortened signal.
%         U     : Position of first "useful" window within shortened signal.
%         V     : Position of last "useful" window within shortened signal.
%         L     : Length of the shortened signal.

% (c) 2015 Hana Bartlová, Pavel Rajmic, Brno University of Technology, Czech Republic
% modified by Ondrej Mokry, 2020

%%
% number of first useful window
S=ceil((s-ceil(w/2))/a)+1;

% middle index of the first window overlapping with the gap
p=1+(S-1)*a;

% correction of the first overlapping window position
switch nargin
    case 6
        offset = 0; neig = 1;
    case 7
        neig = varargin{1};
        offset = 0;
    case 8
        neig = varargin{1};
        offset = varargin{2};
end
offset = mod(offset,a); % changing the offset by the parameter a does not have any effect
p = p + offset;

% if the offset is large enough, the gap can overlap with windows previous
% to the first useful window computed with offset = 0
if p-a+ceil(w/2)-1 >= s
    S = S - 1;
    p = p - a;
end

% first index of the shortened part
q=p-ceil(floor(w/2)/a)*a;

% number of last useful window
F=S+floor((f+floor(w/2)-p)/a);

% middle index of the last useful window
P=p+(F-S)*a;

% find the last sample of the shortened part
Q=P+ceil(ceil(w/2)/a)*a;

% does nothing in case of no groups; otherwise enlarges the signal piece
q = q - (neig-1)*a;
Q = Q + (neig-1)*a;

if P+ceil(w/2)-1>N
    disp('Last "useful" window exceeds length of the signal!')
end

if p-floor(w/2)<1
    disp('First "useful" window exceeds length of the signal!')
end

if Q>N
    disp('Index Q exceeds length of the signal!') %It can be padded by zeros
end

% position of the gap within the shortened signal
u=s-q+1;
v=f-q+1;

% number of the first and last useful window within the shortened signal
U = (p-q)/a+1;
V = U + (F-S);

% final length used
L=Q-q+1;
