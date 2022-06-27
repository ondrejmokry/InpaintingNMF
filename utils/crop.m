function y = crop(x,N)
% CROP crops the input array x to length N along the first dimension.

% read the size and original length along the first dimension
oldsz = size(x);
oldN  = oldsz(1);

% reshape x such that it is 2D
x    = reshape(x,oldN,[]);

% crop
if N < oldN
    y = x(1:N,:);
elseif N > oldN
    y = [x; zeros(N-oldN,size(x,2))];
else
    y = x;
end

% reshape back
sz    = oldsz;
sz(1) = N;
y     = reshape(y,sz);
    
end