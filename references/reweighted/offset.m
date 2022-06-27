function [ off ] = offset(s,f,a,type)
% OFFSET computes the value of offset parameter for the
% function min_sig_supp.m (or, preferably, min_sig_supp_2.m) in a way that
% the processing of the signal window by window is done symmetrically with
% respect to the center of the gap
%
% input:
%   s ..... index of the first missing sample
%   f ..... index of the last missing sample
%   type .. 'full' / 'half' / 'none'
%
% output:
%   off ... value of offset for min_sig_supp.m / min_sig_supp_2.m
%
% Date: 08/07/2020
% By Ondrej Mokry
% Brno University of Technology
% Contact: ondrej.mokry@mensa.cz

switch lower(type)
    case 'full' % middle of the gap coincides with central index of a window
        c = ceil((s+f)/2);
        k = floor((c-1)/a);
        d = 1 + k*a ;
        off = c-d;
    case 'half' % middle of the gap coincides with symmetry axis of two adjacent windows
        c = ceil((s+f)/2);
        k = floor((c-1)/a);
        d = 1 + k*a + ceil(a/2);
        off = c-d;
    otherwise % no offset
        off = 0;
end

end

