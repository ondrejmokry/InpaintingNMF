function proj = proj_time(x, mask, data_gapped)
% PROJ_TIME perfoms projection of vector x onto the set of feasible
% solutions for the inpainting problem in time domain.
%
% Input parameters
%       x               vector of input signal
%       mask            logical vector indicating the reliable samples
%       data_gapped     degraded signal
%
% Date: 08/07/2020
% By Ondrej Mokry, Pavel Zaviska
% Brno University of Technology
% Contact: ondrej.mokry@mensa.cz

proj = x;
proj(mask) = data_gapped(mask);

end