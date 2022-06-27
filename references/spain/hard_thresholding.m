function s = hard_thresholding(a, k)

% Date: 08/07/2020
% By Ondrej Mokry, Pavel Zaviska
% Brno University of Technology
% Contact: ondrej.mokry@mensa.cz

odd = mod(length(a),2);

% taking only half of the spectrum + dc coefficient
a = a(1:floor(length(a)/2)+1);
a(1) = a(1)/2;

% sorting them
[~, ind] = sort(abs(a), 'descend');
s = zeros(length(a),1);
if k < length(a)
    s(ind(1:k)) = a(ind(1:k));
else
    s = a;
end

% compute conjugates of selected coefficients
s(1) = s(1)*2;

if odd
    s_conj = conj(flip(s(2:end)));
else
    s_conj = conj(flip(s(2:end-1)));
end

s = [s; s_conj];

end 