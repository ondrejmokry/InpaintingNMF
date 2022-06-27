function name = fname(folder)
% FNAME sets the filename for saving the results based on the name of the
% script fname is called in.
%
% It returns string formed from the script name and a numeric identifier
% such that no results are overwritten.
%
% Name also includes the folder path 'results/' as a prefix, if not given
% otherwise.
%
% Date: 08/12/2021
% By Ondrej Mokry
% Brno University of Technology
% Contact: ondrej.mokry@mensa.cz

% get the function call stack
stack = dbstack;

% get folder name
if nargin < 1
    folder = 'results/';
end

% read the second one (the first one is fname) and put a prefix
name = [folder, stack(2).name, '_01.mat'];

% iterate to determine the numeric identifier
i = 1;
while isfile(['./', name])
    i = i + 1;
    name = [name(1:end-6), num2str(i,'%02d'), '.mat'];
end

