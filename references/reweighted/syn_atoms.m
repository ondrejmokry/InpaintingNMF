function [ atoms ] = syn_atoms( F )
% auxiliary function for atom weighting
%
% input:    F ...... frame, F = frame('type',{'window',w},a,M)
% output:   atoms .. matrix having "relevant" atoms as columns
%
% Date: 08/07/2020
% By Ondrej Mokry
% Brno University of Technology
% Contact: ondrej.mokry@mensa.cz

%% frame parameters
% m ... number of atoms per one shift
type = F.type;
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
w = length(F.g);
L = framelength(F,w);
Ncoef = frameclength(F,L);

%% atoms synthesis
AT = eye(Ncoef,m);
atoms = frsyn(F,AT);