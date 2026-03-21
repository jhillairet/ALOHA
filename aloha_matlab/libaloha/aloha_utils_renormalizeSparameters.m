function S_ren=aloha_utils_renormalizeSparameters(S, Zref, Zref_new)
%% ALOHA
% Re-normalize a generalized S-matrix to a reference impedance 
% 
% INPUT
%  - S (NxN) : S-matrix to renormalize
%  - Zref (1x1) or (1xN) : reference impedance of the original S-matrix
%  - Zref_new (1x1) or (1xN) : new reference impedance for the renormalized S-matrix
% 
% OUTPUT
%  - S_ren (NxN) : renormalized S-matrix
% 
% Author: JH
% Date: July 2011. Tested OK with both HFSS and Marks1992 formulae.

%% Check the size of the reference impedances. 
% if scalar -> transform into vector
% The length of the vector is the length of the S-matrix 
% (ie the number of ports)
if length(Zref) == 1
    Zref = repmat(Zref, 1, length(S));
end

if length(Zref_new) == 1
    Zref_new = repmat(Zref_new, 1, length(S));
end

Id = eye(length(S)); % identity matrix

%  %% Method 1 (OK)
%  % Ref: HFSS documentation (Renormalized S-Matrices)
%  % Generates impedance matrix
%  sqrt_z = diag(sqrt(Zref));
%  Z = sqrt_z*(Id+S)*inv(Id-S)*sqrt_z;
%  % Tranform the Z-matrix into a new renormalized S-matrix
%  S_ren = inv(diag(sqrt(Zref_new)))*(Z-diag(Zref_new))*inv(Z+diag(Zref_new))*diag(sqrt(Zref_new));

%% Method 2 (also OK). Should this one works also with complex valued impedances??
% Ref: Marks1992 
% creating impedance matrix
U = diag(sqrt(real(Zref))./abs(Zref));
U_new = diag(sqrt(real(Zref_new))./abs(Zref_new));
% Transform the input S-matrix into a Z-matrix
Z = inv(Id - inv(U)*S*U) * (Id + inv(U)*S*U) * diag(Zref);
% Ref: Marks1992 or Wikipedia http://en.wikipedia.org/wiki/Impedance_parameters
S_ren = U*(Z-diag(Zref_new))*inv(Z+diag(Zref_new))*inv(U);