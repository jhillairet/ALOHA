function [ext, platform] = mexexts
%List of Mex files extensions
%  MEXEXTS returns a cell array containing the Mex files platform
%  dependent extensions and another cell array containing the full names
%  of the corresponding platforms.
%
%  See also MEX, MEXEXT

%  Copyright (C) 2003 Guillaume Flandin <Guillaume@artefact.tk>
%  $Revision$Date: 2003/29/04 17:33:43 $

ext = {'.mexglx' '.mexa64' '.mexmac' '.mexmaci' '.mexs64' '.mexw32' '.mexw64'};
 
platform = {'Linux (32-bit)' 'Linux x86-64' 'Macintosh (PPC)' 'Macintosh (Intel)' '64-bit Solaris SPARC' 'Windows (32-bit)' 'Windows x64'};
