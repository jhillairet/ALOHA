%  INSTALLATION :
%  
%  Normally, the code should be able to run after the SVN checkout. 
%  Since the fortran code part used by matlab depends of the
%  computer architecture, some all-ready made
%  binaries exist for i686 32 and 64 bits (Linux). 
%  It is however possible that you need to compile the code on your
%  platform.
% 
%  Before using ALOHA, the directory which contain all the librairies
%  needed by the code has to be included into the Matlab PATH. In order 
%  to do so, one must type into matlab prompt the following command :
%    
%  addpath(genpath([chemin_absolu_du_code_ALOHA_V2/libaloha'])); 
%  
%  This command can be added to the matlab file startup.m
%  in order to automatically run this command at the matlab startup.
%  (cf help matlab)
disp('You must add the "libaloha" directory to the MATLAB PATH');
disp('In order to do so, please change directory (cd) to the aloha root directory, ');
disp('from where you should see the "libaloha" directory.');
disp('Then enter in the MATLAB prompt the command : "addpath(genpath([pwd, ''/libaloha'']));"');
disp('This add the functions located into the directory "libaloha" (and sub-dir) to the MATLAB PATH');
