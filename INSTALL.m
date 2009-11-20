%  INSTALLATION :
%  
%  A priori, le code doit pouvoir tourner apres checkout.
%  
%  Au prealable, le repertoire contenant les librairies
%  d'ALOHA doit etre charge dans le PATH matlab. Pour cela, 
%  on doit effectuer dans matlab : 
%  
%  addpath(genpath([chemin_absolu_du_code_ALOHA_V2/libaloha'])); 
%  
%  Cette commande peut etre rajoutee dans le fichier startup.m de 
%  matlab pour effectuer cette operation a chaque demarrage de matlab
%  (cf help matlab)
disp('You must add the "libaloha" directory to the MATLAB PATH');
disp('In order to do so, please change directory (cd) to the aloha root directory, ');
disp('from where you should see the "libaloha" directory.');
disp('Then enter in the MATLAB prompt the command : "addpath(genpath([pwd, ''/libaloha'']));"');
disp('This add the functions located into the directory "libaloha" (and sub-dir) to the MATLAB PATH');
