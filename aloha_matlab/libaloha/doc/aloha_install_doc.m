function aloha_install_doc(test)
% ALOHA - Installation script for automatic generation of the on-line documentation 
%
% Installation script for automatic generation of the on-line documentation
% (pdflatex and tex4html required for a complete build of the documentation)
%
% INPUT
%  -test: dummy argument. If exists, test mode is activated with reduced documentation
%
% OUTPUT: none
%
% AUTHORS:
%  - Joan Decker and Yves Peysson for the LUKE code
%  - adapted by J.Hillairet for ALOHA
%
dirinit = pwd;

% TODO : faire une fonction pour recuperer de maniere solide
% un path qui va bien !
aloha_root = aloha_utils_getRootPath; 

cd(aloha_root);

m2htmllib_path = fullfile(aloha_root,'libaloha','doc','m2html');
if nargin == 0,
    search_path = {...
                aloha_root, ...
                fullfile(aloha_root,'libaloha'), ...
                fullfile(aloha_root,'libaloha','scenario'), ...
                fullfile(aloha_root,'libaloha','file'), ...
                fullfile(aloha_root,'libaloha','excitation'), ...
                fullfile(aloha_root,'libaloha','itm'), ...
                fullfile(aloha_root,'libaloha','plot'), ...        
                fullfile(aloha_root,'libaloha','doc'), ...  
              };   
else
    disp('----------------------------------------------')
    disp('WARNING: Test mode.')
    disp('----------------------------------------------')
    %
    search_path = {fullfile(aloha_root,'libaloha'), ...
              };
end
%
overview_scripts(1).path = fullfile(aloha_root,'INSTALL.m');overview_scripts(1).format = '';
overview_scripts(2).path = fullfile(aloha_root,'README.m');overview_scripts(2).format = '';
overview_scripts(3).path = fullfile(aloha_root,'VERSIONS.m');overview_scripts(3).format = '';
overview_scripts(4).path = fullfile(aloha_root,'AUTHORS.m');overview_scripts(4).format = '';
%
excludedfoldernames = {};
htmldocbuilder_ALOHA_yp(m2htmllib_path,'off',aloha_root,search_path,excludedfoldernames,overview_scripts);%build html doc

%
% Check that the list of TeX files given below is similar in the html file ~','LUKE/Packages/m2html/templates/LUKE/lukeoverview.html
%  
if nargin == 0,
%      texfilelist = {fullfile(aloha_root,'Project_DKE','Doc','NoticeDKE.tex'),...      
%              fullfile(aloha_root,'Project_DKE','Modules','RT','Doc','raytracing.tex'),...
%              fullfile(aloha_root,'Project_DKE','Modules','FEB','Doc','3Tmodel.tex'),...
%             }; 
else
    disp('----------------------------------------------')
    disp('WARNING: Test mode.')
    disp('----------------------------------------------')
    %
%      texfilelist = {fullfile(aloha_root,'Project_DKE','Modules','FEB','Doc','3Tmodel.tex'),...
%             }; 
end
%  pdfhtlatex_yp(aloha_root,texfilelist);%build html and pdf doc from tex files
%
cd(dirinit);
