function search_path = htmldocbuilder_aloha_yp(m2htmllib_path,recursive_mode,root_path,search_path,excludedfoldernames,overview_scripts)
% ALOHA - HTML documentation builder
%
% Create dynamically all the documentation in html format from *.m, *.c,
% *.f and *.h text files. Recognize binary files. If a TODO comment is
% written somewhere, it is considered in the links.
% Creates also graphes (needs graphviz graph maker which can be downloaded at www.graphviz.org)
%
% INPUT:
%
%   - m2htmllib_path: full path to m2html library (string)
%   - recursive_mode: automatic recursive scan of directories (except directories whose names start with '.' and 'Doc'). 'on' or 'off' (string)
%       (default: 'off')
%   - root_path: full path of the directory to be scanned. (string)
%       (default: local directory)
%   - search_path: list subdirectories to be scanned (relative path to root_path) (cells of strings)
%       (default: all subdirectories are scanned recursively)
%   - excludedfoldernames: list excluded folder names is the build-in search_path (cells of strings)
%       (default: empty)
%   - overview_scripts: list of scripts to be put in the overview html file (structure)
%       (default: none)
%
% OUTPUT:
%
%   - search_path: list subdirectories scanned (cells of strings)
%
% WARNING: 
%
%   - based on a modified version of m2html MatLab library.
%   - overview_scripts.path: absolute path
%   - overview_scripts.format = 'doc', 'html', 'latex', 'ppt', 'xml' (no field or empty for m-files)
%
% by Y. PEYSSON (CEA/DSM/IRFM, yves.peysson@cea.fr) and J. DECKER (CEA/DSM/IRFM, joan.decker@cea.fr)
%
    if nargin == 0,
        error('The modified m2html library absolute path must be specified as input parameter for htmldocbuilder_ALOHA_yp !');
    end
    if nargin == 1,
        recursive_mode = 'on';
        root_path = '';
    end   
    if nargin == 2,
        root_path = '';
    end
    if nargin <=3,
        search_path = {};
        excludedfoldernames = {};
    end
    %
    if isempty(root_path),
        root_path = [pwd,filesep];%Working in the local directories where htmldocbuilder_aloha_yp is launched
    end
    %
    if ~strcmp(root_path(end),filesep),
        root_path = [root_path,filesep];%root_path must be ended by a filesep
    end
    %
    docpath = root_path;
    dochtmlpath = [root_path,'Documentation'];
    if ~exist(dochtmlpath),
        mkdir(docpath,'Documentation');
    else
        rmdir(dochtmlpath,'s');
        mkdir(docpath,'Documentation');
    end
    %
    addpath(m2htmllib_path,'-BEGIN');
    %
    dirinit = pwd;%initial directory
    %
    if isempty(search_path),
        search_path = dirscan_yp(root_path,search_path,excludedfoldernames,recursive_mode);%Search_path set-up
    end
    %
    pathstr = fileparts(root_path);
    ipathstr = find(pathstr==filesep,1,'last');
    %
    for ii = 1:length(search_path),
        search_path{ii} = search_path{ii}(ipathstr+1:end);%remove the absolute path part up to the scanned folder
    end
    %
    cd(root_path);
    cd ..;%for m2html needs
    %
    % todo: initial description in documentation
    %
    try
        m2html('mfiles',search_path, 'htmldir',dochtmlpath,'recursive','off','global','on', 'template','ALOHA','index','menu','graph','on','todo','on');
    catch
        m2html('mfiles',search_path, 'htmldir',dochtmlpath,'recursive','off', 'global','on', 'template','ALOHA','index','menu','graph','off','todo','on');%if graphviz not installed
        %
        disp('WARNING: Graphes have not been generated, since the ''graphviz'' maker is not installed on your computer.');
        disp('WARNING: Download ''graphviz'' maker from www.graphviz.org, and install it first.');
        disp('WARNING: Execute htmldocbuilder_aloha_yp MatLab script do rebuild the documentation.');
    end    
    %
    rmpath(m2htmllib_path);
    %
    mkdir(fullfile(root_path,'Documentation'),'OVERVIEW_files');
    %
    opts.outputDir = fullfile(root_path,'Documentation','OVERVIEW_files');
    opts.evalCode = false;
    %
    if exist('overview_scripts'),
        for ii = 1:length(overview_scripts)
            if ~isempty(overview_scripts(ii).format) & (strcmp(overview_scripts(ii).format,'html') | strcmp(overview_scripts(ii).format,'doc') | strcmp(overview_scripts(ii).format,'ppt') | strcmp(overview_scripts(ii).format,'xml') | strcmp(overview_scripts(ii).format,'rpt') | strcmp(overview_scripts(ii).format,'latex')),            
                opts.format = overview_scripts(ii).format;
            end
            %
            publish(overview_scripts(ii).path,opts);
            %
            if isfield(opts, 'format'),
                opts = rmfield(opts, 'format');
            end
        end
    end
    %
    stat = web(['file://',root_path,'Documentation/index.html'],'-browser');
    if stat == 1 | stat == 2,
        web(['file://',root_path,'Documentation/index.html'],'-helpbrowser');%Open in MatLab help browser    
    end
    disp(['----> The weblink to the documentation is ','file:/',root_path,'Documentation/index.html']);
    %
    cd(dirinit);%back to the initial directory
    %
end



function [search_path] = dirscan_yp(root_path,search_path,excludedfoldernames,recursive_mode)
%ALOHA - Recursive directory scan
%
% Recursive scan of directory tree starting from scanfolder level
%
% INPUT:
%
%   - root_path: absolute directory path name to be scanned (string)
%   - search_path: list of subdirectories (cell of strings)
%   - excludedfoldernames : list of folder name that must ne be scanned (cell of strings)
%   - recursive_mode : 'on' or 'off' (string)
%       (default: 'off')
%
% OUTPUT:
%
%   - search_path: structure of subdirectories scanned (not equal to input when recursive option is 'on') (structure)
%
% by Y. PEYSSON (CEA/DSM/IRFM, yves.peysson@cea.fr) and J. DECKER (CEA/DSM/IRFM, joan.decker@cea.fr)
%
    if nargin == 0,
        error('The root_path of the folder to be scanned must be given as input dirscan_yp.m !');
    end
    if nargin == 1,
        search_path = {};
        excludedfoldernames = {};
        recursive_mode = 'off';
    end
    if nargin == 2,
        excludedfoldernames = {};
        recursive_mode = 'off';
    end
    if nargin == 3,
        recursive_mode = 'off';
    end
    %
    dirinit = pwd;%record the active directory name
    %
    if ~strcmp(root_path(end),filesep),
        root_path = [root_path,filesep];%root_path must be ended by a filesep
    end

    cd(root_path);%move to the root folder to be scanned
    dirlist = dir;
    nsearch_path0 = length(search_path);
    nsearch_path = nsearch_path0;
    %
    for ii = 1:length(dirlist),
        if dirlist(ii).isdir == 1 & ~strcmp(dirlist(ii).name(1),'.') & ~strcmp(dirlist(ii).name,'CVS') & ~strcmp(dirlist(ii).name,'Doc') & ~strcmp(dirlist(ii).name,'Documentation'),
            search_path = [search_path,[root_path,dirlist(ii).name]];
            nsearch_path = nsearch_path + 1;%a valid new directory has been found
        end
    end 
    %
    for ii = nsearch_path0+1:nsearch_path,
        if strcmp(recursive_mode,'on'),
            [search_path] = dirscan_yp(search_path{ii},search_path,excludedfoldernames,recursive_mode);%recursive scan
        end
    end
    %
    cd(dirinit);%back to the initial active directory
end


