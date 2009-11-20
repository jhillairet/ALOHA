function s = what_aloha_yp(folder)
%LUKE - Extension of WHAT MatLab built-in command
%
% This function extends the built-in MatLab command WHAT for listing *.c, *.h
% and *.f type files. The recognized extensions are
%
% - C/C++ type files: *.c, *.cpp, *.c++, *.cxx, *.C
% - Header type files *.h
% - Fortran type files: .for, *.f, *.FOR, *.F
%
% Useful for the automatic documentation builder.
%
% by Y.Peysson <CEA/DSM/IRFM, yves.peysson@cea.fr> J. Decker <CEA/DSM/IRFM, joan.decker@cea.fr>
%
dirinit = pwd;%initial directory
%
s = what(folder);
%
cd(folder);
%
% C/C++ type files (*.c, *.cpp, *.c++, *.cxx, *.C)
%
c_files = struct2cell(dir('*.c'));
cpp_files = struct2cell(dir('*.cpp'));
cplusplus_files = struct2cell(dir('*.c++'));
cxx_files = struct2cell(dir('*.cxx'));
C_files = struct2cell(dir('*.C'));
%
% Header type files (*.h)
%
h_files = struct2cell(dir('*.h'));
%
% Fortran type files (.for, *.f, *.FOR, *.F)
%
for_files = struct2cell(dir('*.for'));
f_files = struct2cell(dir('*.f'));
FOR_files = struct2cell(dir('*.FOR'));
F_files = struct2cell(dir('*.F'));
%
for is = 1:size(s,1),%Case insensitive system (mac OSX, Windows,...)
    s(is).c = {};
    s(is).f = {};
    s(is).h = {};
    %
    ic = 1;
    for ii = 1:size(c_files,2),
        s(is).c{ic,1} = c_files{1,ii};
        ic = ic + 1;
    end
    %
    for ii = 1:size(cpp_files,2),
        s(is).c{ic,1} = cpp_files{1,ii};
        ic = ic + 1;
    end
    %
    for ii = 1:size(cplusplus_files,2),
        s(is).c{ic,1} = cplusplus_files{1,ii};  
        ic = ic + 1;
    end
    %
    for ii = 1:size(cxx_files,2),
        s(is).c{ic,1} = cxx_files{1,ii};  
        ic = ic + 1;
    end
    %
    for ii = 1:size(C_files,2),
        s(is).c{ic,1} = C_files{1,ii};   
        ic = ic + 1;
    end
    %
    ifo = 1;
    for ii = 1:size(for_files,2),
        s(is).f{ifo,1} = for_files{1,ii};     
        ifo = ifo + 1;
    end
    %
    for ii = 1:size(f_files,2),
        s(is).f{ifo,1} = f_files{1,ii};     
        ifo = ifo + 1;
    end
    %
    for ii = 1:size(FOR_files,2),
        s(is).f{ifo,1} = FOR_files{1,ii};     
        ifo = ifo + 1;
    end
    %
    for ii = 1:size(F_files,2),
        s(is).f{ifo,1} = F_files{1,ii};     
        ifo = ifo + 1;
    end
    %
    for ih = 1:size(h_files,2),
        s(is).h{ih,1} = h_files{1,ih};     
    end
end
%
cd(dirinit);



