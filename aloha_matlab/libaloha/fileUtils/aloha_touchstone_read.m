function[S,freq1]=aloha_touchstone(file, header_lines)
% [S,freq1]=aloha_touchstone(file, header_lines)
% 
% ALOHA 
% 
% Read a Touchstone S-parameters file 
% (.snp files, where n is the number of ports of the RF device)
% 
% The data in the the touchstone file are expected to as magnitude / phase[deg]
% 
% INPUT
%  - file : filename (with path if needed)
%  - header_lines : number of header lines to skip (depend of the Touchstone format)
% 
% OUTPUT
%  - S : scattering matrix(ces)
%  - freq1 : frequency(ies)
% 
% AUTHOR: ? 
% Taken from the matlab file exchange
%
% Typical header lines number:
% cst: header lines = 8
% hfss:header lines = 4 + num_of_ports
% comsol: header lines = 1


%function [S,freq1]=read_snp(header_lines)
%

[PathName,FileName,EXT] = fileparts(file);
%FileName,PathName]=uigetfile('*.s*p','Select valid touchstone file');
FileName = [FileName, EXT];
% deduces the number ports from the extension of the file
[t,r]=strtok(FileName,'.');
[t,r]=strtok(r,'.');
s1=sscanf(t,'%c%d%c');
port_num=s1(2);

fname=[PathName, FileName];
fid=fopen(fname,'r');

if(port_num>4)
    line_entries=4;
    positions_per_last_line=rem(port_num,4);
elseif(port_num==2)
    line_entries=4;
else
    line_entries=port_num;
end

freq_p=0;
tt=1;
port=1;

for uu=1:header_lines
    tline=fgets(fid);
end

while 1
    tline=fgets(fid);
    if(~ischar(tline))
        break;
    end
    A=sscanf(tline,'%f');
    A=A.';
    
    if(length(A)==line_entries*2+1) %data line + freq. point
        port=1;
        freq_p=freq_p+1;
        freq1(freq_p)=A(1);
        
        if(port_num~=2)
            ind1=(tt-1)*line_entries+1;
            ind2=tt*line_entries;
            B1(port,ind1:ind2,freq_p)=A(2:2:end);
            B2(port,ind1:ind2,freq_p)=A(3:2:end);
        else
            ind1=1;
            ind2=1;
            B1(1,1:2,freq_p)=A(2:2:5);
            B2(1,1:2,freq_p)=A(3:2:5);
            B1(2,1:2,freq_p)=A(6:2:end);
            B2(2,1:2,freq_p)=A(7:2:end);
        end
        
    elseif(length(A)==line_entries*2) %full data line
        ind1=(tt-1)*line_entries+1;
        ind2=tt*line_entries;
        B1(port,ind1:ind2,freq_p)=A(1:2:end);
        B2(port,ind1:ind2,freq_p)=A(2:2:end);
    elseif(length(A)==0) %HFSS blank lines !!!
        
         
    else %partial data line
        
         ind1=(tt-1)*line_entries+1;
         ind2=(tt-1)*line_entries+positions_per_last_line;
         B1(port,ind1:ind2,freq_p)=A(1:2:end);
         B2(port,ind1:ind2,freq_p)=A(2:2:end);
    end
    
    if(port_num~=2)
        if(ind2==port_num)
            port=port+1;
            tt=1;
        else
            tt=tt+1;
        end
    end
end

for tt=1:length(freq1)
    X=B1(:,:,tt).*cos(pi/180*B2(:,:,tt));
    Y=B1(:,:,tt).*sin(pi/180*B2(:,:,tt));
    S(:,:,tt)=X+j*Y;    
end
fclose(fid);
    
            
