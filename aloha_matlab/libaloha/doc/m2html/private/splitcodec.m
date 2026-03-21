function splitc = splitcodec(code)
%Split a line of C/C++ code in string, comment and other
%  SPLITC = SPLITCODE(CODE) splits line of C/C++ code CODE into a cell
%  array SPLITC where each element is either a character array ('*/'),
%  a comment (//) or (/*) or something else.
%  Note that CODE = [SPLITC{:}]
%
%  See also M2HTML, HIGHLIGHT

%  Copyright (C) 2003 Guillaume Flandin <Guillaume@artefact.tk>
%  $Revision$Date: 2003/29/04 17:33:43 $
%  modified for C3PO/LUKE by Y. Peysson <CEA/DSM/IRFM, yves.peysson@cea.fr>

%- Label quotes in {'transpose', 'beginstring', 'midstring', 'endstring'}
%
jquote = []; 

iquote = findstr(code,'''');
quotetransp = [double('_''.)}]') ...
			   double('A'):double('Z') ...
			   double('0'):double('9') ...
			   double('a'):double('z')];
flagstring = 0;
flagdoublequote = 0;
jquote1 = [];
for i=1:length(iquote)
	if ~flagstring
		if iquote(i) > 1 & any(quotetransp == double(code(iquote(i)-1)))
			% => 'transpose';
		else
			% => 'beginstring';
			jquote1(size(jquote1,1)+1,:) = [iquote(i) length(code)];
			flagstring = 1;
		end
	else % if flagstring
		if flagdoublequote | ...
		   (iquote(i) < length(code) & strcmp(code(iquote(i)+1),''''))
			% => 'midstring';
			flagdoublequote = ~flagdoublequote;
		else
			% => 'endstring';
			jquote1(size(jquote1,1),2) = iquote(i);
			flagstring = 0;
		end
	end
end

iquote = findstr(code,'"');
quotetransp = [double('_''.)}]') ...
			   double('A'):double('Z') ...
			   double('0'):double('9') ...
			   double('a'):double('z')];
flagstring = 0;
flagdoublequote = 0;
jquote2 = [];
for i=1:length(iquote)
	if ~flagstring
		if iquote(i) > 1 & any(quotetransp == double(code(iquote(i)-1)))
			% => 'transpose';
		else
			% => 'beginstring';
			jquote2(size(jquote,1)+1,:) = [iquote(i) length(code)];
			flagstring = 1;
		end
	else % if flagstring
		if flagdoublequote | ...
		   (iquote(i) < length(code) & strcmp(code(iquote(i)+1),''''))
			% => 'midstring';
			flagdoublequote = ~flagdoublequote;
		else
			% => 'endstring';
			jquote2(size(jquote2,1),2) = iquote(i);
			flagstring = 0;
		end
	end
end


jpercent = [];

%- Find if a portion of code is a comment
ipercent = findstr(code,'//');
jpercent1 = [];
for i=1:length(ipercent)
	if isempty(jquote) | ...
	   ~any((ipercent(i) > jquote(:,1)) & (ipercent(i) < jquote(:,2)))
		jpercent1 = [ipercent(i) length(code)];
        jpercent = jpercent1;
		break;
	end
end


%- Find if a portion of code is a comment
ipercent = findstr(code,'/*');
jpercent2 = [];
for i=1:length(ipercent)
	if isempty(jquote) | ...
	   ~any((ipercent(i) > jquote(:,1)) & (ipercent(i) < jquote(:,2)))
		jpercent2 = [ipercent(i) length(code)];
        if ~isempty(jpercent)
            jpercent = [jpercent;jpercent2];
        else
            jpercent = jpercent2; 
        end
		break;
	end
end

%- Find if a portion of code is a comment
ipercent = findstr(code,'*/');
jpercent3 = [];
for i=1:length(ipercent)
	if isempty(jquote) | ...
	   ~any((ipercent(i) > jquote(:,1)) & (ipercent(i) < jquote(:,2)))
		jpercent3 = [ipercent(i) length(code)];
        if ~isempty(jpercent)
            jpercent = [jpercent;jpercent3];
        else
            jpercent = jpercent3; 
        end
		break;
	end
end

if ~isempty(jquote1), jquote = [jquote;jquote1];end
if ~isempty(jquote2), jquote = [jquote;jquote2];end


if ~isempty(jpercent)
    jpercent = sort(jpercent);
    if size(jpercent,1)>=2,jpercent(1,2) = jpercent(2,1)-1;end
    if size(jpercent,1)>=3,jpercent(2,2) = jpercent(3,1)-1;end
end

%- Remove strings inside comments
if ~isempty(jpercent) & ~isempty(jquote)
	jquote(find(jquote(:,1) > jpercent(1)),:) = [];
end

icode = [];
icode = [jquote ; jpercent];
if ~isempty(icode),
    icode = sort(icode);
end


%- Split code in a cell array of strings

splitc = {};
if isempty(icode)
	splitc{1} = code;
elseif icode(1,1) > 1
	splitc{1} = code(1:icode(1,1)-1);
end
for i=1:size(icode,1)
	splitc{end+1} = code(icode(i,1):icode(i,2));
	if i < size(icode,1) & icode(i+1,1) > icode(i,2) + 1
		splitc{end+1} = code((icode(i,2)+1):(icode(i+1,1)-1));
	elseif i == size(icode,1) & icode(i,2) < length(code)
		splitc{end+1} = code(icode(i,2)+1:end);
	end
end
%keyboard
