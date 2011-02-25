function aloha_util_runSWAN(scenario, output_filename)
% Run the previous code SWAN
% 
% Warning : it is assumed that the SWAN binary file is located in the
% current directory
% 
% Warning2 : 

SWAN = 'swan95_deneb';

%% write the SWAN input file into a text file
SWAN_input = aloha_scenario_conversion4swan(scenario);

input_filename = 'ALOHA0'; %filename must 6 char long only !
fid=fopen(input_filename, 'w');
for idx=1:size(SWAN_input,1)
    fprintf(fid,'%s\n', SWAN_input(idx,:));
end
fclose(fid);

%% run SWAN in command line
cmd_line = ['./', SWAN, ' <', input_filename, ' >', output_filename];
[s,w] = system(cmd_line);

%% if it worked
if s == 0
    fid = fopen(output_filename, 'r');
    while feof(fid) == 0
        line = fgetl(fid);
        if ~isempty(findstr(line,'COEFFICIENT DE REFLEXION GLOBAL'))
            disp(line);
        end
    end
    fclose(fid);
end
