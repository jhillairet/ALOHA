function aloha_plot_export(H, figname)
% export a figure into a PDF image file
% for publication purpose
% 
% INPUT
%  h : image Hr
%  figname : figure filename
%  
% OUTPUT
%  none
%  
% AUTHOR: JH
% LAST UPDATE: 
%  - 07/05/2009 : creation

% Backup previous settings
prePaperType = get(H,'PaperType');
prePaperUnits = get(H,'PaperUnits');
preUnits = get(H,'Units');
prePaperPosition = get(H,'PaperPosition');
prePaperSize = get(H,'PaperSize');

% Make changing paper type possible
set(H,'PaperType','<custom>');

% Set units to all be the same
set(H,'PaperUnits','centimeters');
set(H,'Units','centimeters');

% Set the page size and position to match the figure's dimensions
paperPosition = get(H,'PaperPosition');
position = get(H,'Position');
set(H,'PaperPosition',[0,0,position(3:4)]);
set(H,'PaperSize',position(3:4));

% Save the pdf (this is the same method used by "saveas")
%  print(H,'-dpdf',pdfFileName,sprintf('-r%d',dpi))
saveas(H, figname);

% Restore the previous settings
set(H,'PaperType',prePaperType);
set(H,'PaperUnits',prePaperUnits);
set(H,'Units',preUnits);
set(H,'PaperPosition',prePaperPosition);
set(H,'PaperSize',prePaperSize);
