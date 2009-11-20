function vargout=aloha_plot_figure(varargin)
%  Set some properties of a figure in order to have a better exported result. 
%  
%  It consists in : 
%   - increasing the font size (specially usefull for powerpoint slides) 
%   - crop the paper to the figure's size, in order to avoid extra crop operation
%     after image exporting.
%   - increasing the line width to 2
%     
%  INPUT
%   - h : figure handle
%   - figure_name (string) : the name of the figure [optionnal]
%  
%  OUPUT
%   - h [optionnal]: figure handle
% 
%  EXAMPLE
%   aloha_plot_figure(h)
%    Or,
%   aloha_plot_figure(h, 'My figure')
%  
%  AUTHOR:JH
%  LAST UPDATE:
%   - 24/04/2009: now also accept no input argument
%   - 12/11/2008: set page size to A4
%   - 06/11/2008: rename to aloha_plot_figure, set defaut values for line width & font size, return the handle
%   - 22/10/2008: added the name of the figure input (optionnal)
%   - 13/10/2008: creation
%   
    switch nargin
      case 0
        h = figure;
      case 1
        if ishandle(varargin{1})
            h = varargin{1};
        elseif isint(varargin{1});
            h = figure(varargin{1});
        else
            error('Bad argument type. See help.');
        end
      case 2
        h = varargin{1};
        set(h, 'Name', ['-= ', varargin{2}, ' =-']);
      otherwise
        error('Wrong input argument number !');
    end
    
    % enlarge image size
    Pos=get(h, 'Position'); 
    set(h, 'Position', [Pos(1) Pos(2) 860 Pos(4)])

    %% usefull for slides
    % set default value for larger fontsize for better visibilty 
    set(h, 'DefaultAxesFontSize', 15);
    set(h, 'DefaultTextFontSize', 15);
    % set defaut value for larger line width 
    set(h, 'DefaultLineLineWidth', 2);
%      set(get(gca, 'YLabel'), 'FontSize', 15)
%      % adapt the paper to the size of the axes
%      % in order to have a image without white spaces on each sides
%      set(gca, 'Units', 'centimeters')
%      pos=get(gca, 'position');
%      set(gcf, 'PaperUnits', 'centimeters');
%      posTitle=get(get(gca, 'Title'), 'Position');
%      paperSize = [pos(1)+pos(3)+1, pos(2)+pos(4)+posTitle(2)];
%      set(gcf, 'PaperSize', paperSize);% +x cm for title
%      set(gcf, 'PaperUnits', 'centimeters') ;
%      set(gcf, 'PaperPosition', [0, 0, paperSize(1), paperSize(2)]);

    % ouput argument
    if nargout == 1
        vargout(1)=h;
    end