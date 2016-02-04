%  ALOHA (MAJORS) CHANGE-LOG
% 01/2016 : Some work on ALOHA-2D in order to allow any (nz, ny) spectral domain range.  
%           Add calculation of the spectrum and electrical fields.
%
% 12/2011 : Added vacuum gap support with two linear density gradients (version=6). 
%           Corrected bugs for poloidal density mismatch support.
%  
% 12/2011 : Corrected bugs about the ITM antenna description (with on-demand compatibility)
%          and in the EM field evaluation at the antenna mouth.
%  
% 10/2011 : the default ALOHA usage is now to create and process scenarios. The script file aloha.m has thus
%           been removed from SVN. A tutorial directory has been added to introduce the ALOHA usage.  
%
%  06/2011 : the antenna description has been changed to be compliant with ITM LH antenna description 
%           In order to keep the compatibility with old scenarios, the previous antenna description should
%           still work.
%           The code has been also improved to read the Tore Supra database (see the aloha_ondemand... functions)
%
%  02/2011 : The antenna description now follows the ITM standard (antenna_lh CPO). Moreover, the 
%           antenna parameters, such as the waveguide dimensions, are now stored into the scenarios.
%           This allows batch calculation on antenna parameters.   
%
%  05/2010 : the code repository has been moved from CVS to the partenaires zone SVN server.
%           Many little changes have made since last year : mostly improvement of usability or
%           performance improvements (especially for the 2 gradients version). 
%  
%  05/2009 : The C3 Antenna description has been modified in order to take into account its real
%  shape. In particular, the septum width for passive waveguide (between modules) is different than
%  the septum width between all other waveguide.
%  
%  04/2009 : Directivity and Electric fields now works are generated with function
%  calls (input : 'scenario', output : 'scenario'). Idem for plotting directivity 
%  and fields. Electric fields routine had been checked when I was in IPP Prague. As required by 
%  Vladimir Fuch, the electric field computed directly in the spatial domain has now the same spacing
%  between 2 z points. Also, there is a new routine for compute the spectrum using the electric field.
%  TODO : the same for ALOHA2D
%  
%  02/2009 : Spectrum is computed and plotted using some aloha function. This feature ables one 
%  to not compute the spectrum when he makes a aloha run, and compute the spectrum later if needed.
%  
%  01/2009 : ALOHA 2D and ALOHA 1D are now the same code. A new option labelled
%          'version_code' (values : '1D' or '2D') will launch either ALOHA 1D or ALOHA 2D
%          
%          
