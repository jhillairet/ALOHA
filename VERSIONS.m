%  ALOHA (MAJORS) CHANGE-LOG
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