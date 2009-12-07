function varargout = aloha_plot_fftSpectrum(f,z,E)
% plot the FFT-made power spectrum of a field E
% 
% INPUT
%  f : frequency [Hz]
%  z : E field abscisse (1xN)
%  E : E field (1xN)
%  
%  OUPUT
%   none
%   or
%   [n_z, Ez_spectral]
%   
% AUTHOR: J.Hillairet / V.Fuch
% LAST UPDATE:
%  - 04/2009 : creation
 

% load usefull physical constants
aloha_constants;

k0=2*pi*f/celerite;
lambda=celerite/f;

% E field fft
Efft = fftshift(fft(E));

% shifted spatial frequency index
B = length(Efft);
K = -round(B/2) + [1:B];

% parallel index
n_parallel = lambda*K./(abs(max(z)-min(z)));


h=aloha_plot_figure(figure);
    plot(n_parallel, abs(Efft).^2./max(abs(Efft).^2));
    xlabel('n_{//}');
    ylabel('|FFT[E]|^2');
    set(gca, 'XLim', [-20, 20]);
    grid on;

if (nargout == 2)
    varargout{1} = n_parallel;
    varargout{2} = Efft;
end