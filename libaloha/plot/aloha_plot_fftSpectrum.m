function varargout = aloha_plot_fftSpectrum(varargin)
% plot the FFT-made power spectrum of a field E
% 
% INPUT
%  - f : frequency [Hz]
%  - z : E field abscisse (1xN)
%  - E : E field (1xN)
%  Or
%  - scenario : ALOHA scenario file
%  
%  OUPUT
%  - none
%   or
%  - n_z
%  - Ez_spectral
%   
% AUTHOR: J.Hillairet / V.Fuch
% LAST UPDATE:
%  - 04/2009 : creation

% default optionnal input argument value
bool_normalization = true;

% parsing input arguments
if nargin == 3
    f = varargin{1};
    z = varargin{2};
    E = varargin{3};
elseif nargin == 1 % only a scenario
    sc= varargin{1};
    f = sc.antenna.freq;
    z = sc.results.abs_z(1,:);
    E = sc.results.Efield(1,:);
else
    error('bad number of input arguments ! See help.');
end




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

if bool_normalization
    norma = max(abs(Efft).^2);
else
    norma = 1;
end


h=aloha_plot_figure(figure);
    plot(n_parallel, abs(Efft).^2./norma);
    xlabel('n_{//}');
    ylabel('|FFT[E]|^2');
    set(gca, 'XLim', [-20, 20]);
    grid on;

if (nargout == 2)
    varargout{1} = n_parallel;
    varargout{2} = Efft;
end