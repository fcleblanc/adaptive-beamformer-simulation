% Francis Le Blanc
% Part 1: Fixed LCMV Beamformer
%
% Implements a linearly constrained minimum variance (LCMV) beamformer
% for an 18-element linear antenna array. The beamformer is designed
% to pass a desired signal from a known direction while suppressing
% two strong interference sources and noise.
%
% The script computes the optimal weight vector, evaluates the spatial
% response, and compares the raw vs. beamformed output signals.

% Variables
N = 100000; % number of samples for generating 
M = 18; % number of antennas in array
T = 100; % number of samples for beamformer

des_angle = 8; % angle of desired signal
intf_angle_1 = -5; % angle of first interfering signal
intf_angle_2 = 8.2; % angle of second interfering signal

des_amp = 1; % amplitude of desired signal
intf_amp = 120; % amplitude of interfering signals
noise_amp = 0.06; % amplitude of noise signal

% Load audio signal
load ABF.mat;

% Choose signals
des_sound = voice; % desired signal audio
intf_sound_1 = jetnoise; % first interference signal audio
intf_sound_2 = babble; % second interference signal audio
noise = noise_amp*(randn(M,N) + 1j*randn(M,N)); % random noise signal

% Take N samples of each sound as a column vector 
des_samples = des_sound(1:N).';
intf_1_samples = intf_sound_1(1:N).';
intf_2_samples = intf_sound_2(1:N).';

% Create steering vectors
a_des = linarr(M,des_angle);
a_intf_1 = linarr(M, intf_angle_1);
a_intf_2 = linarr(M, intf_angle_2);

% Create total signal
X = des_amp*(a_des*des_samples) + intf_amp*(a_intf_1*intf_1_samples) + intf_amp*(a_intf_2*intf_2_samples) + noise;

% Take first 100 samples of total signal
X_samples = X(:, 1:T); 

% Compute beamformer weights
% Minimize output power subject to constraints A^H g = c
Rxx = (X_samples*X_samples') / T;
A = [a_des, a_intf_1, a_intf_2];
c = [1;0;0]; % enforce as a column vector
Ri = inv(Rxx);
g = Ri*A*inv(A'*Ri*A)*c;
h = conj(g);

% Compute spatial frequency response
[H,w] = dtft(h, 1024);
f = (180/pi)*asin(w/pi);

% Plot spatial frequency response
plot(f, 20*log10(abs(H)));
grid on;
xlabel('Angle (degrees)');
ylabel('|H(e^{j\omega_s})| (dB)');
title('Beamformer Spatial Response');
hold on;

% Mark desired/interferer angles
yL = ylim;
plot([des_angle des_angle], yL, 'k--');
plot([intf_angle_1 intf_angle_1], yL, 'r--');
plot([intf_angle_2 intf_angle_2], yL, 'r--');
hold off;

beam_output = h.' * X;     

sound(real(X(1,:)), 22000);    
pause(5);
sound(real(beam_output), 22000); 

function [A,D,D2,W,T]=linarr(m,ang,sep)
%function [A,D,D2,W,T]=linarr(m,ang,sep)
if nargin < 3
sep=0.5;
end
d=length(ang);
rang=sin(ang(:)*pi/180);
cang=cos(ang(:)*pi/180);
j=sqrt(-1);
r=[0:m-1]';
T=2*pi*r*sep*rang';
A=exp(-j*T);
if nargout > 1
W=2*pi*r*sep*cang';
D=-j*W.*A;
D2=(j*T-W.*W).*A;
end
end

function [H,W] = dtft( h, N )
%DTFT calculate DTFT at N equally spaced frequencies
% usage: H = dtft( h, N )
% h: finite-length input vector, whose length is L
% N: number of frequencies for evaluation over [-pi,pi)
% ==> constraint: N >= L
% if N not included, it is assumed to be equal to length of h
% H: DTFT values (complex)
% W: (2nd output) vector of freqs where DTFT is computed
if nargin < 2
N=length(h);
end
N = fix(N);
L = length(h); h = h(:); %<-- for vectors ONLY !!!
if( N < L )
error('DTFT: # data samples cannot exceed # freq samples')
end
W = (2*pi/N) * [ 0:(N-1) ]';
mid = ceil(N/2) + 1;
W(mid:N) = W(mid:N) - 2*pi; % move [pi,2pi) to [-pi,0)
W = fftshift(W);
H = fftshift( fft( h, N ) ); % move negative freq components
end