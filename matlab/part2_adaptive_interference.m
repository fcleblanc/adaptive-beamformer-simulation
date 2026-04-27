% Francis Le Blanc
% Part 2: Adaptive Beamformer with Emerging Interference
%
% Implements an adaptive beamformer for a linear antenna array using
% an exponentially weighted covariance estimate. The beamformer updates
% its weights over time to respond to changes in the signal environment.
%
% A third interference source is introduced partway through the simulation.
% The algorithm adapts to this new interferer, and the evolution of the
% beamformer response is analyzed at selected time steps.

% Variables
N = 100000; % number of samples for generating 
M = 18; % number of anetennas in array
T = 100; % number of samples for beamformer
e = 0.003; % forgetting factor
intf_start = 50000; % new interference signal introduced at N = 50000

des_angle = 8; % angle of desired signal
intf_angle_1 = -5; % angle of first interfering signal
intf_angle_2 = 25; % angle of second interfering signal
intf_angle_3 = -15; % angle of third interfering signal

des_amp = 1; % amplitude of desired signal
intf_amp = 120; % amplitude of interfering signals
noise_amp = 0.06; % amplitude of noise signal

% Create envelope for new interference signal
intf_env = zeros(1,N); 
intf_env(intf_start:N) = 1;

y = zeros(1,N); % beamformer output 
g_values = struct(); % values of g before and after interferer signal
n_values = [49000, 49250, 49500, 49750, 50000, 50010, 50020, 50030, 50040, 50050, 50100, 50500, 51000];  


% Load audio signal
load ABF.mat;

% Choose signals
des_sound = voice; % desired signal audio
intf_sound_1 = jetnoise; % first interference signal audio
intf_sound_2 = babble; % second interference signal audio
intf_sound_3 = seasons; % third interference signal audio
noise = noise_amp*(randn(M,N) + 1j*randn(M,N)); % random noise signal

% Take N samples of each sound as a column vector 
des_samples = des_sound(1:N);    
des_samples = des_samples(:).';
intf_1_samples = intf_sound_1(1:N); 
intf_1_samples = intf_1_samples(:).';
intf_2_samples = intf_sound_2(1:N); 
intf_2_samples = intf_2_samples(:).';
intf_3_samples = intf_sound_3(1:N); 
intf_3_samples = intf_3_samples(:).';

% Create steering vectors
a_des = linarr(M,des_angle);
a_intf_1 = linarr(M, intf_angle_1);
a_intf_2 = linarr(M, intf_angle_2);
a_intf_3 = linarr(M, intf_angle_3);

% Create total signal
X = des_amp*(a_des*des_samples) + intf_amp*(a_intf_1*intf_1_samples)  ...
                                + intf_amp*(a_intf_2*intf_2_samples) ...
                                + intf_amp*(a_intf_3*(intf_env.*intf_3_samples)) ...
                                + noise;

% Take first 100 samples of total signal
X_samples = X(:, 1:T); 

% Construct vectors for beamformer
Rxx = (X(:,1:T) * X(:,1:T)') / T;
A = [a_des, a_intf_1, a_intf_2];
c = [1;0;0]; % enforce as a column vector

% Update Rxx adaptively
for n = T+1:N
    x = X(:,n);
    Rxx = (1-e)*Rxx + e*(x*x');

    RiA = Rxx \ A;
    g = RiA * ((A' * RiA) \ c);
    h = conj(g);

    y(n) = h.' * x;

    if any(n == n_values)
        g_values.(sprintf('n%d',n)) = h;
    end
end

% Plot spatial frequency responses at n samples
fn = fieldnames(g_values);

% Extract and sort by n
n_list = zeros(size(fn));
for k = 1:numel(fn)
    n_list(k) = sscanf(fn{k}, 'n%d');
end
[~, idx] = sort(n_list);
fn = fn(idx);
n_list = n_list(idx);

for k = 1:numel(fn)
    hk = g_values.(fn{k});
    [H,w] = dtft(hk, 1024);
    angle = (180/pi)*asin(w/pi);

    figure;
    plot(angle, 20*log10(abs(H)));
    grid on;
    xlabel('Angle (degrees)');
    ylabel('|H| (dB)');
    title(sprintf('Beampattern at n = %d', n_list(k)));

    yL = ylim;
    hold on;
    plot([des_angle des_angle], yL, 'k--');
    plot([intf_angle_1 intf_angle_1], yL, 'r--');
    plot([intf_angle_2 intf_angle_2], yL, 'r--');
    plot([intf_angle_3 intf_angle_3], yL, 'm--'); 
    hold off;
end


sound(real(X(1,:)),22000);
pause(5);
sound(real(y),22000);

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
L = length(h); h = h(:); % for vectors ONLY !!!
if( N < L )
error('DTFT: # data samples cannot exceed # freq samples')
end
W = (2*pi/N) * [ 0:(N-1) ]';
mid = ceil(N/2) + 1;
W(mid:N) = W(mid:N) - 2*pi; % move [pi,2pi) to [-pi,0)
W = fftshift(W);
H = fftshift( fft( h, N ) ); % move negative freq components
end
