% Francis Le Blanc
% Part 3: Adaptive Beamformer with Time-Varying Interference
%
% Extends the adaptive beamformer to handle interference sources whose
% directions change over time. Interferer angles vary continuously
% throughout the simulation, creating a dynamic signal environment.
%
% The beamformer updates its covariance estimate using an exponentially
% weighted approach and adjusts its spatial response to track and suppress
% moving interference sources.

% Variables
N = 100000; % number of samples for generating 
M = 18; % number of anetennas in array
T = 100; % number of samples for beamformer
e = 0.001; % forgetting factor

des_angle = 8; % angle of desired signal
intf_angle_1 = -5; % angle of first interfering signal
intf_angle_2 = 25; % angle of second interfering signal
slope_1 = 5/N; % slope of first interfering signal
slope_2 = -5/N; % slope of second interfering signal

des_amp = 1; % amplitude of desired signal
intf_amp = 120; % amplitude of interfering signals
noise_amp = 0.06; % amplitude of noise signal

y = zeros(1,N); % beamformer output 
g_values = struct();
n_values = [1000 20000 40000 60000 80000 99900];

eval_step = 200;                 % evaluate every 200 samples 
eval_n = T+1:eval_step:N;

H_des_dB  = zeros(size(eval_n)); % desired signal magnitude response    
H_i1_dB   = zeros(size(eval_n)); % interference signal 1 magnitude response
H_i2_dB   = zeros(size(eval_n)); % interference signal 2 magnitude response
t_seconds = (eval_n-1)/22000; % magnitude response time step

% Load audio signal
load ABF.mat;

% Choose signals
des_sound = voice; % desired signal audio
intf_sound_1 = jetnoise; % first interference signal audio
intf_sound_2 = babble; % second interference signal audio

% Take N samples of each sound as a column vector 
des_samples = des_sound(1:N);    
des_samples = des_samples(:).';
intf_1_samples = intf_sound_1(1:N); 
intf_1_samples = intf_1_samples(:).';
intf_2_samples = intf_sound_2(1:N); 
intf_2_samples = intf_2_samples(:).';

% Create steering vectors
a_des = linarr(M,des_angle);
a_intf_1 = linarr(M, intf_angle_1);
a_intf_2 = linarr(M, intf_angle_2);


% Create total signal
[X, intf1_ang_t, intf2_ang_t] = X_timevarying( ...
                                M, N, des_samples, intf_1_samples, intf_2_samples, ...
                                des_angle, intf_angle_1, intf_angle_2, ...
                                slope_1, slope_2, des_amp, intf_amp, noise_amp);

% Take first 100 samples of total signal
X_samples = X(:, 1:T); 

% Construct vectors for beamformer, only defined steering vector for
% desired signal
Rxx = (X(:,1:T) * X(:,1:T)') / T;
A = a_des;
c = 1;

% Update Rxx adaptively
for n = T+1:N
    x = X(:,n);
    Rxx = (1-e)*Rxx + e*(x*x');

    RiA = Rxx \ A;
    g = RiA * ((A' * RiA) \ c);
    h = conj(g);

    % Calculate magnitude response of each signal at time step intervals
    if mod(n-(T+1), eval_step) == 0
        idx = (n-(T+1))/eval_step + 1;
        ang_grid = -90:0.25:90;               % degrees 
        Agrid = linarr(M, ang_grid);          % M x length(ang_grid)
        Hgrid = h.' * Agrid;                  % 1 x length(ang_grid)

        H_des_dB(idx) = 20*log10(abs(h.'*linarr(M, des_angle)));

    
        [~,k1] = min(abs(ang_grid - intf1_ang_t(n)));
        [~,k2] = min(abs(ang_grid - intf2_ang_t(n)));

        H_i1_dB(idx) = 20*log10(abs(Hgrid(k1)));
        H_i2_dB(idx) = 20*log10(abs(Hgrid(k2)));
    end

    y(n) = h.' * x;

    if any(n == n_values)
        g_values.(sprintf('n%d',n)) = h;
    end
end

% Plot the magnitude response for each signal over time
figure; grid on; hold on;
plot(t_seconds, H_des_dB, 'DisplayName','Desired (8°)');
plot(t_seconds, H_i1_dB,  'DisplayName','Interferer 1');
plot(t_seconds, H_i2_dB,  'DisplayName','Interferer 2');
xlabel('Time (seconds)');
ylabel('Array response magnitude (dB)');
title('Adaptive beamformer magnitude response vs. time');
legend('Location','best');

% Plot interfering signal angles over time
figure; grid on; hold on;
plot(t_seconds, intf1_ang_t(eval_n), 'DisplayName','θ_i1(t)');
plot(t_seconds, intf2_ang_t(eval_n), 'DisplayName','θ_i2(t)');
xlabel('Time (seconds)'); ylabel('Angle (deg)');
title('Moving interferer angles vs time');
legend;

% Plot spatial frequency response at chosen sample times
fn = fieldnames(g_values);

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
    hold off;
end

% Playback raw signal vs. received signal
%sound(real(X(1,:)),22000); 
%pause(5); 
sound(real(y),22000);

% Functions

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
% constraint: N >= L
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


function [X, intf_angle_1_t, intf_angle_2_t] = X_timevarying( ...
    M, N, des_samples, intf_1_samples, intf_2_samples, ...
    des_angle, intf_angle_1, intf_angle_2, ...
    slope_1, slope_2, des_amp, intf_amp, noise_amp)

    % Enforce row vectors
    des_samples    = des_samples(:).';
    intf_1_samples = intf_1_samples(:).';
    intf_2_samples = intf_2_samples(:).';

    % Signals of length N
    des_samples    = des_samples(1:N);
    intf_1_samples = intf_1_samples(1:N);
    intf_2_samples = intf_2_samples(1:N);

    % Compute time-varying angles
    n = 0:(N-1);
    intf_angle_1_t = intf_angle_1 + slope_1*n;
    intf_angle_2_t = intf_angle_2 + slope_2*n;

    % Desired signal 
    a_des = linarr(M, des_angle);
    X = des_amp * (a_des * des_samples);

    % Add time-varying inteference signals
    for k = 1:N
        a_intf_1_k = linarr(M, intf_angle_1_t(k));
        a_intf_2_k = linarr(M, intf_angle_2_t(k));
        X(:,k) = X(:,k) + intf_amp*a_intf_1_k*intf_1_samples(k) + intf_amp*a_intf_2_k*intf_2_samples(k);
    end

    % Add noise
    X = X + noise_amp*(randn(M,N) + 1j*randn(M,N));
end