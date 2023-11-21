%% Selective and non-selective excitation
% Constants 
gamma   = 267522187.44; % rad/s/T
gamma_hz= gamma/2/pi;
FA      = 45 * pi/180;
ntime   = 400; % number of samples
rf_len  = 4e-3; % second
td      = rf_len/ntime; %  dwell time 
voxel_sz= [3, 3, 3] * 1e-3;
fov     = [-60 45; -30 70; -90 90] * 1e-3; % [-100 100; -130 110; -120 120] * 1e-3;
rf_tbw  = 12;
is_sinc = true;
is_selective = true;

% positions
[prx, pry, prz] = ndgrid(fov(1,1):voxel_sz(1):fov(1,2), fov(2,1):voxel_sz(2):fov(2,2), fov(3,1):voxel_sz(3):fov(3,2));
sz      = size(prz);
pr      = transpose([prx(:), pry(:), prz(:)]);
npos    = size(pr, 2); % number of spatial positions

% off-resonance
b0 = zeros(size(prx));

% sinc pulse
b1 = complex(ones(ntime, 1) * FA/gamma/rf_len); % this is for non-selective
if is_sinc
    t = linspace(-rf_len/2, rf_len/2, ntime); % must be column
    BW = rf_tbw/rf_len;
    x  = pi*t*BW + eps; % didn't use 2pi -> we are interested in full BW not only the positive
    snc = sin(x) ./ x; 
    hamming_window = 0.53836 + 0.46164*cos(2*pi * linspace(-0.5,0.5,ntime));
    rf = transpose(snc .* hamming_window);
    rf = rf / sum(rf); % normalize
    b1 = complex(rf * FA/gamma/td);
end

% gradients
gr = zeros(3, ntime);
if is_selective
    BW = rf_tbw/rf_len; % [Hz]
    gr(3, :) = BW / voxel_sz(3) / gamma_hz;
end

m0 = [zeros(2, npos); ones(1, npos)];

% add phase to RF
b1 = real(b1) + 1i*real(b1);
b1 = b1 / sqrt(2);
b1 = repmat(b1, [1, npos]);

T1 = 100000;
T2 = 100000;

b1= single(b1);
gr= single(gr);
m0= single(m0);
b0= single(b0);
pr= single(pr);
T1 = single(T1);
td = single(td);
T2 = single(T2);

%% run
% B1 (complex), gr, td, b0, pr, t1, T2, m0, save_all
tic
try
clc
result = bloch_mex((b1), (gr), (td), (b0(:)), (pr), (T1), (T2), (m0), false);
catch me_err
me_err
end
toc
%%
close all
spm_viewer(reshape(result(1,end,:), sz), reshape(result(2,end,:), sz), reshape(result(3,end,:), sz));
vin(reshape(result(1,end,:), sz))
vin(reshape(result(2,end,:), sz))
vin(reshape(result(3,end,:), sz))
% for mex file is open error
clear functions

%%
close all
vin(reshape(result(1,100,:), sz))
vin(reshape(result(2,100,:), sz))
vin(reshape(result(3,100,:), sz))


%% off-resonance, first run the example above 
m0 = result;
rot= 90 * pi/180 / td / gamma;
b0 = ones(npos, 1) * rot; % Tesla
b1 = complex(zeros(1,size(b1,2))); % no rf
gr = [0;0;0]; % no gradient

result2    = bloch_mex(b1, gr, td, b0(:), pr, 1000, 1000, sens, m0);
clear functions

close all
vin(reshape(result2(1,:), sz))
vin(reshape(result2(2,:), sz))
vin(reshape(result2(3,:), sz))


