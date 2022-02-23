%% Selective and non-selective excitation
% Constants 
gamma   = 267522187.44;
gamma_hz= gamma/2/pi;
FA      = 45 * pi/180;
ntime   = 400; % number of samples
rf_len  = 4e-3; % second
td      = rf_len/ntime; %  dwell time 
ncoil   = 8;
voxel_sz= [3, 3, 3] * 1e-3;
fov     = [-60 45; -30 70; -90 90] * 1e-3;
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
b1 = complex(ones(ntime, ncoil) * FA/gamma/rf_len/ncoil); % this is for non-selective
if is_sinc
    t = single(linspace(-rf_len/2, rf_len/2, ntime)); % must be column
    BW = rf_tbw/rf_len;
    x  = pi*t*BW + eps; % didn't use 2pi -> we are interested in full BW not only the positive
    snc = sin(x) ./ x; 
    hamming_window = 0.53836 + 0.46164*cos(2*pi * linspace(-0.5,0.5,ntime));
    rf = transpose(snc .* hamming_window);
    rf = repmat(rf / sum(rf), [1 ncoil]); % normalize
    b1 = complex(rf * FA/gamma/td/ncoil);
end

sens = complex(single(ones(ncoil, npos)));

% gradients
gr = zeros(3, ntime);
if is_selective
    BW = rf_tbw/rf_len; % [Hz]
    gr(3, :) = BW / voxel_sz(3) / gamma_hz;
end

m0 = single([zeros(2, npos); ones(1, npos)]);

% add phase to RF
b1 = real(b1) + 1i*real(b1);
b1 = b1 / sqrt(2);

% make single precision 
pr = single(pr);
b1 = single(b1);
gr = single(gr);
b0 = single(b0);
m0 = single(m0);

% run
tic
try
clc
result    = bloch_sim(b1, gr, single(td), b0(:), pr, single(1), single(1), sens, m0);
catch me_err
me_err
end
toc
close all
%   vin(reshape(result(1,:), sz))
% vin(reshape(result(2,:), sz))
% vin(reshape(result(3,:), sz))
% for mex file is open error
clear functions

%% off-resonance, first run the example above 
m0 = result;
rot= 90 * pi/180 / rf_len / gamma;
b0 = single(ones(npos, 1)) * rot; % Tesla
b1 = complex(single(zeros(1,size(b1,2)))); % no rf
gr = single([0;0;0]); % no gradient

result2    = bloch_sim(b1, gr, rf_len, b0, pr, single(10000), single(10000), sens, m0);

close all
vin(reshape(result2(1,:), sz))
vin(reshape(result2(2,:), sz))
vin(reshape(result2(3,:), sz))

