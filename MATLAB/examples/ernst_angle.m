%% Ernst Angle demonstration
% Constants 
gamma   = 267522187.44;
FA      = [0:0.1:15, 16:180];
TR      = [10, 1000] * 1e-3; % second
T1      = [100, 500, 1000] * 1e-3; % second
T2      = 1000;
ntime   = 1;    % number of samples
rf_len  = 1e-6; % second
td      = rf_len/ntime; %  dwell time 


% positions
pr      = [0;0;0];
npos    = size(pr, 2); % number of spatial positions

% off-resonance
b0 = zeros(npos, 1);

% RF pulse
B1 = complex(ones(ntime, npos) * FA * pi/180 /gamma/rf_len); % this is for non-selective

% gradients
gr = zeros(3, ntime);

% Initial magnetization
m0 = [zeros(2, npos); ones(1, npos)];

B1= single(B1);
gr= single(gr);
m0= single(m0);
b0= single(b0);
pr= single(pr);
T1 = single(T1);
td = single(td);
T2 = single(T2);

%%
% B1 (complex), gr, td, b0, pr, t1, T2, m0, save_all
nDummy = 200; % dummy pulses to reach steady-state
for tr = 1:numel(TR)
    for t1 = T1
        mxy = zeros(1, numel(B1));        
        for b1 = 1:numel(B1) 
            result= m0;
            for i=1:nDummy
                % perfect spoiler
                result(1:2) = 0; 
                % apply RF
                result = bloch_mex(complex(B1(b1)), gr, td, b0, pr, t1, T2, result, false);
                % TR fill
                result = bloch_mex(complex(single(0)), gr, TR(tr)-td, b0, pr, t1, T2, result, false);                
            end
            mxy(b1) = norm(result(1:2));
        end
        subplot(1,2,tr); plot(FA, mxy); hold on;
    end
    legend(strcat('T1 = ', arrayfun(@num2str, T1, 'UniformOutput', 0)))
    title(['TR = ' num2str(TR(tr)) ' [Sec]']);
end
clear functions
