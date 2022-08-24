%% Ernst Angle demonstration
% Constants 
gamma   = 267522187.44;
FA      = [0:0.1:15, 16:180];
TR      = [10, 1000] * 1e-3; % second
T1      = [100, 500, 1000] * 1e-3; % second
ntime   = 1;    % number of samples
rf_len  = 1e-6; % second
td      = rf_len/ntime; %  dwell time 
ncoil   = 1;

% positions
pr      = [0;0;0];
npos    = size(pr, 2); % number of spatial positions

% off-resonance
b0 = zeros([npos, 1]);

% RF pulse
B1 = complex(ones(ntime, ncoil) * FA * pi/180 /gamma/rf_len/ncoil); % this is for non-selective
sens = complex(ones(ncoil, npos));

% gradients
gr = zeros(3, ntime);

% Initial magnetization
m0 = [zeros(2, npos); ones(1, npos)];

%%
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
                result = bloch_mex(complex(B1(b1)), gr, td, b0, pr, t1, 1000, sens, result);
                % TR fill
                result = bloch_mex(complex(0), gr, TR(tr)-td, b0, pr, t1, 1000, sens, result);                
            end
            mxy(b1) = norm(result(1:2));
        end
        subplot(1,2,tr); plot(FA, mxy); hold on;
    end
    legend(strcat('T1 = ', arrayfun(@num2str, T1, 'UniformOutput', 0)))
    title(['TR = ' num2str(TR(tr)) ' [Sec]']);
end
