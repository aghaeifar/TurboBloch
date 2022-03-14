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
pr      = single([0;0;0]);
npos    = size(pr, 2); % number of spatial positions

% off-resonance
b0 = single(zeros([npos, 1]));

% RF pulse
B1 = complex(single(ones(ntime, ncoil)) * FA * pi/180 /gamma/rf_len/ncoil); % this is for non-selective
sens = complex(single(ones(ncoil, npos)));

% gradients
gr = single(zeros(3, ntime));

% Initial magnetization
m0 = single([zeros(2, npos); ones(1, npos)]);

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
                result = bloch_mex(complex(B1(b1)), gr, single(td), b0, pr, single(t1), single(1000), sens, single(result));
                % TR fill
                result = bloch_mex(complex(single(0)), gr, single(TR(tr)-td), b0, pr, single(t1), single(1000), sens, result);                
            end
            mxy(b1) = norm(result(1:2));
        end
        subplot(1,2,tr); plot(FA, mxy); hold on;
    end
    legend(strcat('T1 = ', arrayfun(@num2str, T1, 'UniformOutput', 0)))
    title(['TR = ' num2str(TR(tr)) ' [Sec]']);
end
