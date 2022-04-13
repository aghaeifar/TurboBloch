%%
% 2D Excitation with Spiral trajectory 
% Ali Aghaeifar <ali.aghaeifar.mri@gmail.com>
%
%% Timing & Resolution
T   = 2.5e-3;   % Sec, dur of pulse
dt  = 5e-6;     % RF/grad dwell time
t   = 0:dt:T-dt;% time points

fov   = 0.20;   % m
res   = 0.01;   % m, resolution
kmax  = 1/2/res;% cycles/m, max spatial frequency
N     = ceil(kmax*fov); % number of turns -> kmax / 1 / fov where 1/fov is k-space resolution
gamma = 267522187.44 / 2 / pi; % Hz / T

% constant angular rate spiral-in:
k = kmax * (1-t./T) .* exp(1i*2*pi*N*(1-t./T)); % first (1-t./T) indicates spiral-in, second (1-t./T) indicates clock-wise rotation
g = -interp1(t, [0;diff(k)'], t+dt/2, 'spline' ,0) / gamma / dt; % T/m

subplot(2,3,1);
plot(k); xlabel('Kx'); ylabel('Ky');
subplot(2,3,2); 
plot(t, real(g)*1e3, t, imag(g)*1e3); ylabel('Gradient, mT/m')

% build system matrix (2x-oversampled)
xx = -fov/2 : res/2 : fov/2-res/2;
yy = -fov/2 : res/2 : fov/2-res/2;
zz = 0;
[x,y,z] = ndgrid(xx,yy,zz); 
A    = dt*2*pi*gamma*exp( 1i*2*pi*( x(:)*real(k) + y(:)*imag(k) ) );

% build desired pattern
mdes = (2*x.^2 + y.^2 < 0.0025);
mdes = IFFT2D( (hamming(size(x,1))*hamming(size(x,1))').^2 .* FFT2D(mdes) );
subplot(2,3,4);
imagesc(abs(mdes));  caxis([-1 1]); title('Desired Pattern');

% design RF
lambda = 1e8;
b1 = (A'*A + lambda*eye(length(k))) \ (A'*mdes(:)); % = inv(A'*A + lambda*eye(length(k))) * (A'*mdes(:));
% b1 = (A'*A ) \ (A'*mdes(:));

subplot(2,3,3);
plot(t, real(b1)*1e6, t, imag(b1)*1e6); ylabel('RF, (uT)')

% forward projection
nm = reshape (A*b1, size(x));
subplot(2,3,5);
imagesc(abs(nm)); caxis([-1 1]); title('Estimated Pattern');

% Simulate it
gr = transpose([real(g(:)), imag(g(:)), zeros(size(g(:)))]);
b0 = zeros(size(x(:)));
sens = complex(ones(size(b1,2), length(x(:)))); % complex(ones(ncoil, npos));
T1 = 5;
T2 = 2;

pr = transpose([x(:),y(:),z(:)]); 
m0 = zeros(size(pr));
m0(3,:) = 1;

%%
tic
try
result = bloch_mex(b1, gr, dt, b0(:), pr, 100000, 100000, sens, m0, true);
catch me_err
me_err
end
toc
%
result_xy = squeeze(result(1,end,:) + 1j*result(2,end,:));
result_xy = reshape(result_xy, size(x));

subplot(2,3,6);
imagesc(abs(result_xy));  caxis([-1 1]);
