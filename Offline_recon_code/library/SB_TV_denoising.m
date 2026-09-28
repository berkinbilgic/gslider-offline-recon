%% Minimize (over x) ||x-y||^2 + lambda ||Gx||_1 using SB
% x_init    : Initial solution         
% y         : Observation 
% lambda    : Regularization parameter for L1 term
% mu        : Regularization parameter for the slack-consistency term
% maxSBIter : Maximum number of SB iterations

function x = SB_TV_denoising(x_init, y, lambda, mu, maxSbIter)

N = size(x_init);

% Precompute Operators
[k2,k1] = meshgrid(0:N(2)-1, 0:N(1)-1);
Ex = 1 - exp(-2*pi*1i*k1/N(1));
Ey = 1 - exp(-2*pi*1i*k2/N(2));

M2_muE2 = 1 + mu * (abs(Ex).^2 + abs(Ey).^2);

Fy = fft2(y);

% Step 0
x = x_init; 
s = zeros([N,2]); 
Gx =  zeros([N,2]); % Data dimensions = Nx, Ny, #gradient directions

for sbIter = 1:maxSbIter

    % Step 1: soft-thresholding
    Fx = fft2(x);
    Gx(:,:,1) = ifft2(Fx .* Ex);
    Gx(:,:,2) = ifft2(Fx .* Ey);

    Gx_plus_s = Gx + s;

    z = sign(Gx_plus_s) .* max(abs(Gx_plus_s) - 2*lambda/(mu+eps), 0);

    % Step 2: Update each contrast separately
    z_minus_s = z - s;
    temp = Fy + mu*(conj(Ex) .* fft2(z_minus_s(:,:,1)) + conj(Ey) .* fft2(z_minus_s(:,:,2)));
    x = ifft2(temp ./ M2_muE2);

    % Step 3: Update s
    s = s - (z - Gx);

end


