close all; clear; clc;
% Focused lens simulation (ideal aspheric phase)
% Editable parameters
lambda = 532e-9;        % Wavelength (m)
f = 25e-3;              % Focal length (m)
D = 30e-3;               % Aperture diameter (m)
N = 256;               % Grid size (NxN)
L = 50e-3;              % Physical grid size (m) (square extent)
useAngularSpectrum = false; % true: angular spectrum, false: Fresnel
zStart = 0e-3;          % Z start (m)
zEnd = 120e-3;           % Z end (m)
Nz = 801;               % Number of z samples for longitudinal views
plotLog = true;         % true: log-scale intensity for fringe visibility
dynamicRangeDb = 40;    % dB range for log plots

% Radial-model comparison (faster than full 2D for interference patterns)
runRadialModels = true; % true: run 1D radial Fresnel-Hankel and Debye-Wolf
Nr = 1200;              % Radial pupil samples
Nro = 600;              % Radial observation samples (for y axis)
NzRad = 401;            % Z samples for radial models
rObsMax = L/2;          % Max observation radius (m)
Ntheta = 800;           % Angular samples for Debye-Wolf

% Derived parameters
k = 2*pi/lambda;
dx = L/N;
x = (-N/2:N/2-1) * dx;
[X, Y] = meshgrid(x, x);
R = sqrt(X.^2 + Y.^2);

% Aperture mask
aperture = double(R <= D/2);

% Ideal aspheric lens phase for focusing at f
% Optical path difference to a reference plane at z=0
% phi(x,y) = -k * (sqrt(r^2 + f^2) - f)
phi = -k * (sqrt(R.^2 + f^2) - f);

% Lens phase field
lensField = aperture .* exp(1i * phi);

% Propagation to focus
z = f;

% Precompute spectrum for repeated propagation
U0 = fftshift(fft2(ifftshift(lensField)));

if useAngularSpectrum
    % Angular spectrum propagation
    fx = (-N/2:N/2-1) / L;  % spatial frequency (1/m)
    [FX, FY] = meshgrid(fx, fx);
    H = exp(1i * 2*pi*z * sqrt(max(0, (1/lambda)^2 - FX.^2 - FY.^2)));
    Uz = ifftshift(ifft2(fftshift(U0 .* H)));
else
    % Fresnel propagation (paraxial)
    % Quadratic phase transfer function
    fx = (-N/2:N/2-1) / L;
    [FX, FY] = meshgrid(fx, fx);
    H = exp(1i*pi*lambda*z * (FX.^2 + FY.^2));
    Uz = ifftshift(ifft2(fftshift(U0 .* H))) * exp(1i*k*z);
end

% Intensity at focus
I = abs(Uz).^2;
I = I / max(I(:));

% Plot results
figure('Color','w');
subplot(2,2,1);
imagesc(x*1e3, x*1e3, aperture);
axis image; colormap gray; colorbar;
title('Aperture'); xlabel('x (mm)'); ylabel('y (mm)');

subplot(2,2,2);
imagesc(x*1e3, x*1e3, angle(lensField));
axis image; colormap hsv; colorbar;
title('Lens phase (rad)'); xlabel('x (mm)'); ylabel('y (mm)');

subplot(2,2,3);
imagesc(x*1e3, x*1e3, abs(lensField));
axis image; colormap gray; colorbar;
title('Lens amplitude'); xlabel('x (mm)'); ylabel('y (mm)');

subplot(2,2,4);
imagesc(x*1e3, x*1e3, I);
axis image; colormap hot; colorbar;
title('Intensity at focus'); xlabel('x (mm)'); ylabel('y (mm)');

% Optional: line cut through focus
figure('Color','w');
plot(x*1e3, I(N/2+1, :), 'LineWidth', 1.5);
grid on; xlabel('x (mm)'); ylabel('Normalized intensity');
title('Focus line cut (y=0)');

% Z-Y cross-section (x = 0)
zList = linspace(zStart, zEnd, Nz);
Izy = zeros(N, Nz);
Ixz = zeros(N, Nz);
xIdx = N/2 + 1;
yIdx = N/2 + 1;

if useAngularSpectrum
    for zi = 1:Nz
        H = exp(1i * 2*pi*zList(zi) * sqrt(max(0, (1/lambda)^2 - FX.^2 - FY.^2)));
        Uz = ifftshift(ifft2(fftshift(U0 .* H)));
        Izy(:, zi) = abs(Uz(:, xIdx)).^2;
        Ixz(:, zi) = abs(Uz(yIdx, :)).^2;
    end
else
    for zi = 1:Nz
        H = exp(1i*pi*lambda*zList(zi) * (FX.^2 + FY.^2));
        Uz = ifftshift(ifft2(fftshift(U0 .* H))) * exp(1i*k*zList(zi));
        Izy(:, zi) = abs(Uz(:, xIdx)).^2;
        Ixz(:, zi) = abs(Uz(yIdx, :)).^2;
    end
end

Izy = Izy / max(Izy(:));
Ixz = Ixz / max(Ixz(:));

if plotLog
    IzyPlot = 10*log10(Izy + eps);
    IxzPlot = 10*log10(Ixz + eps);
    cmin = -dynamicRangeDb;
    cmax = 0;
else
    IzyPlot = Izy;
    IxzPlot = Ixz;
    cmin = 0;
    cmax = 1;
end

figure('Color','w');
imagesc(zList*1e3, x*1e3, IzyPlot);
axis image; colormap hot; colorbar; caxis([cmin cmax]);
title('Z-Y cross-section (x=0)');
xlabel('z (mm)'); ylabel('y (mm)');

figure('Color','w');
imagesc(zList*1e3, x*1e3, IxzPlot);
axis image; colormap hot; colorbar; caxis([cmin cmax]);
title('Z-X cross-section (y=0)');
xlabel('z (mm)'); ylabel('x (mm)');

% On-axis intensity vs z
IonAxis = squeeze(Izy(xIdx, :));
figure('Color','w');
plot(zList*1e3, IonAxis, 'LineWidth', 1.5);
grid on; xlabel('z (mm)'); ylabel('Normalized intensity');
title('On-axis intensity vs z');

% -------------------------------------------------------------------------
% Radial (cylindrical) models for interference near focus
% -------------------------------------------------------------------------
if runRadialModels
    % Radial grids
    rP = linspace(0, D/2, Nr);              % Pupil radius
    rObs = linspace(0, rObsMax, Nro);       % Observation radius
    zListRad = linspace(zStart, zEnd, NzRad);

    % Pupil phase (ideal asphere)
    phiP = -k * (sqrt(rP.^2 + f^2) - f);
    Up = exp(1i * phiP);                    % Uniform amplitude

    % Fresnel-Hankel propagation (scalar, paraxial)
    UzyH = zeros(Nro, NzRad);
    for zi = 1:NzRad
        zR = zListRad(zi);
        Q = exp(1i * k * rP.^2 / (2*zR));
        integrand = (Up .* Q) .* rP;  % includes r' for Hankel integral
        J0arg = (k/zR) * (rObs(:) * rP);  % Nro x Nr
        J0 = besselj(0, J0arg);
        integral = 2*pi * trapz(rP, J0 .* integrand, 2);
        U = exp(1i*k*zR) ./ (1i*lambda*zR) .* exp(1i*k*rObs(:).^2/(2*zR)) .* integral;
        UzyH(:, zi) = U;
    end
    IzyH = abs(UzyH).^2;
    IzyH = IzyH / max(IzyH(:));

    % Debye-Wolf focal field (aplanatic, scalar)
    NA = (D/2) / sqrt(f^2 + (D/2)^2);
    alpha = asin(NA);
    theta = linspace(0, alpha, Ntheta);
    sinT = sin(theta);
    cosT = cos(theta);
    apod = sinT .* sqrt(cosT);  % aplanatic apodization

    UzyD = zeros(Nro, NzRad);
    for zi = 1:NzRad
        zR = zListRad(zi);
        phase = exp(1i * k * zR .* cosT);
        for ri = 1:Nro
            J0 = besselj(0, k * rObs(ri) .* sinT);
            UzyD(ri, zi) = trapz(theta, apod .* phase .* J0);
        end
    end
    IzyD = abs(UzyD).^2;
    IzyD = IzyD / max(IzyD(:));

    % Mirror to negative y for display
    yFull = [-fliplr(rObs(2:end)), rObs] * 1e3; % mm
    IzyHfull = [flipud(IzyH(2:end, :)); IzyH];
    IzyDfull = [flipud(IzyD(2:end, :)); IzyD];

    if plotLog
        IzyHplot = 10*log10(IzyHfull + eps);
        IzyDplot = 10*log10(IzyDfull + eps);
        cmin = -dynamicRangeDb;
        cmax = 0;
    else
        IzyHplot = IzyHfull;
        IzyDplot = IzyDfull;
        cmin = 0;
        cmax = 1;
    end

    figure('Color','w');
    imagesc(zListRad*1e3, yFull, IzyHplot);
    axis image; colormap hot; colorbar; caxis([cmin cmax]);
    title('Radial Fresnel-Hankel Z-Y cross-section');
    xlabel('z (mm)'); ylabel('y (mm)');

    figure('Color','w');
    imagesc(zListRad*1e3, yFull, IzyDplot);
    axis image; colormap hot; colorbar; caxis([cmin cmax]);
    title('Debye-Wolf Z-Y cross-section');
    xlabel('z (mm)'); ylabel('y (mm)');

    % On-axis comparison
    IonH = IzyH(1, :);
    IonD = IzyD(1, :);
    figure('Color','w');
    plot(zListRad*1e3, IonH, 'LineWidth', 1.5); hold on;
    plot(zListRad*1e3, IonD, 'LineWidth', 1.5);
    grid on; xlabel('z (mm)'); ylabel('Normalized intensity');
    legend('Fresnel-Hankel', 'Debye-Wolf');
    title('On-axis intensity comparison');
end
