function [R,b,T] = lam_1D_2(wl, h1, h2, gp, alp, eps1, eps2, eps3, theta, phi, pol)

%% Initialization
wv = 2*pi./wl; % wavevector
kh1 = wv.*h1; % grating heigth multiplied by the vacuum wavenumber
kh2 = wv.*h2; % grating layer multiplied by the vacuum wavenumber
kg = wl./gp; % wavelength-to-period ratio (grating vector)
no = 20; % # of Fourier modes
ind0 = ceil(no/2); % index of the zero harmonic (0th order diffraction)

%% Incidence wavevector projections:
kx0 = sqrt(eps3).*sin(theta.*pi/180).*cos(phi.*pi/180);

V_inc = zeros(no,2); % matrix of incident field amplitudes
V_inc(ind0,2) = 1; % plane wave coming from the superstrate

%% Matrix diffraction for grating

%% substrate-material transition
SM1 = calc_SMD_interface(no, kx0, kg, eps1, eps2, pol);
SM2 = calc_SMD_layer(no, kx0, kg, kh2, eps2);
SM = mul_SM(SM1, SM2);

%% layer-grating transition
FE = calc_emn_lam(no, alp, eps2, eps3);
SM3 = fmm(no, kx0, kg, kh1, eps2, eps3, FE, pol);

SM = mul_SM(SM, SM3);

% diffraction of a plane wave example
V_dif(:,1) = SM(:,:,1,1)*V_inc(:,1) + SM(:,:,1,2)*V_inc(:,2); % diffraction to the substrate
V_dif(:,2) = SM(:,:,2,1)*V_inc(:,1) + SM(:,:,2,2)*V_inc(:,2); % diffraction to the superstrate

% check the power conservation:
b = fmm_balance(no, V_inc, V_dif, kx0, kg, eps1, eps3, pol);
% fprintf("power balance: %f\n",b); % precision of the power conservation
% calculate the vector of diffraction efficiencies
V_eff = fmm_efficiency(no, V_inc, V_dif, kx0, kg, eps1, eps3, pol);

% calculate the vector of diffraction efficiencies:
R = V_eff(ind0,2);
T = V_eff(ind0,1);

end