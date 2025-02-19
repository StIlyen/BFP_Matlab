clc
clear 
close all
tic

%% Input parameters
np =201; % Number of points for calculation
wl = linspace(0.6,0.69,np); % wavelength of incident light
% E = linspace(2.0,2.4,np);
% wl = 1.2398./E;

% load dispersion_MAPbBr3.mat
% 
% x1 = dispersion_MAPbBr3(:,1); % dispersion wavelength
% y1 = dispersion_MAPbBr3(:,2); % n from dispersion data 
% z1 = dispersion_MAPbBr3(:,3); % k from dispersion data

% % interpolation of "n" & "k"
% n = interp1(x1,y1,wl,'spline');
% k = interp1(x1,z1,wl,'spline');

%% Grating parameters of unit cells
a = 0.420; % period
w = 0.320; % ridge width
alp = w/a; % fill factor

h1 = 0.045; % ridge height (d1)
h2 = 0.075; % grating layer thickness (d2)

n1 = 1.45; % refractive index of substrate
n2 = 1.89; % refractive index of material
n3 = 1.00; % refractive index of superstrate 

eps1 = (n1^2).*ones(length(wl),1);
eps21 = (n2^2).*ones(length(wl),1);
% eps22 = (n + 1i.*k).^2; % permitivity from dispersion data
eps3 = (n3^2).*ones(length(wl),1);

pol='TE'; % polarization
theta = linspace(-10,10,np); % angle of incidence
knorm = sin(theta.*pi/180);
phi = 0;

%% Create empty matrix
R1 = zeros(length(wl),length(theta)); 
b1 = zeros(length(wl),length(theta));
T1 = zeros(length(wl),length(theta));

%% Compute values
for i=1:length(wl)
     for j=1:length(theta)
         [R1(i,j),b1(i,j),T1(i,j)] = lam_1D_2(wl(i), h1, h2, a, alp, eps1(i), eps21(i), eps3(i), theta(j), phi, pol);
    end
end

%% Plotting figures

figure(1)
surfc(theta,wl,R1)
view(2)
colormap turbo
shading interp
title('p = 420 nm w = 320 nm h = 45 nm d = 75 nm')
set(gca,'FontSize',11)
xlabel('k_x/k_0','FontSize',11)
ylabel('Wavelength (um)','FontSize',11)
% xticks([-0.4,-0.2,0.0,0.2,0.4])
% yticks([2.0,2.1,2.2,2.3])
% xlim([-0.2 0.2])
% ylim([2.1 2.4])


toc