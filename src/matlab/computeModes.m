%% Mode Frequencies For Isometric Generalised Plate
%
%   Authors: Michele Ducceschi, Matthew Hamilton, Sebastian Duran
%
%% House Keeping
clearvars
close all

%% physical and elastic parameters
Lx = 1.10; Ly = 0.8; Lz = 5e-3;
ldim    = [Lx Ly Lz];       %-- plate dimensions [x, y, z] in metres
E0      = 9.0e+9 ;        %-- Young's mod [Pa]
rho     = 8765 ;            %-- density [kg/m^3]
nu      = 0.3 ;             %-- poisson's ratio
Nmodes  = 1500;               %-- number of modes to compute
h       = sqrt(Lx*Ly)*0.01; %-- grid spacing
BCs     = ones(4,2) * 1e15; %-- elastic constants around the edges
% BCs(1,:) = 1e15;


%% compute modal frequencies
mode_freqs = zeros(Nmodes,5);
col = 1;

for E = linspace(E0*0.8, E0*1.2, 1)    
    [Q,Om,Nx,Ny,biHarm] = magpie(rho, E, nu, ldim, h, BCs, Nmodes, 'chladni');
    fprintf('--- Modal Freqs (Hz) ---\n');
    fprintf('%.1f\n', Om/2/pi);
    mode_freqs(:, col) = Om/2/pi;
    col = col + 1;
end