%% physical information
length_beam = 96; % in mm
width_beam = 48; % in mm

Young_modulus_1 = 100e9; % in Pa
Young_modulus_2 = 50e9; % in Pa
Young_modulus_3 = 25e9; % in Pa

F = 500; % in N, force at end of beam

%% Gauss quadrature
xi = [-(3/5)^(0.5) 0 (3/5)^(0.5)]; % points
wi = [5/9 8/9 5/9]; % weights

%% 
