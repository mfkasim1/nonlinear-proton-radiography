% to get started, you can do as follow:
% >> cd <to the directory where you clone this project>
% >> addpath('lib')
% >> test_main_forward

N = 100; % number of proton grid (N x N)
NB = 100; % number of 3D grid of the object (NB x NB x NB)

inputs = struct();

inputs.sourceMap = ones(N); % initial distribution of the proton beam (uniform)
inputs.B = zeros([NB, NB, NB, 3]);
inputs.wxBE = 2; % size of the object
inputs.wyBE = 2;
inputs.wzBE = 2;

inputs.beamMass = 1;
inputs.beamCharge = 1;
inputs.beamKEnergy = 1;
% inputs.c = 100;
inputs.l = 1000; % distance from the source to the object
inputs.L = 10000; % distance from the object to the target plane

% options
inputs.opt_withTargetMap = 1;

% magnetic field
% the beam is propagating to z-direction
[x, y, z] = meshgrid(linspace(-1, 1, NB), linspace(-1, 1, NB), linspace(-1, 1, NB));
s = [3.2, 0.2, 0.2];
potential = 1e-3;
potentialProfile = -potential * exp( - (x.*x / 2 / s(1) / s(1)) - (y.*y / 2 / s(2) / s(2)) - (z.*z / 2 / s(3) / s(3)));
inputs.B(:,:,:,3) = - sign(y) .* y.^2 / s(2)^2 .* potentialProfile; % z
inputs.B(:,:,:,2) =   sign(z) .* z.^2 / s(3)^2 .* potentialProfile; % y

outputs = main_forward(inputs);
imagesc(outputs.targetMap); % show the proton radiograph
