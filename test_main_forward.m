N = 100;
NB = 100;

inputs = struct();

inputs.sourceMap = ones(N);
inputs.B = zeros([NB, NB, NB, 3]);
inputs.wxBE = 2;
inputs.wyBE = 2;
inputs.wzBE = 2;

inputs.beamMass = 1;
inputs.beamCharge = 1;
inputs.beamKEnergy = 1;
% inputs.c = 100;
inputs.l = 1000;
inputs.L = 10000;

% options
inputs.opt_withTargetMap = 1;

% magnetic field
[x, y, z] = meshgrid(linspace(-1, 1, NB), linspace(-1, 1, NB), linspace(-1, 1, NB));
s = [0.2, 0.2, 0.2];
potential = 1e-3;
potentialProfile = -potential * exp( - (x.*x / 2 / s(1) / s(1)) - (y.*y / 2 / s(2) / s(2)) - (z.*z / 2 / s(3) / s(3)));
inputs.B(:,:,:,3) = - sign(y) .* y.^2 / s(2)^2 .* potentialProfile; % z
inputs.B(:,:,:,2) =   sign(z) .* z.^2 / s(3)^2 .* potentialProfile; % y

outputs = main_forward(inputs);
imagesc(outputs.targetMap);
