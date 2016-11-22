%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Predict the intensity on the target plane given the deflection potential, Phi, and 
%   the intensity profile on the source plane, sourceMap.
% The beam from source plane is mapped to the target plane as function of its position on the source plane.
% The map is given by:
%   xt(xs,ys) = xs - \partial(Phi(xs,ys))/\partial(xs);
%   yt(xs,ys) = ys - \partial(Phi(xs,ys))/\partial(ys);
% Inputs:
%   * sourceMap: intensity profile on the source plane (Ny x Nx) [opt, if detectorSize, def=ones(detectorSize)]
%   * detectorSize: (Ny, Nx) number of pixels on the detector plane, if not same with the size of sourceMap, it is size of the sourceMap [opt, if sourceMap, def=size(sourceMap)]
%   * B: [Bx; By; Bz] where each component is the magnetic field (NyB x NxB x NzB x 3) [opt]
%   * E: [Ex; Ey; Ez] where each component is the electric field (NyB x NxB x NzB x 3) [opt]
%   * wxBE, wyBE, wzBE: the lengths (length, width, height) (in its own unit) of the objects
%   * beamMass: mass of the beam (in its own unit) [opt, def=1]
%   * beamKEnergy: energy of each particle in the beam (in its own unit) [opt, def=1]
%   * beamCharge: charge of each particle (in its own unit) [opt, def=1]
%   * l, L: length from the object to the beam source and to the target plane, respectively (in its own unit) [opt, def l=99999999999999, L=0]
% Options:
%   * opt_withTargetMap: 1 if we want to get the targetMap, 0 otherwise (targetMap takes quite much computational resources) [def=0]
%   * opt_track: contains (n x 2) or (n x 1)-matrix that shows the (iy, ix) or (iy+(ix-1)*Ny) position of which particles to save the positions and velocity [def=-1 (not tracking)]
% Outputs:
%   * targetMap: intensity profile on the target plane (Ny x Nx)
%   * xTP, yTP: the axis (in its own unit) on the target plane (1 x Nx and 1 x Ny)
%   * ixmap, iymap: where to map each particle on the target plane (in normalised units)
%   * trackedStates: contains (n x 5 x (NzB+2))-matrix that shows the positions and velocity of tracked particles during the interaction
% Assumptions:
%   * Normalised coordinate: the first index is (1,1) and all pixels have size of 1
%   * all the position vector / matrix has the same spacing between its two elements
%   * for most cases, l, L >> wxBE, wyBE, wzBE
%   * the central axis is located at the centre of the object
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function outputs = main_forward(inputs)
    
    %%%%%%%%%%%%%%%% handling the inputs %%%%%%%%%%%%%%%%
    % sourceMap
    if isfield(inputs, 'sourceMap')
        sourceMap = inputs.sourceMap;
    elseif isfield(inputs, 'detectorSize')
        sourceMap = ones(detectorSize);
    end
    [Ny, Nx] = size(sourceMap);
    
    % size of the object
    wxBE = inputs.wxBE;
    wyBE = inputs.wyBE;
    wzBE = inputs.wzBE;
    
    % distance from the object to the beam source (l) and to the target plane (L)
    if isfield(inputs, 'l') l = inputs.l; else l = 99999999999999; end % default: beam travel in parallel directions
    if isfield(inputs, 'L') L = inputs.L; else L = 0; end
    scale = (L/l + 1);
    
    % object B and E-field
    if isfield(inputs, 'B')
        B = inputs.B;
        if ~isfield(inputs, 'E') E = zeros(size(B)); % there is B, no E
        else E = inputs.E; end % there are B and E
    elseif isfield(inputs, 'E') % no B, there is E
        E = inputs.E;
        B = zeros(size(E));
    else % no E and B, the targetMap is just the source map
        targetMap = sourceMap;
        ixmap = [1:size(targetMap,2)];
        iymap = [1:size(targetMap,1)];
        xTP = linspace(-.5, .5, size(targetMap,2)) * scale * wxBE;
        yTP = linspace(-.5, .5, size(targetMap,1)) * scale * wyBE;
        outputs = struct('targetMap', targetMap, 'xTP', xTP, 'yTP', yTP, 'ixmap', ixmap, 'iymap', iymap);
        return;
    end
    [NyB, NxB, NzB, ~] = size(B);
    
    % beam parameters
    if isfield(inputs, 'beamMass') beamMass = inputs.beamMass; else beamMass = 1; end
    if isfield(inputs, 'beamCharge') beamCharge = inputs.beamCharge; else beamCharge = 1; end
    if isfield(inputs, 'beamKEnergy') beamKEnergy = inputs.beamKEnergy; else beamKEnergy = 0; end % default: no beam kinetic energy
    % if isfield(inputs, 'c') c = inputs.c; else c = 99999999999999; end % default: non-relativistic case
    
    % options
    if isfield(inputs, 'opt_withTargetMap') opt_withTargetMap = inputs.opt_withTargetMap; else opt_withTargetMap = 0; end
    if isfield(inputs, 'opt_track') opt_track = inputs.opt_track; else opt_track = -1; end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% pre-processing of the algorithm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % getting the basic parameters for this algorithm
    [ix,iy] = meshgrid([1:Nx], [1:Ny]); % position on the object plane (normalised unit) (Ny x Nx)
    
    % normalising the index of opt_track
    if size(opt_track,2) == 2
        opt_track = opt_track(:,1) + (opt_track(:,2)-1) * Ny;
    end
    isTracking = (sum(opt_track) ~= -1);
    
    % axis on the object (in its own units)
    xBE = linspace(-wxBE/2, wxBE/2, NxB);
    yBE = linspace(-wyBE/2, wyBE/2, NyB);
    zBE = linspace(-wzBE/2, wzBE/2, NzB);
    
    % calculating the beam speed from the energy
    speed = sqrt(2*beamKEnergy/beamMass);
    % speed = c * sqrt(1 - 1./(1 + beamKEnergy/beamMass/c/c)^2);
    
    % beam position on the object plane (in its own unit) (Ny x Nx)
    [xO, yO] = meshgrid(linspace(-wxBE/2, wxBE/2, Nx), linspace(-wyBE/2, wyBE/2, Ny));
    rOsq = xO.^2 + yO.^2;
    rO = sqrt(rOsq);
    sqrtlsqrOsq = sqrt(l.*l - rOsq);
    
    % beam position on the object plane (in normalised units, starting from (1,1))
    ixO = ix;
    iyO = iy;
    
    % initial speed of the beam just before entering the plasma
    vx = speed * xO ./ sqrtlsqrOsq;
    vy = speed * yO ./ sqrtlsqrOsq;
    vz = speed * l ./ sqrtlsqrOsq;
    speed = speed + zeros(size(vx));
    % invgamma = sqrt(1 - speed .* speed / c / c);
    
    % the field position on the beam's plane (used for interpolation in normalised units, starting from (1,1))
    ixB = (xBE+wxBE/2) / wxBE * (Nx-1) + 1;
    iyB = (yBE+wyBE/2) / wyBE * (Ny-1) + 1;
    [ixB, iyB] = meshgrid(ixB, iyB);
    
    % so at this point, ixO and iyO are the position of the beam (in normalised units), vx and vy are the velocity of the beam
    % save the tracked particles
    if isTracking
        trackedStates = zeros(size(opt_track,1), 5, NzB+2); % 5 states
        trackedStates(:,1,1) = ixO(opt_track);
        trackedStates(:,2,1) = iyO(opt_track);
        trackedStates(:,3,1) = vx(opt_track);
        trackedStates(:,4,1) = vy(opt_track);
        trackedStates(:,5,1) = vz(opt_track);
    else
        trackedStates = [];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% simulating the interactions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for izO = [1:NzB]
        % extract the electric and magnetic field component in each position
        Bx = B(:,:,izO,1); Ex = E(:,:,izO,1);
        By = B(:,:,izO,2); Ey = E(:,:,izO,2);
        Bz = B(:,:,izO,3); Ez = E(:,:,izO,3);
        
        % linearly interpolate (and zero extrapolate) the magnetic and electric field in the position
        Bpx = interp2(ixB, iyB, Bx, ixO, iyO, 'linear', 0);
        Bpy = interp2(ixB, iyB, By, ixO, iyO, 'linear', 0);
        Bpz = interp2(ixB, iyB, Bz, ixO, iyO, 'linear', 0);
        Epx = interp2(ixB, iyB, Ex, ixO, iyO, 'linear', 0);
        Epy = interp2(ixB, iyB, Ey, ixO, iyO, 'linear', 0);
        Epz = interp2(ixB, iyB, Ez, ixO, iyO, 'linear', 0);
        
        % get the force
        Fx = beamCharge * (vy .* Bpz - vz .* Bpy + Epx);
        Fy = beamCharge * (vz .* Bpx - vx .* Bpz + Epy);
        Fz = beamCharge * (vx .* Bpy - vy .* Bpx + Epz);
        
        % % now get the force parallel and perpendicular to each beam's velocity (relativistic case)
        % Fparx = Fx .* vx ./ speed;
        % Fpary = Fy .* vy ./ speed;
        % Fparz = Fz .* vz ./ speed;
        % Fperpx = Fx - Fparx;
        % Fperpy = Fy - Fpary;
        % Fperpz = Fz - Fparz;
        
        % % get the acceleration
        % invgamma3 = (invgamma .^ 3);
        % ax = (Fparx .* invgamma3 + Fperpx .* invgamma) / beamMass;
        % ay = (Fpary .* invgamma3 + Fperpy .* invgamma) / beamMass;
        % az = (Fparz .* invgamma3 + Fperpz .* invgamma) / beamMass;
        
        % get the acceleration (non-relativistic case)
        ax = Fx / beamMass;
        ay = Fy / beamMass;
        az = Fz / beamMass;
        
        % get the time required to get to the next xy-grid (the next z-position)
        dt = 1 ./ vz;
        
        % get the next position
        ixO = ixO + vx .* dt;
        iyO = iyO + vy .* dt;
        
        % get the next velocity
        vx = vx + ax .* dt;
        vy = vy + ay .* dt;
        vz = vz + az .* dt;
        speed2 = vx .* vx + vy .* vy + vz .* vz;
        speed = sqrt(speed2);
        % invgamma = sqrt(1 - speed2 / c / c);
        
        % track the states of the selected particles
        if isTracking
            trackedStates(:,1,izO+1) = ixO(opt_track);
            trackedStates(:,2,izO+1) = iyO(opt_track);
            trackedStates(:,3,izO+1) = vx(opt_track);
            trackedStates(:,4,izO+1) = vy(opt_track);
            trackedStates(:,5,izO+1) = vz(opt_track);
        end
    end
    
    % translate from ixO and iyO to xO and yO
    dxO = (xO(1,2) - xO(1,1));
    dyO = (yO(2,1) - yO(1,1));
    xO = (ixO - 1) * dxO + xO(1,1);
    yO = (iyO - 1) * dyO + yO(1,1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% going from the object to the target plane %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % get the axis of the target plane (in its own units)
    xTP = linspace(-.5, .5, Nx) * scale * wxBE;
    yTP = linspace(-.5, .5, Ny) * scale * wyBE;
    
    % get where the beam on the target plane (in its own unit)
    xT = xO + vx ./ vz * L;
    yT = yO + vy ./ vz * L;
    
    % now the position on normalised unit
    ixmap = (xT - xTP(1)) / (xTP(2) - xTP(1)) + 1;
    iymap = (yT - yTP(1)) / (yTP(2) - yTP(1)) + 1;
    
    % ixmap = (ixO + vx ./ vz * L / dxO) / scale + (1 - 1./scale) * (1 + wxBE/2/dxO);
    % iymap = (iyO + vy ./ vz * L / dyO) / scale + (1 - 1./scale) * (1 + wyBE/2/dyO);
    
    % track the states of the selected particles
    if isTracking
        trackedStates(:,1,NzB+2) = ixmap(opt_track);
        trackedStates(:,2,NzB+2) = iymap(opt_track);
        trackedStates(:,3,NzB+2) = vx(opt_track);
        trackedStates(:,4,NzB+2) = vy(opt_track);
        trackedStates(:,5,NzB+2) = vz(opt_track);
    end
    
    % from the mapping function, fill in the intensity on the target plane
    if (opt_withTargetMap)
        targetMap = fill_in_pixels(ix, iy, ixmap, iymap, sourceMap);
    else
        targetMap = 0;
    end
    
    % create the struct of the output
    outputs = struct('targetMap', targetMap, 'xTP', xTP, 'yTP', yTP, 'ixmap', ixmap, 'iymap', iymap, 'trackedStates', trackedStates);
end

