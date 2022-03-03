function s = plot_earth_sphere(XE, YE, ZE, angle)

    % Load and define topographic data
    load('topo.mat','topo','topomap1');
    
    % Rotate data to be consistent with the Earth-Centered-Earth-Fixed
    % coordinate conventions. X axis goes through the prime meridian.
    % http://en.wikipedia.org/wiki/Geodetic_system#Earth_Centred_Earth_Fixed_.28ECEF_or_ECF.29_coordinates
    %
    % Note that if you plot orbit trajectories in the Earth-Centered-
    % Inertial, the orientation of the contintents will be misleading.
    topo2 = [topo(:,181:360) topo(:,1:180)]; %#ok<NODEF>

    % Define surface settings
    props.FaceColor= 'texture';
    props.EdgeColor = 'none';
    props.FaceLighting = 'phong';
    props.Cdata = topo2;
    
    % Create the sphere with Earth topography and adjust colormap
    s = surface(XE,YE,ZE,props,'parent',gca);
    direction = [0 0 1];
    rotate(s, direction, angle);
    colormap(topomap1);

end