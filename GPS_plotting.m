% Plotting

clear all; close all; clc;

load('GPS_data.mat');
load('chief_simdata.mat');

%% ECEF Positions Plot

% Plot each GPS satellite (blue) and spacecraft trajectory (red)
figure(); hold on; grid on;
for ii = 1:32
    plot3(GPS.ECEFs(ii,:,1), GPS.ECEFs(ii,:,2), GPS.ECEFs(ii,:,3), 'b');
end
plot3(simdata.r_ecef(1,:), simdata.r_ecef(2,:), simdata.r_ecef(3,:), 'r', 'LineWidth', 2);

% Earth ellipsoid
rE = 6378.137 * 1000;
[xE, yE, zE] = ellipsoid (0, 0, 0, rE, rE, rE, 20);
s = surf(xE, yE, zE);
s.FaceColor = 'b';
s.EdgeColor = 'k';

% Plot labeling
view(3); axis equal;
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
title('ECEF Positions of Spacecraft and GPS Satellites');

%% ECI Positions Plot

% Plot each GPS satellite (blue) and spacecraft trajectory (red)
figure(); hold on; grid on;
for ii = 1:32
    plot3(GPS.ECIs(ii,:,1), GPS.ECIs(ii,:,2), GPS.ECIs(ii,:,3), 'b');
end
plot3(simdata.r_eci(1,:), simdata.r_eci(2,:), simdata.r_eci(3,:), 'r', 'LineWidth', 2);

% Earth ellipsoid
rE = 6378.137 * 1000;
[xE, yE, zE] = ellipsoid (0, 0, 0, rE, rE, rE, 20);
s = surf(xE, yE, zE);
s.FaceColor = 'b';
s.EdgeColor = 'k';

% Plot labeling
view(3); axis equal;
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
title('ECI Positions of Spacecraft and GPS Satellites');

%% Make a GIF for ECEF

% Various figure specifications
figure(); hold on; grid on; axis equal; %view(3);
view(-45, 10);
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
title('ECEF Positions of Spacecraft and GPS Satellites');
axis([-26623762.708952, 26801305.4091477, ...
      -26494689.6841939, 26437777.7350942, ...
      -22290211.8140644, 22632697.1315625]);

% Plot dummy lines for legend and plot Earth sphere
plot3(0,0,0, 'Color', 'r', 'LineWidth', 2, 'DisplayName', 'S/C Trajectory');
scatter3(0,0,0, 5, 'k', 'filled', 'DisplayName', 'GPS Sat');
scatter3(0,0,0, 5, 'r', 'filled', 'DisplayName', 'Nearby GPS Sat');
plot3([0 0.1], [0 0.1], [0 0.1], 'g:', 'DisplayName', 'Received Measurement');
legend('AutoUpdate', 'off', 'Location', 'Northwest');
s = surf(xE, yE, zE);
s.FaceColor = 'b';
s.EdgeColor = 'k';

% Animated line functions
hs = arrayfun(@(x) animatedline('MaximumNumPoints', 1, 'LineWidth', 5, 'Marker', 'o', 'MarkerFaceColor', 'k'), 1:32);
h = animatedline('Color', 'r', 'LineWidth', 2);
hs2 =  arrayfun(@(x) animatedline('MaximumNumPoints', 2, 'Color', 'g', 'LineStyle', ':'), 1:4);

% Plot the animation!
for kk = 1:length(GPS.closest4)
    
    jj = 1;
    
    for ii = 1:32
        if ismember(ii, GPS.closest4(kk,:))
            hs(ii).Color = 'r';
            hs(ii).MarkerFaceColor = 'r';
            addpoints(hs2(jj), [GPS.ECEFs(ii,kk,1) simdata.r_ecef(1,kk)], [GPS.ECEFs(ii,kk,2) simdata.r_ecef(2,kk)], [GPS.ECEFs(ii,kk,3) simdata.r_ecef(3,kk)]);
            jj = jj + 1;
        else 
            hs(ii).Color = 'k';
            hs(ii).MarkerFaceColor = 'k';
        end
        addpoints(hs(ii), GPS.ECEFs(ii,kk,1), GPS.ECEFs(ii,kk,2), GPS.ECEFs(ii,kk,3));
    end
    addpoints(h, simdata.r_ecef(1,kk), simdata.r_ecef(2,kk), simdata.r_ecef(3,kk));
    drawnow limitrate
    
    % Set to 1 if we want to save as GIF
    if 0
        % Capture the plot as an image 
        frame = getframe(gcf); 
        im = frame2im(frame); 
        [imind,cm] = rgb2ind(im,256); 
        % Write to the GIF File 
        if kk == 1 
          imwrite(imind,cm,'test2.gif','gif', 'Loopcount',inf); 
        else 
          imwrite(imind,cm,'test2.gif','gif','WriteMode','append'); 
        end 
    end
    
    %pause(0.00001);
end
drawnow

%% Make a GIF for ECI

% Various figure specifications
figure(); hold on; grid on; axis equal; %view(3);
view(-45, 10);
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
title('ECI Positions of Spacecraft and GPS Satellites');
axis([-26623762.708952, 26801305.4091477, ...
      -26494689.6841939, 26437777.7350942, ...
      -22290211.8140644, 22632697.1315625]);

% Plot dummy lines for legend and plot Earth sphere
plot3(0,0,0, 'Color', 'r', 'LineWidth', 2, 'DisplayName', 'S/C Trajectory');
scatter3(0,0,0, 5, 'k', 'filled', 'DisplayName', 'GPS Sat');
scatter3(0,0,0, 5, 'r', 'filled', 'DisplayName', 'Nearby GPS Sat');
plot3([0 0.1], [0 0.1], [0 0.1], 'g:', 'DisplayName', 'Received Measurement');
legend('AutoUpdate', 'off', 'Location', 'Northwest');
s = surf(xE, yE, zE);
s.FaceColor = 'b';
s.EdgeColor = 'k';

% Animated line functions
hs = arrayfun(@(x) animatedline('MaximumNumPoints', 1, 'LineWidth', 5, 'Marker', 'o', 'MarkerFaceColor', 'k'), 1:32);
h = animatedline('Color', 'r', 'LineWidth', 2);
hs2 =  arrayfun(@(x) animatedline('MaximumNumPoints', 2, 'Color', 'g', 'LineStyle', ':'), 1:4);

% Plot the animation!
for kk = 1:length(GPS.closest4)
    
    jj = 1;
    
    for ii = 1:32
        if ismember(ii, GPS.closest4(kk,:))
            hs(ii).Color = 'r';
            hs(ii).MarkerFaceColor = 'r';
            addpoints(hs2(jj), [GPS.ECIs(ii,kk,1) simdata.r_eci(1,kk)], [GPS.ECIs(ii,kk,2) simdata.r_eci(2,kk)], [GPS.ECIs(ii,kk,3) simdata.r_eci(3,kk)]);
            jj = jj + 1;
        else 
            hs(ii).Color = 'k';
            hs(ii).MarkerFaceColor = 'k';
        end
        addpoints(hs(ii), GPS.ECIs(ii,kk,1), GPS.ECIs(ii,kk,2), GPS.ECIs(ii,kk,3));
    end
    addpoints(h, simdata.r_eci(1,kk), simdata.r_eci(2,kk), simdata.r_eci(3,kk));
    drawnow limitrate
    
    % Set to 1 if we want to save as GIF
    if 0
        % Capture the plot as an image 
        frame = getframe(gcf); 
        im = frame2im(frame); 
        [imind,cm] = rgb2ind(im,256); 
        % Write to the GIF File 
        if kk == 1 
          imwrite(imind,cm,'test2.gif','gif', 'Loopcount',inf); 
        else 
          imwrite(imind,cm,'test2.gif','gif','WriteMode','append'); 
        end 
    end
    
    %pause(0.00001);
end
drawnow
