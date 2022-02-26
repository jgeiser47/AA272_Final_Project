% Plotting

clear all; close all; clc;

load('GPS_data.mat');
load('chief_simdata.mat');

figure(); hold on; grid on;
for ii = 1:32
    plot3(GPS.ECEFs(ii,:,1), GPS.ECEFs(ii,:,2), GPS.ECEFs(ii,:,3), 'b');
end
plot3(simdata.r_ecef(1,:), simdata.r_ecef(2,:), simdata.r_ecef(3,:), 'r', 'LineWidth', 2);
%plot3(simdata.r_eci(1,:), simdata.r_eci(2,:), simdata.r_eci(3,:), 'r', 'LineWidth', 2);

% Earth ellipsoid
rE = 6378.137 * 1000;
[xE, yE, zE] = ellipsoid (0, 0, 0, rE, rE, rE, 20);
s = surf(xE, yE, zE);
s.FaceColor = 'b';
s.EdgeColor = 'k';

view(3); axis equal;
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
title('ECEF Positions of Spacecraft and GPS Satellites');

%%

figure(); hold on; grid on; axis equal; %view(3);
view(-45, 10);
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
title('ECEF Positions of Spacecraft and GPS Satellites');

plot3(0,0,0, 'Color', 'r', 'LineWidth', 2, 'DisplayName', 'S/C Trajectory');
scatter3(0,0,0, 5, 'k', 'filled', 'DisplayName', 'GPS Sat');
scatter3(0,0,0, 5, 'r', 'filled', 'DisplayName', 'Received Measurement');
legend('AutoUpdate', 'off', 'Location', 'Northwest');

axis([-26623762.708952, 26801305.4091477, ...
      -26494689.6841939, 26437777.7350942, ...
      -22290211.8140644, 22632697.1315625]);
s = surf(xE, yE, zE);
s.FaceColor = 'b';
s.EdgeColor = 'k';

hs = arrayfun(@(x) animatedline('MaximumNumPoints', 1, 'LineWidth', 5, 'Marker', 'o', 'MarkerFaceColor', 'k'), 1:32);
h = animatedline('Color', 'r', 'LineWidth', 2);

for kk = 1:length(GPS.closest4)
    
    for ii = 1:32
        if ismember(ii, GPS.closest4(kk,:))
            hs(ii).Color = 'r';
            hs(ii).MarkerFaceColor = 'r';
        else 
            hs(ii).Color = 'k';
            hs(ii).MarkerFaceColor = 'k';
        end
        addpoints(hs(ii), GPS.ECEFs(ii,kk,1), GPS.ECEFs(ii,kk,2), GPS.ECEFs(ii,kk,3));
    end
    addpoints(h, simdata.r_ecef(1,kk), simdata.r_ecef(2,kk), simdata.r_ecef(3,kk));
    drawnow limitrate
    
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

