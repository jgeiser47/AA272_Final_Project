%==========================================================================
%
% gg_stability_plot  Plot of regions of varying stability under gravity
% gradient torque together with the satellite's stability.
%
%   gg_stability_plot(kR,kT,pp)
%
% Author: Tamas Kis
% Last Update: 2021-09-09
%
% REFERENCES:
%   [1] D'Amico, "Gravity Gradient Torque, Stability, Damping", AA 279C 
%       Lecture 8 Slides (p. 7)
%   [2] Wertz, "Spacecraft Attitude Determination and Control" (pp. 611-
%       612)
%   [3] Kaplan, "Modern Spacecraft Dynamics and Control" (pp. 203-204)
%
%--------------------------------------------------------------------------
%
% ------
% INPUT:
% ------
%   kR      - (1×1 double) R coefficient of principal moments of inertia
%   kT      - (1×1 double) T coefficient of principal moments of inertia
%   pp      - (struct) plot parameters
%
%==========================================================================
function gg_stability_plot(kR,kT,pp)
    
    % equations for implicit line
    f1 = @(x,y) 1+(3.*x)+(x.*y)+(4*sqrt(x.*y));
    f2 = @(x,y) 1+(3.*x)+(x.*y)-(4*sqrt(x.*y));

    % sets up kR domain
    kR1 = (-0.2:0.001:0)';
    kR2 = (-1:0.001:0)';

    % array preallocation
    kT1 = zeros(size(kR1));
    kT2 = zeros(size(kR2));

    % original implicit line
    for i = 1:length(kR1)
        kT1(i)= bisection_method(@(x)f1(x,kR1(i)),-1,-1/3);
    end
    for i = 1:length(kR2)
        kT2(i)= bisection_method(@(x)f2(x,kR2(i)),-1/3,0);
    end

    % finds intersection of original implicit line with kR=kT line
    fint = @(x)interp1(kT2,kR2-kT2,x);
    x_int = bisection_method(fint,-0.3,-0.1);

    % checks for index where intersection occurs
    j = 0;
    for i = 1:length(kT2)
        if kT2(i) > x_int
            j = j+1;
        end
    end

    % portion after intersection
    kT3 = kT2(1:j);
    kR3 = kR2(1:j);

    % portion before intersection
    kT2 = kT2(j:end);
    kR2 = kR2(j:end);

    % initializes figure
    figure('position',pp.plot_position);
    hold on;
    
    % draws regions of stability
    fill([-1,-1,1,1],[-1,1,1,-1],'w','linewidth',pp.pp.line_width,...
        'handlevisibility','off');
    fill([-1,-1,1],[-1,1,1],'y','linewidth',pp.line_width);
    rectangle('position',[-1,0,1,1],'facecolor','g','linewidth',...
        pp.line_width);
    rectangle('position',[-1,-1,1,1],'facecolor','b','linewidth',...
        pp.line_width);
    area([0,1],[-1,-1],'facecolor','b','linewidth',pp.line_width);
    area([-1,0],[-1,0],'facecolor','g','linewidth',pp.line_width);
    area([-1,kT1(1)],[kR1(1),kR1(1)],'facecolor','y','linewidth',...
        pp.line_width,'handlevisibility','off');
    area(kT1,kR1,'facecolor','y','linewidth',pp.line_width,...
        'handlevisibility','off');
    area(kT2,kR2,'facecolor','y','linewidth',pp.line_width,...
        'handlevisibility','off');
    area(kT3,kR3,'facecolor','w','linewidth',pp.line_width,...
        'handlevisibility','off');
    area([kT3(1),0],[-1,-1],'facecolor','w','linewidth',pp.line_width,...
        'handlevisibility','off');
    area([kT2(1),0],[kR2(1),0],'facecolor','y','linewidth',...
        pp.line_width,'handlevisibility','off');
    rectangle('position',[-1,-1,1,1],'facecolor','none','linewidth',...
        pp.line_width);
    rectangle('position',[0,-1,1,1],'facecolor','none','linewidth',...
        pp.line_width);
    
    % satellite
    plot(kT,kR,'r*','markersize',10,'linewidth',pp.line_width);
    hold off;
    
    % axis formatting
    set(gca,'XTickLabel',[]);
    set(gca,'YTickLabel',[]);
    xlabel('$k_{R}$','interpreter','latex','fontsize',pp.axis_font_size);
    ylabel('$k_{T}$','interpreter','latex','fontsize',pp.axis_font_size);
    axis equal;
    axis([-1,1,-1,1]);
    
    % legend
    legend('unstable pitch','unstable yaw and roll',...
        'unstable yaw, roll, and pitch','satellite','interpreter',...
        'latex','fontsize',pp.legend_font_size,'location','northwest');

end