%%simulate_furata()
x =x_sim;  
qube2_rotpen_param
L1 = Lr;
L2 = Lp;
th = x(1,:);
al = x(2,:);
o = [0;0;0];
c = 0.1
f = figure;
for i = 1:length(th)
    p1 = [L1*cos(th(i));L1*sin(th(i));0];
    b1 = [L2*sin(al(i))*sin(th(i));
          -L2*sin(al(i))*cos(th(i));
          L2*cos(al(i))];
    p2 = p1 + b1; 
    hold off
    plot3([o(1),p1(1)],[o(2),p1(2)],[o(3),p1(3)],'k','LineWidth',3);
    grid on
    hold on 
    plot3([p1(1),p2(1)],[p1(2),p2(2)],[p1(3),p2(3)],'k','LineWidth',3);
    view([60 24]);
    %view([90 90]);
    plotcube([c c c],[-c/2 -c/2 -c], 0.5, [1 0 1])
    plot3(p1(1),p1(2),p1(3), 'bo','LineWidth',3);
    plot3(p2(1),p2(2),p2(3), 'bo','LineWidth',3);
    axis([-0.2 0.2 -0.2 0.2 -0.2 0.2])
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    drawnow;
    frame = getframe(f);
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    fn = strcat('anim_gifs/dircol.gif');
    % Write to the GIF File 
    if i == 1 
        imwrite(imind,cm,fn,'gif','DelayTime',0.01, 'Loopcount',inf); 
    else 
        imwrite(imind,cm,fn,'gif','DelayTime',0.01,'WriteMode','append'); 
    end 
    pause(0.01)
end
