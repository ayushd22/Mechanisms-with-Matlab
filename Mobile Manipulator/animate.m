function animate(robot, tip, mobile, q)

filename = 'test.gif';

i = 1;
num = size(q);
while i < num(1)
    if mod(i, 20) == 1
        % Get joint positions
        base = mobile(i, 1:2);
        tb = mobile(i, 3);
        joint1 = [base(1) + robot.L1*cos(tb+q(i, 3)),...
                  base(2) + robot.L1*sin(tb+q(i, 3))];
        joint2 = [joint1(1) + robot.L2*cos(tb+q(i, 3)+q(i, 4)),...
                  joint1(2) + robot.L2*sin(tb+q(i, 3)+q(i, 4))];
        joint3 = tip(i, 1:2);
        pos = [base; joint1; joint2; joint3];
        
        if i == 1      
            % Initialize plot
            animation = plot(pos(:, 1), pos(:, 2), ...
                    '.-','linewidth',4,'markersize',40,'color','black');
            hold on
            
            th = 0:pi/50:2*pi;
            x = cos(th) + base(1);
            y = sin(th) + base(2);
            circle = plot(x,y,'black','linewidth',2);
            
            axis([-2, 12, -2, 10])
           
            xlabel('X')
            ylabel('Y')
            title('Mobile Manipulator Simulation')

            grid on
            hold on
        else
            % Update robot position
            set(animation,'xdata',pos(:, 1),'ydata',pos(:, 2));
            th = 0:pi/50:2*pi;
            x = cos(th) + base(1);
            y = sin(th) + base(2);
            set(circle, 'xdata', x);
            set(circle, 'ydata', y);
            
            % Plot trajectory
            plot(tip(i, 1), tip(i, 2),'.','markersize',2, 'color','blue')
            plot(mobile(i, 1), mobile(i, 2),'.','markersize',2, 'color','red')
%             pause(0.001)
        end
    end
    
    % Save frame to the gif file
    if mod(i, 100) == 1
        CurrFrame = getframe;   
        im = frame2im(CurrFrame);   
        [A, map] = rgb2ind(im,256);
        if i == 1
            imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime', 0);  
        else
            imwrite(A,map,filename,'gif','WriteMode','append','DelayTime', 0);
        end
    end
        
    i = i + 1;
end
    
end
