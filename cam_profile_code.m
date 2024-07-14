%Project Title : Design of Cam with its profile generated in MATLAB

%Lets clear command window before running program
clc
clear all

%Lift and Base Circle Radius
h = input('Lift of Follower (mm) \n');
r = input('Base Radius (mm) \n');

%Ascent and Descent Angles
ascent = input('Ascent Angle \n');
dwell_1 = input('Dwell Angle \n');
descent =  input('Descent Angle \n');
dwell_2 = 360-(ascent + dwell_1 + descent);

%Additional Angles for our convenience
after_ascent = ascent + dwell_1;
after_descent = after_ascent + descent;

%Cam Angle
theta = linspace(0,360,361);

%Let us consider SHM Motion, Uniform Velocity, Cycloidal 
%Adding Motion Selector Window
motion_list = {'Simple Harmonic Motion','Uniform Velocity','Cycloidal'};
[ascent_motion] = listdlg('ListSize',[270,100],'Name','Select Ascent Motion','SelectionMode','single','ListString',motion_list);
[descent_motion] = listdlg('ListSize',[270,100],'Name','Select Descent Motion','SelectionMode','single','ListString',motion_list);

%Ascent Motion Conditions
if ascent_motion == 1 %SHM
    h_ascent = (0.5*h).*(1-cosd((180/ascent).*theta(theta<ascent)));
elseif ascent_motion == 2 %Uniform Velocity
    h_ascent = (h/ascent).*theta(theta<ascent);
elseif ascent_motion == 3 %Cycloidal
    h_ascent = (h/pi)*(((pi/ascent).*theta(theta<ascent)) - 0.5*sind((2*180/ascent).*theta(theta<ascent)));
end

%Descent Motion Conditions
if descent_motion == 1 %SHM
    h_descent = h-((0.5*h).*(1 - cosd((180/descent).*theta(theta<=descent))));
elseif descent_motion == 2 %Uniform Velocity
    h_descent = h- (h/descent).*theta(theta<=descent);
elseif descent_motion == 3 %Cycloidal
    h_descent = h - ((h/pi)*(((pi/descent).*theta(theta<=descent)) - 0.5*sind((2*180/descent).*theta(theta<=descent))));
end

%Plotting Cam Angle Vs Lift
plot(theta(theta<ascent),h_ascent,theta(theta>=after_ascent & theta<=after_descent),h_descent);
set(gca,'XTick',(0:10:after_descent))
set(gca,'YTick',(0:5:h))
set(gcf,'position',[0,0,6000,400])
title('Cam Angle Vs Lift')
xlabel('Cam Angle (degrees)');
ylabel('Lift of follower (mm)');

%Lift during Dwell
h_dwell1 = ones(1,dwell_1).*h;
h_dwell2 = zeros(1,dwell_2);

%Defining radii during different phases of Cam
r1 = r + h_ascent;
r2 = r + h_dwell1;
r3 = r + h_descent;
r4 = r + h_dwell2;

%Joining all radii
r = [r1 r2 r3 r4];

%Convert theta to radians
theta_radians = deg2rad(theta);

%Plotting Cam Profile
figure
polarplot(theta_radians,r);
set(gca,'ThetaZeroLocation','top')

%Converting Polar to Cartesian Coordinate System
[x,y] = pol2cart(theta_radians,r);
x_cord = transpose(x);
y_cord = transpose(y);
z_cord = zeros(361,1);

cam_profile = [x_cord y_cord z_cord];

%Export Cam Profile to Excel as XYZ Coordinates
writematrix(cam_profile,'cam_profile.xls')