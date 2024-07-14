clc; clear all; close all;
options.Interpreter = 'Latex';
% save('four_bar_2021.mat')
% break
tic
L1=45;L2=90;L3=100;L4=120;L5=150;L6=120;L7=150;
%Grashof check
ang_speed=12*2*pi/60;
t=0:0.1:10;
theta=t*ang_speed;
Ax=0; Ay=0;
Dx=L4; Dy=0;
Ex=L4+L5; Ey=0;
for i=1:length(t)
Bx(i)=L1*cos(theta(i));
By(i)=L1*sin(theta(i));
%______________Numerical Solution_________________________
syms beta phi alfa gama
eqn1=L1*sin(theta(i))+L2*sin(beta)-L3*sin(phi)==0;
eqn2=L1*cos(theta(i))+L2*cos(beta)-L3*cos(phi)-L4==0;
sol_ang=solve(eqn1,eqn2,[beta,phi]);
beta_sol(i,:)=eval(sol_ang.beta);
phi_sol(i,:)=eval(sol_ang.phi);
eqn3=L3*sin(phi_sol(i+1))+L7*sin(gama)-L6*sin(alfa)==0;
eqn4=L3*cos(phi_sol(i+1))+L7*cos(gama)-L6*cos(alfa)-L5==0;
sol_ang1=solve(eqn3,eqn4,[gama,alfa]);
gama_sol(i,:)=eval(sol_ang1.gama);
alfa_sol(i,:)=eval(sol_ang1.alfa);
end
%______________Solution Selection__________________________
% Cx(i)=Bx(i)+L2*cos((beta_sol(i,1)));
% Cy(i)=By(i)+L2*sin((beta_sol(i,1)));
% Cx(i)=Bx(i)+L2*cos((beta_sol1(i)));
% Cy(i)=By(i)+L2*sin((beta_sol1(i)));



%______________Solution Selection__________________________
beta_sol1=max(beta_sol');
phi_sol1=max(phi_sol');
Cx=Bx+L2*cos((beta_sol1));
Cy=By+L2*sin((beta_sol1));
gama_sol1=max(gama_sol');
alfa_sol1=max(alfa_sol');
Fx=Cx+L7*cos((gama_sol1));
Fy=Cy+L7*sin((alfa_sol1));
% Cx=Bx'+L2*cos(abs(beta_sol(:,1)));
% Cy=By'+L2*sin(abs(beta_sol(:,1)));
%toc
%Velocities
Vel_F_x=diff(Fx)./diff(t);
Vel_F_y=diff(Fy)./diff(t);
Vel_F=sqrt(Vel_F_x.^2+Vel_F_y.^2);
Ang_vel_F=Vel_F/(L6);
%Accelerations
Acc_F_x=diff(Vel_F_x)./diff(t(1:length(t)-1));
Acc_F_y=diff(Vel_F_y)./diff(t(1:length(t)-1));
Acc_F=sqrt(Acc_F_x.^2+Acc_F_y.^2);
Ang_Acc_F=Acc_F/(L6);
%Jerk
Jerk_F_x=diff(Acc_F_x)./diff(t(1:length(t)-2));
Jerk_F_y=diff(Acc_F_y)./diff(t(1:length(t)-2));
Jerk_F=sqrt(Jerk_F_x.^2+Jerk_F_y.^2);


clc
close all
% for j=1%:3
for i=1:length(Jerk_F)
set(gcf,'units','inches','position',[.5,.5,10,6])
subplot(4,2,[1 3 5 7])
hold on
P1=plot(Bx(1:i),By(1:i),'k','markersize',5);
P2=plot(Cx(1:i),Cy(1:i),'b-','markersize',5);
P3=plot(Fx(1:i),Fy(1:i),'g-','markersize',5);
h1=plot([Ax Bx(i)],[Ay By(i)],'m','linewidth',6);
h2=plot([Bx(i) Cx(i)],[By(i) Cy(i)],'c','linewidth',6);
h3=plot([Cx(i) Dx],[Cy(i) Dy],'color',[255 165 0]/255,'linewidth',6);
h4=plot([Cx(i) Fx(i)],[Cy(i) Fy(i)],'c','linewidth',6);
h5=plot([Fx(i) Ex],[Fy(i) Ey],'color',[255 165 0]/255,'linewidth',6);
plot(0,-5,'b^','markersize',15,'markerfacecolor','b')
plot(Dx,-5,'b^','markersize',15,'markerfacecolor','b')
plot(Ex,-5,'b^','markersize',15,'markerfacecolor','b')
viscircles([0,0],2);
viscircles([Dx,0],2);
viscircles([Ex,0],2);
xlim([-50 400])
ylim([-50 400])
str3='B';
B_text=text(Bx(i),By(i)+1,str3);
str1='F';
F_text=text(Fx(i),Fy(i)+1,str1);
str2='C';
C_text=text(Cx(i),Cy(i)+1,str2);

axis equal
box on; grid on
%__________Angular Position of C______________
subplot(4,2,2)
plot(t(1:i),phi_sol1(1:i)*180/pi,'color',[255 165 0]/255,'linewidth',3)

xlim([0 10])
ylim([90 180])
ylabel('\phi [deg]')
xlabel('time [sec]')
title('\phi angle of point F')
grid on; box on
%__________Velocity of C______________
subplot(4,2,4)
plot(t(1:i),Vel_F(1:i),'b','linewidth',3)
xlim([0 10])
ylim([0 200])
ylabel('v [mm/sec]')
xlabel('time [sec]')
title('Velocity of point F')
grid on
%__________Acceleration of C______________
subplot(4,2,6)
plot(t(1:i),Acc_F(1:i),'r','linewidth',3)
xlim([0 10])
ylim([0 300])
ylabel('a [mm/sec^2]')
xlabel('time [sec]')
title('Acceleration of point F')
grid on
%__________Jerk of C______________
subplot(4,2,8)
plot(t(1:i),Jerk_F(1:i),'k','linewidth',3)
xlim([0 10])
ylim([0 500])
ylabel('Jerk [mm/sec^3]')
xlabel('time [sec]')
title('Jerk of point F')
grid on
%______________________________________
pause(.01)
% delete(P1)
delete(P2)
delete(P3)
% delete(P4)
delete(h1)
delete(h2)
delete(h3)
delete(h4)
delete(h5)
delete(F_text)
delete(C_text)
delete(B_text)
end
% end
set(gcf,'units','inches','position',[.5,.5,10,6])
subplot(2,2,[1 3])
hold on
P1=plot(Bx(1:i),By(1:i),'k','markersize',5);
P2=plot(Cx(1:i),Cy(1:i),'b','markersize',5);
P3=plot(Fx(1:i),Fy(1:i),'g','markersize',5);
h1=plot([0 Bx(i)],[0 By(i)],'m','linewidth',6);
h2=plot([Bx(i) Cx(i)],[By(i) Cy(i)],'c','linewidth',6);
h3=plot([Cx(i) Dx],[Cy(i) Dy],'y','linewidth',6);
h4=plot([Cx(i) Fx(i)],[Cy(i) Fy(i)],'c','linewidth',6);
h5=plot([Fx(i) Ex],[Fy(i) Ey],'color',[255 165 0]/255,'linewidth',6);

plot(0,-5,'b^','markersize',15,'markerfacecolor','b')
plot(Dx,-5,'b^','markersize',15,'markerfacecolor','b')
plot(Ex,-5,'b^','markersize',15,'markerfacecolor','b')
viscircles([0,0],2);
viscircles([Dx,0],2);
viscircles([Ex,0],2);
xlim([-50 400])
ylim([-50 400])
str3='B';
B_text=text(Bx(i),By(i)+1,str3);
str1='F';
F_text=text(Fx(i),Fy(i)+1,str1);
str2='C';
C_text=text(Cx(i),Cy(i)+1,str2);
axis equal
box on; grid on
xlabel('x axis [mm]')
ylabel('y axis [mm]')
title('Simulation')
legend([P1,P2,P3],'Point B trajectory','Point C trajectory','Point F trajectory')
sgtitle('Six Bar Mechanism')
toc
