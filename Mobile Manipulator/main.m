% Assign values to intrinsic parameters
robot = struct('r',0.2,...
               'L',0.5,...
               'L1',2,...
               'L2',2,...
               'L3',2,...
               'd',1);

% Set up default initial position
phi_1 = 0;
phi_2 = 0;
tb = 0;
t1 = 0;
t2 = pi/2;
t3 = pi/2;
q = [phi_1, phi_2, t1, t2, t3];      % the state vector
Fixedwheel = false;     % whether the wheels are locked

%% Test #1, Straight line, wheels not locked
% Comment the following code if you don't run this test

% % Set up desired velocity and time
% V = [1, 0];
% T = 8;   % total time for each trajectory

%% Test #2, Straight line, wheels locked
% Comment the following code if you don't run this test

% % Set up desired velocity and time
% V = [1, 0];
% T = 5;   % total time for each trajectory
% 
% Fixedwheel = true;     % lock the wheels

%% Test #3, Square box, wheels not locked
% Comment the following code if you don't run this test

% Set up desired velocity and time
V = [1, 0;
     0, 1;
     -1, 0;
     0, -1];
T = [7; 7; 7; 7];   % total time for each trajectory

%% Test #4, Square box, wheels locked
% Comment the following code if you don't run this test

% % Set up desired velocity and time
% V = [1, 0;
%      0, 1;
%      -1, 0;
%      0, -1];
% T = [3; 3; 3; 3];   % total time for each trajectory
% Fixedwheel = true;     % lock the wheels

%% Simulation loop
T_step = 0.001;   % one iteration step
loop = T/T_step;
traj = size(T);
num = 1;

% Compute robot positions
x_mobile = [0, 0, 0];
x_mani = [robot.L1*cos(tb+t1)+robot.L2*cos(tb+t1+t2)+robot.L1*cos(tb+t1+t2+t3),...
         robot.L1*sin(tb+t1)+robot.L2*sin(tb+t1+t2)+robot.L1*sin(tb+t1+t2+t3),...
         t1+t2+t3];
x_end = x_mobile + x_mani;

% Initialize data to store
data_q = ones(sum(loop), 5);
data_end = ones(sum(loop), 3);
data_mobile = ones(sum(loop), 3);

for i = 1:traj(1)     
    for k = 1:loop(i)
        % Store the data
        data_q(num, :) = q;
        data_end(num, :) = x_end;
        data_mobile(num, :) = x_mobile;
        
        % Compute Jacobian and its pseudoinverse
        Jacob = jacobian(robot, q, tb);
        if Fixedwheel == true 
            J = Jacob(1:2, :);
            Jprime = Jacob(1:2, 3:5);
            J_inv = Jprime.'*inv(Jprime*Jprime.');    
        else
            J = Jacob(1:2, :);
            J_inv = J.'*inv(J*J.');
        end

        % Compute joint velocities 
        if Fixedwheel == true
            q_dot = [0; 0; J_inv*V(i, :).'];
        else
            q_dot = J_inv*V(i, :).';
        end
        v_end = Jacob*q_dot;  % velocity of the end-effector
        
        % Update joint positions
        q = q + q_dot.'*T_step; 
        
        % Update robot positions
        tb_dot = robot.r/robot.d*(q_dot(2) - q_dot(1));
        tb = tb + tb_dot*T_step; 
        x_end = x_end + v_end.'*T_step;
        x_mani = [robot.L1*cos(tb+q(3))+robot.L2*cos(tb+q(3)+q(4))+robot.L1*cos(tb+q(3)+q(4)+q(5)),...
                  robot.L1*sin(tb+q(3))+robot.L2*sin(tb+q(3)+q(4))+robot.L1*sin(tb+q(3)+q(4)+q(5)),...
                  q(3)+q(4)+q(5)];
        x_mobile = x_end - x_mani;
       
        num = num + 1;
    end
end

%% Plot
figure(1)
plot(data_q(:, 3:5))
legend('theta1', 'theta2', 'theta3')
xlabel('Step')
ylabel('Theta(t)')
title('State of the Mobile Manipulator')

figure(2)
plot(data_end(:,1), data_end(:,2))
hold on
plot(data_mobile(:,1), data_mobile(:,2))
legend('End-effector', 'Base')
axis([0, 10, 0, 10])
xlabel('X')
ylabel('Y')
title('Path of the Mobile Manipulator')

%% Animation
figure(3)
animate(robot, data_end, data_mobile, data_q);  









