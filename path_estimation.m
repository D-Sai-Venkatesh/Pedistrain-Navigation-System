clear all;

%% getting data
data  = load("my_case.mat");
data = data.data;
gyro_bias = [-0.09  -0.010   -0.0077]';

%% init arrays required
timestamp = data.AccR(5:length(data.AccR)-100,1)'; 
acc_s = data.AccR(5:length(data.AccR)-100,2:4)'; 
gyro_s = data.AngVelR(5:length(data.AngVelR)-100,2:4)'; 
data_size = length(timestamp);
acc_mag = sqrt(acc_s(1,:)'.^2 + acc_s(2,:)'.^2 + acc_s(3,:)'.^2);
gyro_mag = sqrt(gyro_s(1,:)'.^2 + gyro_s(2,:)'.^2 + gyro_s(3,:)'.^2);
stance = zeros(data_size,1);
g = 9.8; 


%% init roll pitch yaw
pitch = -asin(acc_s(1,1)/g);
roll = atan(acc_s(2,1)/acc_s(3,1));
yaw = 0;


%% init rotation matrix (it will rotate coordinates from body frame to navigational frame)
C = [cos(pitch)*cos(yaw) (sin(roll)*sin(pitch)*cos(yaw))-(cos(roll)*sin(yaw)) (cos(roll)*sin(pitch)*cos(yaw))+(sin(roll)*sin(yaw));
    cos(pitch)*sin(yaw)  (sin(roll)*sin(pitch)*sin(yaw))+(cos(roll)*cos(yaw))  (cos(roll)*sin(pitch)*sin(yaw))-(sin(roll)*cos(yaw));
    -sin(pitch) sin(roll)*cos(pitch) cos(roll)*cos(pitch)];
C_prev = C;

%% init array to store the outputs
acc_n = nan(3, data_size);
acc_n(:,1) = C*acc_s(:,1);

vel_n = nan(3, data_size);
vel_n(:,1) = [0 0 0]';

pos_n = nan(3, data_size);
pos_n(:,1) = [0 0 0]';

%% init coovarience matrix and other matrices required for the kalman filter
P = zeros(9);
sigma_omega = 0.3; sigma_a = 0.4;

%% init measurement matrix
H = [zeros(3) zeros(3) eye(3)];

%% measurement noise matrix
sigma_v = 1e-2;
R = diag([sigma_v sigma_v sigma_v]).^2;

%% init variables required for stance detction
priv_ZUPT = 0;
c1=0;
c2=0;
c3=0;
count = 0; % number of ZUPTs

for t = 2:data_size
    dt = timestamp(t) - timestamp(t-1);
    
    %% removing gyro bias
    gyro_s1 = gyro_s(:,t) - gyro_bias;
    
    ang_rate_matrix = [0   -gyro_s1(3)   gyro_s1(2);
        gyro_s1(3)  0   -gyro_s1(1);
        -gyro_s1(2)  gyro_s1(1)  0];
    
    %% update equations
    C = C_prev*(2*eye(3)+(ang_rate_matrix*dt))/(2*eye(3)-(ang_rate_matrix*dt));
    acc_n(:,t) = 0.5*(C + C_prev)*acc_s(:,t);
    vel_n(:,t) = vel_n(:,t-1) + ((acc_n(:,t) - [0; 0; g] )+(acc_n(:,t-1) - [0; 0; g]))*dt/2;
    pos_n(:,t) = pos_n(:,t-1) + (vel_n(:,t) + vel_n(:,t-1))*dt/2;
    
    S = [0  -acc_n(3,t)  acc_n(2,t);
        acc_n(3,t)  0  -acc_n(1,t);
        -acc_n(2,t) acc_n(1,t) 0];
    
    F = [eye(3)  zeros(3,3)    zeros(3,3);
        zeros(3,3)   eye(3)  dt*eye(3);
        -dt*S  zeros(3,3)    eye(3) ];
    
    % process noise matrix
    Q = diag([sigma_omega sigma_omega sigma_omega 0 0 0 sigma_a sigma_a sigma_a]*dt).^2;
    
    % propagation of covarince matrix
    P = F*P*F' + Q;
    
   %% stance detection 
   % condition 1
    if(9<acc_mag(t) && acc_mag(t)<11)
         c1 = 1;
     else
         c1 = 0;
     end
     %% condition 2
    window = 15;
    if ((t-window) > 0)
        up = acc_mag((t-window):t);
    else
        up = acc_mag(1:t);
    end
    if((t+window) <= length(data))
        down = acc_mag(t:(t+window));
    else
        down = acc_mag(t:length(data));
    end
    up_down = [up ; down];
    var_a = var(up_down);
    
    if(var_a > 3)
         c2 = 1;
    else
         c2 = 0;
    end
    
    %% condition 3
    if(gyro_mag(t)*(180/pi) < 50)
         c3 = 1;
    else
         c3 = 0;
    end
    
    %% checking if all conditions are satisfied
    if(c1==1)
        if(c2==1)
            if(c3==1)
                 
                if(true)
                    stance(t) = 1;
                    count = count +1 ;
                    priv_ZUPT = timestamp(t);  
                end
            end
        end
    end
    
    
    
    
   %% ZUPT 
    if (stance(t) ==1)
        % kalman gain
        K = (P*(H)')/((H)*P*(H)' + R);
        
        % state correction
        delta_x = K*vel_n(:,t);
        
        % covariace correction
        P = (eye(9) - K*(H)) * P * (eye(9) - K*(H))' + K*R*K'; 
        
        % corrections
        attitude_error = delta_x(1:3);
        pos_error = delta_x(4:6);
        vel_error = delta_x(7:9);
        
        ang_matrix = -[0   -attitude_error(3,1)   attitude_error(2,1);
            attitude_error(3,1)  0   -attitude_error(1,1);
            -attitude_error(2,1)  attitude_error(1,1)  0];
        
        % rotation matrix correction
        C = (2*eye(3)+(ang_matrix))/(2*eye(3)-(ang_matrix))*C;
        
        vel_n(:,t)=vel_n(:,t)-vel_error;
        pos_n(:,t)=pos_n(:,t)-pos_error;
    end
     
    C_prev = C; 
    
end

scatter(timestamp, stance);

figure;
box on;
hold on;
angle = 90+90+90; 
rotation_matrix = [cosd(angle) -sind(angle);
    sind(angle) cosd(angle)];
pos_r = zeros(2,data_size);
for idx = 1:data_size
    pos_r(:,idx) = rotation_matrix*[pos_n(1,idx) pos_n(2,idx)]';
end
plot(pos_r(1,:),pos_r(2,:),'LineWidth',2,'Color','r');
start = plot(pos_r(1,1),pos_r(2,1),'Marker','X','LineWidth',2);
stop = plot(pos_r(1,end),pos_r(2,end),'Marker','o','LineWidth',2);

xlabel('x (m)');
ylabel('y (m)');
title('Estimated 2D path');
legend([start;stop],'Start','End');
axis equal;
grid;
hold off;
