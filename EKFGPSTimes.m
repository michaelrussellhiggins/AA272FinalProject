function [UTCsec, mu_pos, ub, lb] = EKFGPSTimes(fileheader)

gpsdata = readtable(strcat('Data/GPSPositionClip2/', fileheader, '.csv'));
acceldata = readtable(strcat('Data/IMUReadings/', fileheader, '.xlsx'), 'Sheet', 'accel');
gyrodata = readtable(strcat('Data/IMUReadings/', fileheader, '.xlsx'), 'Sheet', 'gyro');
orientdata = readtable(strcat('Data/IMUReadings/', fileheader, '.xlsx'), 'Sheet', 'orientation');

gpsdata = gpsdata(~(isnan(gpsdata.utc)),:);

numGPS = size(gpsdata, 1);

UTCsec = [];

UTCsec(1) = gpsdata{1,1}/(10^3);

Q_ekf = diag([0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01]);
R_ekf = diag([1000, 1000, 1000, 1, 1, 1, 1, 1, 1, 1, 1, 1]);

mu_t_t(:,1) = zeros(15,1);

P_t_t = diag([5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5]);

ub(1,1) = mu_t_t(1,1) + 1.96*P_t_t(1,1);
ub(2,1) = mu_t_t(2,1) + 1.96*P_t_t(2,2);
ub(3,1) = mu_t_t(3,1) + 1.96*P_t_t(3,3);
lb(1,1) = mu_t_t(1,1) - 1.96*P_t_t(1,1);
lb(2,1) = mu_t_t(2,1) - 1.96*P_t_t(2,2);
lb(3,1) = mu_t_t(3,1) - 1.96*P_t_t(3,3);

for i = 2:numGPS

    UTCsec(i) = gpsdata{i,1}/(10^3);

    [accelUTC accelidx] = min(abs(acceldata{:,1} - gpsdata{i,1}));
    [orientUTC orientidx] = min(abs(orientdata{:,1} - gpsdata{i,1}));
    [gyroUTC gyroidx] = min(abs(gyrodata{:,1} - gyrodata{i,1}));

    z = [gpsdata{i, 8:10}, acceldata{accelidx, 9:11}, orientdata{orientidx, 3:5}, gyrodata{gyroidx, 9:11}]';

    deltat = (UTCsec(i) - UTCsec(i-1));

    % Predict
    mu_t_t_minus = get_F(mu_t_t(:,i-1), deltat)*mu_t_t(:,i-1);
    P_t_t_minus = get_F(mu_t_t(:,i-1), deltat)*P_t_t*get_F(mu_t_t(:,i-1), deltat)' + Q_ekf;

    % Update
    y_tilde_t = z - get_H(mu_t_t(:,i-1))*mu_t_t_minus;
    K_t = P_t_t_minus*get_H(mu_t_t(:,i-1))'*inv(R_ekf + get_H(mu_t_t(:,i-1))*P_t_t_minus*get_H(mu_t_t(:,i-1))');
    mu_t_t(:,i) = mu_t_t_minus + K_t*y_tilde_t;
    P_t_t = (eye(15) - K_t*get_H(mu_t_t(:,i-1)))*P_t_t_minus;

    ub(1,i) = mu_t_t(1,i) + 1.96*P_t_t(1,1);
    ub(2,i) = mu_t_t(2,i) + 1.96*P_t_t(2,2);
    ub(3,i) = mu_t_t(3,i) + 1.96*P_t_t(3,3);
    lb(1,i) = mu_t_t(1,i) - 1.96*P_t_t(1,1);
    lb(2,i) = mu_t_t(2,i) - 1.96*P_t_t(2,2);
    lb(3,i) = mu_t_t(3,i) - 1.96*P_t_t(3,3);

end

mu_pos = mu_t_t(1:3,:);

end

function F = get_F(mu_cur, deltat)

    F = zeros(15,15);

    F(1,1) = 1;
    F(2,2) = 1;
    F(3,3) = 1;
    F(4,4) = 1;
    F(5,5) = 1;
    F(6,6) = 1;
    F(7,7) = 1;
    F(8,8) = 1;
    F(9,9) = 1;
    F(10,10) = 1;
    F(11,11) = 1;
    F(12,12) = 1;
    F(13,13) = 1;
    F(14,14) = 1;
    F(15,15) = 1;
    
    F(1,4) = deltat;
    F(2,5) = deltat;
    F(3,6) = deltat;
    F(4,7) = deltat;
    F(5,8) = deltat;
    F(6,9) = deltat;
    F(10,13) = deltat;
    F(11,14) = deltat;
    F(12,15) = deltat;

end

function H = get_H(mu_cur)

    H = zeros(12,15);

    theta = mu_cur(10);
    psi = mu_cur(11);
    phi = mu_cur(12);
    
    xdd_e = mu_cur(7);
    xdd_n = mu_cur(8);
    xdd_u = mu_cur(9);

    g = 9.81;

    H(1,1) = 1;
    H(2,2) = 1;
    H(3,3) = 1;

    H(4,7) = cosd(psi)*cosd(phi);
    H(4,8) = -sind(psi)*cosd(phi);
    H(4,9) = -sind(phi);
    H(4,10) = 0;
    H(4,11) = -cosd(psi)*sind(phi)*xdd_e + sind(psi)*sind(phi)*xdd_n - cosd(phi)*xdd_u - cosd(phi)*g;
    H(4,12) = -sind(psi)*cosd(phi)*xdd_e - cosd(psi)*cosd(phi)*xdd_n;

    H(5,7) = sind(psi)*cosd(theta) - cosd(psi)*sind(phi)*sind(theta);
    H(5,8) = cosd(psi)*cosd(theta) + sind(psi)*sind(phi)*sind(theta);
    H(5,9) = -cosd(phi)*sind(theta);
    H(5,10) = (-sind(psi)*sind(theta)-cosd(psi)*sind(phi)*cosd(theta))*xdd_e + (-cosd(psi)*sind(theta) + sind(psi)*sind(phi)*cosd(theta))*xdd_n + (-cosd(phi)*cosd(theta))*xdd_u - cosd(phi)*cosd(theta)*g;
    H(5,11) = (-cosd(psi)*cosd(phi)*sind(theta))*xdd_e + (sind(psi)*cosd(phi)*sind(theta))*xdd_n + (sind(phi)*sind(theta))*xdd_u + (sind(phi)*sind(theta))*g;
    H(5,12) = (cosd(psi)*cosd(theta) + sind(psi)*sind(phi)*sind(theta))*xdd_e + (-sind(psi)*cosd(theta) + cosd(psi)*sind(phi)*sind(theta))*xdd_n;

    H(6,7) = sind(psi)*sind(theta) + cosd(psi)*sind(phi)*cosd(theta);
    H(6,8) = cosd(psi)*sind(theta) - sind(psi)*sind(phi)*cosd(theta);
    H(6,9) = cosd(phi)*cosd(theta);
    H(6,10) = (sind(psi)*cosd(theta) - cosd(psi)*sind(phi)*sind(theta))*xdd_e + (cosd(psi)*cosd(theta) + sind(psi)*sind(phi)*sind(theta))*xdd_n + (-cosd(phi)*sind(theta))*xdd_u + (-cosd(phi)*sind(theta))*g;
    H(6,11) = (cosd(psi)*cosd(phi)*cosd(theta))*xdd_e + (-sind(psi)*cosd(phi)*cosd(theta))*xdd_n + (-sind(phi)*cosd(theta))*xdd_u + (-sind(phi)*cosd(theta))*g;
    H(6,12) = (cosd(psi)*sind(theta) - sind(psi)*sind(phi)*cosd(theta))*xdd_e + (-sind(psi)*sind(theta) - cosd(psi)*sind(phi)*cosd(theta))*xdd_n;

    H(7,10) = 1;
    H(8,11) = 1;
    H(9,12) = 1;
    H(10,13) = 1;
    H(11,14) = 1;
    H(12,15) = 1;

end