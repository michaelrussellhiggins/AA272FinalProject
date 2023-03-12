function FuseGPSIMU2(fileheader)

gpsdata = readtable(strcat('Data/GPSPositionDirect/', fileheader, '.csv'));
acceldata = readtable(strcat('Data/IMUReadings/', fileheader, '.xlsx'), 'Sheet', 'accel');
orientdata = readtable(strcat('Data/IMUReadings/', fileheader, '.xlsx'), 'Sheet', 'orientation');

numGPS = size(gpsdata, 1);

UTCsec = [];

UTCsec(1) = gpsdata{1,1}/(10^3);

Q_ekf = diag([1, 1, 1, 1, 1, 1, 1, 1, 1]);
R_ekf = diag([5, 5, 5, 0.1, 0.1, 0.1]);

mu_t_t(:,1) = zeros(9,1);

P_t_t = diag([5, 5, 5, 1, 1, 1, 1, 1, 1]);

for i = 2:numGPS

    UTCsec(i) = gpsdata{i,1}/(10^3);

    [accelUTC accelidx] = min(abs(acceldata{:,1} - gpsdata{i,1}));
    [orientUTC orientidx] = min(abs(orientdata{:,1} - gpsdata{i,1}));

    a_XYZ = acceldata{accelidx, 9:11}';

    psi = orientdata{orientidx, 3};
    phi = orientdata{orientidx, 4};
    theta = orientdata{orientidx, 5};

    R_psi = [cosd(psi), sind(psi), 0; -sind(psi), cosd(psi), 0; 0, 0, 1];
    R_phi = [cosd(phi), 0, sind(phi); 0, 1, 0; -sind(phi), 0, cosd(phi)];
    R_theta = [1, 0, 0; 0, cosd(theta), sind(theta); 0, -sind(theta), cosd(theta)];

    a_ENU = R_psi*R_phi*R_theta*a_XYZ - [0; 0; 9.81];

    z = [gpsdata{i, 8:10}, a_ENU']';

    deltat = (UTCsec(i) - UTCsec(i-1));

    % Predict
    mu_t_t_minus = mu_t_t(:,i-1);
    P_t_t_minus = get_F(mu_t_t(:,i-1), deltat)*P_t_t*get_F(mu_t_t(:,i-1), deltat)' + Q_ekf;

    % Update
    y_tilde_t = (z - meas_mdl(mu_t_t(:,i-1))');
    K_t = P_t_t_minus*get_H(mu_t_t(:,i-1))'*inv(R_ekf + get_H(mu_t_t(:,i-1))*P_t_t_minus*get_H(mu_t_t(:,i-1))');
    mu_t_t(:,i) = mu_t_t_minus + K_t*y_tilde_t;
    P_t_t = (eye(9) - K_t*get_H(mu_t_t(:,i-1)))*P_t_t_minus;

end

figure(1)

plot(UTCsec, mu_t_t(1,:))
hold on
plot(UTCsec, mu_t_t(2,:))

figure(2)

plot(mu_t_t(1,:), mu_t_t(2,:))

end

function noisy_meas = meas_mdl(mu_cur)

    noisy_meas = [];

    noisy_meas(1) = mu_cur(1);
    noisy_meas(2) = mu_cur(2);
    noisy_meas(3) = mu_cur(3);
    noisy_meas(4) = mu_cur(7);
    noisy_meas(5) = mu_cur(8);
    noisy_meas(6) = mu_cur(9);

end

function F = get_F(mu_cur, deltat)

    F = zeros(9,9);
    F(1,1) = 1;
    F(2,2) = 1;
    F(3,3) = 1;
    F(4,4) = 1;
    F(5,5) = 1;
    F(6,6) = 1;
    F(7,7) = 1;
    F(8,8) = 1;
    F(9,9) = 1;
    F(1,4) = deltat;
    F(2,5) = deltat;
    F(3,6) = deltat;
    F(4,7) = deltat;
    F(5,8) = deltat;
    F(5,9) = deltat;

end

function H = get_H(mu_cur)

    H = zeros(6,9);
    H(1,1) = 1;
    H(2,2) = 1;
    H(3,3) = 1;
    H(4,7) = 1;
    H(5,8) = 1;
    H(6,9) = 1;

end