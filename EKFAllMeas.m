function FuseGPSIMU4(fileheader)

gpsdata = readtable(strcat('Data/GPSPositionDirect/', fileheader, '.csv'));
acceldata = readtable(strcat('Data/IMUReadings/', fileheader, '.xlsx'), 'Sheet', 'accel');
gyrodata = readtable(strcat('Data/IMUReadings/', fileheader, '.xlsx'), 'Sheet', 'gyro');
orientdata = readtable(strcat('Data/IMUReadings/', fileheader, '.xlsx'), 'Sheet', 'orientation');

numGPS = size(gpsdata, 1);
numaccel = size(acceldata, 1);
numgyro = size(gyrodata, 1);
numorient = size(orientdata, 1);

latestGPS = gpsdata{1, :};
latestaccel = acceldata{1, :};
latestgyro = gyrodata{1, :};
latestorient = orientdata{1, :};

GPSUTC = latestGPS(1)/(10^3);
AccelUTC = latestaccel(1)/(10^3);
GyroUTC = latestaccel(1)/(10^3);
OrientUTC = latestorient(1)/(10^3);

GPSidx = 1;
Accelidx = 1;
Gyroidx = 1;
Orientidx = 1;

UTCsec = [];

[mintime, idx] = min([GPSUTC, AccelUTC, OrientUTC, GyroUTC]);
UTCsec(1) = mintime;

if idx == 1
    GPSidx = GPSidx + 1;
    latestGPS = gpsdata{GPSidx, :};
elseif idx == 2
    Accelidx = Accelidx + 1;
    latestaccel = acceldata{Accelidx, :};
elseif idx == 3
    Orientidx = Orientidx + 1;
    latestorient = orientdata{Orientidx, :};
elseif idx == 4
    Gyroidx = Gyroidx + 1;
    latestgyro = orientdata{Gyroidx, :};
end

i = 1;

Q_ekf = diag([0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01]);
R_ekf = diag([1000, 1000, 1000, 1, 1, 1, 1, 1, 1, 1, 1, 1]);

mu_t_t(:,1) = zeros(15,1);

P_t_t = diag([5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5]);

while (GPSidx < numGPS) && (Accelidx < numaccel) && (Orientidx < numorient) && (Gyroidx < numgyro)

    i = i+1;

    GPSUTC = latestGPS(1)/(10^3);
    AccelUTC = latestaccel(1)/(10^3);
    OrientUTC = latestorient(1)/(10^3);
    GyroUTC = latestorient(1)/(10^3);
       
    [mintime, idx] = min([GPSUTC, AccelUTC, OrientUTC, GyroUTC]);
    UTCsec(i) = mintime;

    if idx == 1
        GPSidx = GPSidx + 1;
        latestGPS = gpsdata{GPSidx, :};
    elseif idx == 2
        Accelidx = Accelidx + 1;
        latestaccel = acceldata{Accelidx, :};
    elseif idx == 3
        Orientidx = Orientidx + 1;
        latestorient = orientdata{Orientidx, :};
    elseif idx == 4
        Gyroidx = Gyroidx + 1;
        latestorient = gyrodata{Gyroidx, :};
    end

    z = [latestGPS(8:10), latestaccel(9:11), latestorient(3:5), latestgyro(9:11)]';

    deltat = (UTCsec(i) - UTCsec(i-1));

    % Predict
    mu_t_t_minus = get_F(mu_t_t(:,i-1), deltat)*mu_t_t(:,i-1);
    P_t_t_minus = get_F(mu_t_t(:,i-1), deltat)*P_t_t*get_F(mu_t_t(:,i-1), deltat)' + Q_ekf;

    % Update
    y_tilde_t = z - get_H(mu_t_t(:,i-1))*mu_t_t_minus;
    K_t = P_t_t_minus*get_H(mu_t_t(:,i-1))'*inv(R_ekf + get_H(mu_t_t(:,i-1))*P_t_t_minus*get_H(mu_t_t(:,i-1))');
    mu_t_t(:,i) = mu_t_t_minus + K_t*y_tilde_t;
    P_t_t = (eye(15) - K_t*get_H(mu_t_t(:,i-1)))*P_t_t_minus;

end

figure(1)

plot(UTCsec, mu_t_t(1,:))
hold on
plot(UTCsec, mu_t_t(2,:))

figure(2)

plot(mu_t_t(1,:), mu_t_t(2,:))

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