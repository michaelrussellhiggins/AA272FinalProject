function [UTCsec, mu_pos, ub, lb] = KFAllTimes(fileheader)

gpsdata = readtable(strcat('Data/GPSPositionClip2/', fileheader, '.csv'));
acceldata = readtable(strcat('Data/IMUReadingsClip/', fileheader, '.xlsx'), 'Sheet', 'accel');
orientdata = readtable(strcat('Data/IMUReadingsClip/', fileheader, '.xlsx'), 'Sheet', 'orientation');

gpsdata = gpsdata(~(isnan(gpsdata.utc)),:);

numGPS = size(gpsdata, 1);
numaccel = size(acceldata, 1);
numorient = size(orientdata, 1);

latestGPS = gpsdata{1, :};
latestaccel = acceldata{1, :};
latestorient = orientdata{1, :};

GPSUTC = latestGPS(1)/(10^3);
AccelUTC = latestaccel(1)/(10^3);
OrientUTC = latestorient(1)/(10^3);

UTCsec = [];

[maxtime, idx] = max([GPSUTC, AccelUTC, OrientUTC]);
UTCsec(1) = maxtime;

[GPSUTC firstgps] = min(abs(gpsdata{:,1} - UTCsec(1)*(10^3)));
[accelUTC firstaccel] = min(abs(acceldata{:,1} - UTCsec(1)*(10^3)));
[orientUTC firstorient] = min(abs(orientdata{:,1} - UTCsec(1)*(10^3)));

if (gpsdata{firstgps,1} - UTCsec(1)*(10^3)) < 0
    GPSidx = firstgps + 2;
else
    GPSidx = firstgps + 1;
end
if (acceldata{firstaccel,1} - UTCsec(1)*(10^3)) < 0
    Accelidx = firstaccel + 2;
else
    Accelidx = firstaccel + 1;
end
if (orientdata{firstorient,1} - UTCsec(1)*(10^3)) < 0
    Orientidx = firstorient + 2;
else
    Orientidx = firstorient + 1;
end

latestGPS = gpsdata{GPSidx - 1, :};
latestaccel = acceldata{Accelidx - 1, :};
latestorient = orientdata{Orientidx - 1, :};

i = 1;

Q_ekf = diag([0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]);
R_ekf = diag([1000, 1000, 1000, 1, 1, 1]);

mu_t_t(:,1) = zeros(9,1);

P_t_t = diag([5, 5, 5, 1, 1, 1, 1, 1, 1]);

ub(1,1) = mu_t_t(1,1) + 1.96*P_t_t(1,1);
ub(2,1) = mu_t_t(2,1) + 1.96*P_t_t(2,2);
ub(3,1) = mu_t_t(3,1) + 1.96*P_t_t(3,3);
lb(1,1) = mu_t_t(1,1) - 1.96*P_t_t(1,1);
lb(2,1) = mu_t_t(2,1) - 1.96*P_t_t(2,2);
lb(3,1) = mu_t_t(3,1) - 1.96*P_t_t(3,3);

while (GPSidx <= numGPS) || (Accelidx <= numaccel) || (Orientidx <= numorient)

    i = i+1;

    if GPSidx > numGPS
        GPSUTC = NaN;
    else
        GPSUTC = latestGPS(1)/(10^3);
    end
    if Accelidx > numaccel
        AccelUTC = NaN;
    else
        AccelUTC = latestaccel(1)/(10^3);
    end
    if Orientidx > numorient
        OrientUTC = NaN;
    else
        OrientUTC = latestorient(1)/(10^3);
    end
       
    [mintime, idx] = min([GPSUTC, AccelUTC, OrientUTC]);

    if idx == 1
        latestGPS = gpsdata{GPSidx, :};
        GPSidx = GPSidx + 1;
    elseif idx == 2
        latestaccel = acceldata{Accelidx, :};
        Accelidx = Accelidx + 1;
    elseif idx == 3
        latestorient = orientdata{Orientidx, :};
        Orientidx = Orientidx + 1;
    end

    UTCsec(i) = mintime;

    a_XYZ = latestaccel(9:11)';

    psi = latestorient(1, 3);
    phi = latestorient(1, 4);
    theta = latestorient(1, 5);

    R_psi = [cosd(psi), sind(psi), 0; -sind(psi), cosd(psi), 0; 0, 0, 1];
    R_phi = [cosd(phi), 0, sind(phi); 0, 1, 0; -sind(phi), 0, cosd(phi)];
    R_theta = [1, 0, 0; 0, cosd(theta), sind(theta); 0, -sind(theta), cosd(theta)];

    a_ENU = R_psi*R_phi*R_theta*a_XYZ - [0; 0; 9.81];

    z = [latestGPS(8:10), a_ENU']';

    deltat = (UTCsec(i) - UTCsec(i-1));

    % Predict
    mu_t_t_minus = get_F(mu_t_t(:,i-1), deltat)*mu_t_t(:,i-1);
    P_t_t_minus = get_F(mu_t_t(:,i-1), deltat)*P_t_t*get_F(mu_t_t(:,i-1), deltat)' + Q_ekf;

    % Update
    y_tilde_t = z - meas_mdl(mu_t_t_minus)';
    K_t = P_t_t_minus*get_H(mu_t_t_minus)'*inv(R_ekf + get_H(mu_t_t_minus)*P_t_t_minus*get_H(mu_t_t_minus)');
    mu_t_t(:,i) = mu_t_t_minus + K_t*y_tilde_t;
    P_t_t = (eye(9) - K_t*get_H(mu_t_t_minus))*P_t_t_minus;

    ub(1,i) = mu_t_t(1,i) + 1.96*P_t_t(1,1);
    ub(2,i) = mu_t_t(2,i) + 1.96*P_t_t(2,2);
    ub(3,i) = mu_t_t(3,i) + 1.96*P_t_t(3,3);
    lb(1,i) = mu_t_t(1,i) - 1.96*P_t_t(1,1);
    lb(2,i) = mu_t_t(2,i) - 1.96*P_t_t(2,2);
    lb(3,i) = mu_t_t(3,i) - 1.96*P_t_t(3,3);

end

mu_pos = mu_t_t(1:3,:);

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
    F(6,9) = deltat;

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