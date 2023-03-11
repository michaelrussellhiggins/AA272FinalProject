function FuseGPSIMU3(fileheader)

gpsdata = readtable(strcat('Data/GPSPositionDirect/', fileheader, '.csv'));
acceldata = readtable(strcat('Data/IMUReadings/', fileheader, '.xlsx'), 'Sheet', 'accel');
magdata = readtable(strcat('Data/IMUReadings/', fileheader, '.xlsx'), 'Sheet', 'mag');
gyrodata = readtable(strcat('Data/IMUReadings/', fileheader, '.xlsx'), 'Sheet', 'gyro');
orientdata = readtable(strcat('Data/IMUReadings/', fileheader, '.xlsx'), 'Sheet', 'orientation');

numGPS = size(gpsdata, 1);
numaccel = size(acceldata, 1);
numgyro = size(gyrodata, 1);
nummag = size(magdata, 1);
numorient = size(orientdata, 1);

filt = insfilterAsync;

refloc = [0, 0, 0];

filt.ReferenceLocation = refloc;
filt.State = [1; zeros(27,1)];

Rmag = 80;
Rvel = 0.0464;
Racc = 800;
Rgyro = 1e-4;
Rpos = 34;

p = zeros(numaccel, 3);
q = zeros(numaccel, 1, 'quaternion');

GPSidx = 1;
magidx = 1;

imuFs = 125;

for i = 1:numaccel

    predict(filt, 1./imuFs);

    fuseaccel(filt, acceldata{i, 9:11}, 800);
    fusegyro(filt, gyrodata{i, 9:11}, 1e-4);

    if ~mod(i, fix(5*imuFs/4))
        fusemag(filt, magdata{magidx, 9:11}, 80);
        magidx = magidx + 1;
    end

    if ~mod(i, fix(imuFs/100))
        fusegps(filt, gpsdata{GPSidx, 5:7}, 34);
        GPSidx = GPSidx + 1;
    end

    [p(i,:),q(i)] = pose(filt);

end