%% Import GPS and IMU Positions from csv files

close all
clear
clc

%%
gpsposfiles = dir('Data/GPSPositionDirect/*.csv');

gpsfilenames = {gpsposfiles.name};
fileidx = 25;
gpsfile = gpsfilenames{fileidx};

gpsdata = readtable(['Data/GPSPositionDirect/' gpsfile]);

%%

imuposfiles = dir('Data/IMUReadings/*.xlsx');

imufilenames = {imuposfiles.name};
imufile = imufilenames{fileidx};

acceldata = readtable(['Data/IMUReadings/' imufile],'Sheet','accel');
gyrodata = readtable(['Data/IMUReadings/' imufile],'Sheet','gyro');
magdata = readtable(['Data/IMUReadings/' imufile],'Sheet','mag');
orientdata = readtable(['Data/IMUReadings/' imufile],'Sheet','orientation');

timeA = (acceldata.utcTimeMillis-acceldata.utcTimeMillis(1))/1000;
accX = acceldata.CalAccelXMps2;
accY = acceldata.CalAccelYMps2;
accZ = acceldata.CalAccelZMps2;

timeO = (orientdata.utcTimeMillis-orientdata.utcTimeMillis(1))/1000;
yaw = orientdata.yawDeg;
roll = orientdata.rollDeg;
pitch = orientdata.pitchDeg;
%%

r_ENU = [gpsdata.E gpsdata.N gpsdata.U];

loc = 0;
loc_char = gpsfile(9:11);
if (strcmp(loc_char,'Out') == 1) || (strcmp(loc_char,'Dra') == 1)
    loc = 1;
elseif (strcmp(loc_char,'Cur') == 1) || (strcmp(loc_char,'Fly') == 1)
    loc = 2;
elseif (strcmp(loc_char,'Cor') == 1) || (strcmp(loc_char,'Pos') == 1)
    loc = 3;
end

%%
figure;
plot(r_ENU(:,1),r_ENU(:,2),'-o')
axis equal
grid on
xlabel('E-axis(m)')
ylabel('N-axis(m)')

%%
%M = movmean(accX,100);

figure;
plot(accX)
hold on
%plot(M)
plot(accY)
plot(accZ)
grid on
ylabel('Acceleration (m/s^2)')
xlabel('index')
%%


figure;
plot(yaw)
hold on
plot(roll)
plot(pitch)
grid on
ylabel('Orientation (deg)')
xlabel('index')

%%

figure;
yyaxis right
plot(timeO,yaw)
hold on
plot(timeO,roll)
plot(timeO,pitch)
ylabel('Orientation (deg)')
yyaxis left
plot(timeA,accX)
plot(timeA,accY)
plot(timeA,accZ)
grid on
ylabel('Acceleration (m/s^2)')
xlabel('time from log start (s)')
legend('X acc','Y acc','Z acc','yaw','roll','pitch')

%%
load 'ClipData/accelstartstop.mat'
accelcliptimes = zeros(72,2);
for i = 1:length(gpsposfiles)
    imufile = imufilenames{i};
    acceldata = readtable(['Data/IMUReadings/' imufile],'Sheet','accel');
    accelcliptimes(i,1) = acceldata.utcTimeMillis(accelstartstop(i,1));
    accelcliptimes(i,2) = acceldata.utcTimeMillis(accelstartstop(i,2));
end
%% orient clip

orientstartstop = zeros(72,2);
orientcliptimes = zeros(72,2);

for i = 1:72
    imufile = imufilenames{i};
    orientdata = readtable(['Data/IMUReadings/' imufile],'Sheet','orientation');
    orientutc = orientdata.utcTimeMillis;

    minOstop = 1e10;
    minOstart = 1e10;

    for j = 1:length(orientutc)
        odiffstart = abs(orientutc(j) - accelcliptimes(i,1));
        odiffstop = abs(orientutc(j) - accelcliptimes(i,2));
        if odiffstart < minOstart
            orientcliptimes(i,1) = orientutc(j);
            orientstartstop(i,1) = j;
            minOstart = odiffstart;
        end
        if odiffstop < minOstop
            orientcliptimes(i,2) = orientutc(j);
            orientstartstop(i,2) = j;
            minOstop = odiffstop;
        end
    end
end
%% gyro clip

gyrostartstop = zeros(72,2);
gyrocliptimes = zeros(72,2);

for i = 1:72
    imufile = imufilenames{i};
    gyrodata = readtable(['Data/IMUReadings/' imufile],'Sheet','gyro');
    gyroutc = gyrodata.utcTimeMillis;

    minYstop = 1e10;
    minYstart = 1e10;

    for j = 1:length(gyroutc)
        ydiffstart = abs(gyroutc(j) - accelcliptimes(i,1));
        ydiffstop = abs(gyroutc(j) - accelcliptimes(i,2));
        if ydiffstart < minYstart
            gyrocliptimes(i,1) = gyroutc(j);
            gyrostartstop(i,1) = j;
            minYstart = ydiffstart;
        end
        if ydiffstop < minYstop
            gyrocliptimes(i,2) = gyroutc(j);
            gyrostartstop(i,2) = j;
            minYstop = ydiffstop;
        end
    end
end
%% gps clip

gpsstartstop = zeros(72,2);
gpscliptimes = zeros(72,2);

for i = 1:72
    gpsfile = gpsfilenames{i};
    gpsdata = readtable(['Data/GPSPositionDirect/' gpsfile]);  
    gpsutc = gpsdata.utc;

    minGstop = 1e10;
    minGstart = 1e10;
    
    for j = 1:length(gpsutc)
        gdiffstart = abs(gpsutc(j) - accelcliptimes(i,1));
        gdiffstop = abs(gpsutc(j) - accelcliptimes(i,2));
        if gdiffstart < minGstart
            gpscliptimes(i,1) = gpsutc(j);
            gpsstartstop(i,1) = j;
            minGstart = gdiffstart;
        end
        if gdiffstop < minGstop
            gpscliptimes(i,2) = gpsutc(j);
            gpsstartstop(i,2) = j;
            minGstop = gdiffstop;
        end
    end
end

%%
% figure;
% plot(orientdata.utcTimeMillis)
% hold on
% plot(accelcliptimes(idx,1)*ones(size(orientdata,1),1),'g')
% plot(accelcliptimes(idx,2)*ones(size(orientdata,1),1),'r')
% 
% figure;
% plot(gpsdata.utc)
% hold on
% plot(accelcliptimes(idx,1)*ones(size(gpsdata,1),1),'g')
% plot(accelcliptimes(idx,2)*ones(size(gpsdata,1),1),'r')




