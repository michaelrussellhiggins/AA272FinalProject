%% Import GPS and IMU Positions from csv files

close all
clear
clc

%%
gpsposfiles = dir('Data/GPSPosition2/*.csv');

gpsfilenames = {gpsposfiles.name};
fileidx = 40;
gpsfile = gpsfilenames{fileidx};

gpsdata = readtable(['Data/GPSPosition2/' gpsfile]);

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