%% Compare to ground truth
close all
clear
clc

lldata = readtable('routecheckpoints.xlsx');
%% get ground truth enu

for i = 1:3
    [E, N, U] = geodetic2enu(lldata.dragLat(i),lldata.dragLon(i),16.5,lldata.dragLat(1),lldata.dragLon(1),16.5,wgs84Ellipsoid);
    dragENU(i,:) = [E N U];
    [E, N, U] = geodetic2enu(lldata.outLat(i),lldata.outLon(i),16.5,lldata.outLat(1),lldata.outLon(1),16.5,wgs84Ellipsoid);
    outENU(i,:) = [E N U];
    [E, N, U] = geodetic2enu(lldata.curlLat(i),lldata.curlLon(i),16.5,lldata.curlLat(1),lldata.curlLon(1),16.5,wgs84Ellipsoid);
    curlENU(i,:) = [E N U];
    [E, N, U] = geodetic2enu(lldata.flyLat(i),lldata.flyLon(i),16.5,lldata.flyLat(1),lldata.flyLon(1),16.5,wgs84Ellipsoid);
    flyENU(i,:) = [E N U];
    [E, N, U] = geodetic2enu(lldata.postLat(i),lldata.postLon(i),16.5,lldata.postLat(1),lldata.postLon(1),16.5,wgs84Ellipsoid);
    postENU(i,:) = [E N U];
    [E, N, U] = geodetic2enu(lldata.cornerLat(i),lldata.cornerLon(i),16.5,lldata.cornerLat(1),lldata.cornerLon(1),16.5,wgs84Ellipsoid);
    cornerENU(i,:) = [E N U];
end

groundtruth = table(dragENU, outENU, curlENU, flyENU, postENU, cornerENU);
%save groundtruth

%% compare

imufiles = dir('Data/IMUReadingsClip/*.xlsx');
imufilenames = {imufiles.name};
idx = 35;
imufile = imufilenames{idx};

rte = 0;
loc_char = imufile(9:11);
if (strcmp(loc_char,'Dra') == 1)
    rte = 1;
elseif (strcmp(loc_char,'Out') == 1)
    rte = 2;
elseif (strcmp(loc_char,'Cur') == 1)
    rte = 3;
elseif (strcmp(loc_char,'Fly') == 1)
    rte = 4;
elseif (strcmp(loc_char,'Pos') == 1)
    rte = 5;
elseif (strcmp(loc_char,'Cor') == 1)
    rte = 6;
end

gtENU = table2array(groundtruth(:,rte));

[UTCsec, mu_pos, ub, lb] = EKFGPSTimes(imufilenames{idx}(1:end-5));

figure(2);
hold on
plot(gtENU(:,1),gtENU(:,2))
axis equal





