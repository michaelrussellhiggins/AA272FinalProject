%% Compare to ground truth
close all
clear
clc

lldata = readtable('routecheckpoints.xlsx');
%% get ground truth enu

% for i = 1:3
%     [E, N, U] = geodetic2enu(lldata.dragLat(i),lldata.dragLon(i),16.5,lldata.dragLat(1),lldata.dragLon(1),16.5,wgs84Ellipsoid);
%     dragENU(i,:) = [E N U];
%     [E, N, U] = geodetic2enu(lldata.outLat(i),lldata.outLon(i),16.5,lldata.outLat(1),lldata.outLon(1),16.5,wgs84Ellipsoid);
%     outENU(i,:) = [E N U];
%     [E, N, U] = geodetic2enu(lldata.curlLat(i),lldata.curlLon(i),16.5,lldata.curlLat(1),lldata.curlLon(1),16.5,wgs84Ellipsoid);
%     curlENU(i,:) = [E N U];
%     [E, N, U] = geodetic2enu(lldata.flyLat(i),lldata.flyLon(i),16.5,lldata.flyLat(1),lldata.flyLon(1),16.5,wgs84Ellipsoid);
%     flyENU(i,:) = [E N U];
%     [E, N, U] = geodetic2enu(lldata.postLat(i),lldata.postLon(i),16.5,lldata.postLat(1),lldata.postLon(1),16.5,wgs84Ellipsoid);
%     postENU(i,:) = [E N U];
%     [E, N, U] = geodetic2enu(lldata.cornerLat(i),lldata.cornerLon(i),16.5,lldata.cornerLat(1),lldata.cornerLon(1),16.5,wgs84Ellipsoid);
%     cornerENU(i,:) = [E N U];
% end
% 
% groundtruth = table(dragENU, outENU, curlENU, flyENU, postENU, cornerENU);
%save groundtruth
load 'ClipData/groundtruth.mat'
%% simple compare

% imufiles = dir('Data/IMUReadingsClip/*.xlsx');
% imufilenames = {imufiles.name};
% idx = 29;
% imufile = imufilenames{idx};
% 
% rte = 0;
% loc_char = imufile(9:11);
% if (strcmp(loc_char,'Dra') == 1)
%     rte = 1;
% elseif (strcmp(loc_char,'Out') == 1)
%     rte = 2;
% elseif (strcmp(loc_char,'Cur') == 1)
%     rte = 3;
% elseif (strcmp(loc_char,'Fly') == 1)
%     rte = 4;
% elseif (strcmp(loc_char,'Pos') == 1)
%     rte = 5;
% elseif (strcmp(loc_char,'Cor') == 1)
%     rte = 6;
% end
% 
% gtENU = table2array(groundtruth(:,rte));
% 
% [UTCsec, mu_pos, ub, lb] = KFGPSTimes(imufilenames{idx}(1:end-5));

% figure;
% hold on
% plot(gtENU(:,1),gtENU(:,2))
% plot(mu_pos(1,:),mu_pos(2,:))
% axis equal

% figure(1);
% hold on
% plot(UTCsec,mu_pos(1,:))

%% 6x6 plot

figure;
imufiles = dir('Data/IMUReadingsClip/*.xlsx');
imufilenames = {imufiles.name};

for i=1:6
    for j = 1:6
        subidx = j + (i-1)*6;

        rte = 0;
        if i == 1
            rte = 1;
            fileidx = 12 + j; % compensates for alphabetical ordering of files
            plotname = "Drag Trial " + string(j)
        elseif (i == 2)
            rte = 2;
            fileidx = 24 + j;
            plotname = "Out Trial " + string(j)
        elseif (i == 3)
            rte = 3;
            fileidx = 6 + j;
            plotname = "Curl Trial " + string(j)
        elseif (i == 4)
            rte = 4;
            fileidx = 18 + j;
            plotname = "Fly Trial " + string(j)
        elseif (i == 5)
            rte = 5;
            fileidx = 30 + j;
            plotname = "Post Trial " + string(j)
        elseif (i == 6)
            rte = 6;
            fileidx = 0 + j;
            plotname = "Corner Trial " + string(j)
        end

        imufile = imufilenames{fileidx};

        gtENU = table2array(groundtruth(:,rte));
        [UTCsecO, mu_posO, ~, ~] = KFGPSTimes(imufilenames{fileidx}(1:end-5));
        [UTCsecD, mu_posD, ~, ~] = KFGPSTimes(imufilenames{fileidx+36}(1:end-5));

        subplot(6,6,subidx)
        hold on
        plot(gtENU(:,1),gtENU(:,2),'-ko')
        plot(mu_posO(1,:),mu_posO(2,:),'-b')
        plot(mu_posD(1,:),mu_posD(2,:),'-r')
        axis equal
        grid on
        xlabel('E-axis (m)')
        ylabel('N-axis (m)')
        title(plotname)

    end
end

%% diff kf implementations

compidx = 20;
[UTCsec1, mu_pos1, ub1, lb1] = KFGPSTimes(imufilenames{compidx}(1:end-5));
[UTCsec2, mu_pos2, ub2, lb2] = KFAllTimes(imufilenames{compidx}(1:end-5));
[UTCsec3, mu_pos3, ub3, lb3] = EKFGPSTimes(imufilenames{compidx}(1:end-5));
[UTCsec4, mu_pos4, ub4, lb4] = EKFAllMeas(imufilenames{compidx}(1:end-5));
gtENUcomp = table2array(groundtruth(:,4));

figure;
subplot(2,2,1)
plot(gtENUcomp(:,1),gtENUcomp(:,2),'-ko')
hold on
plot(mu_pos1(1,:),mu_pos1(2,:),'-b')
axis equal
grid on
xlabel('E-axis (m)')
ylabel('N-axis (m)')
title('Kalman Filter, GPS Times Only')

subplot(2,2,2)
plot(gtENUcomp(:,1),gtENUcomp(:,2),'-ko')
hold on
plot(mu_pos2(1,:),mu_pos2(2,:),'-b')
axis equal
grid on
xlabel('E-axis (m)')
ylabel('N-axis (m)')
title('Kalman Filter, All Times')

subplot(2,2,3)
plot(gtENUcomp(:,1),gtENUcomp(:,2),'-ko')
hold on
plot(mu_pos3(1,:),mu_pos3(2,:),'-b')
axis equal
grid on
xlabel('E-axis (m)')
ylabel('N-axis (m)')
title('Extended Kalman Filter, GPS Times Only')

subplot(2,2,4)
plot(gtENUcomp(:,1),gtENUcomp(:,2),'-ko')
hold on
plot(mu_pos4(1,:),mu_pos4(2,:),'-b')
axis equal
grid on
xlabel('E-axis (m)')
ylabel('N-axis (m)')
title('Extended Kalman Filter, All Measurements')

%% calc error

pos_error = [];
for i = 1:72
    imufile = imufilenames{i};
    [UTCsec, mu_pos, ub, lb] = EKFGPSTimes(imufilenames{i}(1:end-5));

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


        t_GT = [UTCsec(1) UTCsec(floor(end/2)) UTCsec(end)];
        gtE = interp1(t_GT,gtENU(:,1),UTCsec);
        gtN = interp1(t_GT,gtENU(:,2),UTCsec);

        pos_error{i,1} = {sqrt((gtE - mu_pos(1,:)).^2 +(gtN - mu_pos(2,:)).^2)};



end

%%
% save pos_errorEKFall.mat pos_errorEKFall


%% error plots

load pos_errorEKFall.mat

mean_error = zeros(72,1);

for i=1:72
    mean_error(i) = mean(cell2mat(pos_errorEKFall{i}));
end

figure;
plot(mean_error)
xlabel('Trial #')
ylabel('Mean Position Error (m)')
title('Position Error for each Trial, EKF All Measurements')


