%% Compare to ground truth
close all
clear
clc
bad = [15, 25, 27, 29, 36, 39, 57, 64];

lldata = readtable('routecheckpoints.xlsx');
%% get ground truth enu
% 
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
% save groundtruth.mat groundtruth
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

        UTCsecD = zeros(2,1);
        mu_posD = zeros(2,1);
        UTCsecO = zeros(2,1);
        mu_posO = zeros(2,1);

        gtENU = table2array(groundtruth(:,rte));
        if sum(fileidx == bad) == 0
            [UTCsecO, mu_posO, ~, ~] = KFAllMeas(imufilenames{fileidx}(1:end-5));
        end
        if sum((fileidx+36) == bad) == 0
            [UTCsecD, mu_posD, ~, ~] = KFAllMeas(imufilenames{fileidx+36}(1:end-5));
        end

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
imufilenames{compidx}
[UTCsec1, mu_pos1, ub1, lb1] = KFGPSTimes(imufilenames{compidx}(1:end-5));
[UTCsec2, mu_pos2, ub2, lb2] = KFAllMeas(imufilenames{compidx}(1:end-5));
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
title('Kalman Filter, All Measurements')

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

pos_error = cell(72,1);

for i = 1:72
    imufile = imufilenames{i};
    
    UTCsec = 0;
    mu_pos = 0;
    ub = 0;
    lb = 0;

     if sum(i == bad) == 0
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
end

%%
pos_errorEKFgps = pos_error;
save pos_errorEKFgps.mat pos_errorEKFgps


%% error plots

load pos_errorEKFall.mat
load pos_errorEKFgps.mat
load pos_errorKFall.mat
load pos_errorKFgps.mat

mean_error = zeros(72,1);

for i=1:72
    mean_error1(i) = mean(cell2mat(pos_errorKFgps{i}));
    mean_error2(i) = mean(cell2mat(pos_errorKFall{i}));
    mean_error3(i) = mean(cell2mat(pos_errorEKFgps{i}));
    mean_error4(i) = mean(cell2mat(pos_errorEKFall{i}));
end

figure;
plot(mean_error1,'-o')
hold on
plot(mean_error2,'-o')
plot(mean_error3,'-o')
plot(mean_error4,'-o')
grid on
xlabel('Trial #')
ylabel('Mean Position Error (m)')
title('Position Error for each Trial')
legend('KF GPS Only','KF All Measurements','EKF GPS Only','EKF All Measurements')

%% means of means

av_mean_error1 = mean(mean_error1,'omitnan')
av_mean_error2 = mean(mean_error2,'omitnan')
av_mean_error3 = mean(mean_error3,'omitnan')
av_mean_error4 = mean(mean_error4,'omitnan')

save('finalerror.mat',"av_mean_error1","av_mean_error2","av_mean_error3","av_mean_error4")

