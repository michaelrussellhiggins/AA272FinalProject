function [r_enu] = ecefTOenu(r_ecef,loc)
%% finds position vector in ENU coordinates using geodetic Earth
% Find ENU coordinates (km) from ECEF coordinates and station position
% r_ecef = ECEF postion vector (km)
% lat_station = ground station latitude (deg)
% lon_station = ground station longitude (deg)
    
    switch loc
        case 1
            lat_station = 37.430780;
            lon_station = -122.155026;
        case 2
            lat_station = 37.431294;
            lon_station = -122.155059;
        case 3
            lat_station = 37.431108;
            lon_station = -122.155158;
    end

    r_E = 6.3781e6; % m
    e = .0818; % Earth eccentricity
    den1 = 1 - (e*sind(lat_station))^2;
    N = r_E/sqrt(den1); % modified radius of curvature

%     r_station = N * [cosd(lat_station)*cosd(lon_station);
%                        cosd(lat_station)*sind(lon_station);
%                        sind(lat_station)]; % km

    r_station = lla2ecef([lat_station lon_station 16.5])';

    eE = [-sind(lon_station); 
           cosd(lon_station);
           0];
    eN = [-sind(lat_station)*cosd(lon_station);
          -sind(lat_station)*sind(lon_station);
          cosd(lat_station)];
    eU = [cosd(lat_station)*cosd(lon_station);
          cosd(lat_station)*sind(lon_station);
          sind(lat_station)];
      
    rotENU = [eE eN eU]';

    r_enu = rotENU*[r_ecef - r_station];
end