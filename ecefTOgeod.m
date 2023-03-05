function [lat_gd, lon_gd, h_gd] = ecefTOgeod(r_ecef)
%% finds position vector in geodetic coordinates ( phi', lambda', h')
% r_ecef = radius in geocentric coordinates (km)
% lon_gc = geocentric longitude
% lat_gc = geocentric latitude
% lon_gd = geocentric longitude
% lat_gd = geocentric latitude
% all angles in deg
    
    x = r_ecef(1);
    y = r_ecef(2);
    z = r_ecef(3);
    
    r_E = 6378e3; % km
    r_geoc = sqrt(x^2 + y^2 + z^2);
    h_gc = r_geoc - r_E;
    lon_gc = atan2d(y,x);
    lat_gc = asind(z/r_geoc);
    

    lon_gd = lon_gc;
    e = .0818; % Earth eccentricity
    
    xy = (r_E + h_gc)*cosd(lat_gc);
    z = (r_E + h_gc)*sind(lat_gc);
    tol = 1e-8;
    lat0 = lat_gc; % first guess for lat_gd 
    i = 1;
    while i < 20
        
        den1 = 1 - (e*sind(lat0))^2;
        N = r_E/sqrt(den1); % modified radius of curvature

        num1 = z + N*(e^2)*sind(lat0);
        lat1 = atan2d(num1,xy);
        
        i = i + 1;
        
        if abs(lat1 - lat0) < tol
            break
        else
           lat0 = lat1;
        end
    end
    
    if abs(lat1 - lat0) < tol
        lat_gd = lat1;
        h_gd = (xy/cosd(lat_gd)) - N;
    else
        lat_gd = -1000;
        h_gd = -1000;
    end

end