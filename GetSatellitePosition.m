function GetSatellitePosition(gpsfile)

% Loads data

data = readtable(strcat('Data/GPSMeasurements/', gpsfile));
dirName = 'Ephem';

numData = size(data, 1);

data(:,44) = num2cell(zeros(numData,1));
data.Properties.VariableNames{44} = 'X';

data(:,45) = num2cell(zeros(numData,1));
data.Properties.VariableNames{45} = 'Y';

data(:,46) = num2cell(zeros(numData,1));
data.Properties.VariableNames{46} = 'Z';

data(:,47) = num2cell(zeros(numData,1));
data.Properties.VariableNames{47} = 'B';

allRxMillis = int64(data{end,2} - data{end,5})*1e-6;

fctSeconds = 1e-3*double(allRxMillis);
utcTime = Gps2Utc([],fctSeconds);
allGpsEph = GetNasaHourlyEphemeris(utcTime,dirName);

for i = 1:numData

    allRxMillis = int64(data{i,2} - data{i,5})*1e-6;
    fctSeconds = 1e-3*double(allRxMillis);
    [gpsEph,iSv] = ClosestGpsEph(allGpsEph, data{i,11}, fctSeconds);
    tx_time = data{i,14}/1e9;
    sat_pos = get_sat_ECEF(gpsEph, tx_time);
    B = get_clock_error(gpsEph, tx_time);

    data(i,44) = {sat_pos(1)};
    data(i,45) = {sat_pos(2)};
    data(i,46) = {sat_pos(3)};
    data(i,47) = {B};

end

writetable(data, strcat('Data/GPSOnly/', gpsfile))


% Compute clock error

function clock_error = get_clock_error(ephem, tx_time)

    % Constants
    
    F = -4.442807633e-10;
    mu = 3.986005e14;

    % Ephemeris data extraction

    af0 = ephem.af0;
    af1 = ephem.af1;
    af2 = ephem.af2;
    deltan = ephem.Delta_n;
    M_0 = ephem.M0;
    e = ephem.e;
    sqrt_a = ephem.Asqrt;
    t_oe = ephem.Toe;
    TGD = ephem.TGD;
    t_oc = ephem.Toc;

    t = tx_time;

    a = sqrt_a^2;
    n = sqrt(mu/(a^3)) + deltan;
    t_k = t - t_oe;
    M_k = M_0 + n*t_k;
    E_k = NewtonRaphson_Eccen(M_k, e);

    deltat_r = F*e*(sqrt_a)*sin(E_k);

    deltat_sv = af0 + af1*(t - t_oc) + af2*(t - t_oc)^2 + deltat_r - TGD;

    clock_error = deltat_sv;

end



% Put one row of ephem and one time into get_sat_ECEF

function sat_pos = get_sat_ECEF(ephem, tx_time)

    % Constants

    mu = 3.986005e14;
    Omegadot_E = 7.2921151467e-5;

    % Ephemeris data extraction

    C_rs = ephem.Crs;
    deltan = ephem.Delta_n;
    M_0 = ephem.M0;
    C_uc = ephem.Cuc;
    e = ephem.e;
    C_us = ephem.Cus;
    sqrt_a = ephem.Asqrt;
    t_oe = ephem.Toe;
    C_ic = ephem.Cic;
    Omega_0 = ephem.OMEGA;
    C_is =ephem.Cis;
    i_0 = ephem.i0;
    C_rc = ephem.Crc;
    omega = ephem.omega;
    Omegadot = ephem.OMEGA_DOT;
    idot = ephem.IDOT;

    t = tx_time;

    % Step 1
    a = sqrt_a^2;

    % Step 2
    n = sqrt(mu/(a^3)) + deltan;

    % Step 3
    t_k = t - t_oe;

    % Step 4
    M_k = M_0 + n*t_k;

    % Step 5
    E_k = NewtonRaphson_Eccen(M_k, e);

    % Step 6
    sin_nu_k = (sqrt(1 - e^2)*sin(E_k))/(1 - e*cos(E_k));
    cos_nu_k = (cos(E_k) - e)/(1 - e*cos(E_k));
    nu_k = atan2(sin_nu_k, cos_nu_k);

    % Step 7
    phi_k = nu_k + omega;

    % Step 8
    delta_phi_k = C_us*sin(2*phi_k) + C_uc*cos(2*phi_k);

    % Step 9
    u_k = phi_k + delta_phi_k;

    % Step 10
    deltar_k = C_rs*sin(2*phi_k) + C_rc*cos(2*phi_k);

    % Step 11
    deltai_k = C_is*sin(2*phi_k) + C_ic*cos(2*phi_k);

    % Step 12
    Omega_k = Omega_0 - Omegadot_E*t + Omegadot*t_k;

    % Step 13
    r_k = a*(1 - e*cos(E_k)) + deltar_k;

    % Step 14
    i_k = i_0 + idot*t_k + deltai_k;

    % Step 15
    x_p = r_k*cos(u_k);

    % Step 16
    y_p = r_k*sin(u_k);

    % Step 17
    x_ECEF = x_p*cos(Omega_k) - y_p*cos(i_k)*sin(Omega_k);

    % Step 18
    y_ECEF = x_p*sin(Omega_k) + y_p*cos(i_k)*cos(Omega_k);

    % Step 19
    z_ECEF = y_p*sin(i_k);

    sat_pos = [x_ECEF, y_ECEF, z_ECEF];

end

% Calculates the eccentric anomaly using the Newton Raphson method

function E = NewtonRaphson_Eccen(M, e)

    epsilon = 1e-8;
    
    E_i = pi;
    E_iplus = E_i - (E_i - e*sin(E_i) - M)/(1 - e*cos(E_i));

    while abs(E_iplus - E_i) > epsilon

        E_i = E_iplus;

        E_iplus = E_i - (E_i - e*sin(E_i) - M)/(1 - e*cos(E_i));

    end

    E = E_i;


end

end