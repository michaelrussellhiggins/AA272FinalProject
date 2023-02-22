clc
close all
clear all

% Loads data

GPS_data_table = readtable('gnss_log.csv');
GPS_data = table2array(GPS_data_table);

% Loads data

ephem_data_table = readtable('ephem.csv', 'Range', 'A:AD');
ephem_data = table2array(ephem_data_table);

% Problem 1b

% GPS position for for time instane

sat_pos = [];
positions = [];

for i = 1:7
    
    Svid = GPS_data(i, 6);
    tx_time = GPS_data(i, 9)/1e9;

    sat_pos(i,:) = get_sat_ECEF(ephem_data(Svid,:), tx_time);

    positions(i,:) = [Svid, sat_pos(i,:)];

end


% Problem 1c

% GPS position for 15 time instances

sat_pos = [];
error = [];

for i = 1:238
    
    Svid = GPS_data(i, 6);
    tx_time = GPS_data(i, 9)/1e9;

    sat_pos(i,:) = get_sat_ECEF(ephem_data(Svid,:), tx_time);

    error(i,:) = [sat_pos(i, 1) - GPS_data(i, 28), sat_pos(i, 2) - GPS_data(i, 29), sat_pos(i, 3) - GPS_data(i, 30)];

end

figure(1)

plot(1:238, error(:,1))
hold on
plot(1:238, error(:,2))
hold on
plot(1:238, error(:,3))
title('Error between compute ECEF position and GNSS log position')
xlabel('Measurement number')
ylabel('Error [m]')
legend('X [m]', 'Y [m]', 'Z [m]')
grid on


% Problem 1d

% GPS position for first time instance with clock bias

clock_baises = [];

for i = 1:7
    
    Svid = GPS_data(i, 6);
    tx_time = GPS_data(i, 9)/1e9;
    
    clock_error = get_clock_error(ephem_data(Svid,:), tx_time);

    c = 2.99792458e8;

    clock_baises(i,:) = [Svid, clock_error*c];

end


% Problem 1e

clock_errors = [];

for i = 1:238
    
    Svid = GPS_data(i, 6);
    tx_time = GPS_data(i, 9)/1e9;
    
    clock_error = get_clock_error(ephem_data(Svid,:), tx_time);

    c = 2.99792458e8;

    clock_errors(i) = [clock_error*c - GPS_data(i, 31)];

end

figure(2)

plot(1:238, clock_errors)
title('Error between computed clock bias and GNSS clock bias')
xlabel('Measurement number')
ylabel('Error [m]')
grid on


% Bonus 2

% GPS position for 15 time instances

sat_pos = [];
clock_corrected_error = [];

for i = 1:238
    
    Svid = GPS_data(i, 6);
    tx_time = GPS_data(i, 9)/1e9;

    clock_error = get_clock_error(ephem_data(Svid,:), tx_time);

    sat_pos(i,:) = get_sat_ECEF(ephem_data(Svid,:), tx_time - clock_error);

    clock_corrected_error(i,:) = [sat_pos(i, 1) - GPS_data(i, 28), sat_pos(i, 2) - GPS_data(i, 29), sat_pos(i, 3) - GPS_data(i, 30)];

end

figure(3)

plot(1:238, clock_corrected_error(:,1))
hold on
plot(1:238, clock_corrected_error(:,2))
hold on
plot(1:238, clock_corrected_error(:,3))
title('Error between compute ECEF position and GNSS log position with clock error corrected')
xlabel('Measurement number')
ylabel('Error [m]')
legend('X [m]', 'Y [m]', 'Z [m]')
grid on


% Compute clock error

function clock_error = get_clock_error(ephem, tx_time)

    % Constants
    
    F = -4.442807633e-10;
    mu = 3.986005e14;

    % Ephemeris data extraction

    af0 = ephem(2);
    af1 = ephem(3);
    af2 = ephem(4);
    deltan = ephem(7);
    M_0 = ephem(8);
    e = ephem(10);
    sqrt_a = ephem(12);
    t_oe = ephem(13);
    TGD = ephem(27);
    t_oc = ephem(29);

    t = tx_time;

    a = sqrt_a^2;
    n = sqrt(mu/(a^3)) + deltan;
    t_k = t - t_oe;
    M_k = M_0 + n*t_k;
    E_k = NewtonRaphson_Eccen(M_k, e);

    deltat_r = F*e^(sqrt_a)*sin(E_k);

    deltat_sv = af0 + af1*(t - t_oc) + af2*(t - t_oc)^2 + deltat_r - TGD;

    clock_error = deltat_sv;

end



% Put one row of ephem and one time into get_sat_ECEF

function sat_pos = get_sat_ECEF(ephem, tx_time)

    % Constants

    mu = 3.986005e14;
    Omegadot_E = 7.2921151467e-5;

    % Ephemeris data extraction

%     Sv_clock_bias = ephem(2);
%     Sv_clock_drift = ephem(3);
%     Sv_clock_drift_rate = ephem(4);
%     IODE = ephem(5);
    C_rs = ephem(6);
    deltan = ephem(7);
    M_0 = ephem(8);
    C_uc = ephem(9);
    e = ephem(10);
    C_us = ephem(11);
    sqrt_a = ephem(12);
    t_oe = ephem(13);
    C_ic = ephem(14);
    Omega_0 = ephem(15);
    C_is = ephem(16);
    i_0 = ephem(17);
    C_rc = ephem(18);
    omega = ephem(19);
    Omegadot = ephem(20);
    idot = ephem(21);
%     codesL2 = ephem(22);
%     GPSWeek = ephem(23);
%     L2Pflag = ephem(24);
%     SVacc = ephem(25);
%     health = ephem(26);
%     TGD = ephem(27);
%     iodc = ephem(28);
%     trans_time = ephem(29);
%     fit_intvl = ephem(30);
%     Svid = ephem(31);

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

%     while delta_phi_k > 1e-8
% 
%         delta_phi_k = C_us*sin(2*u_k) + C_uc*cos(2*u_k);
%         u_k = phi_k + delta_phi_k;
% 
%     end

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