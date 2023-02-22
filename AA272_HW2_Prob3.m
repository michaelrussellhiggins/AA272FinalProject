clc
clear all
close all

% Loads data

GPS_data_table = readtable('gnss_log.csv');
GPS_data = table2array(GPS_data_table);

% Initial guess

x0(1,:) = [0, 0, 0, 0];

% Indices

k = 2;

start_index = 1;
end_index = 2;

% Computes lestimated position using each group of time-synced GPS signals

while GPS_data(end_index, 24) ~= 1.590629989281280e+05

    % Finds indices of synced data

    while GPS_data(start_index, 24) == GPS_data(end_index, 24)
        end_index = end_index + 1;
    end

    % Extracts time-synced data

    end_index = end_index - 1;
    X = GPS_data(start_index:end_index, 28);
    Y = GPS_data(start_index:end_index, 29);
    Z = GPS_data(start_index:end_index, 30);
    B = GPS_data(start_index:end_index, 31);
    prange = GPS_data(start_index:end_index, 27);

    % Computes position estimate

    x0(k,:) = solve_pos(x0(k-1,:), X, Y, Z, B, prange);

    % Adjusts indices

    k = k + 1;
    start_index = end_index + 1;
    end_index = start_index + 1;
    
end

% Plots data

plot3(x0(2:end,1),x0(2:end,2),x0(2:end,3))
xlabel('X (m)')
ylabel('Y (m)')
zlabel('Z (m)')
title('Newton-Raphson solution for GPS data (ECEF Frame)')

% Part 3b

x0_1 = solve_pos([0, 0, 0, 0], X(1:7), Y(1:7), Z(1:7), B(1:7), prange(1:7));
x0_2 = solve_pos([0, 0, 6378000, 0], X(1:7), Y(1:7), Z(1:7), B(1:7), prange(1:7));
x0_3 = solve_pos([-2703952.62995156, -4263104.28204005, 3885071.81339034, 0], X(1:7), Y(1:7), Z(1:7), B(1:7), prange(1:7));
x0_4 = solve_pos([100000, 100000, 100000, 0], X(1:7), Y(1:7), Z(1:7), B(1:7), prange(1:7));
x0_5 = solve_pos([-100000, -100000, -100000, 0], X(1:7), Y(1:7), Z(1:7), B(1:7), prange(1:7));

error_1 = norm(x0_1(1:3) - x0_1(1:3));
error_2 = norm(x0_2(1:3) - x0_1(1:3));
error_3 = norm(x0_3(1:3) - x0_1(1:3));
error_4 = norm(x0_4(1:3) - x0_1(1:3));
error_5 = norm(x0_5(1:3) - x0_1(1:3));

function pseudo_theor = get_expected_pseudoranges(x_est, X, Y, Z, B)

    pseudo_theor = [];

    for i = 1:length(X)
        pseudo_theor(i) = norm([X(i), Y(i), Z(i)] - x_est(1:3)) + x_est(4) - B(i);
    end

end

function G = get_geometry_matrix(x_est, X, Y, Z)

    G = zeros(1, 4);

    for i = 1:length(X)
        G(i, 1:3) = -([X(i), Y(i), Z(i)] - x_est(1:3))/norm([X(i), Y(i), Z(i)] - x_est(1:3));
        G(i, 4) = 1;
    end

end

function x0 = solve_pos(x0, X, Y, Z, B, prange)

    dx0 = inf;

    while norm(dx0) > 1
        G = get_geometry_matrix(x0, X, Y, Z);
        pseudo_theor = get_expected_pseudoranges(x0, X, Y, Z, B);
        dp = prange - pseudo_theor';
        dx0 = inv(G'*G)*G'*dp;
        x0 = x0 + dx0';
    end

end