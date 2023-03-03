function GetUserPosition(gpsfile)

% Loads data

data = readtable(strcat('Data/GPSOnly/', gpsfile));

data = data(~(data.ConstellationType ~= 1),:);

numData = size(data, 1);

% Initial guess

utc(1) = [0];
x0(1,:) = [0, 0, 0, 0];

% Indices

k = 2;

start_index = 1;
end_index = 2;

% Computes estimated position using each group of time-synced GPS signals

while start_index <= numData

    % Finds indices of synced data

    while (data{start_index, 41} == data{end_index, 41})
        end_index = end_index + 1;
        if end_index == numData + 1
            break
        end
    end

    % Extracts time-synced data

    end_index = end_index - 1;
    X = data{start_index:end_index, 44};
    Y = data{start_index:end_index, 45};
    Z = data{start_index:end_index, 46};
    B = data{start_index:end_index, 47};
    prange = data{start_index:end_index, 43};

    % Computes position estimate

    utc(k) = data{start_index, 41};
    x0(k,:) = solve_pos(x0(k-1,:), X, Y, Z, B, prange);

    % Adjusts indices

    k = k + 1;
    start_index = end_index + 1;
    end_index = start_index + 1;
    
end

utc = utc(2:end)';
X = x0(2:end,1);
Y = x0(2:end,2);
Z = x0(2:end,3);

PositionData = table(utc, X, Y, Z);

writetable(PositionData, strcat('Data/GPSPositionData/', gpsfile))

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

end