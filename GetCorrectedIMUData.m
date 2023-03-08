function GetCorrectedIMUData(imufile)

    % Accelerometer

    data = readtable(strcat('Data/Parsed/', imufile), 'Sheet', 'accel');
    
    numData = size(data, 1);

    data(:,9) = num2cell(zeros(numData,1));
    data.Properties.VariableNames{9} = 'CalAccelXMps2';

    data(:,10) = num2cell(zeros(numData,1));
    data.Properties.VariableNames{10} = 'CalAccelYMps2';

    data(:,11) = num2cell(zeros(numData,1));
    data.Properties.VariableNames{11} = 'CalAccelZMps2';

    for i = 1:numData

        data(i,9) = {data{i,3} - data{i,6}};
        data(i,10) = {data{i,4} - data{i,7}};
        data(i,11) = {data{i,5} - data{i,8}};

    end

    writetable(data, strcat('Data/IMUReadings/', imufile), 'Sheet', 'accel')

    % Gyroscope

    data = readtable(strcat('Data/Parsed/', imufile), 'Sheet', 'gyro');
    
    numData = size(data, 1);

    data(:,9) = num2cell(zeros(numData,1));
    data.Properties.VariableNames{9} = 'CalGyroXRadPerSec';

    data(:,10) = num2cell(zeros(numData,1));
    data.Properties.VariableNames{10} = 'CalGyroYRadPerSec';

    data(:,11) = num2cell(zeros(numData,1));
    data.Properties.VariableNames{11} = 'CalGyroZRadPerSec';

    for i = 1:numData

        data(i,9) = {data{i,3} - data{i,6}};
        data(i,10) = {data{i,4} - data{i,7}};
        data(i,11) = {data{i,5} - data{i,8}};

    end

    writetable(data, strcat('Data/IMUReadings/', imufile), 'Sheet', 'gyro')

    % Magnetometer

    data = readtable(strcat('Data/Parsed/', imufile), 'Sheet', 'mag');
    
    numData = size(data, 1);

    data(:,9) = num2cell(zeros(numData,1));
    data.Properties.VariableNames{9} = 'CalMagXMicroT';

    data(:,10) = num2cell(zeros(numData,1));
    data.Properties.VariableNames{10} = 'CalMagYMicroT';

    data(:,11) = num2cell(zeros(numData,1));
    data.Properties.VariableNames{11} = 'CalMagZMicroT';

    for i = 1:numData

        data(i,9) = {data{i,3} - data{i,6}};
        data(i,10) = {data{i,4} - data{i,7}};
        data(i,11) = {data{i,5} - data{i,8}};

    end

    writetable(data, strcat('Data/IMUReadings/', imufile), 'Sheet', 'mag')

    % Orientation

    data = readtable(strcat('Data/Parsed/', imufile), 'Sheet', 'orientation');

    writetable(data, strcat('Data/IMUReadings/', imufile), 'Sheet', 'orientation')

end