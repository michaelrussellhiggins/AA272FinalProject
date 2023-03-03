function GetPseudoranges(gpsfile)

data = readtable(strcat('Data/GPSOnly/', gpsfile));

% Add a new columns

numData = size(data, 1);

data(:,37) = num2cell(zeros(numData,1));
data.Properties.VariableNames{37} = 'WeekNumber';

data(:,38) = num2cell(zeros(numData,1));
data.Properties.VariableNames{38} = 'CorrectedRxTime';

data(:,39) = num2cell(zeros(numData,1));
data.Properties.VariableNames{39} = 'AdjustedBiasRxClock';

data(:,40) = num2cell(zeros(numData,1));
data.Properties.VariableNames{40} = 'RxTimeGPSFrame';

data(:,41) = num2cell(zeros(numData,1));
data.Properties.VariableNames{41} = 'RxTimeWeekFrame';

data(:,42) = num2cell(zeros(numData,1));
data.Properties.VariableNames{42} = 'Pseudorangens';

data(:,43) = num2cell(zeros(numData,1));
data.Properties.VariableNames{43} = 'Pseudorangem';

% Compute pseudorange

secondsinweek = 604800;
speedoflight = 299792458;

for i = 1:numData

    FullBiasNanos = data{i,5};
    TimeNanos = data{i,2};
    TimeOffsetNanos = data{i,12};
    BiasNanos = data{i,6};
    ReceivedSvTimeNanos = data{i,14};

    WeekNumber = floor(-1*FullBiasNanos/(secondsinweek*(10^9)));
    CorrectedRxTime = TimeNanos + TimeOffsetNanos;
    AdjustedBiasRxClock = FullBiasNanos+BiasNanos;
    RxTimeGPSFrame = CorrectedRxTime - AdjustedBiasRxClock;
    RxTimeWeekFrame = RxTimeGPSFrame - (WeekNumber*secondsinweek*(10^9));
    Pseudorangens = RxTimeWeekFrame - ReceivedSvTimeNanos;
    Pseudorangem = Pseudorangens*speedoflight/(10^9);

    data(i,37) = {WeekNumber};
    data(i,38) = {CorrectedRxTime};
    data(i,39) = {AdjustedBiasRxClock};
    data(i,40) = {RxTimeGPSFrame};
    data(i,41) = {RxTimeWeekFrame};
    data(i,42) = {Pseudorangens};
    data(i,43) = {Pseudorangem};

end

writetable(data, strcat('Data/GPSOnly/', gpsfile))

end