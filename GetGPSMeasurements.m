function GetGPSMeasurements(gnssfile)

data = readtable(gnssfile);

data = data(~(data.ConstellationType ~= 1),:);

writetable(data, strcat('Data/GPSOnly/', gnssfile))

end