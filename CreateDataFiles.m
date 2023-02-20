function CreateDataFiles(filename)

textfilename = strcat(filename,'.txt');
GPSfilename = strcat(filename, 'GPS.csv');
accelfilename = strcat(filename, 'Accel.csv');
magfilename = strcat(filename, 'Mag.csv');
gyrofilename = strcat(filename, 'Gyro.csv');
orientationfilename = strcat(filename, 'Orientation.csv');
combinedfilename = strcat(filename, '.xlsx');

textfile = fopen(textfilename,'r');
GPSfile = fopen(GPSfilename,'w');
accelfile = fopen(accelfilename,'w');
magfile = fopen(magfilename,'w');
gyrofile = fopen(gyrofilename,'w');
orientationfile = fopen(orientationfilename,'w');

while true
    thisline = fgetl(textfile);
    if ~ischar(thisline)
        break; 
    elseif isempty(strfind(thisline,'Raw,')) & isempty(strfind(thisline,'UncalAccel,')) & isempty(strfind(thisline,'UncalMag,')) & isempty(strfind(thisline,'UncalGyro,')) & isempty(strfind(thisline,'OrientationDeg,'))
        continue;
    elseif ~isempty(strfind(thisline,'UncalAccel,'))
        thisline = strrep(thisline,'UncalAccel,',''); 
        thisline = strrep(thisline,'#','');
        thisline = strrep(thisline,' ','');
        fprintf(accelfile,'%s\n',thisline);
    elseif ~isempty(strfind(thisline,'UncalMag,'))
        thisline = strrep(thisline,'UncalMag,',''); 
        thisline = strrep(thisline,'#','');
        thisline = strrep(thisline,' ','');
        fprintf(magfile,'%s\n',thisline);
    elseif ~isempty(strfind(thisline,'UncalGyro,'))
        thisline = strrep(thisline,'UncalGyro,',''); 
        thisline = strrep(thisline,'#','');
        thisline = strrep(thisline,' ','');
        fprintf(gyrofile,'%s\n',thisline);
    elseif ~isempty(strfind(thisline,'OrientationDeg,'))
        thisline = strrep(thisline,'OrientationDeg,',''); 
        thisline = strrep(thisline,'#','');
        thisline = strrep(thisline,' ','');
        fprintf(orientationfile,'%s\n',thisline);
    elseif ~isempty(strfind(thisline,'Raw,'))
        thisline = strrep(thisline,'Raw,',''); 
        thisline = strrep(thisline,'#','');
        thisline = strrep(thisline,' ','');
        fprintf(GPSfile,'%s\n',thisline);
    else
        error('Unknown data in data file')
    end
end

data = readcell(accelfilename);
writecell(data,combinedfilename,'Sheet','accel');
data = readcell(gyrofilename);
writecell(data,combinedfilename,'Sheet','gyro');
data = readcell(magfilename);
writecell(data,combinedfilename,'Sheet','mag');
data = readcell(orientationfilename);
writecell(data,combinedfilename,'Sheet','orientation');

fclose(textfile);
fclose(accelfile);
fclose(magfile);
fclose(gyrofile);
fclose(orientationfile);
fclose(GPSfile);

delete(accelfilename);
delete(magfilename);
delete(gyrofilename);
delete(orientationfilename);

end