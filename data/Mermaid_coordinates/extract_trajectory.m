% Read the filenames from the external file
fid = fopen('filenames.txt', 'r');
rawFilenames = textscan(fid, '%s');
fclose(fid);

% Extract the cell array from the result
rawFilenames = rawFilenames{1};

% Loop through the filenames
for idx = 1:length(rawFilenames)
    % Remove newline characters
    geoCsvFilename = strtrim(rawFilenames{idx});

    disp(geoCsvFilename);
    % Read GeoCSV file
    G = readGeoCSV(geoCsvFilename);

    % Combine vectors into a table
    data = table(G.StartTime', G.Latitude', G.Longitude', G.WaterPressure', 'VariableNames', {'Datetime', 'Latitude', 'Longitude', 'WaterPressure'});

    % Find rows with NaN values in Latitude, Longitude or WaterPressure
    nanRows = any(ismissing(data(:, {'Latitude', 'Longitude','WaterPressure'})), 2);

    % Remove rows with NaN values
    data = data(~nanRows, :);

    % Construct the output file path
    [~, filename, ~] = fileparts(geoCsvFilename);
    %outputPath = [filename, '.txt'];

    prefix = regexp(filename, 'P\d+', 'match');
    outputPath = [prefix{1}, '_traj.txt']; % Concatenate the prefix with '_traj.txt'

    % Write the table to a .txt file
    writetable(data, outputPath, 'Delimiter', '\t');
end


