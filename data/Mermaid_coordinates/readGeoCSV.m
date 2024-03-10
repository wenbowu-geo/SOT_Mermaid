function [G, C] = readGeoCSV(filename)
% [G, C] = READGEOCSV(filename)
%
% GeoCSV reader to parse files formatted (roughly) per IRIS/EarthCube GeoCSV v2.0.0.4:
% http://geows.ds.iris.edu/documents/GeoCSV.pdf
%
% Input:
% filename     GeoCSV file name
%
% Output:
% G            Structure organizing data, with fieldnames generated from header
% C            Structure organizing comments
%
% READGEOCSV requires these three "#keyword: value" comments, which must come
% before a single uncommented header line:
%    "#delimiter: <delim>"
%    "#field_type: <ftype>"
%    "#field_unit: <funit>"
%
% Any combination of comment lines and uncommented data lines may then follow,
% but comments must always follow "#keyword: value" format.
%
% NB, the requirement of specific known keyword-comments differs slightly from
% the v2.0.4 specification above, specifically the "Minimal IRIS Station example."
%
% Datetime field types must be IS08601 format ending in "Z" with millisecond
% precision: 'uuuu-MM-ddTHH:mm:ss.SSSZ'.
%
% Ex: (FYI - https://github.com/joelsimon/omnia/blob/master/plotbits/longitude360.m)
%    G = READGEOCSV('P0006.GeoCSV')
%    figure; scatter(G.StartTime, -G.WaterPressure/100)
%    title('dives'); xlabel('date'); ylabel('depth (m)')
%    figure; scatter(G.Latitude, longitude360(G.Longitude))
%
% Author: Joel D. Simon
% Contact: jdsimon@alumni.princeton.edu | joeldsimon@gmail.com
% Last modified: 27-Apr-2023, Version 9.3.0.948333 (R2017b) Update 9 on MACI64
% Last tested: GeoCSV v2.2.0-1, automaid v3.6.0-I

% Future:
% *Recognize/identify/act on other known keyword comments?

%% ___________________________________________________________________________ %%
%% Read GeoCSV file and determine which lines are COMMENTS, HEADER, and DATA
%% ___________________________________________________________________________ %%

% Open GeoCSV file.
fid = fopen(filename, 'r');
if fid == -1
    error('Cannot read %s. Check path and permissions.', filename)

end

% Read file one line at a time.
hdr = {};
data = {};
fline = 0;
dline = 0;
while true
    l = fgetl(fid);
    fline = fline  + 1;
    if ~ischar(l)
        % End of file.
        break

    end

    % Strip off any quotes that protect non-delimiting delimiters (e.g.,
    % "#delimiter: ','"). If a line starts with a double quote it must also end
    % with a double quote.
    if strcmp(l(1), '"')
        if ~strcmp(l(end), '"')
            error(['Line beginning with double quote ''"'', e.g. to protect ' ...
                   '"#delimiter: ," keyword-comment, must also end with double quote'])

        end
        % If the line starts with a double quote it also ends
        l(1) = [];
        l(end) = [];

    end

    %% COMMENT line(s) start with "#" and come before HEADER line.
    if strcmp(l(1), '#')

        % Strip hash symbol.
        l(1) = [];

        % Parse comment line to check for known keywords, which are formatted as
        % "#keyword: value"
        spl = regexp(l, ':', 'split', 'once');
        if length(spl) ~= 2
            error(['Bad comment: Comment lines must be formatted as ' ...
                  '"#keyword: value"\nCurrently (L%i): "#%s"'], fline, l)

        end
        key = lower(strtrim(spl{1}));
        val = strtrim(spl{2});

        % Strip any quotes around special (necessarily protected) keyword-value pairs.
        if contains(key, {'delimiter', 'lineterminator'})
            val = strrep(val, '"', '');
            val = strrep(val, '''', '');

        end

        % Organize comment structure.
        C.(key) = val;

    %% HEADER line is first uncommented line after any COMMENT line(s)
    elseif isempty(hdr)
        hdr = l;

        % Check for three required keyword comments, which fully describe how
        % to organize the output.
        if ~isfield(C, 'delimiter')
            error('Did not find required "#delimiter: <char>" comment line')

        end
        hdr = strsplit(hdr, C.delimiter);

        % Field type, e.g., string or datetime (for computers)
        if ~isfield(C, 'field_type')
            error('Did not find required "#field_type: <char>" comment line')

        end
        ftype = strsplit(C.field_type, C.delimiter);

        % Field unit, e.g., iso8601 or degreees_north (for humans)
        if ~isfield(C, 'field_unit')
            error('Did not find required "#field_unit: <char>" comment line')

        end
        funit = strsplit(C.field_unit, C.delimiter);

        % Check metadata and data lines are self consistent.
        if length(hdr) ~= length(ftype) || length(ftype) ~= length(funit)
            error(['Number of columns differ between field_unit, field_type, '...
                  'and header lines'])

        end

        % Initialize output by dynamically generating fieldnames from the header.
        for i = 1:length(hdr)
            G.(hdr{i}) = {};
            G.(sprintf('%sUnit', hdr{i})) = funit{i};

        end

    %% DATA lines are all uncommented lines after HEADER
    else
        % Split the data line along the delimiter.
        data = strsplit(l, C.delimiter);
        dline = dline + 1;

        % The number of data columns matches the header, which where used to
        % initialze the output structre.  Loop over the length of the split
        % data line and place each in its respective output strcutre fieldname.
        for i = 1:length(hdr)
            G.(hdr{i}){dline} = data{i};

        end
    end
end
fclose(fid);

%% ___________________________________________________________________________ %%
%% Cast strings to format descibed in the field type comment
%% ___________________________________________________________________________ %%

% The ouput struct is complete at this point and the GeoCSV file is closed, but
% the fields are all cell arrays of string; convert e.g., latitudes, to double.

% Datetime in IS0 8601 format with millisecond precision.
dtime_fmt = 'uuuu-MM-dd''T''HH:mm:ss.SSS''Z''';

% Loop over each fieldname (header/column of GeocSV) and cast string to
% specified format.
for i = 1:length(hdr)
    switch lower(ftype{i})
      case {'string', 'unitless'}
        % Do nothing

      case 'datetime'
        G.(hdr{i}) = datetime(G.(hdr{i}), 'Format', dtime_fmt, 'TimeZone', 'UTC');

      case {'float' 'integer'}
        % Float and integer are cast as double, but likely printed as single.
        G.(hdr{i}) = str2double(G.(hdr{i}));

      otherwise
        error('Unrecognized field type: %s', ftype{i})

    end
end
