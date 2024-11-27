function KRAKEN_mode_debug(frequency)
    % KRAKEN_mode_debug.m
    % frequency: the frequency to use in the script (passed as argument)

    % Replace the frequency value in the configuration file with the command-line input
    system(sprintf('sed -e "2s/.*/%f #frequency/" < Twave_example.env > Twave.env', frequency));

    % Run the Krakenc function
    krakenc('Twave');

    % Clear any existing variable or data
    clear read_modes_bin;

    % Read mode data from the updated 'Twave.mod'
    [ Modes ] = read_modes('Twave.mod', [1 2]);

    % Plot the mode (optional, uncomment if needed)
    % figure;
    % plotmode('Twave.mod', [1 2]);

    % Save the first mode data to a file
    mode_save = [Modes.phi(:,1) Modes.z];
    dlmwrite(sprintf('mode_1_%.2fHz.txt', frequency), mode_save, 'delimiter', '\t');

    % Optionally, save other modes as needed (uncomment to save the second mode)
    % mode_save = [Modes.phi(:,2) Modes.z];
    % dlmwrite(sprintf('mode_2_%.2fHz.txt', frequency), mode_save, 'delimiter', '\t');
end
