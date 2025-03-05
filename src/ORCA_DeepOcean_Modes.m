function ORCA_DeepOcean_Modes(nfreq_input, freq_min, freq_max, ...
    thickness_sedi, vp_sedi, vs_sedi, rho_sedi, vp_lowerhalf, ...
    vs_lowerhalf, rho_lowerhalf)
   % This script is modified from the script ORCA_MEX_example.m of the MATLAB wrapper by YT Lin (UCSD) for ORCA package.
   % Detailed introduction of ORCA is available at:
   % https://oalib-acoustics.org/models-and-software/normal-modes/
   system(['install_name_tool -add_rpath ...' ...
    '/opt/homebrew/anaconda3/envs/x86_env/lib']);
   % Add ORCA_MEX path so MATLAB can find it
   addpath('/Users/wenbowu/work/mytools/ORCA_mex_linux');
   %clear
   % Plot if enabled
   do_plot = true;

   % Sound speed profile in the upper half-space
   svp_in.uphalf_cp = 343.0;
   svp_in.uphalf_cs = 0.0;
   svp_in.uphalf_rho = 0.00121;
   svp_in.uphalf_ap = 0.0;
   svp_in.uphalf_as = 0.0;
   
   svp_in.ctol = 0;
   % Read sound speed file
   filename = 'average1D_VpRhoTS.txt';
   data = readmatrix(filename, 'NumHeaderLines', 2); % Skip first 2 lines
   % Extract first two columns, remove NaN rows
   data = data(:, 1:2);
   data = data(~isnan(data(:, 2)), :);
   % Assign water column sound speed profile
   svp_in.wssp = data;
   WaterDepth = svp_in.wssp(end, 1)
   svp_in.nsvp = size(svp_in.wssp, 1);

   % Water density
   svp_in.wrho = 1.0;
   svp_in.walphs = 0;
   
   % Bottom properties
   svp_in.nlayb = 1;
   svp_in.btm_env = [1, thickness_sedi, vp_sedi, vp_sedi, vs_sedi, vs_sedi, ...
                     rho_sedi, rho_sedi, -0.2, -0.2, -0.4, -0.4, 0, 0, 0, 0];
   %Example of two layers
   %svp_in.btm_env = [1  10 1600.0 1600.0 0 0 1.5 1.5 -.2 -.2 0 0 0 0 0 0 ];
   %                 1  30 1600.0 1600.0 0 0 1.5 1.5 -.2 -.2 0 0 0 0 0 0 ];
   % Lower half-space properties
   svp_in.lowhalf_cp = vp_lowerhalf;
   svp_in.lowhalf_cs = vs_lowerhalf;
   svp_in.lowhalf_rho = rho_lowerhalf;
   svp_in.lowhalf_ap = -0.2;
   svp_in.lowhalf_as = -0.4;
   
   svp_in.ntop = 0;
   svp_in.above_sea = [1, 50, 1650, 1650, 0, 0, 1.5, 1.5, -0.1, -0.1, 0, 0, 0, 0, 0, 0];

   % Mode calculation parameters
   opt_in.nmode = 10; % Number of modes to be computed
   opt_in.cphmax = 4000; % Phase velocity threshold,  above which are not computed.
   opt_in.rmin = 1.0;
   opt_in.rmax = 5000.0;
   opt_in.phfac = 4;
   opt_in.dbcut = 50;
   opt_in.Aih_l = -1;
   opt_in.Aih_u = -1;

   % Frequency settings
   opt_in.nf = -nfreq_input;
   opt_in.fcw_n = 2;
   %disp(['freq_min: ', num2str(freq_min), ' freq_max: ', num2str(freq_max)]);
   opt_in.fcw = [freq_min, freq_max];

   % Number of layers and depth range to be saved and plotted
   opt_in.nzm = -1250;
   opt_in.zm_n = 2;
   %opt_in.zm = [0, 1.5 * WaterDepth];
   opt_in.zm = [0,  WaterDepth];

   % Output mode functions: 0=no, 1=p-wave, 2=s-wave, 3=both p & s
   iimf = 1;

   % Compute modes using ORCA_MEX
   [nMODES, kn, freq, phi, phi_z, vg] = ORCA_MEX(svp_in, opt_in, iimf);
   vg(vg == 0) = NaN;

   % Save the eigenfunctions
   for ifreq = 1:length(freq)
       mode_save = [real(phi(:, 1, ifreq)), phi_z];
       filename = sprintf('mode_1_%.2fHz.txt', freq(ifreq));
       dlmwrite(filename, mode_save, 'delimiter', '\t');
   end

   % Plot if enabled
   if do_plot
      figure;
      hold on;
      imode_plot = 1;

      for ifreq = 1:length(freq)
          plot(real(phi(:, imode_plot, 1)), phi_z, 'LineWidth', 3);
      end

      % Plot ocean bottom line
      line([min(real(phi(:, imode_plot, :)), [], 'all'), max(real(phi(:, imode_plot, :)), [], 'all')], ...
           [WaterDepth, WaterDepth], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 2);

      % Label ocean bottom, moving text slightly below
      text(max(real(phi(:, imode_plot, :)), [], 'all') / 3, WaterDepth * 1.05, ...
          'Ocean bottom', 'FontSize', 12, 'Color', 'b');

      % Title with frequency range
      title(sprintf("Mode %d for %.2f - %.2f Hz", imode_plot, freq(1), freq(end)), ...
          'FontSize', 12, 'Color', 'b');

      ylabel('Depth (m)');
      xlabel('Mode functions');
      axis tight;
      grid on;
      set(gca, 'YDir', 'reverse');
      hold off;
   end
end
