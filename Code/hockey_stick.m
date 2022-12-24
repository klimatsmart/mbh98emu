% This script performs the regression analysis on the proxy and instrumental data,
% and calculates cross-validation statistics for the temperature reconstruction.

% Reconstruction steps.
start_step = [1400, 1450, 1500, 1600, 1650, 1700, 1730, 1750, 1760, 1780, 1800, 1820];
end_step = [start_step(2 : end) - 1, 1980];
length_step = end_step - start_step + 1;
length_recon = end_step(end) - start_step(1) + 1;

% Instrumental period.
start_instr = 1854;
end_instr = 1993;
year_count = end_instr - start_instr + 1;
month_count = year_count * 12;

% Calibration period.
start_cal = 1902;
end_cal = 1980;
length_cal = end_cal - start_cal + 1;

% Verification period.
start_ver = 1854;
end_ver = 1901;
length_ver = end_ver - start_ver + 1;

% Useful time indices.
idx_cal_monthly = 12 * (start_cal - start_instr) + 1 : 12 * (end_cal - start_instr + 1);       % 1902–1980
idx_cal_end_monthly = 12 * (start_cal - start_instr) + 1 : 12 * (end_instr - start_instr + 1); % 1902–1993
idx_cal = start_cal - start_instr + 1 : end_cal - start_instr + 1;                             % 1902–1980
idx_ver = start_ver - start_instr + 1 : end_ver - start_instr + 1;                             % 1854–1901

% Load dense temperature data.
idx_dense = importdata("../Data/Instrumental/idx_dense.txt");
[fid, msg] = fopen("../Data/Instrumental/anomalies_dense.bin", "r");
if fid == -1
  error(msg);
endif
T_dense_monthly = fread(fid, [month_count, numel(idx_dense)], "int16");
T_dense_monthly = single(T_dense_monthly);
fclose(fid);

% Load sparse temperature data.
idx_sparse = importdata("../Data/Instrumental/idx_sparse.txt");
[fid, msg] = fopen("../Data/Instrumental/anomalies_sparse.bin", "r");
if fid == -1
  error(msg);
endif
T_sparse = fread(fid, [year_count, numel(idx_sparse)], "int16");
T_sparse = single(T_sparse);
fclose(fid);

% Convert to degrees.
T_dense_monthly /= 100;
T_sparse        /= 100;

% Verification grid cells as a subset of calibration grid cells.
idx_sparse_subset_of_dense = zeros(size(idx_sparse));
for i = 1 : numel(idx_sparse)
  idx_sparse_subset_of_dense(i) = find(idx_dense == idx_sparse(i), 1);
endfor

% Coordinates of grid cell centers.
grid_size = single(5);
lon_count = 360 / grid_size;
lat_count = 180 / grid_size;
lon = linspace(0.5 * grid_size, 360 - 0.5 * grid_size, lon_count);
lon(lon > 180) -= 360;
lat = linspace(90 - 0.5 * grid_size, -90 + 0.5 * grid_size, lat_count);
[lon, lat] = meshgrid(lon, lat);
lon = reshape(lon, [], 1);
lat = reshape(lat, [], 1);
lon_dense  = lon(idx_dense);
lat_dense  = lat(idx_dense);
lat_sparse = lat(idx_sparse);

% Grid cell statistics.
mu_dense    = mean(T_dense_monthly(idx_cal_monthly, :), 1);
sigma_dense = std(detrend(T_dense_monthly(idx_cal_end_monthly, :)), 1, 1);

% Center the data.
T_dense_monthly = T_dense_monthly - mu_dense;
T_sparse        = T_sparse - mu_dense(idx_sparse_subset_of_dense);

% Compute singular value decomposition.
T_dense_monthly_standardized = T_dense_monthly ./ sigma_dense;
[U_monthly, S, V] = svd(T_dense_monthly_standardized(idx_cal_end_monthly, :) .* cosd(lat_dense)');

% Calculate annual averages of the left singular vectors.
U = zeros(size(U_monthly, 1) / 12, size(U_monthly, 2), "single");
for i = 1 : size(U, 1)
  U(i, :) = mean(U_monthly(12 * i - 11 : 12 * i, :), 1);
endfor

% Calculate annual averages of dense subset.
T_dense = zeros(size(T_dense_monthly, 1) / 12, size(T_dense_monthly, 2), "single");
for i = 1 : size(T_dense, 1)
  T_dense(i, :) = mean(T_dense_monthly(12 * i - 11 : 12 * i, :), 1);
endfor

% Calibration and verification data.
T_dense_cal = T_dense(idx_cal, :);
T_sparse_ver = T_sparse(idx_ver, :);

% Calculate global mean temperature based on dense subset.
w_glob_dense = cosd(lat_dense);
w_glob_dense /= sum(w_glob_dense);
T_glob_dense = T_dense * w_glob_dense;
T_glob_dense_cal = T_glob_dense(idx_cal);

% Calculate Northern Hemisphere mean temperature based on dense subset.
w_nhem_dense = cosd(lat_dense) .* (lat_dense > 0);
w_nhem_dense /= sum(w_nhem_dense);
T_nhem_dense = T_dense * w_nhem_dense;
T_nhem_dense_cal = T_nhem_dense(idx_cal);

% Calculate detrended Northern Hemisphere mean temperature based on dense subset.
T_detr_dense_cal = detrend(T_nhem_dense_cal);

% Calculate Nino mean temperature based on dense subset.
w_nino_dense = cosd(lat_dense) .* (lat_dense > -5 & lat_dense < 5 & lon_dense > -150 & lon_dense < -90);
w_nino_dense /= sum(w_nino_dense);
T_nino_dense = T_dense * w_nino_dense;
T_nino_dense_cal = T_nino_dense(idx_cal);

% Calculate global mean temperature based on sparse subset.
w_glob_sparse = cosd(lat_sparse);
w_glob_sparse /= sum(w_glob_sparse);
T_glob_sparse = T_sparse * w_glob_sparse;
T_glob_sparse_ver = T_glob_sparse(idx_ver);

% Calculate Northern Hemisphere mean temperature based on sparse subset.
w_nhem_sparse = cosd(lat_sparse) .* (lat_sparse > 0);
w_nhem_sparse /= sum(w_nhem_sparse);
T_nhem_sparse = T_sparse * w_nhem_sparse;
T_nhem_sparse_ver = T_nhem_sparse(idx_ver);

% Calculate RE statistics for eigenvector-filtered instrumental data.
group = {'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j'};
EOFLIST = {
[1 : 40],
[1 : 20],
[1, 2, 3, 4, 5, 7, 9, 11, 14, 15, 16],
[1, 2, 3, 4, 5, 7, 9, 11, 15],
[1, 2, 3, 5, 6, 8, 11, 15],
[1 : 5],
[1, 2, 5, 11, 15],
[1, 2, 11, 15],
[1, 2],
1};
RE_glob_filt_cal = nan(1, numel(group));
RE_nhem_filt_cal = nan(1, numel(group));
RE_detr_filt_cal = nan(1, numel(group));
RE_nino_filt_cal = nan(1, numel(group));
RE_mult_filt_cal = nan(1, numel(group));
for i = 1 : numel(group)
  % Filter gridded temperature.
  T_dense_filt_standardized = U(:, EOFLIST{i}) * S(EOFLIST{i}, EOFLIST{i}) * V(:, EOFLIST{i})' ./ cosd(lat_dense)';
  T_dense_filt = T_dense_filt_standardized .* sigma_dense;
  
  % Calculate average global and Northern Hemisphere filtered temperature.
  T_glob_dense_filt = T_dense_filt * w_glob_dense;
  T_nhem_dense_filt = T_dense_filt * w_nhem_dense;
  T_nino_dense_filt = T_dense_filt * w_nino_dense;
  
  % Calibration and verification statistics.
  T_glob_dense_filt_cal = T_glob_dense_filt(1 : length_cal);
  T_nhem_dense_filt_cal = T_nhem_dense_filt(1 : length_cal);
  T_detr_dense_filt_cal = detrend(T_nhem_dense_filt_cal);
  T_nino_dense_filt_cal = T_nino_dense_filt(1 : length_cal);
  T_dense_filt_cal = T_dense_filt(1 : length_cal, :);
  RE_glob_filt_cal(i) = 1 - sumsq(T_glob_dense_cal - T_glob_dense_filt_cal) / sumsq(T_glob_dense_cal);
  RE_nhem_filt_cal(i) = 1 - sumsq(T_nhem_dense_cal - T_nhem_dense_filt_cal) / sumsq(T_nhem_dense_cal);
  RE_detr_filt_cal(i) = 1 - sumsq(T_detr_dense_cal - T_detr_dense_filt_cal) / sumsq(T_detr_dense_cal);
  RE_nino_filt_cal(i) = 1 - sumsq(T_nino_dense_cal - T_nino_dense_filt_cal) / sumsq(T_nino_dense_cal);
  RE_mult_filt_cal(i) = 1 - sum(cumsum((T_dense_cal - T_dense_filt_cal).^2)(:)) / sum(cumsum(T_dense_cal.^2)(:));
endfor

% Calculate temperature reconstruction for all time steps.
t = [];
T_nhem_recon = [];
RE_glob_cal = nan(1, numel(start_step));
RE_glob_ver = nan(1, numel(start_step));
RE_nhem_cal = nan(1, numel(start_step));
RE_nhem_ver = nan(1, numel(start_step));
RE_detr_cal = nan(1, numel(start_step));
RE_nino_cal = nan(1, numel(start_step));
RE_mult_cal = nan(1, numel(start_step));
RE_mult_ver = nan(1, numel(start_step));
for i = 1 : numel(start_step)
  % Import proxy data.
  P = importdata(sprintf("../Data/Proxy/Data/data%d.txt", start_step(i)));
  t_step = P(:, 1);
  idx_cal_step = t_step >= start_cal & t_step <= end_cal;
  idx_ver_step = t_step >= start_ver & t_step <= end_ver;
  P(:, 1) = [];

  % Extrapolate proxy data.
  for j = 1 : size(P, 2)
    idx = find(~isnan(P(:, j)), 1, "last");
    P(idx + 1 : end, j) = P(idx, j);
  endfor

  % Standardize proxy data over calibration period.
  P -= mean(P(idx_cal_step, :), 1);
  P ./= std(detrend(P(idx_cal_step, :)), 1, 1);
  
  % Import MBH selection of EOFs.
  eoflist = importdata(sprintf("../Data/EOFs/eoflist%d.txt", start_step(i)));

  % Reconstruct singular vectors.
  G = U(1 : length_cal, eoflist) \ P(idx_cal_step, :);
  U_recon = P / G;
  
  % Rescale singular vectors.
  U_recon = U_recon ./ std(U_recon(idx_cal_step, :), 1, 1) .* std(U(1 : length_cal, eoflist), 1, 1);
  
  % Reconstruct gridded temperature.
  T_dense_recon_standardized = U_recon * S(eoflist, eoflist) * V(:, eoflist)' ./ cosd(lat_dense)';
  T_dense_recon = T_dense_recon_standardized .* sigma_dense;
  T_sparse_recon = T_dense_recon(:, idx_sparse_subset_of_dense);
  
  % Calculate average global reconstructed temperature.
  T_glob_dense_recon = T_dense_recon * w_glob_dense;
  T_glob_sparse_recon = T_sparse_recon * w_glob_sparse;

  % Calculate average Northern Hemisphere reconstructed temperature.
  T_nhem_dense_recon = T_dense_recon * w_nhem_dense;
  T_nhem_sparse_recon = T_sparse_recon * w_nhem_sparse;
  
  % Calculate average nino reconstructed temperature.
  T_nino_dense_recon = T_dense_recon * w_nino_dense;
  
  % Calibration and verification statistics.
  T_dense_recon_cal       = T_dense_recon(idx_cal_step, :);
  T_sparse_recon_ver      = T_sparse_recon(idx_ver_step, :);
  T_glob_dense_recon_cal  = T_glob_dense_recon(idx_cal_step);
  T_glob_sparse_recon_ver = T_glob_sparse_recon(idx_ver_step);
  T_nhem_dense_recon_cal  = T_nhem_dense_recon(idx_cal_step);
  T_nhem_sparse_recon_ver = T_nhem_sparse_recon(idx_ver_step);
  T_detr_dense_recon_cal  = detrend(T_nhem_dense_recon_cal);
  T_nino_dense_recon_cal  = T_nino_dense_recon(idx_cal_step);
  RE_glob_cal(i) = 1 - sumsq(T_glob_dense_cal - T_glob_dense_recon_cal) / sumsq(T_glob_dense_cal);
  RE_glob_ver(i) = 1 - sumsq(T_glob_sparse_ver - T_glob_sparse_recon_ver) / sumsq(T_glob_sparse_ver);
  RE_nhem_cal(i) = 1 - sumsq(T_nhem_dense_cal - T_nhem_dense_recon_cal) / sumsq(T_nhem_dense_cal);
  RE_nhem_ver(i) = 1 - sumsq(T_nhem_sparse_ver - T_nhem_sparse_recon_ver) / sumsq(T_nhem_sparse_ver);
  RE_detr_cal(i) = 1 - sumsq(T_detr_dense_cal - T_detr_dense_recon_cal) / sumsq(T_detr_dense_cal);
  RE_nino_cal(i) = 1 - sumsq(T_nino_dense_cal - T_nino_dense_recon_cal) / sumsq(T_nino_dense_cal);
  RE_mult_cal(i) = 1 - sum(cumsum((T_dense_cal - T_dense_recon_cal).^2)(:)) / sum(cumsum(T_dense_cal.^2)(:));
  RE_mult_ver(i) = 1 - sum(cumsum((T_sparse_ver - T_sparse_recon_ver).^2)(:)) / sum(cumsum(T_sparse_ver.^2)(:));
  
  % Splice.
  t = [t; t_step(1 : length_step(i))];
  T_nhem_recon = [T_nhem_recon; T_nhem_dense_recon(1 : length_step(i))];
endfor

% Write reconstruction to file.
mkdir "../Data/Emulation";
[fid, msg] = fopen("../Data/Emulation/nhem_recon.txt", "w");
if fid == -1
  error(msg);
endif
fputs(fid, "Year   NH recon\n");
for i = 1 : numel(t)
  fprintf(fid, "%4d %11.7f", t(i), T_nhem_recon(i));
  if i ~= numel(t)
    fprintf(fid, "\n");
  endif
endfor
fclose(fid);

% Write cross-validation statistics to file.
[fid, msg] = fopen("../Data/Emulation/cross_validation.txt", "w");
if fid == -1
  error(msg);
endif
fputs(fid, "Filtered instrumental data RE statistics:\n\n");
fputs(fid, "Group  GLB   NH    DET   NIN   MULT\n");
fputs(fid, "-----------------------------------\n");
for i = 1 : numel(group)
    fprintf(fid, "%c     %5.2f %5.2f %5.2f %5.2f %5.2f\n",
            group{i}, RE_glob_filt_cal(i), RE_nhem_filt_cal(i),
                      RE_detr_filt_cal(i), RE_nino_filt_cal(i), RE_mult_filt_cal(i));
endfor
fputs(fid, "\nReconstruction RE statistics:\n\n");
fputs(fid, "       Calibration                    Verification\n");
fputs(fid, "Step   GLB   NH    DET   NIN   MLT    GLB   NH    MLTA\n"); 
fputs(fid, "------------------------------------------------------\n");
for i = numel(start_step) : -1 : 1
  fprintf(fid, "%4d  %5.2f %5.2f %5.2f %5.2f %5.2f  %5.2f %5.2f %5.2f\n",
          start_step(i), RE_glob_cal(i), RE_nhem_cal(i), RE_detr_cal(i), RE_nino_cal(i),
                         RE_mult_cal(i), RE_glob_ver(i), RE_nhem_ver(i), RE_mult_ver(i));
endfor
fclose(fid);
