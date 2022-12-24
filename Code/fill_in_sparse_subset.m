% This script fills in missing data in the sparse subset of the instrumental record
% and saves the data as binary.

% Time parameters.
start_instr = 1854;
end_instr = 1993;
start_ver = 1854;
end_ver = 1980;
idx_ver = (start_ver - start_instr) * 12 + 1 : (end_ver - start_instr + 1) * 12;
month_count = (end_instr - start_instr + 1) * 12;

% Grid parameters.
grid_size = 5;
lat_count = 180 / grid_size;
lon_count = 360 / grid_size;
cell_count = lat_count * lon_count;

% Load temperature data.
[fid, msg] = fopen("../Data/Instrumental/anomalies_all.bin");
if fid == -1
  error(msg);
endif
T_raw = fread(fid, [month_count, cell_count], "int16");
fclose(fid);

% Rearrange to agree with MBH grid cell order.
M_idx = reshape(1 : cell_count, lon_count, lat_count)';
v_idx = reshape(M_idx, 1, cell_count);
T_raw = T_raw(:, v_idx);

% Select sparse subset.
T_raw_ver = T_raw(idx_ver, :);
idx_sparse = [];
for i = 1 : size(T_raw_ver, 2)
  idx     = find(T_raw_ver(:, i) ~= -9999);
  idx_nan = find(T_raw_ver(:, i) == -9999);
  if numel(idx_nan) > 120
    continue
  elseif idx(1) > 24 || idx(end) <= size(T_raw_ver, 1) - 24
    continue
  elseif max(diff(idx)) > 24
    continue
  else
    idx_sparse = [idx_sparse, i];
  endif
endfor
T_raw = T_raw(:, idx_sparse);

% Fill in missing data.
T_monthly = fill_in_matrix(T_raw);

% Calculate annual averages.
T = zeros(size(T_monthly, 1) / 12, size(T_monthly, 2), "single");
for i = 1 : size(T, 1)
  for j = 1 : 12
    T(i, :) += T_monthly(12 * (i - 1) + j, :);
  endfor
  T(i, :) /= 12;
endfor

% Round to two decimal places.
for i = 1 : size(T, 1)
  for j = 1 : size(T, 2)
    T(i, j) = str2double(sprintf("%.2f", T(i, j)));
  endfor
endfor

% Convert back to hundredths of a degree.
T = round(100 * T);

% Save infilled temperature data as binary.
[fid, msg] = fopen("../Data/Instrumental/anomalies_sparse.bin", "w");
if fid == -1
  error(msg);
endif
fwrite(fid, T, "int16");
fclose(fid);

% Save grid cell indices.
dlmwrite("../Data/Instrumental/idx_sparse.txt", idx_sparse');
