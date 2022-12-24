% This script fills in missing data in the dense subset of the instrumental
% record and saves the data as binary.

% Time parameters.
start_instr = 1854;
end_instr = 1993;
start_cal = 1902;
end_cal = 1980;
idx_cal = (start_cal - start_instr) * 12 + 1 : (end_cal - start_instr + 1) * 12;
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

% Select dense subset.
T_raw_cal = T_raw(idx_cal, :);
idx_dense = [];
for i = 1 : size(T_raw_cal, 2)
  idx = find(T_raw_cal(:, i) ~= -9999);
  if isempty(idx)
    continue
  elseif idx(1) > 36 || idx(end) <= size(T_raw_cal, 1) - 36
    continue
  elseif max(diff(idx)) > 36
    continue
  else
    idx_dense = [idx_dense, i];
  endif
endfor
T_raw = T_raw(:, idx_dense);

% Sign error (comment out this line to fix the error).
T_raw(T_raw < -999 & T_raw ~= -9999) *= -1;

% Fill in missing data.
T = fill_in_matrix(T_raw);

% Round to two decimal places.
for i = 1 : size(T, 1)
  for j = 1 : size(T, 2)
    if T_raw(i, j) == -9999
      T(i, j) = str2double(sprintf("%.2f", T(i, j)));
    endif
  endfor
endfor

% Convert back to hundredths of a degree.
T = round(100 * T);

% Save infilled temperature data as binary.
[fid, msg] = fopen("../Data/Instrumental/anomalies_dense.bin", "w");
if fid == -1
  error(msg);
endif
fwrite(fid, T, "int16");
fclose(fid);

% Save grid cell indices.
dlmwrite("../Data/Instrumental/idx_dense.txt", idx_dense');
