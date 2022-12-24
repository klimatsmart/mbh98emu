% This script reads the file anomalies-new and saves the data as binary.

% Skip if already done.
if exist("../Data/Instrumental/anomalies_all.bin")
  return
endif

% Start and end of instrumental record.
start_instr = 1854;
end_instr = 1993;
month_count = (end_instr - start_instr + 1) * 12;

% Grid parameters.
grid_size = 5; % 5x5 degrees.
lat_count = 180 / grid_size;
lon_count = 360 / grid_size;
cell_count = lat_count * lon_count;

% Data matrix.
T = zeros(month_count, cell_count);

% Read gridded temperature data.
[fid, msg] = fopen("../Data/Instrumental/anomalies-new", "r");
if fid == -1
  error(msg);
endif
line_count = 0;
col = 0;
month = 0;
cells_per_line = 18;
lines_per_month = cell_count / cells_per_line + 1;
while ischar(tline = fgetl(fid))
  if mod(line_count, lines_per_month) == 0
    month += 1;
    col = 0;
  else
    Tline = sscanf(tline, "%d");
    T(month, col + (1 : cells_per_line)) = Tline;
    col += cells_per_line;
  endif
  line_count += 1;
endwhile
fclose(fid);

% Save data as binary for ease of use.
[fid, msg] = fopen("../Data/Instrumental/anomalies_all.bin", "w");
if fid == -1
  error(msg);
endif
fwrite(fid, T, "int16");
fclose(fid);
