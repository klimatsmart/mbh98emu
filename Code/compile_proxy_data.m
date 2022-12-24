% This script extracts all necessary proxy data from the mbh98.tar archive.

% Make directory for data files.
mkdir "../Data/Proxy/Data";

% Unpack proxy data.
if ~exist("../Data/Proxy/mbh98")
  untar("../Data/Proxy/mbh98.tar", "../Data/Proxy/mbh98");
endif

% List of proxy rosters.
L = dir("../Data/Proxy/Lists");

% For each proxy list, write a data file.
for i = 1 : numel(L)
  if ~L(i).isdir
    filename_in = L(i).name;
    filename_out = filename_in([1 : 4, 9 : end]);

    start_year = str2double(filename_in(9 : 12));
    end_year = 1980;
    t = (start_year : end_year)';
    
    % Read proxy list.
    datalist = importdata(fullfile("../Data/Proxy/Lists", filename_in));
    
    % Proxy matrix.
    P = nan(numel(t), numel(datalist) + 1);

    % First column is time.
    P(:, 1) = t;

    % Add proxy series to matrix.
    for j = 1 : numel(datalist)
      % Import proxy data.
      p = importdata(fullfile("../Data/Proxy/mbh98", datalist{j}));
      
      % Truncate.
      p = p(p(:, 1) >= t(1) & p(:, 1) <= t(end), :);
      
      % Add to matrix.
      idx = p(:, 1) - t(1) + 1;
      P(idx, j + 1) = p(:, 2);
    endfor
    
    % Write matrix to file.
    [fid, msg] = fopen(fullfile("../Data/Proxy/Data", filename_out), "w");
    if fid == -1
      error(msg);
    endif
    for j = 1 : size(P, 1)
      fprintf(fid, "%4d", P(j, 1));
      for k = 2 : size(P, 2)
        fprintf(fid, "%14.6e", P(j, k));
      endfor
      if j < size(P, 1)
        fprintf(fid, "\n");
      endif
    endfor
    fclose(fid);
  endif
endfor
