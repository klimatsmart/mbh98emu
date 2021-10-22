% Fill in the columns of X by linear interpolation.
function Y = infill_matrix(X)
  Y = single(X);
  Y(X == -9999) = nan;
  Y /= 100;
  n = size(Y, 1);
  for i = 1 : size(Y, 2)
    % Extrapolate backwards.
    idx = find(~isnan(Y(:, i)), 1);
    Y(1 : idx - 1, i) = Y(idx, i);
    
    % Interpolate linearly.
    idx = find(~isnan(Y(:, i)));
    gap = diff(idx);
    for j = 1 : numel(gap)
      if gap(j) > 1
        t1 = idx(j);
        t2 = idx(j + 1);
        T1 = Y(t1, i);
        T2 = Y(t2, i);
        dT = T2 - T1;
        dt = t2 - t1;
        for k = 1 : dt - 1
          Y(t1 + k, i) = T1 + (k * dT) / dt;
        endfor
      endif
    endfor
    
    % Extrapolate forwards.
    idx = find(~isnan(Y(:, i)), 1, "last");
    if idx < n - 1
      % Extrapolate towards zero.
      t1 = idx;
      t2 = n + 1;
      T1 = Y(t1, i);
      T2 = 0;
      dT = T2 - T1;
      dt = t2 - t1;
      for k = 1 : dt - 1
        Y(t1 + k, i) = T1 + (k * dT) / dt;
      endfor
    elseif idx == n - 1
      % Edge case.
      Y(n, i) = 0.5 * (Y(n - 2, i) + Y(n - 1, i));
    endif
  endfor
endfunction
