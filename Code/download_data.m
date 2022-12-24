% Download MBH98 reconstruction.
mkdir "../Data/Original reconstruction";
if ~exist("../Data/Original reconstruction/nhmean.dat")
  [~, success, msg] = urlwrite("http://www.meteo.psu.edu/holocene/public_html/shared/research/MANNETAL98/nhmean.dat", "../Data/Original reconstruction/nhmean.dat");
  if ~success
    error(msg);
  endif
endif

% Download instrumental data.
mkdir "../Data/Instrumental";
if ~exist("../Data/Instrumental/anomalies-new")
  [~, success, msg] = urlwrite("http://www.meteo.psu.edu/holocene/public_html/shared/research/MANNETAL98/INSTRUMENTAL/anomalies-new", "../Data/Instrumental/anomalies-new");
  if ~success
    error(msg);
  endif
endif

% Download proxy data.
if ~exist("../Data/Proxy/mbh98.tar")
  [~, success, msg] = urlwrite("https://www.sealevel.info/FOIA/2009/FOIA/documents/mbh98-osborn/mbh98.tar", "../Data/Proxy/mbh98.tar");
  if ~success
    error(msg);
  endif
endif
