% Load emulation and original Northern Hemisphere temperature reconstruction.
nhem_recon = importdata("../Data/Emulation/nhem_recon.txt", " ", 1).data;
nhem_MBH98 = importdata("../Data/Original reconstruction/nhmean.dat", " ", 1).data;

% Plot reconstruction.
length_recon = size(nhem_recon, 1);
plot(nhem_MBH98(1 : length_recon, 1), nhem_MBH98(1 : length_recon, 2),
     "r-o", "MarkerFaceColor", "r", "MarkerSize", 5,
     nhem_recon(:, 1), nhem_recon(:, 2),
     "b-o", "MarkerFaceColor", "b", "MarkerSize", 3,
     [1400, 2000], [0, 0],
     "k--");
title("Northern Hemisphere temperature reconstruction");
xlabel("Year");
ylabel("Temperature anomaly (°C)\nwith respect to 1902–1980");
legend({"MBH98", "Emulation"}, "location", "northwest");

% Display cross-validation statistics.
type "../Data/Emulation/cross_validation.txt";
