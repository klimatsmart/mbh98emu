% Download data.
printf("Downloading data...\n");
download_data;

% Read instrumental data and save as binary.
printf("Reading instrumental data...\n");
read_instrumental_data;

% Fill in missing data.
printf("Infilling dense subset...\n");
infill_dense_subset;
printf("Infilling sparse subset...\n");
infill_sparse_subset;

% Write proxy data to files.
printf("Writing proxy data to files...\n");
compile_proxy_data;

% Calculate reconstruction.
printf("Calculating reconstruction...\n");
hockey_stick;

printf("Done.\n");
