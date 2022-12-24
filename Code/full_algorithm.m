printf("Progress:\n");

printf("Downloading data...\n");
download_data;

printf("Reading instrumental data...\n");
read_instrumental_data;

printf("Filling in dense subset...\n");
fill_in_dense_subset;
printf("Filling in sparse subset...\n");
fill_in_sparse_subset;

printf("Writing proxy data to files...\n");
compile_proxy_data;

printf("Calculating reconstruction...\n");
hockey_stick;

printf("Done.\n\n");

display_results;
