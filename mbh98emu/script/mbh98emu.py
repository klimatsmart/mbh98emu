#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script emulates the MBH98 temperature reconstruction.

Tested with Python 3.9.2 and the following packages installed:
    
    * requests  2.25.1
    * NumPy     1.19.5
    * SciPy     1.6.0
    * pandas    1.1.5
"""
import pathlib
import shutil
import tarfile
import csv
import requests
import numpy as np
import scipy.signal as sps
import pandas as pd

# Prepare directories.
ROOT_PATH = pathlib.Path(__file__).resolve().parents[1]
CONFIG_PATH = ROOT_PATH.joinpath("config")
DOWNLOAD_PATH = ROOT_PATH.joinpath("downloads")
INSTRUMENTAL_PATH = ROOT_PATH.joinpath("instrumental")
PROXY_PATH = ROOT_PATH.joinpath("proxy")
RECONSTRUCTION_PATH = ROOT_PATH.joinpath("reconstruction")
VALIDATION_PATH = ROOT_PATH.joinpath("validation")
for path in [RECONSTRUCTION_PATH, VALIDATION_PATH]:
    if path.is_dir():
        shutil.rmtree(path)
for path in [DOWNLOAD_PATH, INSTRUMENTAL_PATH, PROXY_PATH,
             RECONSTRUCTION_PATH, VALIDATION_PATH]:
    path.mkdir(exist_ok=True)

# Instrumental, calibration and verification periods.
INSTR_START = 1854
INSTR_END = 1993
CAL_START = 1902
CAL_END = 1980
VER_START = 1854
VER_END = 1901


def get_instrumental_data():
    url = "http://www.meteo.psu.edu/holocene/public_html/shared/research/MANNETAL98/INSTRUMENTAL/anomalies-new"
    path = DOWNLOAD_PATH.joinpath("anomalies-new")
    if not path.is_file():
        print(f"Downloading instrumental data from {url}...")
        get_data_from_url(url, path)


def get_proxy_data():
    url = "https://www.sealevel.info/FOIA/2009/FOIA/documents/mbh98-osborn/mbh98.tar"
    path = DOWNLOAD_PATH.joinpath("mbh98.tar")
    if not path.is_file():
        print(f"Downloading proxy data from {url}...")
        get_data_from_url(url, path)


def get_data_from_url(url, path):
    r = requests.get(url, stream=True)
    r.raise_for_status()
    with open(path, "wb") as f:
        for chunk in r.iter_content(chunk_size=1024**2):
            f.write(chunk)


def prepare_instrumental_data():
    print("Preparing instrumental data...")
    create_instrumental_dataframe()
    prepare_dense_subset()
    prepare_sparse_subset()
    compute_gridbox_statistics()
    center_annual_means()
    compute_instrumental_svd()


def read_instrumental_data():
    month_count = 12 * (INSTR_END - INSTR_START + 1)
    box_count = 72 * 36
    date = []
    temperature = np.full((month_count, box_count), np.nan)
    lines_per_month = box_count/18 + 1
    path = DOWNLOAD_PATH.joinpath("anomalies-new")
    with open(path, "r") as f:
        month = -1
        for count, line in enumerate(f):
            if count % lines_per_month == 0:
                date.append(f"{int(line[0:5])}-{int(line[5:10])}")
                month += 1
                col = 0
            else:
                temperature[month, col:col+18] = [int(line[i:i+5])
                                                  for i in range(0, 90, 5)]
                col += 18
    return date, temperature


def rearrange_instrumental_data(t):
    # Rearrange data to match MBH98 gridbox order.
    lon_count = 72
    lat_count = 36
    box_count = lon_count * lat_count
    index = np.arange(box_count)
    index = np.reshape(index, (lat_count, lon_count), order="C")
    index = index.flatten(order="F")
    t[:] = t[:, index]


def instrumental_gridpoints():
    lat = np.arange(90, -90, -5) - 2.5
    lon = np.arange(0, 360, 5) + 2.5
    lon[lon > 180] -= 360
    lon, lat = np.meshgrid(lon, lat)
    lat = lat.flatten(order="F")
    lon = lon.flatten(order="F")
    return lat, lon


def create_instrumental_dataframe():
    date, temperature = read_instrumental_data()
    rearrange_instrumental_data(temperature)
    lat, lon = instrumental_gridpoints()
    index = pd.to_datetime(date)
    columns = pd.MultiIndex.from_arrays((lat, lon), names=("lat", "lon"))
    df = pd.DataFrame(data=temperature, index=index, columns=columns)
    path = INSTRUMENTAL_PATH.joinpath("anomalies-new.pkl")
    df.to_pickle(path)


def dense_subset():
    path = INSTRUMENTAL_PATH.joinpath("anomalies-new.pkl")
    df = pd.read_pickle(path)
    t = df.loc[str(CAL_START):str(CAL_END), :].to_numpy()
    columns = []
    for i, col in enumerate(t.T):
        val_index = (col != -9999).nonzero()[0]
        val_index = np.insert(val_index, 0, -1)
        val_index = np.append(val_index, col.size)
        if np.diff(val_index).max() <= 36:
            columns.append(df.columns[i])
    return df.loc[:, columns]


def sparse_subset():
    path = INSTRUMENTAL_PATH.joinpath("anomalies-new.pkl")
    df = pd.read_pickle(path)
    t = df.loc[str(VER_START):str(CAL_END), :].to_numpy()
    columns = []
    for i, col in enumerate(t.T):
        nan_index = (col == -9999).nonzero()[0]
        val_index = (col != -9999).nonzero()[0]
        val_index = np.insert(val_index, 0, -1)
        val_index = np.append(val_index, col.size)
        if nan_index.size <= 120 and np.diff(val_index).max() <= 24:
            columns.append(df.columns[i])
    return df.loc[:, columns]


def prepare_dense_subset():
    df = dense_subset()
    convert_to_degrees(df)
    invert_anomalies_below_minus_10(df)
    fill_in_instrumental_data(df)
    round_to_2_decimal_places(df)
    df = df.loc[str(CAL_START):]
    path = INSTRUMENTAL_PATH.joinpath("dense_subset_monthly.pkl")
    df.to_pickle(path)


def prepare_sparse_subset():
    df = sparse_subset()
    convert_to_degrees(df)
    fill_in_instrumental_data(df)
    df = annual_means(df)
    round_to_2_decimal_places(df)
    path = INSTRUMENTAL_PATH.joinpath("sparse_subset.pkl")
    df.to_pickle(path)


def fill_in_instrumental_data(df):
    t = df.to_numpy()
    n = t.shape[0]
    for i in range(t.shape[1]):
        # Indices with data.
        j = (~np.isnan(t[:, i])).nonzero()[0]
        if j.size > 0:
            # Constant extrapolation to the left.
            t[:j[0], i] = t[j[0], i]
            
            # Linear interpolation.
            interpolate32(t[:, i])
            
            # Extrapolation to the right.
            if j[-1] < n - 2:
                # Linear extrapolation towards zero.
                gap = n - j[-1]
                t[j[-1]+1:, i] = linspace32(t[j[-1], i], 0, gap+1)[1:-1]
            elif j[-1] == n - 2:
                # Special rule when only the last value is missing.
                t[-1, i] = 0.5 * add32(t[-3, i], t[-2, i])
    df[:] = t


def interpolate32(t):
    index = (~np.isnan(t)).nonzero()[0]
    gaps = np.diff(index)
    for i, gap in enumerate(gaps):
        if gap > 1:
            t0 = t[index[i]]
            t1 = t[index[i+1]]
            t[index[i]+1:index[i+1]] = linspace32(t0, t1, gap+1)[1:-1]


def linspace32(t0, t1, n):
    k = np.arange(n)
    t = add32(t0, div32(mul32(k, sub32(t1, t0)), n-1))
    t[-1] = np.float32(t1)
    return t


def add32(a, b):
    return np.add(a, b, dtype=np.float32)


def sub32(a, b):
    return np.subtract(a, b, dtype=np.float32)


def mul32(a, b):
    return np.multiply(a, b, dtype=np.float32)


def div32(a, b):
    return np.divide(a, b, dtype=np.float32)


def array_sum32(a):
    s = np.float32(0)
    for num in a:
        s = add32(s, num)
    return s


def column_mean32(a):
    s = np.zeros(a.shape[1], dtype=np.float32)
    for row in a:
        s = add32(s, row)
    mu = div32(s, a.shape[0])
    return mu


def convert_to_degrees(df):
    t = df.to_numpy()
    t[t == -9999] = np.nan
    t /= 100


def invert_anomalies_below_minus_10(df):
    # This function emulates a programming error.
    t = df.to_numpy()
    t[t <= -10] *= -1


def round_to_2_decimal_places(df):
    t = df.to_numpy()
    t[:] = t.round(2)


def annual_means(df):
    years = sorted(set(df.index.year))
    df_annual = pd.DataFrame(index=years, columns=df.columns,
                              dtype=np.float32)
    for year in years:
        t = df.loc[df.index.year == year].to_numpy()
        df_annual.loc[year, :] = column_mean32(t)
    return df_annual


def compute_gridbox_statistics():
    # Load instrumental data.
    path = INSTRUMENTAL_PATH.joinpath("dense_subset_monthly.pkl")
    df = pd.read_pickle(path)
    t = df.to_numpy()

    # Means.
    month_count = 12 * (CAL_END - CAL_START + 1)
    mu = t[:month_count, :].mean(axis=0)
    
    # Detrended standard deviations.
    sigma = sps.detrend(t, axis=0).std(axis=0, ddof=0)
    
    # Save statistics.
    path = INSTRUMENTAL_PATH.joinpath("statistics")
    path.mkdir(exist_ok=True)
    mu_path = path.joinpath("mu.npy")
    sigma_path = path.joinpath("sigma.npy")
    np.save(mu_path, mu)
    np.save(sigma_path, sigma)


def center_annual_means():
    # Load instrumental data.
    dense_path = INSTRUMENTAL_PATH.joinpath("dense_subset_monthly.pkl")
    sparse_path = INSTRUMENTAL_PATH.joinpath("sparse_subset.pkl")
    df_dense_monthly = pd.read_pickle(dense_path)
    df_sparse = pd.read_pickle(sparse_path)
    
    # Load means for dense subset.
    mu_path = INSTRUMENTAL_PATH.joinpath("statistics", "mu.npy")
    mu = np.load(mu_path)
    
    # Center instrumental data.
    mu = pd.Series(data=mu, index=df_dense_monthly.columns)
    df_dense_monthly -= mu
    df_sparse -= mu[df_sparse.columns] # Not properly centered.
    
    # Annual means for dense subset.
    df_dense = df_dense_monthly.groupby(df_dense_monthly.index.year).mean()
    
    # Save centered data.
    dense_path = INSTRUMENTAL_PATH.joinpath("dense_subset_centered.pkl")
    sparse_path = INSTRUMENTAL_PATH.joinpath("sparse_subset_centered.pkl")
    pd.to_pickle(df_dense, dense_path)
    pd.to_pickle(df_sparse, sparse_path)


def compute_instrumental_svd():
    # Load dense subset of instrumental data.
    path = INSTRUMENTAL_PATH.joinpath("dense_subset_monthly.pkl")
    df_monthly = pd.read_pickle(path)

    # Statistics for standardization.
    mu_path = INSTRUMENTAL_PATH.joinpath("statistics", "mu.npy")
    sigma_path = INSTRUMENTAL_PATH.joinpath("statistics", "sigma.npy")
    mu = np.load(mu_path)
    sigma = np.load(sigma_path)

    # Standardize data.
    t = df_monthly.to_numpy()
    t_zscores = (t - mu) / sigma

    # Compute singular value decomposition of standardized and area
    # weighted instrumental data.
    lat = np.array(df_monthly.columns.get_level_values(0))
    w = np.cos(np.deg2rad(lat))
    t_weighted = t_zscores * w
    u_monthly, s, vt = np.linalg.svd(t_weighted)
    u_monthly = pd.DataFrame(data=u_monthly, index=df_monthly.index)
    u = u_monthly.groupby(u_monthly.index.year).mean()
    vt = pd.DataFrame(data=vt, columns=df_monthly.columns)
    
    # Save svd.
    svd_path = INSTRUMENTAL_PATH.joinpath("svd")
    svd_path.mkdir(exist_ok=True)
    u_path = svd_path.joinpath("u.pkl")
    s_path = svd_path.joinpath("s.npy")
    vt_path = svd_path.joinpath("vt.pkl")
    u.to_pickle(u_path)
    np.save(s_path, s)
    vt.to_pickle(vt_path)


def prepare_proxy_data():
    print("Preparing proxy data...")
    extract_proxy_archive()
    create_proxy_matrices()


def extract_proxy_archive():
    tar_path = DOWNLOAD_PATH.joinpath("mbh98.tar")
    untar_path = PROXY_PATH.joinpath("mbh98")
    untar(tar_path, untar_path)


def untar(tar, path):
    if not path.is_dir():
        with tarfile.TarFile(tar, "r") as f:
            f.extractall(path)


def create_proxy_matrices():
    lists_path = CONFIG_PATH.joinpath("proxy")
    data_path = PROXY_PATH.joinpath("networks")
    data_path.mkdir(exist_ok=True)
    for step in reconstruction_steps():
        datalist_path = lists_path.joinpath(f"datalist{step:04d}.dat")
        proxy = proxy_matrix(datalist_path)
        file_path = data_path.joinpath(f"data{step:04d}.dat")
        save_proxy_matrix(proxy, file_path)


def proxy_matrix(datalist_path):
    # Create a proxy data matrix with time in the first column.
    mbh98_path = PROXY_PATH.joinpath("mbh98")
    with open(datalist_path, "r") as f:
        relative_paths = f.read().splitlines()
    year = int(datalist_path.stem.strip("datalist"))
    n = CAL_END - year + 1
    m = len(relative_paths) + 1
    proxy = np.full((n, m), np.nan)
    proxy[:, 0] = np.arange(year, CAL_END + 1)
    for i, relative_path in enumerate(relative_paths):
        data_path = mbh98_path.joinpath(*relative_path.split("/"))
        p = pd.read_table(data_path, header=None)
        p = p[0].str.split(expand=True).astype(float).to_numpy()
        p = submatrix(p, year, CAL_END)
        index = p[:, 0].astype(int) - year
        proxy[index, i+1] = p[:, 1]
    return proxy


def save_proxy_matrix(proxy, path):
    m = proxy.shape[1]
    np.savetxt(path, proxy, fmt="%04d" + "%14.6e" * (m-1))


def fill_in_proxy_matrix(p):
    for i in range(p.shape[1]):
        if np.isnan(p[-1, i]):
            index = (~np.isnan(p[:, i])).nonzero()[0]
            p[index[-1]:, i] = p[index[-1], i]


def standardize_proxy_matrix(p, t0, t1):
    # Standardize proxy matrix p over time period [t0, t1].
    p_cal = submatrix(p, t0, t1)
    mu = p_cal[:, 1:].mean(axis=0)
    sigma = sps.detrend(p_cal[:, 1:], axis=0).std(axis=0, ddof=0)
    p[:, 1:] -= mu
    p[:, 1:] /= sigma


def submatrix(x, t0, t1):
    # Submatrix of x where the first column is in [t0, t1].
    return x[(x[:, 0] >= t0) & (x[:, 0] <= t1), :]


def reconstruct_temperature():
    print("Generating reconstruction...")
    for step in reconstruction_steps():
        recon = reconstructed_temperature_field(step)
        analyze_regions(recon)
    
    # Compile results.
    concatenate_nhem_reconstructions()
    make_re_table()


def reconstruction_steps():
    path = CONFIG_PATH.joinpath("steps.csv")
    with open(path, "r", newline="") as f:
        reader = csv.reader(f)
        steps = next(reader)
    steps = [int(n) for n in steps]
    return steps


def reconstructed_temperature_field(step):
    u_recon, s_recon, vt_recon = reconstructed_svd(step)
    
    # Undo weighting and scaling.
    sigma_path = INSTRUMENTAL_PATH.joinpath("statistics", "sigma.npy")
    sigma = np.load(sigma_path)
    lat = np.array(vt_recon.columns.get_level_values(0))
    w = np.cos(np.deg2rad(lat))
    recon_weighted = (u_recon.to_numpy() * s_recon) @ vt_recon.to_numpy()
    recon_zscores = recon_weighted / w
    recon = recon_zscores * sigma
    recon = pd.DataFrame(data=recon, index=u_recon.index,
                         columns=vt_recon.columns)
    return recon


def reconstructed_svd(step):
    # Load proxy matrix.
    filename = f"data{step:04d}.dat"
    proxy_path = PROXY_PATH.joinpath("networks", filename)
    p = np.genfromtxt(proxy_path)
    fill_in_proxy_matrix(p)
    standardize_proxy_matrix(p, CAL_START, CAL_END)
    p_cal = submatrix(p, CAL_START, CAL_END)
    index_cal = (p[:, 0] >= CAL_START) & (p[:, 0] <= CAL_END)
    
    # Load instrumental PC retention.
    filename = f"eoflist{step:04d}.csv"
    eofs_path = CONFIG_PATH.joinpath("instrumental", filename)
    with open(eofs_path, "r", newline="") as f:
        reader = csv.reader(f)
        eofs = next(reader)
    eofs = [int(n) - 1 for n in eofs]
    
    # Load instrumental svd.
    svd_path = INSTRUMENTAL_PATH.joinpath("svd")
    u_path = svd_path.joinpath("u.pkl")
    s_path = svd_path.joinpath("s.npy")
    vt_path = svd_path.joinpath("vt.pkl")
    u = pd.read_pickle(u_path)
    s = np.load(s_path)
    vt = pd.read_pickle(vt_path)
    
    # Calibration period.
    u = u.loc[CAL_START:CAL_END, eofs]
    u = u.to_numpy()
    
    # Calibrate proxy data against instrumental PCs.
    g = np.linalg.lstsq(u, p_cal[:, 1:], rcond=None)[0]
    
    # Reconstruct PCs.
    u_recon = np.linalg.lstsq(g.T, p[:, 1:].T, rcond=None)[0].T
    
    # Rescale PCs.
    sigma_recon = np.std(u_recon[index_cal, :], axis=0, ddof=0)
    sigma_instr = np.std(u, axis=0, ddof=0)
    u_recon *= sigma_instr / sigma_recon
    
    # Create PC dataframe.
    u_recon = pd.DataFrame(data=u_recon, index=p[:, 0].astype(int))
    
    return u_recon, s[eofs], vt.iloc[eofs]


def analyze_regions(recon):
    # Load instrumental data.
    dense_path = INSTRUMENTAL_PATH.joinpath("dense_subset_centered.pkl")
    sparse_path = INSTRUMENTAL_PATH.joinpath("sparse_subset_centered.pkl")
    dense_instr = pd.read_pickle(dense_path)
    sparse_instr = pd.read_pickle(sparse_path)
    
    # Instrumental means.
    dense_glob_instr = regional_mean(dense_instr, "glob")
    dense_nhem_instr = regional_mean(dense_instr, "nhem")
    dense_detr_instr = sps.detrend(dense_nhem_instr.loc[CAL_START:CAL_END])
    dense_detr_instr = pd.DataFrame(data=dense_detr_instr,
                                    index=np.arange(CAL_START, CAL_END + 1))
    dense_nino_instr = regional_mean(dense_instr, "nino")
    sparse_glob_instr = regional_mean(sparse_instr, "glob")
    sparse_nhem_instr = regional_mean(sparse_instr, "nhem")
    
    # Reconstructed means.
    dense_recon = recon
    dense_glob_recon = regional_mean(dense_recon, "glob")
    dense_nhem_recon = regional_mean(dense_recon, "nhem")
    dense_detr_recon = sps.detrend(dense_nhem_recon.loc[CAL_START:CAL_END])
    dense_detr_recon = pd.DataFrame(data=dense_detr_recon,
                                    index=np.arange(CAL_START, CAL_END + 1))
    dense_nino_recon = regional_mean(dense_recon, "nino")
    sparse_recon = recon.loc[:, sparse_instr.columns]
    sparse_glob_recon = regional_mean(sparse_recon, "glob")
    sparse_nhem_recon = regional_mean(sparse_recon, "nhem")
    
    # Reduction of Error (RE) statistics.
    glob_cal_re = re_statistic(dense_glob_recon, dense_glob_instr,
                               CAL_START, CAL_END)
    nhem_cal_re = re_statistic(dense_nhem_recon, dense_nhem_instr,
                               CAL_START, CAL_END)
    detr_cal_re = re_statistic(dense_detr_recon, dense_detr_instr,
                               CAL_START, CAL_END)
    nino_cal_re = re_statistic(dense_nino_recon, dense_nino_instr,
                               CAL_START, CAL_END)
    mult_cal_re = re_mult_statistic(dense_recon, dense_instr,
                                    CAL_START, CAL_END)
    glob_ver_re = re_statistic(sparse_glob_recon, sparse_glob_instr,
                               VER_START, VER_END)
    nhem_ver_re = re_statistic(sparse_nhem_recon, sparse_nhem_instr,
                               VER_START, VER_END)
    mult_ver_re = re_mult_statistic(sparse_recon, sparse_instr,
                                    VER_START, VER_END)
    
    # Save nhem reconstruction.
    step = recon.index[0]
    step_path = RECONSTRUCTION_PATH.joinpath("steps")
    step_path.mkdir(exist_ok=True)
    nhem_path = step_path.joinpath(f"nhem{step:04d}.txt")
    save_series(nhem_path, dense_nhem_recon)
    
    # Save calibration RE statistics.
    step_path = VALIDATION_PATH.joinpath("steps")
    step_path.mkdir(exist_ok=True)
    re_path = step_path.joinpath(f"cal_re{step:04d}.csv")
    header = ["GLB", "NH", "DET", "NIN", "MLT"]
    data = [glob_cal_re, nhem_cal_re, detr_cal_re, nino_cal_re, mult_cal_re]
    with open(re_path, "w") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerow(data)
    
    # Save verification RE statistics.
    re_path = step_path.joinpath(f"ver_re{step:04d}.csv")
    header = ["GLB", "NH", "MLTA"]
    data = [glob_ver_re, nhem_ver_re, mult_ver_re]
    with open(re_path, "w") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerow(data)


def save_series(path, series):
    # Save pandas series.
    array = np.empty((series.size, 2))
    array[:, 0] = series.index
    array[:, 1] = series.to_numpy()
    save_array(path, array)


def save_array(path, array):
    np.savetxt(path, array, header="Year   Temperature", fmt="%04d %11.7f",
               comments="")


def regional_mean(df, region):
    lat = np.array(df.columns.get_level_values(0))
    lon = np.array(df.columns.get_level_values(1))
    if region == "glob":
        mask = np.full(lat.shape, True)
    elif region == "nhem":
        mask = lat > 0
    elif region == "nino":
        mask = (lat > -5) & (lat < 5) & (lon > -150) & (lon < -90)
    else:
        raise ValueError("Invalid region.")
    weights = np.cos(np.deg2rad(lat)) * mask
    weights /= array_sum32(weights)
    mu = df @ weights
    return mu


def re_statistic(recon, target, start_year, end_year):
    recon = recon.loc[start_year:end_year].to_numpy()
    target = target.loc[start_year:end_year].to_numpy()
    res = recon - target
    ssq_res = np.sum(res**2)
    ssq_target = np.sum(target**2)
    return 1 - ssq_res/ssq_target


def re_mult_statistic(recon, target, start_year, end_year):
    recon = recon.loc[start_year:end_year].to_numpy()
    target = target.loc[start_year:end_year].to_numpy()
    res = recon - target
    ssq_res = np.cumsum(res**2, axis=0).sum() # Incorrectly calculated.
    ssq_target = np.cumsum(target**2, axis=0).sum() # Ditto.
    return 1 - ssq_res/ssq_target


def concatenate_nhem_reconstructions():
    steps = sorted(reconstruction_steps())
    years = np.arange(steps[0], CAL_END + 1)
    recon = np.empty((years.size, 2))
    recon[:, 0] = years
    for step in steps:
        path = RECONSTRUCTION_PATH.joinpath("steps", f"nhem{step:04d}.txt")
        recon_step = np.genfromtxt(path, skip_header=1)
        index = recon_step[:, 0].astype(int) - years[0]
        recon[index, 1] = recon_step[:, 1]
    path = RECONSTRUCTION_PATH.joinpath("nhem_recon.txt")
    save_array(path, recon)


def make_re_table():
    table_path = VALIDATION_PATH.joinpath("validation.txt")
    with open(table_path, "w") as f:
        f.write("Reconstruction RE statistics:\n\n")
        f.write("       Calibration                    Verification\n")
        f.write("Step   GLB   NH    DET   NIN   MLT    GLB   NH    MLTA\n")
        f.write("------------------------------------------------------\n")
        for step in reversed(sorted(reconstruction_steps())):
            line = f"{step:04d}"
            # Calibration RE.
            path = VALIDATION_PATH.joinpath("steps", f"cal_re{step:04d}.csv")
            df = pd.read_csv(path)
            data = df.to_numpy()
            line = line + " "
            for value in data[0, :]:
                line = line + f"{value:6.2f}"
            # Verification RE.
            path = VALIDATION_PATH.joinpath("steps", f"ver_re{step:04d}.csv")
            df = pd.read_csv(path)
            data = df.to_numpy()
            line = line + " "
            for value in data[0, :]:
                line = line + f"{value:6.2f}"
            f.write(line + "\n")


def main():
    get_instrumental_data()
    get_proxy_data()
    prepare_instrumental_data()
    prepare_proxy_data()
    reconstruct_temperature()


if __name__ == "__main__":
    main()
