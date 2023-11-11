#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""MBH98 hockey stick emulation."""

import argparse
import pathlib
import shutil
import tarfile
import csv

import requests
import numpy as np
import pandas as pd

import svdalg

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

# Other constants.
PI = 3.14159265359
EOF_MAX = 40

# Fast emulation option.
parser = argparse.ArgumentParser()
parser.add_argument("--fast", action="store_true",
                    help="perform faster but less accurate emulation")
args = parser.parse_args()
FAST = args.fast


def svd(a):
    if FAST:
        u, s, vt = np.linalg.svd(a, full_matrices=False)
        v = vt.T
    else:
        u, s, v = svdalg.svd(a)
    return u, s, v


def lstsq(a, b):
    # Return the least-squares solution to a*x = b.
    if FAST:
        x = np.linalg.lstsq(a, b, rcond=None)[0]
    else:
        u, s, v = svdalg.svd(a)
        y = np.zeros((a.shape[1], b.shape[1]))
        for i in range(y.shape[0]):
            for j in range(y.shape[1]):
                for k in range(a.shape[0]):
                    y[i, j] += np.float32(u[k, i]) * b[k, j] / s[i]
        x = v.astype(np.float32) @ y
    return x


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
    process_instrumental_svd()


def read_instrumental_data():
    month_count = 12 * (INSTR_END - INSTR_START + 1)
    box_count = 72 * 36
    date = []
    temperature = np.full((month_count, box_count), np.nan)
    lines_per_month = box_count/18 + 1
    with open(DOWNLOAD_PATH.joinpath("anomalies-new"), "r") as f:
        month = -1
        for count, line in enumerate(f):
            if count % lines_per_month == 0:
                date.append(f"{int(line[5:10])}-{int(line[0:5])}")
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
    index = pd.to_datetime(date, format="%Y-%m")
    columns = pd.MultiIndex.from_arrays((lat, lon), names=("lat", "lon"))
    df = pd.DataFrame(data=temperature, index=index, columns=columns)
    df.to_pickle(INSTRUMENTAL_PATH.joinpath("anomalies-new.pkl"))


def dense_subset():
    df = pd.read_pickle(INSTRUMENTAL_PATH.joinpath("anomalies-new.pkl"))
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
    df = pd.read_pickle(INSTRUMENTAL_PATH.joinpath("anomalies-new.pkl"))
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
    df.to_pickle(INSTRUMENTAL_PATH.joinpath("dense_subset_monthly.pkl"))


def prepare_sparse_subset():
    df = sparse_subset()
    convert_to_degrees(df)
    fill_in_instrumental_data(df)
    df = annual_means(df)
    round_to_2_decimal_places(df)
    df.to_pickle(INSTRUMENTAL_PATH.joinpath("sparse_subset.pkl"))


def fill_in_instrumental_data(df):
    t = df.to_numpy()
    n = t.shape[0]
    for i in range(t.shape[1]):
        # Indices with data.
        j = (~np.isnan(t[:, i])).nonzero()[0]
        
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
            t[-1, i] = div32(add32(t[-3, i], t[-2, i]), 2)


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
    for element in a:
        s = add32(s, element)
    return s


def column_mean32(a):
    s = np.zeros(a.shape[1], dtype=np.float32)
    for row in a:
        s = add32(s, row)
    mu = div32(s, a.shape[0])
    return mu


def linear_fit32(y, x0=0):
    # Fit lines to the columns of y.
    n, m = y.shape
    sum_x = np.zeros(m, dtype=np.float32)
    sum_xx = np.zeros(m, dtype=np.float32)
    sum_y = np.zeros(m, dtype=np.float32)
    sum_xy = np.zeros(m, dtype=np.float32)
    for i in range(n):
        x = x0 + i + 1
        sum_x += np.float32(x)
        sum_xx += np.float32(x**2)
        sum_y += y[i, :]
        sum_xy += y[i, :] * np.float32(x)
    n = np.float32(n)
    alpha = (sum_xx*sum_y - sum_x*sum_xy) / (n*sum_xx - sum_x**2)
    beta = (n*sum_xy - sum_x*sum_y) / (n*sum_xx - sum_x**2)
    return alpha, beta


def convert_to_degrees(df):
    t = df.to_numpy()
    t[t == -9999] = np.nan
    t /= 100


def invert_anomalies_below_minus_10(df):
    t = df.to_numpy()
    t[t <= -10] *= -1


def round_to_2_decimal_places(df):
    t = df.to_numpy()
    t[:] = t.round(2)


def annual_means(df):
    years = sorted(set(df.index.year))
    df_annual = pd.DataFrame(index=years, columns=df.columns,
                             dtype=np.float64)
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
    mu = np.zeros(t.shape[1], dtype=np.float32)
    for i in range(month_count):
        mu += t[i, :]
    mu = (mu / month_count).astype(np.float64)
    t -= mu
    
    # Detrended standard deviations.
    month_count = 12 * (INSTR_END - CAL_START + 1)
    start_month = 12 * (CAL_START - INSTR_START)
    sigma = np.zeros(t.shape[1])
    alpha, beta = linear_fit32(t, start_month)
    for i in range(t.shape[0]):
        month = start_month + i + 1
        detrended = t[i, :] - beta*np.float32(month) - alpha
        sigma += detrended**2
    sigma = np.sqrt(sigma / month_count)
    
    # Save statistics.
    np.save(INSTRUMENTAL_PATH.joinpath("mu.npy"), mu)
    np.save(INSTRUMENTAL_PATH.joinpath("sigma.npy"), sigma)


def center_annual_means():
    # Load instrumental data.
    dense_path = INSTRUMENTAL_PATH.joinpath("dense_subset_monthly.pkl")
    sparse_path = INSTRUMENTAL_PATH.joinpath("sparse_subset.pkl")
    df_dense_monthly = pd.read_pickle(dense_path)
    df_sparse = pd.read_pickle(sparse_path)
    
    # Load means for dense subset.
    mu = np.load(INSTRUMENTAL_PATH.joinpath("mu.npy"))
    
    # Center instrumental data.
    mu = pd.Series(data=mu, index=df_dense_monthly.columns)
    df_dense_monthly -= mu
    df_sparse -= mu[df_sparse.columns] # Not properly centered.
    
    # Annual means for dense subset.
    df_dense = df_dense_monthly.groupby(df_dense_monthly.index.year).mean()
    
    # Save centered data.
    dense_path = INSTRUMENTAL_PATH.joinpath("dense_subset_centered.pkl")
    sparse_path = INSTRUMENTAL_PATH.joinpath("sparse_subset_centered.pkl")
    df_dense.to_pickle(dense_path)
    df_sparse.to_pickle(sparse_path)


def compute_instrumental_svd():
    # Load dense subset of instrumental data.
    path = INSTRUMENTAL_PATH.joinpath("dense_subset_monthly.pkl")
    df_monthly = pd.read_pickle(path)
    
    # Load statistics for standardization.
    mu = np.load(INSTRUMENTAL_PATH.joinpath("mu.npy"))
    sigma = np.load(INSTRUMENTAL_PATH.joinpath("sigma.npy"))
    
    # Standardize data.
    t = df_monthly.to_numpy()
    t_zscores = (t - mu) / sigma
    
    # Compute singular value decomposition of standardized and area
    # weighted instrumental data.
    lat = np.array(df_monthly.columns.get_level_values(0))
    w = np.cos(lat*PI/180)
    t_weighted = t_zscores * w
    u_monthly, s, v = svd(t_weighted)
    
    # Make svd dataframes.
    u_monthly = pd.DataFrame(data=u_monthly, index=df_monthly.index)
    v = pd.DataFrame(data=v, index=df_monthly.columns)
    
    # Save svd.
    svd_path = INSTRUMENTAL_PATH.joinpath("svd")
    svd_path.mkdir(exist_ok=True)
    u_monthly.to_pickle(svd_path.joinpath("u_monthly.pkl"))
    v.to_pickle(svd_path.joinpath("v.pkl"))
    np.save(svd_path.joinpath("s.npy"), s)


def process_instrumental_svd():
    # Load full precision svd.
    svd_path = INSTRUMENTAL_PATH.joinpath("svd")
    u_monthly = pd.read_pickle(svd_path.joinpath("u_monthly.pkl"))
    v = pd.read_pickle(svd_path.joinpath("v.pkl"))
    s = np.load(svd_path.joinpath("s.npy"))
    
    # Save svd rounded to 6 and 14 significant figures.
    u_path = svd_path.joinpath("u_monthly_rounded.txt")
    v_path = svd_path.joinpath("v_rounded.txt")
    s_path = svd_path.joinpath("s.txt")
    np.savetxt(u_path, u_monthly.to_numpy(dtype=np.float32), fmt="%.5e")
    np.savetxt(v_path, v.to_numpy(dtype=np.float32), fmt="%.5e")
    np.savetxt(s_path, s, fmt="%.13e")
    
    # Load rounded singular vectors and save as binary.
    v_rounded = np.genfromtxt(v_path)
    v_rounded = pd.DataFrame(data=v_rounded, index=v.index)
    v_rounded.to_pickle(svd_path.joinpath("v_rounded.pkl"))
    u_monthly_rounded = np.genfromtxt(u_path)
    u_monthly_rounded = pd.DataFrame(data=u_monthly_rounded,
                                     index=u_monthly.index)
    
    # Standardize annual means of left singular vectors.
    u = annual_means(u_monthly_rounded)
    data = u.to_numpy()
    mu = np.zeros(data.shape[1])
    cal_length = np.float32(CAL_END - CAL_START + 1)
    for i in range(int(cal_length)):
        mu += data[i, :]
    mu /= cal_length
    data -= mu
    
    # Standardize by detrended standard deviations.
    sigma = np.zeros(data.shape[1])
    alpha, beta = linear_fit32(data[:int(cal_length), :])
    for i in range(int(cal_length)):
        detrended = data[i, :] - beta*np.float32(i+1) - alpha
        sigma += detrended**2
    sigma = np.sqrt(sigma / cal_length)
    data /= sigma
    
    # Save standardized singular vectors and statistics.
    u.to_pickle(svd_path.joinpath("u_zscores.pkl"))
    np.save(svd_path.joinpath("mu_u.npy"), mu)
    np.save(svd_path.joinpath("sigma_u.npy"), sigma)
    
    # Undo detrending and save again rounded to 14 significant figures.
    data *= sigma
    data += mu
    np.savetxt(svd_path.joinpath("u.txt"), data, fmt="%.13e")
    
    # Weights.
    s = np.genfromtxt(s_path)
    ssum = array_sum32(s[:EOF_MAX]) / np.float32(EOF_MAX)
    w = s / ssum
    np.save(svd_path.joinpath("weights.npy"), w)


def prepare_proxy_data():
    print("Preparing proxy data...")
    extract_proxy_archive()
    create_proxy_matrices()


def extract_proxy_archive():
    tar_path = DOWNLOAD_PATH.joinpath("mbh98.tar")
    untar_path = PROXY_PATH.joinpath("mbh98")
    untar(tar_path, untar_path)


def untar(tar_path, untar_path):
    if not untar_path.is_dir():
        with tarfile.TarFile(tar_path, "r") as f:
            f.extractall(untar_path)


def create_proxy_matrices():
    datalists_path = CONFIG_PATH.joinpath("proxy")
    data_path = PROXY_PATH.joinpath("networks")
    data_path.mkdir(exist_ok=True)
    for step in reconstruction_steps():
        p = proxy_matrix(datalists_path.joinpath(f"datalist{step}.dat"))
        save_proxy_matrix(p, data_path.joinpath(f"data{step}.dat"))


def proxy_matrix(datalist_path):
    # Return a proxy data matrix with year in the first column.
    mbh98_path = PROXY_PATH.joinpath("mbh98")
    with open(datalist_path, "r") as f:
        relative_paths = f.read().splitlines()
    step = int(datalist_path.stem.strip("datalist"))
    n = CAL_END - step + 1
    m = len(relative_paths) + 1
    p = np.full((n, m), np.nan)
    p[:, 0] = np.arange(step, CAL_END + 1)
    for i, relative_path in enumerate(relative_paths):
        data_path = mbh98_path.joinpath(*relative_path.split("/"))
        data = read_proxy_data(data_path)
        data = submatrix(data, step, CAL_END)
        index = data[:, 0].astype(int) - step
        p[index, i+1] = data[:, 1]
    return p


def read_proxy_data(path):
    p = pd.read_table(path, header=None)
    p = p[0].str.split(expand=True).astype(float).to_numpy()
    return p


def save_proxy_matrix(p, path):
    m = p.shape[1]
    np.savetxt(path, p, fmt="%4d" + "%14.6e" * (m-1))


def fill_in_proxy_matrix(p):
    for i in range(p.shape[1]):
        if np.isnan(p[-1, i]):
            index = (~np.isnan(p[:, i])).nonzero()[0]
            p[index[-1]:, i] = p[index[-1], i]


def standardize_proxy_matrix(p):
    # Standardize proxy matrix p over the calibration period.
    p_cal = submatrix(p, CAL_START, CAL_END).copy()
    mu = np.zeros(p.shape[1]-1)
    cal_length = np.float32(CAL_END - CAL_START + 1)
    for i in range(int(cal_length)):
        mu += p_cal[i, 1:]
    mu /= cal_length
    p_cal[:, 1:] -= mu
    p[:, 1:] -= mu
    
    # Standardize by standard deviations.
    sigma = np.zeros(p.shape[1]-1)
    for i in range(int(cal_length)):
        sigma += p_cal[i, 1:]**2
    sigma = np.sqrt(sigma / cal_length)
    p_cal[:, 1:] /= sigma
    p[:, 1:] /= sigma
    
    # Standardize again by detrended standard deviations.
    sigma = np.zeros(p.shape[1]-1)
    alpha, beta = linear_fit32(p_cal[:, 1:])
    for i in range(int(cal_length)):
        detrended = p_cal[i, 1:] - beta*np.float32(i+1) - alpha
        sigma += np.square(detrended, dtype=np.float32)
    sigma = np.sqrt(sigma / cal_length)
    p[:, 1:] /= sigma


def submatrix(x, t0, t1):
    # Return a submatrix of x where the first column is in [t0, t1].
    return x[(x[:, 0] >= t0) & (x[:, 0] <= t1), :]


def reconstruct_temperature():
    print("Generating reconstruction...")
    for step in reconstruction_steps():
        recon = reconstructed_temperature_field(step)
        analyze_regions(recon)


def reconstruction_steps():
    with open(CONFIG_PATH.joinpath("steps.csv"), "r", newline="") as f:
        reader = csv.reader(f)
        steps = next(reader)
    steps = [int(n) for n in steps]
    return steps


def reconstructed_temperature_field(step):
    u_recon, s, v = reconstructed_svd(step)
    
    # Undo weighting and scaling.
    sigma = np.load(INSTRUMENTAL_PATH.joinpath("sigma.npy"))
    lat = np.array(v.index.get_level_values(0))
    w = np.cos(lat*PI/180)
    recon_weighted = (u_recon.to_numpy() * s) @ v.to_numpy().T
    recon_zscores = recon_weighted / w
    recon = recon_zscores * sigma
    recon = pd.DataFrame(data=recon, index=u_recon.index, columns=v.index)
    return recon


def reconstructed_svd(step):
    # Load proxy matrix.
    p = np.genfromtxt(PROXY_PATH.joinpath("networks", f"data{step}.dat"))
    p = p.astype(np.float32).astype(np.float64)
    fill_in_proxy_matrix(p)
    standardize_proxy_matrix(p)
    p_cal = submatrix(p, CAL_START, CAL_END)
    index_cal = (p[:, 0] >= CAL_START) & (p[:, 0] <= CAL_END)
    
    # Load instrumental PC selection.
    eofs = pc_selection(step)
    
    # Load instrumental svd.
    svd_path = INSTRUMENTAL_PATH.joinpath("svd")
    u_z = pd.read_pickle(svd_path.joinpath("u_zscores.pkl"))
    v = pd.read_pickle(svd_path.joinpath("v_rounded.pkl"))
    s = np.genfromtxt(svd_path.joinpath("s.txt"))
    mu = np.load(svd_path.joinpath("mu_u.npy"))
    sigma = np.load(svd_path.joinpath("sigma_u.npy"))
    w = np.load(svd_path.joinpath("weights.npy"))
    
    # Calibration period.
    u_z_cal = u_z.loc[CAL_START:CAL_END, eofs]
    u_z_cal = u_z_cal.to_numpy()
    
    # Calibrate proxy data against instrumental PCs.
    g = lstsq(w[eofs] * u_z_cal, p_cal[:, 1:])
    
    # Reconstruct PCs.
    u_z_recon = lstsq(g.T, p[:, 1:].T).T
    u_recon = u_z_recon * sigma[eofs] / w[eofs] + mu[eofs]
    
    # Rescale PCs.
    sigma_recon = np.zeros(len(eofs), dtype=np.float32)
    sigma_instr = np.zeros(len(eofs), dtype=np.float32)
    u_recon_cal = u_recon[index_cal, :]
    u_cal = u_z_cal * sigma[eofs] + mu[eofs]
    cal_length = np.float32(CAL_END - CAL_START + 1)
    for i in range(int(cal_length)):
        sigma_recon += u_recon_cal[i, :]**2
        sigma_instr += u_cal[i, :]**2
    sigma_recon = np.sqrt(sigma_recon / cal_length)
    sigma_instr = np.sqrt(sigma_instr / cal_length)
    u_recon = u_recon * sigma_instr / sigma_recon
    
    # Save reconstructed singular vectors.
    step_path = RECONSTRUCTION_PATH.joinpath("steps")
    step_path.mkdir(exist_ok=True)
    np.save(step_path.joinpath(f"rpc{step}.npy"), u_recon)
    
    # Create PC dataframe.
    u_recon = pd.DataFrame(data=u_recon, index=p[:, 0].astype(int))
    return u_recon, s[eofs], v.iloc[:, eofs]


def pc_selection(step):
    # Load hard-coded PC selection.
    path = CONFIG_PATH.joinpath("instrumental", f"eoflist{step}.csv")
    with open(path, "r", newline="") as f:
        reader = csv.reader(f)
        eofs = next(reader)
    eofs = [int(n) - 1 for n in eofs]
    return eofs


def analyze_regions(recon):
    # Load instrumental data.
    dense_path = INSTRUMENTAL_PATH.joinpath("dense_subset_centered.pkl")
    sparse_path = INSTRUMENTAL_PATH.joinpath("sparse_subset_centered.pkl")
    dense_instr = pd.read_pickle(dense_path)
    sparse_instr = pd.read_pickle(sparse_path)
    
    # Instrumental means.
    dense_glob_instr = regional_mean(dense_instr, "glob")
    dense_nhem_instr = regional_mean(dense_instr, "nhem")
    dense_detr_instr = detrend_series(dense_nhem_instr.loc[CAL_START:CAL_END])
    dense_nino_instr = regional_mean(dense_instr, "nino")
    sparse_glob_instr = regional_mean(sparse_instr, "glob")
    sparse_nhem_instr = regional_mean(sparse_instr, "nhem")
    
    # Reconstructed means.
    dense_recon = recon
    dense_glob_recon = regional_mean(dense_recon, "glob")
    dense_nhem_recon = regional_mean(dense_recon, "nhem")
    dense_detr_recon = detrend_series(dense_nhem_recon.loc[CAL_START:CAL_END])
    dense_nino_recon = regional_mean(dense_recon, "nino")
    sparse_recon = recon.loc[:, sparse_instr.columns]
    sparse_glob_recon = regional_mean(sparse_recon, "glob")
    sparse_nhem_recon = regional_mean(sparse_recon, "nhem")
    
    # Reduction of Error (RE) statistics.
    glob_cal_re = re_statistic(dense_glob_recon, dense_glob_instr, "cal")
    nhem_cal_re = re_statistic(dense_nhem_recon, dense_nhem_instr, "cal")
    detr_cal_re = re_statistic(dense_detr_recon, dense_detr_instr, "cal")
    nino_cal_re = re_statistic(dense_nino_recon, dense_nino_instr, "cal")
    mult_cal_re = mult_re_statistic(dense_recon, dense_instr, "cal")
    glob_ver_re = re_statistic(sparse_glob_recon, sparse_glob_instr, "ver")
    nhem_ver_re = re_statistic(sparse_nhem_recon, sparse_nhem_instr, "ver")
    mult_ver_re = mult_re_statistic(sparse_recon, sparse_instr, "ver")
    
    # Save Northern Hemisphere reconstruction.
    step = recon.index[0]
    step_path = RECONSTRUCTION_PATH.joinpath("steps")
    step_path.mkdir(exist_ok=True)
    np.savetxt(step_path.joinpath(f"nhem{step}.txt"),
               series_to_array(dense_nhem_recon),
               header="Year   Temperature", fmt="%4d %11.7f", comments="")
    
    # Save calibration RE statistics.
    step_path = VALIDATION_PATH.joinpath("steps")
    step_path.mkdir(exist_ok=True)
    header = ["GLB", "NH", "DET", "NIN", "MLT"]
    data = [glob_cal_re, nhem_cal_re, detr_cal_re, nino_cal_re, mult_cal_re]
    with open(step_path.joinpath(f"cal_re{step}.csv"), "w") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerow(data)
    
    # Save verification RE statistics.
    header = ["GLB", "NH", "MLTA"]
    data = [glob_ver_re, nhem_ver_re, mult_ver_re]
    with open(step_path.joinpath(f"ver_re{step}.csv"), "w") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerow(data)


def series_to_array(s):
    a = np.empty((s.size, 2))
    a[:, 0] = s.index
    a[:, 1] = s.to_numpy()
    return a


def detrend_series(s):
    # Return a detrended pandas series.
    data = s.to_numpy().copy()
    alpha, beta = linear_fit32(data[:, np.newaxis])
    data = data - beta*np.arange(1, data.size+1) - alpha
    return pd.Series(data=data, index=s.index)


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
    weights = np.cos(lat*PI/180) * mask
    mu = (df @ weights) / array_sum32(weights)
    return mu


def re_statistic(recon, target, period):
    if period == "cal":
        start_year, end_year = CAL_START, CAL_END
    elif period == "ver":
        start_year, end_year = VER_START, VER_END
    else:
        raise ValueError("Invalid period.")
    recon = recon.loc[start_year:end_year].to_numpy()
    target = target.loc[start_year:end_year].to_numpy()
    res = recon - target
    ssq_res = np.sum(res**2)
    ssq_target = np.sum(target**2)
    return 1 - ssq_res/ssq_target


def mult_re_statistic(recon, target, period):
    if period == "cal":
        start_year, end_year = CAL_START, CAL_END
    elif period == "ver":
        start_year, end_year = VER_START, VER_END
    else:
        raise ValueError("Invalid period.")
    recon = recon.loc[start_year:end_year].to_numpy()
    target = target.loc[start_year:end_year].to_numpy()
    res = recon - target
    ssq_res = np.cumsum(res**2, axis=0).sum() # Incorrectly calculated.
    ssq_target = np.cumsum(target**2, axis=0).sum() # Ditto.
    return 1 - ssq_res/ssq_target


def summarize_results():
    print("Summarizing results...")
    splice_nhem_reconstructions()
    splice_pc_reconstructions()
    make_re_table()


def splice_nhem_reconstructions():
    steps = sorted(reconstruction_steps())
    years = np.arange(steps[0], CAL_END + 1)
    recon = np.empty((years.size, 2))
    recon[:, 0] = years
    for step in steps:
        path = RECONSTRUCTION_PATH.joinpath("steps", f"nhem{step}.txt")
        recon_step = np.genfromtxt(path, skip_header=1)
        index = years >= step
        recon[index, 1] = recon_step[:, 1]
    path = RECONSTRUCTION_PATH.joinpath("nhem_recon.txt")
    np.savetxt(path, recon, header="Year   Temperature", fmt="%4d %11.7f",
               comments="")


def splice_pc_reconstructions():
    steps = sorted(reconstruction_steps())
    years = np.arange(steps[0], CAL_END + 1)
    for i in range(5):
        recon = np.full((years.size, 2), np.nan)
        recon[:, 0] = years
        for step in steps:
            eofs = pc_selection(step)
            if i in eofs:
                path = RECONSTRUCTION_PATH.joinpath("steps", f"rpc{step}.npy")
                recon_step = np.load(path)
                index = years >= step
                recon[index, 1] = recon_step[:, eofs.index(i)]
        recon = recon[~np.isnan(recon[:, 1]), :]
        path = RECONSTRUCTION_PATH.joinpath(f"rpc{i+1:02d}.txt")
        np.savetxt(path, recon, header=f"Year   PC#{i+1}", fmt="%4d %12.8f",
                   comments="")


def make_re_table():
    with open(VALIDATION_PATH.joinpath("validation.txt"), "w") as f:
        f.write("Reconstruction RE statistics:\n\n")
        f.write("       Calibration                    Verification\n")
        f.write("Step   GLB   NH    DET   NIN   MLT    GLB   NH    MLTA\n")
        f.write("------------------------------------------------------\n")
        for step in reversed(sorted(reconstruction_steps())):
            line = f"{step:4d}"
            # Calibration RE.
            path = VALIDATION_PATH.joinpath("steps", f"cal_re{step}.csv")
            df = pd.read_csv(path)
            data = df.to_numpy()
            line = line + " "
            for value in data[0, :]:
                line = line + f"{value:6.2f}"
            # Verification RE.
            path = VALIDATION_PATH.joinpath("steps", f"ver_re{step}.csv")
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
    summarize_results()


if __name__ == "__main__":
    main()
