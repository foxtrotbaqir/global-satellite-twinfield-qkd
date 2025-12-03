#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TF-QKD Satellite Constellation Simulation (GPU Accelerated, SGP4)


import sys, time, socket
import numpy as np
import cupy as cp
import pandas as pd
import matplotlib.pyplot as plt
from astropy.time import Time
from astropy import units as u
from sgp4.api import Satrec, jday
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra
import numpy as _np
from functools import lru_cache as _lru

# ---------------------------- switches --------------------------------
DEBUG = 0
DIAG  = 0

# ----------------------------- TLE ------------------------------------
MICIUS_TLE = [
    "1 41732U 16051A   20175.61167824  .00000206  00000-0  11664-4 0  9993",
    "2 41732  97.3688 229.5573 0011311 112.8982  53.6457 15.12909341207844"
]

# --------------------------- helpers ----------------------------------
def myLla2ecef(lat_deg, lon_deg, alt_m=0.0):
    a = 6378137.0
    f = 1/298.257223563
    e2 = f*(2-f)
    lat = np.deg2rad(lat_deg); lon = np.deg2rad(lon_deg)
    N = a / np.sqrt(1 - e2*np.sin(lat)**2)
    x = (N+alt_m) * np.cos(lat) * np.cos(lon)
    y = (N+alt_m) * np.cos(lat) * np.sin(lon)
    z = ((1-e2)*N+alt_m) * np.sin(lat)
    return np.array([x,y,z], dtype=np.float64)

def backgroundNoise_series(times_astro, lat, lon):
    out = np.empty(len(times_astro), dtype=np.float32)
    for i, t in enumerate(times_astro):
        jd = t.jd
        d = jd - 2451545.0
        g = np.deg2rad((357.529 + 0.98560028*d) % 360.0)
        q = np.deg2rad((280.459 + 0.98564736*d) % 360.0)
        L = q + np.deg2rad(1.915)*np.sin(g) + np.deg2rad(0.020)*np.sin(2*g)
        eps = np.deg2rad(23.439 - 0.00000036*d)
        RA = np.arctan2(np.cos(eps)*np.sin(L), np.cos(L))
        dec = np.arcsin(np.sin(eps)*np.sin(L))
        T = (jd - np.floor(jd)) * 24.0
        theta = np.deg2rad((280.16 + 360.9856235*(d+T/24.0)) % 360.0)
        H = theta + np.deg2rad(lon) - RA
        latr = np.deg2rad(lat)
        el = np.arcsin(np.sin(latr)*np.sin(dec) + np.cos(latr)*np.cos(dec)*np.cos(H))
        if el > 0:
            out[i] = 800.0
        else:
            day = t.to_datetime(timezone=None).day
            phase = (day % 29) / 29.0
            out[i] = 100.0 if phase > 0.5 else 10.0
    return out

# --- CPU scalar link efficiencies ------
def _eta_link_ul_cpu(L):
    lam=1.55e-6; Dtx=1.0; Drx=0.3; alpha=3e-3; sigma=3e-7
    wL = (lam/(np.pi*Dtx))*L
    return (Drx/(2*wL))**2 * np.exp(-alpha*(L/1e3)) * np.exp(-2*(sigma*L)**2 / wL**2) * 0.8*0.3

def _eta_link_dl_cpu(L):
    lam=1.55e-6; Dtx=0.3; Drx=0.3; alpha=3e-3; sigma=1.5e-7
    wL = (lam/(np.pi*Dtx))*L
    return (Drx/(2*wL))**2 * np.exp(-alpha*(L/1e3)) * np.exp(-2*(sigma*L)**2 / wL**2) * 0.8*0.3

def _eta_link_isl_cpu(L):
    lam=1.55e-6; Dtx=0.3; Drx=0.3; sigma=1e-7
    wL = (lam/(np.pi*Dtx))*L
    return (Drx/(2*wL))**2 * np.exp(-2*(sigma*L)**2 / wL**2) * 0.8*0.8

# --------------------------- SKR calculation --------------------------
def _tfqkd_skr_asym_cpu(etaA, etaB, bg_photons):
    Clock_rate=_np.float32(1e9); d=_np.float32(0.5); M=_np.float32(16.0); eta_sync=_np.float32(0.95)
    f=_np.float32(1.15); e_opt=_np.float32(0.03); E_M=_np.float32(0.01275)
    P_dc = _np.float32(10.0/Clock_rate)
    P_bg = _np.float32(bg_photons/Clock_rate)
    P_tot = _np.float32(P_dc + P_bg)

    etaA=_np.float32(etaA); etaB=_np.float32(etaB)
    muB=_np.float32(0.5); nuB=_np.float32(0.1); wB=_np.float32(1e-4)
    muA=_np.float32(muB*(etaB/etaA)); nuA=_np.float32(nuB*(etaB/etaA)); wA=wB

    exp=_np.exp
    Qmu=_np.float32(1-(1-P_tot)**2*exp(-(muA*etaA+muB*etaB)))
    Qnu=_np.float32(1-(1-P_tot)**2*exp(-(nuA*etaA+nuB*etaB)))
    Qw =_np.float32(1-(1-P_tot)**2*exp(-(wA*etaA+wB*etaB)))

    Emu=_np.float32(0.5+(1/(2*Qmu))*(1-P_tot)*(
         exp(-(muA*etaA+muB*etaB)*(1-(e_opt+E_M)))
        -exp(-(muA*etaA+muB*etaB)*(e_opt+E_M))))

    y0=_np.float32((nuB*Qw*exp(wB)-wB*Qnu*exp(nuB))/(nuB-wB))
    y1=_np.float32(((muB**2)*Qnu*exp(nuB)-(muB**2)*Qw*exp(wB)
        -(nuB**2-wB**2)*(Qmu*exp(muB)-y0)) /(muB*(muB*nuB-muB*wB-nuB**2+wB**2)))

    e1=_np.float32((Emu*Qmu*exp(muB)-_np.float32(0.5)*y0)/(y1*muB))
    Q1=_np.float32(exp(-muB)*muB*y1)

    def h(x):
        x=_np.clip(x, _np.float32(1e-12), _np.float32(1-1e-12))
        return -x*_np.log2(x)-(1-x)*_np.log2(1-x)

    N=_np.float32(1e9)
    delta=_np.float32(5*_np.sqrt(_np.log(2/_np.float32(1e-10))/N))
    R_QKD=_np.float32(Q1*(1-h(e1-delta))-f*Qmu*h(Emu+delta))
    return float(_np.maximum(_np.float32(0), _np.float32(eta_sync)*(d*Clock_rate/M)*R_QKD))

@_lru(maxsize=4096)
def _skr_cached(q_etaA, q_etaB, q_bg):
    etaA = _np.float32(q_etaA)*_np.float32(1e-8)
    etaB = _np.float32(q_etaB)*_np.float32(1e-8)
    bg   = _np.float32(q_bg  )*_np.float32(1.0)
    return _tfqkd_skr_asym_cpu(etaA, etaB, bg)

def tfqkd_skr_asym_cpu_cached(etaA, etaB, bg_photons):
    q_etaA = int(_np.round(float(etaA)/1e-8))
    q_etaB = int(_np.round(float(etaB)/1e-8))
    q_bg   = int(_np.round(float(bg_photons)))
    return _skr_cached(q_etaA, q_etaB, q_bg)

# ----------------------- constellation build --------------------------
def build_walker_sso(base_tle1, base_tle2, n_orbits, sats_per_orbit):
    sats=[]
    for i in range(n_orbits):
        for j in range(sats_per_orbit):
            f=base_tle2.split()
            raan0=float(f[3]); ma0=float(f[6])
            f[3]=f"{(raan0+i*(360/n_orbits))%360.0:8.4f}"
            f[6]=f"{(ma0  +j*(360/sats_per_orbit))%360.0:8.4f}"
            sats.append(Satrec.twoline2rv(base_tle1," ".join(f)))
    return sats

# ------------- ISL candidate pairs at t0 (<= 1500 km) -----------------
def isl_candidates(sat_ecef_t0):
    V = sat_ecef_t0.astype(np.float32)          # (3,N)
    N = V.shape[1]
    Vg = cp.asarray(V)
    mask_i=[]; mask_j=[]
    chunk = 1024
    for a in range(0,N,chunk):
        va = Vg[:, a:a+chunk]                    # (3,ca)
        diff = va[:,:,None] - Vg[:,None,:]
        dij = cp.linalg.norm(diff, axis=0)       # (ca,N)
        for ii in range(dij.shape[0]):
            i_global = a+ii
            sel = (dij[ii] < 1500e3)
            sel[:i_global+1] = False
            js = cp.where(sel)[0]
            if js.size:
                mask_i.append(np.full(js.size, i_global, dtype=np.int32))
                mask_j.append(js.get().astype(np.int32))
    if mask_i:
        I = np.concatenate(mask_i); J = np.concatenate(mask_j)
    else:
        I = np.empty(0, dtype=np.int32); J = np.empty(0, dtype=np.int32)
    return I, J

# ---------------------- path reconstruction ---------------------------
def reconstruct_path(parents, src, dst):
    path=[]; u=dst
    while u != -9999 and u != src:
        path.append(u); u = parents[u]
    if u == -9999: return None
    path.append(src); path.reverse()
    return path

# ---------- bottleneck hop selection (distance-only on GPU) -----------
def path_min_skr_bottleneck_fast(path_ids, nodes_dev, bg, N):
    if len(path_ids) < 3:
        return 0.0
    p = np.asarray(path_ids, dtype=np.int32)
    A_ids = cp.asarray(p[:-2], dtype=cp.int32)
    B_ids = cp.asarray(p[1:-1], dtype=cp.int32)
    C_ids = cp.asarray(p[2:  ], dtype=cp.int32)

    A = nodes_dev[A_ids]; B = nodes_dev[B_ids]; C = nodes_dev[C_ids]
    d1 = B - A; d2 = B - C
    d1sq = cp.sum(d1*d1, axis=1); d2sq = cp.sum(d2*d2, axis=1)
    idx = int(cp.argmax(cp.maximum(d1sq, d2sq)).get())

    L1 = float(cp.sqrt(d1sq[idx]).get())
    L2 = float(cp.sqrt(d2sq[idx]).get())

    a_id = int(A_ids[idx].get()); b_id = int(B_ids[idx].get()); c_id = int(C_ids[idx].get())
    a_is_gs = (a_id == 0) or (a_id == N+1)
    b_is_gs = (b_id == 0) or (b_id == N+1)
    c_is_gs = (c_id == 0) or (c_id == N+1)

    # compute efficiencies
    eta1 = (_eta_link_ul_cpu(L1) if a_is_gs and not b_is_gs else
            _eta_link_dl_cpu(L1) if (not a_is_gs) and b_is_gs else
            _eta_link_isl_cpu(L1))
    eta2 = (_eta_link_ul_cpu(L2) if c_is_gs and not b_is_gs else
            _eta_link_dl_cpu(L2) if (not c_is_gs) and b_is_gs else
            _eta_link_isl_cpu(L2))

    # Call SKR + expose Qmu, Emu
    Clock_rate=1e9; P_dc=10/Clock_rate; P_bg=bg/Clock_rate; P_tot=P_dc+P_bg
    muB=0.5; nuB=0.1; wB=1e-4
    muA=muB*(eta2/eta1); nuA=nuB*(eta2/eta1); wA=wB

    exp=np.exp
    Qmu=1-(1-P_tot)**2*exp(-(muA*eta1+muB*eta2))
    Emu=0.5+(1/(2*Qmu))*(1-P_tot)*(exp(-(muA*eta1+muB*eta2)*(1-(0.03+0.01275))) -
                                   exp(-(muA*eta1+muB*eta2)*(0.03+0.01275)))

    skr = tfqkd_skr_asym_cpu_cached(np.float32(eta1), np.float32(eta2), float(bg))
    return skr, eta1, eta2, Qmu, Emu
# ------------------------------- main ---------------------------------
def main():
    # ===== Config =====
    n_orbits=25; sats_per_orbit=44
    UL_MAX = 1200e3; ISL_MAX = 1500e3
    UL_MAX2 = UL_MAX*UL_MAX; ISL_MAX2 = ISL_MAX*ISL_MAX
    stepSec=1; sim_factor = 30

    # ===== Ground stations =====
    GS1_NAME, GS2_NAME = "Munich", "NewYork"
    gs1_lat,gs1_lon = 48.1351, 11.5820
    gs2_lat,gs2_lon = 40.7128,-74.0060
    gs1_ecef = myLla2ecef(gs1_lat,gs1_lon,0.0).astype(np.float32)
    gs2_ecef = myLla2ecef(gs2_lat,gs2_lon,0.0).astype(np.float32)

    # ===== Time grid =====
    startTime=Time("2025-07-31 12:00:00",scale="utc")
    stopTime = startTime + sim_factor*u.day
    times = pd.date_range(start=startTime.to_datetime(timezone=None),
                          end=stopTime.to_datetime(timezone=None),
                          freq=f"{stepSec}s", tz="UTC")
    T = len(times)

    # ===== Constellation =====
    sats = build_walker_sso(MICIUS_TLE[0], MICIUS_TLE[1], n_orbits, sats_per_orbit)

    # ===== Pre-propagate all positions (SGP4 once) =====
    pos = np.empty((T, 3, len(sats)), dtype=np.float32)
    for si, sat in enumerate(sats):
        xyz = np.empty((T,3), dtype=np.float32)
        for tii, ts in enumerate(times):
            jd, fr = jday(ts.year, ts.month, ts.day, ts.hour, ts.minute,
                          ts.second + ts.microsecond*1e-6)
            e, r, v = sat.sgp4(jd, fr)
            xyz[tii,:] = np.nan if e!=0 else (np.array(r, dtype=np.float64) * 1000.0).astype(np.float32)
        pos[:, :, si] = xyz
    valid_mask = ~np.isnan(pos).any(axis=(0,1))
    pos = pos[:, :, valid_mask]
    N = pos.shape[2]

    # ===== ISL candidates at t0 =====
    I_isl, J_isl = isl_candidates(pos[0])
    M = len(I_isl)

    # ===== Background noise =====
    times_astropy = [Time(ts.to_pydatetime(), scale="utc") for ts in times]
    bg_series = backgroundNoise_series(times_astropy, gs1_lat, gs1_lon)

    # ===== CSV header (enriched) =====
    CSV_FILE="tfqkd_paths_opt_dataset.csv"
    with open(CSV_FILE,"w") as f:
        f.write(",".join([
            "time",
            # Path 1
            "path1","skr1_bps","etaA1","etaB1","Qmu1","Emu1",
            # Path 2
            "path2","skr2_bps","etaA2","etaB2","Qmu2","Emu2",
            # Path 3
            "path3","skr3_bps","etaA3","etaB3","Qmu3","Emu3",
            # Aggregates
            "skr_xor_bps",
            "skr_sum_bps",
            "skr_min_bps",
            "k_paths",
            "connected",
            "handover",
            "window_id"
        ]) + "\n")


    # ===== Node indexing =====
    SRC = 0; DST = N+1; NN = N + 2

    # ===== Device preallocs =====
    gs1_g = cp.asarray(gs1_ecef)
    gs2_g = cp.asarray(gs2_ecef)
    nodes_dev = cp.empty((N + 2, 3), dtype=cp.float32)

    host = socket.gethostname()
    wall_start = time.time()

    # Progress cadence
    PROG_N = 200  # update every 200 steps

    # --- KPI state across steps ---
    prev_active_names = set()
    curr_window_id = 0
    prev_connected = 0

    # ===== Main loop =====
    for ti, ts in enumerate(times):
        sat_gpu = cp.asarray(pos[ti])  # (3,N)
        nodes_dev[0] = gs1_g
        nodes_dev[1:N+1] = sat_gpu.T
        nodes_dev[N+1] = gs2_g

        # 1) ACTIVE EDGE DISCOVERY (squared distances; sqrt only for kept)
        diff1 = nodes_dev[1:N+1] - gs1_g
        diff2 = nodes_dev[1:N+1] - gs2_g
        d1sq = cp.sum(diff1*diff1, axis=1)
        d2sq = cp.sum(diff2*diff2, axis=1)
        ul1_idx = cp.where(d1sq < UL_MAX2)[0].get().astype(np.int32)
        ul2_idx = cp.where(d2sq < UL_MAX2)[0].get().astype(np.int32)

        if M:
            vi = sat_gpu[:, I_isl]; vj = sat_gpu[:, J_isl]
            dksq = cp.sum((vi - vj)*(vi - vj), axis=0)
            keep = cp.where(dksq < ISL_MAX2)[0].get()
        else:
            keep = np.empty(0, dtype=np.int32)

        # 2) BUILD BASE EDGE LIST (rows, cols, data) ONCE PER STEP
        m1 = ul1_idx.size; m2 = ul2_idx.size; mk = keep.size
        nedges = 2*(m1 + m2 + mk)

        if nedges == 0:
            A_base = csr_matrix((NN, NN), dtype=np.float32)
            edge_rows = edge_cols = None
        else:
            rows = np.empty(nedges, dtype=np.int32)
            cols = np.empty(nedges, dtype=np.int32)
            data = np.empty(nedges, dtype=np.float32)
            cursor = 0

            if m1:
                d1 = (cp.sqrt(d1sq[ul1_idx]).get() / 1e3).astype(np.float32)
                s = np.zeros(m1, dtype=np.int32)
                t = ul1_idx + 1
                rows[cursor:cursor+m1] = s; cols[cursor:cursor+m1] = t; data[cursor:cursor+m1] = d1; cursor += m1
                rows[cursor:cursor+m1] = t; cols[cursor:cursor+m1] = s; data[cursor:cursor+m1] = d1; cursor += m1

            if m2:
                d2 = (cp.sqrt(d2sq[ul2_idx]).get() / 1e3).astype(np.float32)
                s = ul2_idx + 1
                t = np.full(m2, DST, dtype=np.int32)
                rows[cursor:cursor+m2] = s; cols[cursor:cursor+m2] = t; data[cursor:cursor+m2] = d2; cursor += m2
                rows[cursor:cursor+m2] = t; cols[cursor:cursor+m2] = s; data[cursor:cursor+m2] = d2; cursor += m2

            if mk:
                ii = I_isl[keep] + 1
                jj = J_isl[keep] + 1
                dk = (cp.sqrt(dksq[keep]).get() / 1e3).astype(np.float32)
                rows[cursor:cursor+mk] = ii; cols[cursor:cursor+mk] = jj; data[cursor:cursor+mk] = dk; cursor += mk
                rows[cursor:cursor+mk] = jj; cols[cursor:cursor+mk] = ii; data[cursor:cursor+mk] = dk; cursor += mk

            A_base = csr_matrix((data, (rows, cols)), shape=(NN, NN), dtype=np.float32)
            edge_rows, edge_cols, edge_data = rows, cols, data  # keep for masking

        # 3) K=3 NODE-DISJOINT SHORTEST PATHS
        paths = []
        if nedges > 0:
            dist, preds = dijkstra(A_base, directed=False, indices=SRC, return_predecessors=True)
            if np.isfinite(dist[DST]):
                p1 = reconstruct_path(preds.astype(np.int64), SRC, DST)
                if p1 is not None:
                    paths.append(p1)

            removed = set()
            for _ in range(2):
                if not paths: break
                removed.update(paths[-1][1:-1])
                if not removed: break
                rm = np.fromiter(removed, dtype=np.int32)
                touch = np.isin(edge_rows, rm) | np.isin(edge_cols, rm)
                keep_e = ~touch
                if not np.any(keep_e): break
                A_masked = csr_matrix((edge_data[keep_e], (edge_rows[keep_e], edge_cols[keep_e])), shape=(NN, NN))
                dist, preds = dijkstra(A_masked, directed=False, indices=SRC, return_predecessors=True)
                if not np.isfinite(dist[DST]): break
                pN = reconstruct_path(preds.astype(np.int64), SRC, DST)
                if pN is None: break
                paths.append(pN)
                if len(paths) == 3: break

        # 4) SKR via bottleneck hop
        bg = float(bg_series[ti])
        perPathSKR=[0.0,0.0,0.0]; pathStrs=["(none)","(none)","(none)"]
        perPathEtas = [(0.0,0.0)]*3
        perPathQmu  = [0.0]*3
        perPathEmu  = [0.0]*3
        for pi, p in enumerate(paths[:3]):
            names=[]
            for nid in p:
                names.append("Munich" if nid==SRC else ("NewYork" if nid==DST else f"Sat-{nid-1}"))
            skr_val, etaA, etaB, Qmu, Emu = path_min_skr_bottleneck_fast(p, nodes_dev, bg, N)
            perPathSKR[pi] = skr_val
            pathStrs[pi]  = "->".join(names)

            # store diagnostics
            perPathEtas[pi] = (etaA, etaB)
            perPathQmu[pi]  = Qmu
            perPathEmu[pi]  = Emu


        # Aggregate parallel KPIs
        active_path_names = []
        for pi, p in enumerate(paths[:3]):
            if perPathSKR[pi] > 0.0:
                active_path_names.append(pathStrs[pi])

        k_paths = len([r for r in perPathSKR if r > 0.0])
        connected = 1 if k_paths > 0 else 0
        skr_xor = float(min([r for r in perPathSKR if r > 0.0])) if connected else 0.0
        skr_sum = float(sum([r for r in perPathSKR if r > 0.0]))
        skr_min = skr_xor if connected else 0.0

        # Handover detection 
        active_set = set(active_path_names)
        handover = 1 if (active_set != prev_active_names) and connected else 0
        prev_active_names = active_set

        # Window labeling (contiguous connected runs)
        if connected and not prev_connected:
            curr_window_id += 1
        prev_connected = connected

        with open("tfqkd_paths_opt_dataset.csv","a") as f:
            f.write(",".join(map(str, [
                str(ts),
                # Path 1
                pathStrs[0], perPathSKR[0], perPathEtas[0][0], perPathEtas[0][1], perPathQmu[0], perPathEmu[0],
                # Path 2
                pathStrs[1], perPathSKR[1], perPathEtas[1][0], perPathEtas[1][1], perPathQmu[1], perPathEmu[1],
                # Path 3
                pathStrs[2], perPathSKR[2], perPathEtas[2][0], perPathEtas[2][1], perPathQmu[2], perPathEmu[2],
                # Aggregates
                skr_xor,
                skr_sum,
                skr_min,
                k_paths,
                connected,
                handover,
                curr_window_id if connected else 0
            ])) + "\n")

        # ===== live progress line =====
        if ((ti + 1) % PROG_N == 0) or (ti + 1 == T):
            elapsed = time.time() - wall_start
            done = ti + 1
            rate = done / elapsed if elapsed > 0 else 0.0
            remain = (T - done) / rate if rate > 0 else 0.0
            msg = (f"[{time.strftime('%H:%M:%S')}] "
                   f"{done:,}/{T:,}  ({100.0*done/T:5.1f}%)  "
                   f"{rate:7.1f} steps/s  ETA {remain/3600:6.2f} h")
            sys.stdout.write("\r" + msg)
            sys.stdout.flush()
            with open("progress.txt", "w") as pf:
                pf.write(msg + "\n")

    print()

    # ===== Summary & KPI Plots =====
    df = pd.read_csv("tfqkd_paths_opt_dataset.csv", parse_dates=["time"])

    # --- Core tallies ---
    total_bits_bottleneck = (df["skr_xor_bps"]).sum()
    total_bits_parallel   = (df["skr_sum_bps"]).sum()

    # Hourly/Daily using bottleneck (conservative supply)
    hourly = df.groupby(df["time"].dt.floor("h"))["skr_xor_bps"].sum()
    daily  = df.groupby(df["time"].dt.date)["skr_xor_bps"].sum()

    # Availability & percentiles (instantaneous, multi-path sum)
    inst = df["skr_sum_bps"].values
    inst_sorted = np.sort(inst)
    def pct(x): return np.percentile(inst, x)
    p50, p90, p95 = pct(50), pct(90), pct(95)
    outage_prob = float((inst <= 1.0).mean())  # Pr{sum-SKR <= 1 bps}
    avail_1k  = float((inst >= 1e3).mean())
    avail_5k  = float((inst >= 5e3).mean())
    avail_10k = float((inst >= 1e4).mean())
    avail_50k = float((inst >= 5e4).mean())

    # Connectivity & disjointness time shares
    k_counts = df["k_paths"].values
    share_k0 = float((k_counts==0).mean())
    share_k1 = float((k_counts==1).mean())
    share_k2 = float((k_counts==2).mean())
    share_k3 = float((k_counts==3).mean())

    # Handover rate (per hour)
    handovers = df["handover"].sum()
    dur_hours = (df["time"].iloc[-1] - df["time"].iloc[0]).total_seconds()/3600.0
    handover_rate_per_hour = handovers / max(dur_hours,1e-9)

    # Access/outage durations
    win_durs = (df[df["window_id"]>0]
                .groupby("window_id")["time"]
                .agg(lambda s: (s.iloc[-1]-s.iloc[0]).total_seconds()+1.0))
    conn = (df["connected"].values.astype(int))
    idx_change = np.flatnonzero(np.diff(np.r_[conn[0], conn]))
    runs = np.split(np.arange(len(conn)), idx_change+1)
    outage_secs = []
    for r in runs:
        if conn[r[0]]==0:
            t0 = df["time"].iloc[r[0]]
            t1 = df["time"].iloc[r[-1]]
            outage_secs.append((t1 - t0).total_seconds()+1.0)
    outage_durs = pd.Series(outage_secs, dtype="float64")

    # Time-to-key (TTK) ECDF using parallel sum
    dt = 1.0  # simulation time step
    cum_bits = np.cumsum(df["skr_sum_bps"].values * dt)
    def ttk_seconds(K_bits):
        idx = np.searchsorted(cum_bits, K_bits, side="left")
        return np.nan if idx>=len(cum_bits) else (idx+1) * dt
    ttk_1M  = ttk_seconds(1e6)
    ttk_10M = ttk_seconds(1e7)

    # ------- Summary log -------
    def scale_val(val):
        if val<1e3: return f"{val:.2f} bps"
        elif val<1e6: return f"{val/1e3:.2f} kbps"
        elif val<1e9: return f"{val/1e6:.2f} Mbps"
        elif val<1e12: return f"{val/1e9:.2f} Gbps"
        else: return f"{val/1e12:.2f} Tbps"

    with open("tfqkdopt_summary.log","w") as f:
        f.write("=== TF-QKD Simulation Summary ===\n\n")
        f.write(f"Run start: {df['time'].min()}\n")
        f.write(f"Run end:   {df['time'].max()}\n")
        f.write(f"Timesteps: {len(df)}\n\n")
        f.write("--- Totals (delivered) ---\n")
        f.write(f"Bottleneck XOR total: {scale_val(total_bits_bottleneck)}\n")
        f.write(f"Parallel-sum total : {scale_val(total_bits_parallel)}\n\n")
        f.write("--- Instantaneous SKR (parallel) ---\n")
        f.write(f"P50: {scale_val(p50)}  P90: {scale_val(p90)}  P95: {scale_val(p95)}\n")
        f.write(f"Outage prob (<=1 bps): {100*outage_prob:.2f}%\n")
        f.write(f"Availability ≥1/5/10/50 kbps: {100*avail_1k:.1f}% / {100*avail_5k:.1f}% / {100*avail_10k:.1f}% / {100*avail_50k:.1f}%\n\n")
        f.write("--- Disjoint-path availability ---\n")
        f.write(f"k=0: {100*share_k0:.1f}%  k=1: {100*share_k1:.1f}%  k=2: {100*share_k2:.1f}%  k=3: {100*share_k3:.1f}%\n\n")
        f.write(f"Handover rate: {handover_rate_per_hour:.2f} per hour\n\n")
        if len(win_durs):
            f.write(f"Access duration (median / P90): {np.median(win_durs):.1f}s / {np.percentile(win_durs,90):.1f}s\n")
        if len(outage_durs):
            f.write(f"Outage duration (median / P90): {np.median(outage_durs):.1f}s / {np.percentile(outage_durs,90):.1f}s\n")
        f.write("\n--- Time-to-key (parallel sum) ---\n")
        f.write(f"TTK(1 Mbit):  {('NA' if np.isnan(ttk_1M) else f'{ttk_1M:.1f} s')}\n")
        f.write(f"TTK(10 Mbit): {('NA' if np.isnan(ttk_10M) else f'{ttk_10M:.1f} s')}\n\n")
        f.write("--- Hourly totals (XOR bottleneck) ---\n")
        for ts,val in hourly.items():
            f.write(f"{ts}: {scale_val(val)}\n")
        f.write("\n--- Daily totals (XOR bottleneck) ---\n")
        for ts,val in daily.items():
            f.write(f"{ts}: {scale_val(val)}\n")

    # ------- Plot (XOR) -------
    total_seconds=(df["time"].max()-df["time"].min()).total_seconds()
    if total_seconds>2*24*3600:
        df_plot=daily.reset_index()
        x_vals,y_vals=df_plot["time"],df_plot["skr_xor_bps"]
        xlabel,title="Date (UTC)","Daily Accumulated TF-QKD SKR (XOR)"
    elif total_seconds>6*3600:
        df_plot=hourly.reset_index()
        x_vals,y_vals=df_plot["time"],df_plot["skr_xor_bps"]
        xlabel,title="Hour (UTC)","Hourly Accumulated TF-QKD SKR (XOR)"
    else:
        df_plot=df[["time","skr_xor_bps"]]
        x_vals,y_vals=df_plot["time"],df_plot["skr_xor_bps"]
        xlabel,title="Time (UTC)","TF-QKD SKR Time Series (XOR)"

    max_val=y_vals.max()
    if max_val<1e3: scale,unit=1,"bps"
    elif max_val<1e6: scale,unit=1e3,"kbps"
    elif max_val<1e9: scale,unit=1e6,"Mbps"
    elif max_val<1e12: scale,unit=1e9,"Gbps"
    else: scale,unit=1e12,"Tbps"

    plt.figure(figsize=(12,6))
    plt.plot(x_vals, y_vals/scale, linewidth=2.8)
    plt.xlabel(xlabel); plt.ylabel(f"Accumulated XOR SKR [{unit}]")
    plt.title(title); plt.xticks(rotation=45)
    plt.grid(True,axis="both",linestyle="--",alpha=0.5)
    plt.tight_layout()
    plt.savefig("tfqkd_skropt_accumulated.png", dpi=150)
    plt.close()

    # (1) Continuity curve (survival) for instantaneous parallel SKR
    xs = inst_sorted
    ccdf = 1.0 - np.arange(1, len(xs)+1)/len(xs)
    plt.figure(figsize=(10,6))
    plt.semilogx(np.maximum(xs,1e-3), ccdf, linewidth=2.2)
    plt.xlabel("Instantaneous multi-path SKR [bps] (log scale)")
    plt.ylabel("Pr{SKR ≥ r}")
    plt.title("Continuity Curve (Survival Function) of SKR")
    plt.grid(True, which="both", linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.savefig("tfqkd_continuity_curve.png", dpi=150)
    plt.close()

    # (2) Access vs Outage duration histograms (seconds)
    plt.figure(figsize=(10,6))
    if len(win_durs):
        plt.hist(win_durs.values, bins=50, alpha=0.6, label="Access durations")
    if len(outage_durs):
        plt.hist(outage_durs.values, bins=50, alpha=0.6, label="Outage durations")
    plt.xlabel("Duration [s]"); plt.ylabel("Count")
    plt.title("Access / Outage Duration Distributions")
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.legend()
    plt.tight_layout()
    plt.savefig("tfqkd_access_outage_hist.png", dpi=150)
    plt.close()

    # (3) k-disjoint availability bars
    plt.figure(figsize=(8,5))
    shares = [share_k0, share_k1, share_k2, share_k3]
    plt.bar([0,1,2,3], [100*s for s in shares])
    plt.xticks([0,1,2,3], ["k=0","k=1","k=2","k=3"])
    plt.ylabel("Time share [%]")
    plt.title("k-Disjoint Path Availability")
    plt.grid(True, axis="y", linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.savefig("tfqkd_kdisjoint_availability.png", dpi=150)
    plt.close()

    # (4) Cumulative delivered bits over time (parallel sum)
    plt.figure(figsize=(10,6))
    plt.plot(df["time"], cum_bits/1e6, linewidth=2.0)
    plt.xlabel("Time (UTC)")
    plt.ylabel("Cumulative delivered bits [Mbit]")
    plt.title("Cumulative Secret Bits (Parallel Sum)")
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.savefig("tfqkd_cumulative_bits.png", dpi=150)
    plt.close()

    print("Simulation complete. CSV, log and KPI plots saved.")

# ----------------------------------------------------------------------
if __name__=="__main__":
    main()
