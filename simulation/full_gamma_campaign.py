#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time
import traceback
import zipfile

import numpy as np
import pandas as pd
try:
    import pyarrow.parquet as pq
except ModuleNotFoundError:
    pq = None
import yaml

EXPECTED_COMMIT = '601fe3036725af876d6ab46a7cfe15f0c0987be5'
EXPECTED_SOURCE_SHA = '3214d6e3886669e70788d864b780490c3af0fce05ab6f532ef26aad876604f08'
EXPECTED_BINARY_SHA = '7eeba558e14dd8e2c94ab3f03bc29fdbabe102e9f51f4c212f89fbcb869ae0e2'
EXPECTED_STAGE2A_PARTICLE_SHA = 'c18c4e1298691689ac9a9c0e7ba7689390d49e2f61d70b58d6b2f1ad9e456518'
EXPECTED_STAGE2A_RESPONSE_CSV_SHA = '9238f4f902e613961984ba8357a045cf4028bf624cc51ad902b0066cf903f398'
KNOWN_URQMD_FAILURE = 'UrQMD terminating without collision'
START_RESERVE_BYTES = 14 * 2**30
HARD_RESERVE_BYTES = 10 * 2**30
MIN_START_EACH_BYTES = 11 * 2**30
FOLD_MAX_ROWS = 8_000_000
FOLD_MAX_BYTES = 1_610_612_736

ROOT = Path(__file__).resolve().parent
SOFTWARE_WORK = Path(os.environ.get('PENTAGON_CORSIKA_SOFTWARE_WORK', str(Path.home() / 'Desktop' / 'Pentagon_CORSIKA_Response_Stage1'))).expanduser().resolve()
WORK = Path(os.environ.get('PENTAGON_REPRO_SIM_WORK', str(Path.home() / 'Desktop' / 'Pentagon_CORSIKA_FINAL_REPRODUCIBILITY'))).expanduser().resolve()
SRC = SOFTWARE_WORK / 'corsika8-source'
BINARY = SOFTWARE_WORK / 'corsika8-build' / 'applications' / 'c8_air_shower'
PROV = SOFTWARE_WORK / 'BUILD_PROVENANCE.json'
FOLD_SCRIPT = ROOT / 'support' / 'pentagon_unthinned_effective_area.py'
GRID = WORK / 'geant4_response' / 'particle_grid_v1'
HIGH_EM = WORK / 'geant4_response' / 'high_energy_em_extension_v1'
LOGDIR = WORK / 'overnight_stage2_R3_logs'
ANALYSIS_DIR = WORK / 'overnight_stage2_R3_analysis'
MANIFEST_DIR = WORK / 'overnight_stage2_R3_manifests'
RESULT_DIR = WORK / 'OVERNIGHT_STAGE2_R3_COMPACT_RESULT'
RESULT_ZIP = WORK / 'PENTAGON_CORSIKA_OVERNIGHT_STAGE2_R3_RESULT.zip'
FAIL_ZIP = WORK / 'PENTAGON_CORSIKA_OVERNIGHT_STAGE2_R3_FAILURE_DIAGNOSTIC.zip'
STATE_JSON = WORK / 'OVERNIGHT_STAGE2_STATE_R3.json'

START_TIME = time.time()
state = {
    'protocol': 'PENTAGON-CORSIKA-PAPER-FINAL-REPRODUCIBILITY-V1',
    'started_unix': START_TIME,
    'time_window_policy': 'No time cutoff; runtime is a soft planning consideration only',
    'tasks': {},
    'classification_100tev': None,
    'hard_disk_reserve_gib': 10,
    'geomagnetic_scope': 'paired-seed azimuth quartets at 300 TeV for zenith 20 and 40 deg',
    'notes': [],
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def free_bytes() -> int:
    return shutil.disk_usage(WORK).free


def gib(n: int | float) -> float:
    return float(n) / 2**30


def elapsed() -> float:
    return time.time() - START_TIME



def save_state() -> None:
    state['updated_unix'] = time.time()
    STATE_JSON.write_text(json.dumps(state, indent=2, sort_keys=True) + '\n')


def run_quiet(args, **kwargs):
    return subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, **kwargs)


def ensure_file(path: Path, label: str) -> None:
    if not path.is_file():
        raise RuntimeError(f'Missing {label}: {path}')


def provenance_checks() -> None:
    if pq is None:
        raise RuntimeError('pyarrow is missing from the active pentagon-corsika8 environment')
    ensure_file(BINARY, 'validated CORSIKA executable')
    ensure_file(PROV, 'Stage-1 provenance')
    ensure_file(FOLD_SCRIPT, 'detector fold script')
    ensure_file(GRID / 'GEANT4_PARTICLE_RESPONSE_GRID_MANIFEST.json', 'GEANT4 response manifest')
    ensure_file(HIGH_EM / 'GEANT4_HIGH_ENERGY_EM_EXTENSION_MANIFEST.json', 'high-energy EM extension manifest')
    if not (SRC / '.git').is_dir():
        raise RuntimeError(f'Missing CORSIKA source checkout: {SRC}')
    commit = run_quiet(['git', '-C', str(SRC), 'rev-parse', 'HEAD'], check=True).stdout.strip()
    if commit != EXPECTED_COMMIT:
        raise RuntimeError(f'CORSIKA commit mismatch: {commit}')
    source_sha = sha256(SRC / 'applications' / 'c8_air_shower.cpp')
    if source_sha != EXPECTED_SOURCE_SHA:
        raise RuntimeError(f'Minimal application source hash mismatch: {source_sha}')
    binary_sha = sha256(BINARY)
    print(f'[PROVENANCE] executable SHA-256: {binary_sha}')
    prov = json.loads(PROV.read_text())
    if prov.get('corsika_commit') != EXPECTED_COMMIT or prov.get('minimal_source_sha256') != EXPECTED_SOURCE_SHA:
        raise RuntimeError('Stage-1 BUILD_PROVENANCE does not match the validated minimal executable')
    help_run = subprocess.run([str(BINARY), '--help'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if help_run.returncode != 0:
        raise RuntimeError(f'Existing minimal executable no longer starts: --help rc={help_run.returncode}')
    print('[PROVENANCE] Validated Stage-1 R18 executable/source/runtime PASSED.')


def stage2a_checks() -> Path:
    """Generate/reuse the reference 300-TeV, 20-deg, az=0 shower from its fixed seed.

    Unlike the historical R3 wrapper, the reproducibility edition does not require
    a previously archived pilot hash; it validates the newly generated physical output.
    """
    raw = WORK / 'stage2a_gamma_0p30PeV_z20_seed_2026082101'
    meta = validate_shower(raw, 300000, 20, 0) if raw.exists() else None
    if meta is not None:
        print('[PROVENANCE] Existing reference 300-TeV shower validated and reused.')
        return raw
    if raw.exists():
        shutil.rmtree(raw)
    ensure_disk_for_new_shower()
    cmd=[str(BINARY),'--energy','300000','--zenith','20','--azimuth','0',
         '--observation-level','1200','--geomagnetic-year','2017.0','--site-latitude','35.704','--site-longitude','51.351',
         '--emcut','0.0005','--hadcut','0.3','--mucut','0.3','--emthin','1e-6','--max-weight','1',
         '--nevent','1','--seed','2026082101','--verbosity','info','--filename',str(raw)]
    log=LOGDIR/'reference_300tev_z20_seed_2026082101.log'; log.parent.mkdir(parents=True,exist_ok=True)
    print('[SHOWER pilot] Starting reference 300 TeV gamma shower from scratch.')
    with log.open('wb') as h:
        rc=subprocess.run(cmd,stdout=h,stderr=subprocess.STDOUT).returncode
    if rc!=0:
        raise RuntimeError(f'Reference 300-TeV shower failed; see {log}')
    meta=validate_shower(raw,300000,20,0)
    if meta is None:
        raise RuntimeError('Reference 300-TeV shower output failed validation')
    print('[PROVENANCE] Fresh reference 300-TeV shower validated.')
    return raw


def read_primary(out: Path) -> tuple[float, float, float]:
    p = yaml.safe_load((out / 'primary' / 'summary.yaml').read_text())
    if isinstance(p, dict) and 'shower_0' in p:
        p = p['shower_0']
    energy = float(p['total_energy'])
    nx = float(p['nx']); ny = float(p['ny']); nz = float(p['nz'])
    zen = math.degrees(math.acos(max(-1.0, min(1.0, -nz))))
    if abs(math.sin(math.radians(zen))) < 1e-10:
        az = 0.0
    else:
        az = math.degrees(math.atan2(ny, nx)) % 360.0
    return energy, zen, az

def angle_distance_deg(a: float, b: float) -> float:
    return abs((a - b + 180.0) % 360.0 - 180.0)


def validate_shower(out: Path, expected_energy: float, expected_zenith: float, expected_azimuth: float = 0.0) -> dict | None:
    try:
        summary = yaml.safe_load((out / 'summary.yaml').read_text())
        if int(summary.get('showers', -1)) != 1:
            return None
        energy, zen, az = read_primary(out)
        if abs(energy - expected_energy) > max(1e-6, 1e-9 * expected_energy):
            return None
        if abs(zen - expected_zenith) > 1e-6:
            return None
        if expected_zenith > 1e-6 and angle_distance_deg(az, expected_azimuth) > 1e-6:
            return None
        q = out / 'particles' / 'particles.parquet'
        meta = pq.read_metadata(q)
        rows = int(meta.num_rows)
        weights = pq.read_table(q, columns=['weight']).column('weight').to_numpy(zero_copy_only=False)
        if len(weights) != rows or (len(weights) and not np.allclose(weights, 1.0)):
            return None
        return {
            'ground_rows': rows,
            'particle_parquet_bytes': q.stat().st_size,
            'particle_parquet_sha256': sha256(q),
            'all_particle_weights_one': True,
            'primary_azimuth_deg': az,
        }
    except Exception:
        return None

def output_name(label: str, energy_gev: int, zenith: int, azimuth: int, seed: int) -> Path:
    e = {100000:'0p10', 200000:'0p20', 300000:'0p30'}[energy_gev]
    return WORK / f'stage2overnightR3_{label}_gamma_{e}PeV_z{zenith}_az{azimuth}_seed_{seed}'


def legacy_candidates(label: str, energy_gev: int, zenith: int, azimuth: int, seed: int) -> list[Path]:
    out = [output_name(label, energy_gev, zenith, azimuth, seed)]
    if energy_gev == 100000 and zenith == 20 and azimuth == 0:
        if seed == 2026082111:
            out.append(WORK / 'stage2b_gamma_0p10PeV_z20_slot1_seed_2026082111')
        if seed == 2026082112:
            out.append(WORK / 'stage2b_gamma_0p10PeV_z20_slot2_seed_2026082112')
    return out



def task_can_start(label: str, predicted_seconds: int) -> bool:
    # Runtime is deliberately not a stopping criterion. Disk space is the hard guard.
    free = free_bytes()
    if free < MIN_START_EACH_BYTES:
        state['tasks'][label] = {'status':'skipped_low_disk_before_start','free_gib':gib(free)}
        save_state()
        print(f'[DISK] Skipping {label}: less than {gib(MIN_START_EACH_BYTES):.1f} GiB free.')
        return False
    return True


def write_manifest(label: str, out: Path, seed: int, energy_gev: int, zenith: int, azimuth: int, runtime: float, meta: dict, reused: bool) -> Path:
    MANIFEST_DIR.mkdir(parents=True, exist_ok=True)
    dest = MANIFEST_DIR / f'{label}_manifest.json'
    record = {
        'protocol':'PENTAGON-CORSIKA-PAPER-FINAL-REPRODUCIBILITY-V1',
        'purpose':'bounded manuscript-supporting gamma response map; not event-by-event energy reconstruction',
        'label':label,
        'energy_gev':float(energy_gev),
        'energy_pev':float(energy_gev)/1e6,
        'zenith_deg':float(zenith),
        'azimuth_deg':float(azimuth),
        'seed':int(seed),
        'runtime_seconds':float(runtime),
        'reused_existing_output':bool(reused),
        'corsika_output':str(out),
        **meta,
        'binary_sha256':sha256(BINARY),
        'minimal_source_sha256':sha256(SRC / 'applications' / 'c8_air_shower.cpp'),
        'created_unix':time.time(),
    }
    dest.write_text(json.dumps(record, indent=2, sort_keys=True) + '\n')
    return dest


def run_one_attempt(label: str, energy_gev: int, zenith: int, azimuth: int, seed: int) -> tuple[str, Path | None, float, dict | None]:
    for cand in legacy_candidates(label, energy_gev, zenith, azimuth, seed):
        if cand.is_dir():
            meta = validate_shower(cand, energy_gev, zenith, azimuth)
            if meta:
                print(f'[SHOWER {label}] Reusing valid existing output: {cand}')
                write_manifest(label, cand, seed, energy_gev, zenith, azimuth, 0.0, meta, True)
                return 'ok', cand, 0.0, meta
    out = output_name(label, energy_gev, zenith, azimuth, seed)
    if out.exists():
        moved = out.with_name(out.name + f'.incomplete_{int(time.time())}')
        out.rename(moved)
    LOGDIR.mkdir(parents=True, exist_ok=True)
    log = LOGDIR / f'{label}_seed_{seed}.log'
    args = [
        str(BINARY), '--energy', str(energy_gev), '--zenith', str(zenith), '--azimuth', str(azimuth),
        '--observation-level', '1200',
        '--geomagnetic-year', '2017.0', '--site-latitude', '35.704', '--site-longitude', '51.351',
        '--emcut', '0.0005', '--hadcut', '0.3', '--mucut', '0.3',
        '--emthin', '1e-6', '--max-weight', '1',
        '--nevent', '1', '--seed', str(seed), '--verbosity', 'info', '--filename', str(out),
    ]
    print(f'\n[SHOWER {label}] Starting gamma E={energy_gev/1e6:.2f} PeV, zenith={zenith} deg, azimuth={azimuth} deg, seed={seed}.')
    print(f'[SHOWER {label}] One CORSIKA process; hard disk reserve {gib(HARD_RESERVE_BYTES):.0f} GiB.')
    t0 = time.time()
    disk_abort = False
    with log.open('w') as fh:
        proc = subprocess.Popen(args, cwd=str(WORK), stdout=fh, stderr=subprocess.STDOUT, env=os.environ.copy())
        while proc.poll() is None:
            time.sleep(30)
            if free_bytes() < HARD_RESERVE_BYTES:
                disk_abort = True
                fh.write(f'\n[DISK GUARD] Free disk below {gib(HARD_RESERVE_BYTES):.1f} GiB; terminating pid {proc.pid}.\n')
                fh.flush()
                proc.terminate()
                try:
                    proc.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    proc.kill(); proc.wait()
                break
        rc = proc.wait()
    runtime = time.time() - t0
    if disk_abort:
        return 'disk_abort', None, runtime, None
    meta = validate_shower(out, energy_gev, zenith, azimuth) if rc == 0 else None
    if rc == 0 and meta:
        write_manifest(label, out, seed, energy_gev, zenith, azimuth, runtime, meta, False)
        print(f'[SHOWER {label}] VALID: rows={meta["ground_rows"]:,}, raw={gib(meta["particle_parquet_bytes"]):.3f} GiB, runtime={runtime/3600:.2f} h, weights=1.')
        return 'ok', out, runtime, meta
    text = log.read_text(errors='replace') if log.exists() else ''
    if KNOWN_URQMD_FAILURE in text:
        print(f'[SHOWER {label}] Exact known UrQMD technical failure; replacement seed is permitted.')
        return 'urqmd', None, runtime, None
    print(f'[SHOWER {label}] Unrecognized failure rc={rc}; automatic arbitrary retry is forbidden.')
    return 'unknown_failure', None, runtime, None

def run_task(label: str, energy_gev: int, zenith: int, azimuth: int, preferred: int, replacement: int, predicted_seconds: int) -> Path | None:
    if label in state['tasks'] and state['tasks'][label].get('status') == 'completed':
        p = Path(state['tasks'][label]['output'])
        if validate_shower(p, energy_gev, zenith, azimuth):
            print(f'[RESUME] State already records {label} complete; reusing.')
            return p
    if not task_can_start(label, predicted_seconds):
        return None
    before = free_bytes()
    status, out, runtime, meta = run_one_attempt(label, energy_gev, zenith, azimuth, preferred)
    used_seed = preferred
    if status == 'urqmd':
        status, out, runtime2, meta = run_one_attempt(label, energy_gev, zenith, azimuth, replacement)
        runtime += runtime2
        used_seed = replacement
    if status == 'disk_abort':
        raise RuntimeError(f'{label}: hard disk reserve reached; campaign stopped')
    if status == 'unknown_failure':
        raise RuntimeError(f'{label}: unrecognized CORSIKA failure; campaign stopped without arbitrary retry')
    if status != 'ok' or out is None or meta is None:
        state['tasks'][label] = {'status':'not_completed','last_status':status}
        save_state()
        return None
    state['tasks'][label] = {
        'status':'completed', 'energy_gev':energy_gev, 'zenith_deg':zenith, 'azimuth_deg':azimuth, 'seed':used_seed,
        'output':str(out), 'runtime_seconds':runtime, 'ground_rows':meta['ground_rows'],
        'particle_parquet_bytes':meta['particle_parquet_bytes'], 'free_gib_before':gib(before),
        'free_gib_after':gib(free_bytes()), 'planning_runtime_seconds':predicted_seconds,
    }
    save_state()
    return out

def fold_boundary_touched(json_path: Path) -> bool:
    try:
        j = json.loads(json_path.read_text())
        return any(bool(a['candidate_grid']['candidate_touches_grid_boundary']) for a in j.get('audits', []))
    except Exception:
        return False


def fold_one(label: str, shower_dir: Path, base_seed: int, replicas: int = 3, windows=(5,10,20)) -> bool:
    ANALYSIS_DIR.mkdir(parents=True, exist_ok=True)
    LOGDIR.mkdir(parents=True, exist_ok=True)
    prefix = ANALYSIS_DIR / f'{label}_response'
    csv_path = prefix.with_suffix('.csv')
    json_path = prefix.with_suffix('.json')
    if csv_path.exists() and json_path.exists():
        try:
            j = json.loads(json_path.read_text())
            requested_extent_ok = float(j.get('extent_m', 0)) >= 150.0 and not fold_boundary_touched(json_path)
            if int(j.get('replicas', -1)) == replicas and sorted(float(x) for x in j.get('pulse_windows_ns', [])) == sorted(float(x) for x in windows):
                particle_sha = sha256(shower_dir / 'particles' / 'particles.parquet')
                if j.get('particles_sha256') == particle_sha and requested_extent_ok:
                    print(f'[FOLD {label}] Reusing completed contained {replicas}-replica, {list(windows)} ns fold.')
                    return True
        except Exception:
            pass
    e,z,a = read_primary(shower_dir)
    meta = validate_shower(shower_dir, e, z, a)
    if not meta:
        print(f'[FOLD {label}] Cannot validate raw shower; fold skipped.')
        return False
    if meta['ground_rows'] > FOLD_MAX_ROWS or meta['particle_parquet_bytes'] > FOLD_MAX_BYTES:
        print(f'[FOLD {label}] Resource gate skipped fold: rows={meta["ground_rows"]:,}, bytes={meta["particle_parquet_bytes"]}.')
        return False

    def do_fold(extent: int) -> bool:
        log = LOGDIR / f'{label}_fold_extent{extent}.log'
        args = [
            sys.executable, str(FOLD_SCRIPT), '--shower-dir', str(shower_dir),
            '--response-grid', str(GRID), '--high-energy-em-extension', str(HIGH_EM),
            '--extent-m', str(extent), '--step-m', '0.5', '--replicas', str(replicas), '--base-seed', str(base_seed),
            '--pulse-windows-ns', *[str(x) for x in windows], '--query-workers', '1', '--output-prefix', str(prefix),
        ]
        print(f'[FOLD {label}] Running {replicas} response replicas at pulse windows {list(windows)} ns, extent +/-{extent} m.')
        with log.open('w') as fh:
            rc = subprocess.run(args, stdout=fh, stderr=subprocess.STDOUT, cwd=str(WORK), env=os.environ.copy()).returncode
        return rc == 0 and csv_path.exists() and json_path.exists()

    if not do_fold(150):
        print(f'[FOLD {label}] Failed at +/-150 m; raw CORSIKA shower remains valid and no shower rerun is authorized for a fold failure.')
        return False
    if fold_boundary_touched(json_path):
        print(f'[FOLD {label}] Candidate region touches +/-150 m boundary; expanding POST-PROCESSING only to +/-200 m.')
        if not do_fold(200):
            print(f'[FOLD {label}] Expanded fold failed; retaining valid raw shower and the prior diagnostic logs.')
            return False
        if fold_boundary_touched(json_path):
            print(f'[FOLD {label}] NOTE: candidate region still touches +/-200 m boundary; record as a limitation, do not rerun CORSIKA automatically.')
    print(f'[FOLD {label}] Completed.')
    return True

def classify_100(labels: list[str]) -> str:
    details = {}
    all_established = True
    all_zero = True
    all_contained = True
    for label in labels:
        p = ANALYSIS_DIR / f'{label}_response.csv'
        jpath = ANALYSIS_DIR / f'{label}_response.json'
        if not p.exists() or not jpath.exists():
            details[label] = {'status':'missing_fold'}
            all_established = False; all_zero = False
            continue
        d = pd.read_csv(p)
        d = d[np.isclose(d['pulse_window_ns'].astype(float), 10.0)]
        j = json.loads(jpath.read_text())
        contained = not any(bool(a['candidate_grid']['candidate_touches_grid_boundary']) for a in j['audits'])
        all_contained &= contained
        rec = {'contained':contained}
        for thr in (0.5,1.0):
            g = d[np.isclose(d['threshold_mip'].astype(float), thr)]
            vals = g['accepted_core_points'].astype(float).to_numpy()
            rec[str(thr)] = {'values':vals.tolist(), 'median':float(np.median(vals)) if len(vals) else None}
            if not len(vals) or float(np.median(vals)) < 100:
                all_established = False
            if not len(vals) or not np.all(vals == 0):
                all_zero = False
        details[label] = rec
    if not all_contained:
        classification = 'transitional'
    elif all_established:
        classification = 'established'
    elif all_zero:
        classification = 'zero'
    else:
        classification = 'transitional'
    state['classification_100tev'] = {'class':classification, 'details':details, 'rule':'10 ns; 0.5 and 1 MIP; established requires median accepted-core count >=100 at both thresholds in both showers'}
    save_state()
    print(f'\n[DECISION] 100 TeV / 20 deg classification: {classification.upper()}')
    return classification



IGRF_THRAN_2017_UT = np.array([27.922276339253266, -2.2868709563595777, -39.333909160217246], dtype=float)


def primary_vector(out: Path) -> np.ndarray:
    p = yaml.safe_load((out / 'primary' / 'summary.yaml').read_text())
    if isinstance(p, dict) and 'shower_0' in p:
        p = p['shower_0']
    v = np.asarray([float(p['nx']), float(p['ny']), float(p['nz'])], dtype=float)
    n = np.linalg.norm(v)
    if not np.isfinite(n) or n <= 0:
        raise RuntimeError(f'Invalid primary direction in {out}')
    return v / n


def geomagnetic_ground_metrics(label: str, out: Path) -> dict:
    q = out / 'particles' / 'particles.parquet'
    cols = ['pdg', 'x', 'y', 'time', 'kinetic_energy', 'weight']
    d = pd.read_parquet(q, columns=cols)
    if len(d) == 0:
        raise RuntimeError(f'No ground particles for geomagnetic audit: {label}')
    if not np.allclose(d['weight'].to_numpy(float), 1.0):
        raise RuntimeError(f'Non-unit particle weights in geomagnetic audit: {label}')
    v = primary_vector(out)
    B = IGRF_THRAN_2017_UT.copy()
    Bmag = float(np.linalg.norm(B))
    g3 = np.cross(v, B)
    bperp = float(np.linalg.norm(g3))
    alpha = float(np.degrees(np.arcsin(np.clip(bperp / Bmag, 0.0, 1.0))))
    gh = g3[:2]
    ghn = float(np.linalg.norm(gh))
    if ghn > 1e-12:
        u = gh / ghn
        w = np.array([-u[1], u[0]], dtype=float)
    else:
        u = np.array([1.0, 0.0], dtype=float)
        w = np.array([0.0, 1.0], dtype=float)

    x = d['x'].to_numpy(float)
    y = d['y'].to_numpy(float)
    xy = np.column_stack([x, y])
    r = np.hypot(x, y)
    pdg = d['pdg'].to_numpy(int)
    t = d['time'].to_numpy(float)
    ke = d['kinetic_energy'].to_numpy(float)

    def centroid(mask: np.ndarray, weights: np.ndarray | None = None):
        if int(mask.sum()) == 0:
            return [float('nan'), float('nan')]
        a = xy[mask]
        if weights is None:
            return [float(a[:,0].mean()), float(a[:,1].mean())]
        ww = np.asarray(weights[mask], dtype=float)
        ww = np.where(np.isfinite(ww) & (ww > 0), ww, 0.0)
        if float(ww.sum()) <= 0:
            return [float(a[:,0].mean()), float(a[:,1].mean())]
        c = np.average(a, axis=0, weights=ww)
        return [float(c[0]), float(c[1])]

    me = pdg == 11
    pe = pdg == -11
    em_ch = me | pe
    gamma = pdg == 22
    mu = np.abs(pdg) == 13
    proton_neutron = np.isin(np.abs(pdg), [2112, 2212])
    pion = np.isin(np.abs(pdg), [111, 211])

    cem = centroid(me); cep = centroid(pe)
    delta = np.asarray(cep) - np.asarray(cem)
    delta_lorentz = float(np.dot(delta, u)) if np.all(np.isfinite(delta)) else float('nan')
    delta_cross = float(np.dot(delta, w)) if np.all(np.isfinite(delta)) else float('nan')

    if int(em_ch.sum()) >= 2:
        a = xy[em_ch]
        c = a.mean(axis=0)
        rel = a - c
        pu = rel @ u; pw = rel @ w
        rms_u = float(np.sqrt(np.mean(pu**2)))
        rms_w = float(np.sqrt(np.mean(pw**2)))
        ellipticity = float(rms_u / rms_w) if rms_w > 0 else float('nan')
    else:
        rms_u = rms_w = ellipticity = float('nan')

    energy_cem = centroid(me, ke); energy_cep = centroid(pe, ke)
    edelta = np.asarray(energy_cep) - np.asarray(energy_cem)
    edelta_lorentz = float(np.dot(edelta, u)) if np.all(np.isfinite(edelta)) else float('nan')

    e_mask = em_ch | gamma
    radial = r[e_mask]
    time_em = t[e_mask]
    record = {
        'simulation_label': label,
        'shower_dir': str(out),
        'particles_sha256': sha256(q),
        'ground_particles': int(len(d)),
        'primary_direction_north_west_up': [float(z) for z in v],
        'igrf13_field_north_west_up_uT': [float(z) for z in B],
        'field_magnitude_uT': Bmag,
        'B_perp_uT': bperp,
        'geomagnetic_angle_deg': alpha,
        'lorentz_positive_horizontal_unit_north_west': [float(z) for z in u],
        'n_electron': int(me.sum()),
        'n_positron': int(pe.sum()),
        'n_gamma': int(gamma.sum()),
        'n_muon': int(mu.sum()),
        'n_proton_neutron': int(proton_neutron.sum()),
        'n_pion': int(pion.sum()),
        'electron_centroid_xy_m': cem,
        'positron_centroid_xy_m': cep,
        'positron_minus_electron_centroid_along_v_cross_B_m': delta_lorentz,
        'positron_minus_electron_centroid_cross_v_cross_B_m': delta_cross,
        'energy_weighted_eplus_minus_eminus_centroid_along_v_cross_B_m': edelta_lorentz,
        'em_charged_rms_along_v_cross_B_m': rms_u,
        'em_charged_rms_cross_v_cross_B_m': rms_w,
        'em_charged_ellipticity_along_over_cross': ellipticity,
        'em_radial_p50_m': float(np.percentile(radial, 50)) if len(radial) else float('nan'),
        'em_radial_p68_m': float(np.percentile(radial, 68)) if len(radial) else float('nan'),
        'em_radial_p90_m': float(np.percentile(radial, 90)) if len(radial) else float('nan'),
        'em_time_p10_s': float(np.percentile(time_em, 10)) if len(time_em) else float('nan'),
        'em_time_p90_s': float(np.percentile(time_em, 90)) if len(time_em) else float('nan'),
        'note': 'Exploratory real-field diagnostic; one paired-seed quartet per zenith is not a standalone precision geomagnetic measurement.',
    }
    return record


def build_geomagnetic_audits(pilot_raw: Path) -> None:
    print('\n[GEOMAG] Building exploratory ground-footprint and detector-response audits.')
    configs = [
        ('geomag_300tev_z20_az0_pilot', pilot_raw, 300000, 20, 0),
    ]
    z40 = state['tasks'].get('300tev_z40_a', {})
    if z40.get('status') == 'completed':
        configs.append(('geomag_300tev_z40_az0', Path(z40['output']), 300000, 40, 0))
    for label, zen, az in [
        ('geomag_300tev_z20_az90',20,90), ('geomag_300tev_z20_az180',20,180), ('geomag_300tev_z20_az270',20,270),
        ('geomag_300tev_z40_az90',40,90), ('geomag_300tev_z40_az180',40,180), ('geomag_300tev_z40_az270',40,270),
    ]:
        t = state['tasks'].get(label, {})
        if t.get('status') == 'completed':
            configs.append((label, Path(t['output']), 300000, zen, az))

    rows = []
    for label, out, energy, zen, az in configs:
        try:
            rec = geomagnetic_ground_metrics(label, out)
            rec.update({'energy_gev':energy,'energy_pev':energy/1e6,'zenith_deg':zen,'azimuth_deg':az})
            rows.append(rec)
        except Exception as exc:
            state['notes'].append(f'Geomagnetic ground audit failed for {label}: {exc}')
            print(f'[GEOMAG] Ground audit failed for {label}: {exc}')
    if rows:
        ANALYSIS_DIR.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(rows).to_csv(ANALYSIS_DIR/'GEOMAGNETIC_GROUND_FOOTPRINT_AUDIT.csv', index=False)
        (ANALYSIS_DIR/'GEOMAGNETIC_GROUND_FOOTPRINT_AUDIT.json').write_text(json.dumps(rows, indent=2, sort_keys=True) + '\n')

    # Summarize detector response as a function of B_perp. The pilot refold and z40/az0 fold
    # use pre-existing labels, while the other six are explicit geomag labels.
    label_map = {
        'geomag_300tev_z20_az0_pilot': 'pilot_300tev_z20_refold',
        'geomag_300tev_z40_az0': '300tev_z40_a',
        'geomag_300tev_z20_az90': 'geomag_300tev_z20_az90',
        'geomag_300tev_z20_az180': 'geomag_300tev_z20_az180',
        'geomag_300tev_z20_az270': 'geomag_300tev_z20_az270',
        'geomag_300tev_z40_az90': 'geomag_300tev_z40_az90',
        'geomag_300tev_z40_az180': 'geomag_300tev_z40_az180',
        'geomag_300tev_z40_az270': 'geomag_300tev_z40_az270',
    }
    blookup = {r['simulation_label']: r for r in rows}
    dro = []
    for glabel, fold_label in label_map.items():
        csvp = ANALYSIS_DIR / f'{fold_label}_response.csv'
        if not csvp.exists() or glabel not in blookup:
            continue
        d = pd.read_csv(csvp)
        for (window, thr), g in d.groupby(['pulse_window_ns','threshold_mip']):
            dro.append({
                'geomagnetic_label':glabel,
                'fold_label':fold_label,
                'zenith_deg':blookup[glabel]['zenith_deg'],
                'azimuth_deg':blookup[glabel]['azimuth_deg'],
                'B_perp_uT':blookup[glabel]['B_perp_uT'],
                'geomagnetic_angle_deg':blookup[glabel]['geomagnetic_angle_deg'],
                'pulse_window_ns':float(window),
                'threshold_mip':float(thr),
                'n_response_replicas':int(len(g)),
                'accepted_core_points_median':float(g['accepted_core_points'].astype(float).median()),
                'area_proxy_m2_median':float(g['effective_area_m2'].astype(float).median()),
                'p68_angular_error_deg_median':float(g['p68_angular_error_deg'].astype(float).median()),
                'scope_note':'Exploratory B_perp/azimuth diagnostic; shower fluctuation is not fully averaged with one paired-seed quartet.',
            })
    if dro:
        pd.DataFrame(dro).to_csv(ANALYSIS_DIR/'GEOMAGNETIC_DETECTOR_RESPONSE_AUDIT.csv', index=False)
        state['geomagnetic_audit_records'] = len(dro)
    state['geomagnetic_ground_records'] = len(rows)
    save_state()
    print(f'[GEOMAG] Ground configurations audited: {len(rows)}; detector-response records: {len(dro)}.')

def collect_results() -> None:
    if RESULT_DIR.exists(): shutil.rmtree(RESULT_DIR)
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    (RESULT_DIR / 'logs').mkdir(exist_ok=True)
    (RESULT_DIR / 'analysis').mkdir(exist_ok=True)
    (RESULT_DIR / 'manifests').mkdir(exist_ok=True)
    # Copy package documentation.
    for name in ['README_FIRST.txt','OVERNIGHT_PROTOCOL.md','OVERNIGHT_SCIENTIFIC_AUDIT.md','LITERATURE_AND_REFEREE_SCOPE.md','GEOMAGNETIC_EXPLORATORY_DESIGN.md','FOCUSED_SOURCE_TRANSIT_ZENITH_AUDIT.csv']:
        p = ROOT / name
        if p.exists(): shutil.copy2(p, RESULT_DIR / name)
    # Existing pilot provenance.
    for p in [WORK/'STAGE2A_PILOT_MANIFEST.json', WORK/'stage2a_pilot_analysis'/'pilot_0p30PeV_z20_coarse_response.csv', WORK/'stage2a_pilot_analysis'/'pilot_0p30PeV_z20_coarse_response.json']:
        if p.exists(): shutil.copy2(p, RESULT_DIR / ('ORIGINAL_STAGE2A_' + p.name))
    for p in MANIFEST_DIR.glob('*.json'): shutil.copy2(p, RESULT_DIR/'manifests'/p.name)
    for p in LOGDIR.glob('*'): 
        if p.is_file(): shutil.copy2(p, RESULT_DIR/'logs'/p.name)
    for p in ANALYSIS_DIR.glob('*'):
        if p.is_file(): shutil.copy2(p, RESULT_DIR/'analysis'/p.name)
    save_state()
    shutil.copy2(STATE_JSON, RESULT_DIR / STATE_JSON.name)

    # Build combined fold records with simulation labels/energies/zeniths.
    frames=[]
    for p in sorted(ANALYSIS_DIR.glob('*_response.csv')):
        label=p.name[:-len('_response.csv')]
        try:
            d=pd.read_csv(p)
        except Exception:
            continue
        if label == 'pilot_300tev_z20_refold':
            energy,zen,az,seed=300000,20,0,2026082101
        else:
            t=state['tasks'].get(label,{})
            energy=int(t.get('energy_gev',0)); zen=int(t.get('zenith_deg',0)); az=int(t.get('azimuth_deg',0)); seed=t.get('seed')
        d.insert(0,'simulation_label',label); d.insert(1,'energy_gev',energy); d.insert(2,'zenith_deg',zen); d.insert(3,'azimuth_deg',az); d.insert(4,'corsika_seed',seed)
        frames.append(d)
    if frames:
        allr=pd.concat(frames,ignore_index=True)
        allr.to_csv(RESULT_DIR/'OVERNIGHT_RESPONSE_ALL_RECORDS.csv',index=False)
        # Summarize response replicas without claiming calibrated effective area.
        summary=[]
        keys=['energy_gev','zenith_deg','azimuth_deg','pulse_window_ns','threshold_mip']
        for key,g in allr.groupby(keys,dropna=False):
            e,z,az,w,t=key
            areas=g['effective_area_m2'].astype(float)
            acc=g['accepted_core_points'].astype(float)
            p68=g['p68_angular_error_deg'].astype(float)
            summary.append({
                'energy_gev':float(e),'energy_pev':float(e)/1e6,'zenith_deg':float(z),'azimuth_deg':float(az),'pulse_window_ns':float(w),'threshold_mip':float(t),
                'n_fold_records':int(len(g)),'n_simulation_labels':int(g['simulation_label'].nunique()),
                'accepted_core_points_min':float(acc.min()),'accepted_core_points_median':float(acc.median()),'accepted_core_points_max':float(acc.max()),
                'area_proxy_m2_min':float(areas.min()),'area_proxy_m2_median':float(areas.median()),'area_proxy_m2_max':float(areas.max()),
                'p68_deg_median':float(p68.median()),
            })
        pd.DataFrame(summary).to_csv(RESULT_DIR/'OVERNIGHT_RESPONSE_SUMMARY.csv',index=False)

    # Particle-class detector-deposit audit. This is diagnostic only and does not alter triggers.
    class_rows=[]
    for jp in sorted(ANALYSIS_DIR.glob('*_response.json')):
        label=jp.name[:-len('_response.json')]
        try:
            j=json.loads(jp.read_text())
        except Exception:
            continue
        if label == 'pilot_300tev_z20_refold':
            energy,zen,az=300000,20,0
        else:
            t=state['tasks'].get(label,{})
            energy=int(t.get('energy_gev',0)); zen=int(t.get('zenith_deg',0)); az=int(t.get('azimuth_deg',0))
        for audit in j.get('audits',[]):
            by_class=audit.get('sampling',{}).get('by_class',{})
            total=sum(float(v.get('sampled_detector_deposit_mev',0.0)) for v in by_class.values())
            for particle,rec in by_class.items():
                dep=float(rec.get('sampled_detector_deposit_mev',0.0))
                class_rows.append({
                    'simulation_label':label,'energy_gev':energy,'zenith_deg':zen,'azimuth_deg':az,
                    'response_seed':audit.get('response_seed'),'particle_class':particle,
                    'ground_particles_in_class':int(rec.get('particles',0)),
                    'downward_particles_in_class':int(rec.get('downward_particles',rec.get('particles',0))),
                    'nonzero_detector_deposits':int(rec.get('nonzero_deposits',0)),
                    'sampled_detector_deposit_mev':dep,
                    'fraction_of_sampled_detector_deposit':dep/total if total>0 else float('nan'),
                    'energy_below_response_grid':int(rec.get('energy_below_grid',0)),
                    'energy_above_response_grid':int(rec.get('energy_above_grid',0)),
                    'angle_below_response_grid':int(rec.get('angle_below_grid',0)),
                    'angle_above_response_grid':int(rec.get('angle_above_grid',0)),
                })
    if class_rows:
        pd.DataFrame(class_rows).to_csv(RESULT_DIR/'PARTICLE_CLASS_DETECTOR_RESPONSE_AUDIT.csv',index=False)
    # Raw output inventory.
    with (RESULT_DIR/'RAW_OUTPUTS_RETAIN_LOCALLY.txt').open('w') as f:
        f.write('Keep these raw CORSIKA directories on the Mac until review. The compact ZIP excludes particles.parquet.\n\n')
        f.write(str(WORK/'stage2a_gamma_0p30PeV_z20_seed_2026082101')+'\n')
        f.write('  particles.parquet SHA256: '+EXPECTED_STAGE2A_PARTICLE_SHA+'\n')
        for label,t in state['tasks'].items():
            if t.get('status')=='completed':
                f.write(str(t['output'])+'\n')
                mp=MANIFEST_DIR/f'{label}_manifest.json'
                if mp.exists():
                    m=json.loads(mp.read_text())
                    f.write('  particles.parquet SHA256: '+m['particle_parquet_sha256']+'\n')
                    f.write('  particles.parquet bytes: '+str(m['particle_parquet_bytes'])+'\n')

    if RESULT_ZIP.exists(): RESULT_ZIP.unlink()
    with zipfile.ZipFile(RESULT_ZIP,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in RESULT_DIR.rglob('*'):
            if p.is_file(): z.write(p, arcname=str(Path(RESULT_DIR.name)/p.relative_to(RESULT_DIR)))
    print(f'\n[PACKAGE] Compact result ZIP: {RESULT_ZIP}')
    print(f'[PACKAGE] SHA-256: {sha256(RESULT_ZIP)}')


def make_failure_zip(exc: BaseException) -> None:
    try:
        diag = WORK / 'OVERNIGHT_STAGE2_R3_FAILURE_DIAGNOSTIC'
        if diag.exists(): shutil.rmtree(diag)
        diag.mkdir(parents=True)
        (diag/'FAILURE.txt').write_text(''.join(traceback.format_exception(type(exc),exc,exc.__traceback__)))
        save_state()
        if STATE_JSON.exists(): shutil.copy2(STATE_JSON,diag/STATE_JSON.name)
        for srcdir,name in [(LOGDIR,'logs'),(MANIFEST_DIR,'manifests'),(ANALYSIS_DIR,'analysis')]:
            if srcdir.exists(): shutil.copytree(srcdir,diag/name,dirs_exist_ok=True)
        if FAIL_ZIP.exists(): FAIL_ZIP.unlink()
        with zipfile.ZipFile(FAIL_ZIP,'w',compression=zipfile.ZIP_DEFLATED) as z:
            for p in diag.rglob('*'):
                if p.is_file(): z.write(p,arcname=str(Path(diag.name)/p.relative_to(diag)))
        print(f'\n[FAILED] Overnight Stage 2 stopped: {exc}')
        print(f'[FAILED] Diagnostic ZIP: {FAIL_ZIP}')
    except Exception as e:
        print(f'[FAILED] Also failed to create diagnostic ZIP: {e}')


def plan_only() -> None:
    print('FINAL REPRODUCIBILITY GAMMA-CAMPAIGN PLAN (no simulation executed)')
    print('Time: no automatic cutoff; runtime estimates are planning information only.')
    print(f'Hard disk reserve: {gib(HARD_RESERVE_BYTES):.0f} GiB; new shower starts require >= {gib(MIN_START_EACH_BYTES):.0f} GiB free.')
    print('Core response map: 100, 200, 300 TeV x zenith 0, 20, 40 deg; baseline 20-deg points are replicated. Geomagnetic extension: 300 TeV azimuth quartets at zenith 20 and 40 deg.')
    print('Priority: energy baseline/fluctuations -> 3x3 energy-zenith map -> paired 300-TeV azimuth quartets at 20 and 40 deg -> geomagnetic footprint/response audit.')
    print('Every valid shower: 3 detector-response replicas x 5,10,20 ns x 0.25,0.5,1,2 MIP.')
    print('There is no manual stop-flag mechanism in R3; the campaign runs the predeclared plan unless the hard disk guard intervenes.')


def main() -> int:
    ap=argparse.ArgumentParser()
    ap.add_argument('--plan-only',action='store_true')
    args=ap.parse_args()
    if args.plan_only:
        plan_only(); return 0
    LOGDIR.mkdir(parents=True,exist_ok=True); ANALYSIS_DIR.mkdir(parents=True,exist_ok=True); MANIFEST_DIR.mkdir(parents=True,exist_ok=True)
    print('============================================================')
    print(' PENTAGON PAPER - FINAL GAMMA RESPONSE REPRODUCIBILITY')
    print('============================================================')
    print(f'Workspace: {WORK}')
    print('Time policy: NO automatic time stop; runtime is a soft planning consideration only.')
    print(f'Hard disk reserve: {gib(HARD_RESERVE_BYTES):.0f} GiB')
    print('[NETWORK] No internet access is used by this program.')
    try:
        provenance_checks()
        pilot_raw=stage2a_checks()
        if free_bytes() < START_RESERVE_BYTES:
            raise RuntimeError(f'Less than {gib(START_RESERVE_BYTES):.0f} GiB free at start: {gib(free_bytes()):.2f} GiB')
        print(f'[DISK] Free at start: {gib(free_bytes()):.2f} GiB')
        print('[SCIENCE] Core target: 100/200/300 TeV x zenith 0/20/40 deg response map plus paired 300-TeV geomagnetic azimuth quartets at 20/40 deg.')
        save_state()

        # Existing pilot: enrich detector-systematics sampling, no new CORSIKA.
        fold_one('pilot_300tev_z20_refold',pilot_raw,2026082691,replicas=3,windows=(5,10,20))

        # Priority 1: baseline energy trend + shower-to-shower fluctuations at representative 20 deg.
        p100a=run_task('100tev_z20_a',100000,20,0,2026082111,2026082113,4500)
        if p100a: fold_one('100tev_z20_a',p100a,2026082291)
        p100b=run_task('100tev_z20_b',100000,20,0,2026082112,2026082114,4500)
        if p100b: fold_one('100tev_z20_b',p100b,2026082391)
        classify_100(['100tev_z20_a','100tev_z20_b'])

        p200a=run_task('200tev_z20_a',200000,20,0,2026082121,2026082123,7200)
        if p200a: fold_one('200tev_z20_a',p200a,2026082491)
        p200b=run_task('200tev_z20_b',200000,20,0,2026082122,2026082124,7200)
        if p200b: fold_one('200tev_z20_b',p200b,2026082494)

        p300b=run_task('300tev_z20_b',300000,20,0,2026082102,2026082105,9720)
        if p300b: fold_one('300tev_z20_b',p300b,2026082591)
        p300c=run_task('300tev_z20_c',300000,20,0,2026082103,2026082106,9720)
        if p300c: fold_one('300tev_z20_c',p300c,2026082692)

        # Priority 2: central-energy zenith span across source-transit/analysis range.
        p300z0=run_task('300tev_z0_a',300000,0,0,2026082131,2026082133,9000)
        if p300z0: fold_one('300tev_z0_a',p300z0,2026082695)
        p300z40=run_task('300tev_z40_a',300000,40,0,2026082141,2026082143,11520)
        if p300z40: fold_one('300tev_z40_a',p300z40,2026082791)

        # Priority 3: complete the 3x3 energy-zenith response map with single qualitative corner showers.
        p100z0=run_task('100tev_z0_a',100000,0,0,2026082151,2026082153,4200)
        if p100z0: fold_one('100tev_z0_a',p100z0,2026082794)
        p100z40=run_task('100tev_z40_a',100000,40,0,2026082161,2026082163,5400)
        if p100z40: fold_one('100tev_z40_a',p100z40,2026082891)
        p200z0=run_task('200tev_z0_a',200000,0,0,2026082171,2026082173,6600)
        if p200z0: fold_one('200tev_z0_a',p200z0,2026082894)
        p200z40=run_task('200tev_z40_a',200000,40,0,2026082181,2026082183,8700)
        if p200z40: fold_one('200tev_z40_a',p200z40,2026082991)

        # Priority 4: exploratory geomagnetic sub-study.
        # Use paired nominal seeds within each zenith quartet so that azimuth, and therefore B_perp,
        # is the intended changing variable. The existing Stage-2A pilot supplies z20/az0.
        p300z20az90=run_task('geomag_300tev_z20_az90',300000,20,90,2026082101,2026082107,9720)
        if p300z20az90: fold_one('geomag_300tev_z20_az90',p300z20az90,2026083191)
        p300z20az180=run_task('geomag_300tev_z20_az180',300000,20,180,2026082101,2026082108,9720)
        if p300z20az180: fold_one('geomag_300tev_z20_az180',p300z20az180,2026083194)
        p300z20az270=run_task('geomag_300tev_z20_az270',300000,20,270,2026082101,2026082109,9720)
        if p300z20az270: fold_one('geomag_300tev_z20_az270',p300z20az270,2026083197)

        # The z40/az0 member is the 300tev_z40_a shower above. Complete the matched azimuth quartet.
        p300z40az90=run_task('geomag_300tev_z40_az90',300000,40,90,2026082141,2026082144,11520)
        if p300z40az90: fold_one('geomag_300tev_z40_az90',p300z40az90,2026083291)
        p300z40az180=run_task('geomag_300tev_z40_az180',300000,40,180,2026082141,2026082145,11520)
        if p300z40az180: fold_one('geomag_300tev_z40_az180',p300z40az180,2026083294)
        p300z40az270=run_task('geomag_300tev_z40_az270',300000,40,270,2026082141,2026082146,11520)
        if p300z40az270: fold_one('geomag_300tev_z40_az270',p300z40az270,2026083297)

        build_geomagnetic_audits(pilot_raw)

        state['completed_unix']=time.time(); state['elapsed_hours']=elapsed()/3600; state['free_gib_end']=gib(free_bytes()); save_state()
        collect_results()
        print('\n============================================================')
        print(' FINAL GAMMA RESPONSE REPRODUCIBILITY COMPLETE')
        print('============================================================')
        print(f'Elapsed: {elapsed()/3600:.2f} h')
        print(f'Free disk: {gib(free_bytes()):.2f} GiB')
        print(f'100 TeV diagnostic classification: {state["classification_100tev"]["class"] if isinstance(state.get("classification_100tev"),dict) else state.get("classification_100tev")}')
        print(f'Compact result ZIP: {RESULT_ZIP}')
        print(f'Compact result ZIP SHA-256: {sha256(RESULT_ZIP)}')
        print('Do NOT start another CORSIKA campaign. Upload the compact ZIP to ChatGPT for review.')
        return 0
    except BaseException as exc:
        state['fatal_error']=str(exc); state['elapsed_hours']=elapsed()/3600; save_state()
        make_failure_zip(exc)
        return 2


if __name__=='__main__':
    raise SystemExit(main())
