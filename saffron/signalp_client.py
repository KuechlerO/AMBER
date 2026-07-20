"""Local SignalP 6.0 CLI wrapper (fast mode) with parse helpers and mock fallback."""

from __future__ import annotations

import csv
import os
import re
import shutil
import subprocess
import tempfile
import uuid
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from django.conf import settings


CS_RE = re.compile(
    r'(?:between\s+)?(?P<a>\d+)\s*(?:and|-)\s*(?:CS\s+)?(?P<b>\d+).*?(?:Pr\s*=?\s*|prob(?:ability)?\s*=?\s*)?(?P<p>0?\.\d+|1(?:\.0+)?)?',
    re.IGNORECASE,
)


@dataclass
class SignalPPrediction:
    seq_id: str
    prediction: str
    probabilities: dict[str, float] = field(default_factory=dict)
    cs_before: int | None = None
    cs_after: int | None = None
    cs_prob: float | None = None
    sp_start: int | None = None
    sp_end: int | None = None
    regions: list[dict[str, Any]] = field(default_factory=list)
    plot_path: str | None = None
    plot_url: str | None = None

    @property
    def sp_prob(self) -> float | None:
        if not self.prediction:
            return None
        return self.probabilities.get(self.prediction)

    def class_prob(self, class_name: str) -> float | None:
        return self.probabilities.get(class_name)


def signalp_bin() -> str:
    return os.environ.get('SIGNALP6_BIN') or getattr(settings, 'SIGNALP6_BIN', 'signalp6')


def signalp_model_dir() -> str | None:
    env = os.environ.get('SIGNALP6_MODEL_DIR') or getattr(settings, 'SIGNALP6_MODEL_DIR', None)
    return env or None


def signalp_cache_dir() -> Path:
    raw = os.environ.get('SIGNALP_CACHE_DIR') or getattr(
        settings, 'SIGNALP_CACHE_DIR', None
    )
    if raw:
        path = Path(raw)
    else:
        path = Path(getattr(settings, 'CACHE_DIR', Path(tempfile.gettempdir()))) / 'saffron_signalp'
    path.mkdir(parents=True, exist_ok=True)
    return path


def _probe_signalp_cli() -> dict[str, Any]:
    """
    Probe the SignalP binary.

    The predector conda package ships a stub that exits 1 until the licensed
    tarball is registered via `signalp6-register signalp-6.0h.fast.tar.gz`.
    """
    bin_path = signalp_bin()
    resolved = shutil.which(bin_path) or (bin_path if Path(bin_path).is_file() else None)
    if not resolved:
        return {
            'available': False,
            'registered': False,
            'bin': bin_path,
            'detail': 'signalp6 binary not found on PATH',
        }
    try:
        proc = subprocess.run(
            [resolved, '-h'],
            capture_output=True,
            text=True,
            timeout=20,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        return {
            'available': False,
            'registered': False,
            'bin': resolved,
            'detail': f'could not run signalp6: {exc}',
        }
    out = f'{proc.stdout or ""}{proc.stderr or ""}'
    if 'has not been installed yet' in out or 'signalp6-register' in out:
        return {
            'available': False,
            'registered': False,
            'bin': resolved,
            'detail': (
                'predector::signalp6 stub is present but not registered. '
                'Download signalp-6.0h.fast.tar.gz from DTU and run: '
                'signalp6-register /path/to/signalp-6.0h.fast.tar.gz'
            ),
        }
    # Registered CLI: -h may return 0 or print usage with another code
    ok = proc.returncode == 0 or 'usage' in out.lower() or '--fastafile' in out or 'SignalP' in out
    return {
        'available': bool(ok),
        'registered': bool(ok),
        'bin': resolved,
        'detail': 'SignalP 6.0 CLI ready' if ok else (out.strip()[:400] or f'exit {proc.returncode}'),
    }


def signalp_available() -> bool:
    if os.environ.get('SIGNALP_MOCK', '').lower() in ('1', 'true', 'yes'):
        return False
    return bool(_probe_signalp_cli()['available'])


def signalp_status() -> dict[str, Any]:
    """UI-facing SignalP backend status."""
    if os.environ.get('SIGNALP_MOCK', '').lower() in ('1', 'true', 'yes'):
        return {
            'available': False,
            'registered': False,
            'backend': 'mock',
            'bin': signalp_bin(),
            'model_dir': signalp_model_dir(),
            'detail': 'SIGNALP_MOCK is set',
        }
    probe = _probe_signalp_cli()
    return {
        'available': probe['available'],
        'registered': probe.get('registered', False),
        'backend': 'cli' if probe['available'] else 'mock',
        'bin': probe.get('bin') or signalp_bin(),
        'model_dir': signalp_model_dir(),
        'detail': probe.get('detail') or '',
    }


def parse_cs_position(text: str | None) -> tuple[int | None, int | None, float | None]:
    if not text or not str(text).strip() or str(text).strip() == '-':
        return None, None, None
    m = CS_RE.search(str(text))
    if not m:
        digits = re.findall(r'\d+', str(text))
        if len(digits) >= 2:
            return int(digits[0]), int(digits[1]), None
        return None, None, None
    a = int(m.group('a'))
    b = int(m.group('b'))
    p = float(m.group('p')) if m.group('p') else None
    return a, b, p


def parse_prediction_results(path: Path) -> dict[str, SignalPPrediction]:
    """Parse SignalP prediction_results.txt into id -> SignalPPrediction."""
    text = path.read_text(encoding='utf-8', errors='replace')
    raw_lines = [ln for ln in text.splitlines() if ln.strip()]
    if not raw_lines:
        return {}

    # SignalP 6 writes the column header as a commented line: "# ID\tPrediction\t..."
    header_line = None
    data_lines: list[str] = []
    for ln in raw_lines:
        if ln.startswith('#'):
            body = ln.lstrip('#').strip()
            if body.upper().startswith('ID') and 'Prediction' in body:
                header_line = body
            continue
        data_lines.append(ln)
    if not data_lines:
        return {}
    if header_line is None:
        # Fallback: first non-comment line is the header (older layouts)
        header_line, data_lines = data_lines[0], data_lines[1:]

    reader = csv.DictReader([header_line, *data_lines], delimiter='\t')
    out: dict[str, SignalPPrediction] = {}
    for row in reader:
        # Normalize header keys
        norm = {(k or '').strip(): (v or '').strip() for k, v in row.items() if k}
        seq_id = norm.get('ID') or norm.get('Name') or next(iter(norm.values()), '')
        if not seq_id:
            continue
        prediction = norm.get('Prediction') or 'OTHER'
        # Map verbose labels → short class names used elsewhere
        pred_map = {
            'Signal Peptide (Sec/SPI)': 'SP',
            'Lipoprotein signal peptide (Sec/SPII)': 'LIPO',
            'Tat signal peptide (Tat/SPI)': 'TAT',
            'Tat lipoprotein signal peptide (Tat/SPII)': 'TATLIPO',
            'Pilin-like signal peptide (Sec/SPIII)': 'PILIN',
            'Other': 'OTHER',
            'OTHER': 'OTHER',
            'SP': 'SP',
            'LIPO': 'LIPO',
            'TAT': 'TAT',
            'TATLIPO': 'TATLIPO',
            'PILIN': 'PILIN',
        }
        prediction = pred_map.get(prediction, prediction)

        probs: dict[str, float] = {}
        alias = {
            'OTHER': 'OTHER',
            'Other': 'OTHER',
            'SP': 'SP',
            'SP(Sec/SPI)': 'SP',
            'Sec/SPI': 'SP',
            'LIPO': 'LIPO',
            'SP(Sec/SPII)': 'LIPO',
            'TAT': 'TAT',
            'TAT(Tat/SPI)': 'TAT',
            'TATLIPO': 'TATLIPO',
            'PILIN': 'PILIN',
        }
        for key, val in norm.items():
            short = alias.get(key)
            if not short and key.upper().startswith('SP('):
                short = 'SP'
            if not short:
                continue
            if val in (None, ''):
                continue
            try:
                probs[short] = float(val)
            except ValueError:
                pass
        cs_before, cs_after, cs_prob = parse_cs_position(norm.get('CS Position') or norm.get('CS'))
        out[seq_id] = SignalPPrediction(
            seq_id=seq_id,
            prediction=prediction,
            probabilities=probs,
            cs_before=cs_before,
            cs_after=cs_after,
            cs_prob=cs_prob,
        )
    return out


def parse_region_gff3(path: Path) -> dict[str, list[dict[str, Any]]]:
    """Parse region_output.gff3 keyed by sequence id."""
    if not path.is_file():
        return {}
    by_id: dict[str, list[dict[str, Any]]] = {}
    for line in path.read_text(encoding='utf-8', errors='replace').splitlines():
        if not line.strip() or line.startswith('#'):
            continue
        parts = line.split('\t')
        if len(parts) < 5:
            continue
        seq_id, _src, feature, start, end = parts[0], parts[1], parts[2], parts[3], parts[4]
        try:
            start_i, end_i = int(start), int(end)
        except ValueError:
            continue
        by_id.setdefault(seq_id, []).append({
            'feature': feature,
            'start': start_i,
            'end': end_i,
        })
    return by_id


def _attach_regions(preds: dict[str, SignalPPrediction], regions: dict[str, list[dict]]) -> None:
    for seq_id, pred in preds.items():
        regs = regions.get(seq_id) or []
        pred.regions = regs
        sp_feats = [r for r in regs if 'signal' in r['feature'].lower() or r['feature'] in ('SP', 'signal_peptide')]
        if sp_feats:
            pred.sp_start = min(r['start'] for r in sp_feats)
            pred.sp_end = max(r['end'] for r in sp_feats)
        elif pred.cs_before:
            pred.sp_start = 1
            pred.sp_end = pred.cs_before
        elif regs:
            pred.sp_start = min(r['start'] for r in regs)
            pred.sp_end = max(r['end'] for r in regs)


def _mock_predict(sequences: dict[str, str]) -> dict[str, SignalPPrediction]:
    """Deterministic stand-in when SignalP CLI is unavailable (tests / local stub)."""
    out: dict[str, SignalPPrediction] = {}
    hydrophobic = set('AILMFVW')
    for seq_id, seq in sequences.items():
        n = min(len(seq), 40)
        head = seq[:n].upper()
        hydro_frac = sum(1 for c in head if c in hydrophobic) / max(n, 1)
        has_sp = hydro_frac >= 0.35 and len(seq) >= 15
        if has_sp:
            cs = max(15, min(30, int(12 + hydro_frac * 20)))
            cs = min(cs, len(seq) - 1)
            pred = SignalPPrediction(
                seq_id=seq_id,
                prediction='SP',
                probabilities={
                    'OTHER': round(1 - hydro_frac, 4),
                    'SP': round(min(0.99, hydro_frac + 0.4), 4),
                    'LIPO': 0.0,
                    'TAT': 0.0,
                    'TATLIPO': 0.0,
                    'PILIN': 0.0,
                },
                cs_before=cs,
                cs_after=cs + 1,
                cs_prob=round(0.7 + hydro_frac * 0.2, 4),
                sp_start=1,
                sp_end=cs,
                regions=[{'feature': 'signal_peptide', 'start': 1, 'end': cs}],
            )
        else:
            pred = SignalPPrediction(
                seq_id=seq_id,
                prediction='OTHER',
                probabilities={
                    'OTHER': round(min(0.99, 0.6 + (1 - hydro_frac) * 0.3), 4),
                    'SP': round(hydro_frac * 0.3, 4),
                    'LIPO': 0.0,
                    'TAT': 0.0,
                    'TATLIPO': 0.0,
                    'PILIN': 0.0,
                },
                regions=[],
            )
        out[seq_id] = pred
    return out


def _seq_id_from_plot_stem(stem: str) -> str:
    """Map SignalP plot basenames to FASTA ids (6.0i uses output_<ID>_plot.*)."""
    name = stem
    if name.endswith('_plot'):
        name = name[: -len('_plot')]
    if name.startswith('output_'):
        name = name[len('output_') :]
    return name


def _collect_plots(output_dir: Path, preds: dict[str, SignalPPrediction], media_subdir: Path) -> None:
    """Prefer SignalP's own PNGs; otherwise render from *_plot.txt (probability vs position)."""
    media_subdir.mkdir(parents=True, exist_ok=True)
    for png in output_dir.glob('*plot.png'):
        seq_id = _seq_id_from_plot_stem(png.stem)
        dest = media_subdir / f'{seq_id}_plot.png'
        shutil.copy2(png, dest)
        if seq_id in preds:
            preds[seq_id].plot_path = str(dest)

    for txt in output_dir.glob('*plot.txt'):
        seq_id = _seq_id_from_plot_stem(txt.stem)
        if seq_id not in preds:
            continue
        if preds[seq_id].plot_path and Path(preds[seq_id].plot_path).is_file():
            continue
        series = parse_plot_txt(txt)
        if not series:
            continue
        dest = media_subdir / f'{seq_id}_plot.png'
        _render_probability_plot(
            dest,
            title=f'SignalP 6.0 · {preds[seq_id].seq_id} · {preds[seq_id].prediction}',
            series=series,
            cs_before=preds[seq_id].cs_before,
        )
        preds[seq_id].plot_path = str(dest)


def parse_plot_txt(path: Path) -> dict[str, list[float]]:
    """
    Parse SignalP SEQUENCE_plot.txt into named probability series keyed by column.

    Typical columns include position plus region / type probabilities
    (e.g. OTHER, SP, n, h, c, Sec/SPI n, …). Exact headers vary by version.
    """
    if not path.is_file():
        return {}
    with path.open(encoding='utf-8', errors='replace') as fh:
        rows = []
        header = None
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith('#'):
                # "# pos\taa\tlabel\tOther\tSec/SPI n\t..."
                body = line.lstrip('#').strip()
                if body.lower().startswith('pos') or body.lower().startswith('name='):
                    if body.lower().startswith('pos'):
                        header = body.split('\t')
                continue
            if header is None:
                header = line.split('\t')
                continue
            rows.append(line.split('\t'))
    if not header or not rows:
        return {}

    cols = {h.strip(): i for i, h in enumerate(header)}
    pos_key = next((k for k in cols if k.lower() in ('pos', 'position', '#', 'aa', 'idx', 'index')), None)

    # Prefer region marginals used on the DTU web plot
    rename = {
        'Other': 'OTHER',
        'OTHER': 'OTHER',
        'Sec/SPI n': 'n',
        'Sec/SPI h': 'h',
        'Sec/SPI c': 'c',
    }
    preferred = []
    for h in header:
        if h == pos_key or h.lower() in ('aa', 'amino', 'residue', 'label', 'pred', 'prediction'):
            continue
        preferred.append(h)

    series: dict[str, list[float]] = {}
    for h in preferred:
        key = rename.get(h, h)
        series.setdefault(key, [])
        i = cols[h]
        for row in rows:
            try:
                series[key].append(float(row[i]))
            except (IndexError, ValueError):
                series[key].append(0.0)

    # Keep the most informative eukarya curves if present
    keep_order = ['n', 'h', 'c', 'OTHER']
    if any(k in series for k in keep_order):
        series = {k: series[k] for k in keep_order if k in series and any(abs(x) > 1e-6 for x in series[k])}
    else:
        series = {k: v for k, v in series.items() if any(abs(x) > 1e-6 for x in v)}
    return series


def _synthetic_region_series(pred: SignalPPrediction, length: int) -> dict[str, list[float]]:
    """Build approximate n/h/c/OTHER curves for mock mode (not a SignalP model)."""
    length = max(length, pred.sp_end or pred.cs_before or 40, 20)
    xs = list(range(1, length + 1))
    sp_end = int(pred.sp_end or pred.cs_before or max(15, length // 4))
    sp_end = min(sp_end, length)
    n_end = max(1, sp_end // 3)
    h_end = max(n_end + 1, (2 * sp_end) // 3)
    c_end = sp_end

    def _bell(center: float, width: float, peak: float, x: int) -> float:
        return peak * max(0.0, 1.0 - abs(x - center) / max(width, 1e-6))

    n_vals, h_vals, c_vals, other = [], [], [], []
    sp_peak = float(pred.probabilities.get('SP') or 0.7)
    for x in xs:
        n = _bell((1 + n_end) / 2, max(2, n_end), sp_peak * 0.85, x) if x <= c_end else 0.0
        h = _bell((n_end + h_end) / 2, max(2, h_end - n_end), sp_peak, x) if x <= c_end else 0.0
        c = _bell((h_end + c_end) / 2, max(2, c_end - h_end), sp_peak * 0.75, x) if x <= c_end else 0.0
        o = max(0.0, 1.0 - (n + h + c))
        if x > c_end:
            o = min(0.99, 0.85 + 0.1 * (x - c_end) / max(1, length - c_end))
            n = h = c = 0.0
        n_vals.append(round(n, 4))
        h_vals.append(round(h, 4))
        c_vals.append(round(c, 4))
        other.append(round(o, 4))
    return {'n': n_vals, 'h': h_vals, 'c': c_vals, 'OTHER': other}


def _render_probability_plot(
    dest: Path,
    *,
    title: str,
    series: dict[str, list[float]],
    sequence: str | None = None,
    cs_before: int | None = None,
    subtitle: str | None = None,
) -> None:
    """Draw probability (y) vs amino-acid position (x), SignalP-style."""
    import matplotlib

    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    if not series:
        raise ValueError('No series to plot')

    length = max(len(v) for v in series.values())
    positions = list(range(1, length + 1))
    # Stable colour map similar to SignalP web plots
    colors = {
        'n': '#1f77b4',
        'h': '#ff7f0e',
        'c': '#2ca02c',
        'OTHER': '#7f7f7f',
        'SP': '#d62728',
        'LIPO': '#9467bd',
        'TAT': '#8c564b',
        'TATLIPO': '#e377c2',
        'PILIN': '#17becf',
    }

    fig, ax = plt.subplots(figsize=(9.5, 4.2), dpi=120)
    for name, values in series.items():
        ys = list(values) + [0.0] * (length - len(values))
        ax.plot(positions, ys[:length], label=name, color=colors.get(name, None), linewidth=1.8)

    if cs_before and 1 <= cs_before <= length:
        ax.axvline(cs_before, color='#111827', linestyle='--', linewidth=1.2, label=f'CS {cs_before}')

    ax.set_xlabel('Amino acid position')
    ax.set_ylabel('Probability')
    ax.set_ylim(-0.02, 1.05)
    ax.set_xlim(1, length)
    ax.set_title(title, fontsize=11)
    if subtitle:
        ax.text(0.01, 1.02, subtitle, transform=ax.transAxes, fontsize=8, color='#b45309')
    ax.grid(True, alpha=0.25)
    ax.legend(loc='upper right', fontsize=8, framealpha=0.9)

    if sequence and len(sequence) == length and length <= 80:
        ax2 = ax.twiny()
        ax2.set_xlim(ax.get_xlim())
        ticks = list(range(1, length + 1, max(1, length // 20)))
        ax2.set_xticks(ticks)
        ax2.set_xticklabels([sequence[i - 1] for i in ticks], fontsize=7)
        ax2.set_xlabel('Sequence')

    dest.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(dest, format='png')
    plt.close(fig)


def _attach_mock_plots(
    preds: dict[str, SignalPPrediction],
    *,
    job_id: str | None,
    sequences: dict[str, str] | None = None,
) -> None:
    job = job_id or uuid.uuid4().hex
    plot_dir = signalp_cache_dir() / job / 'plots'
    plot_dir.mkdir(parents=True, exist_ok=True)
    sequences = sequences or {}
    for seq_id, pred in preds.items():
        safe = re.sub(r'[^A-Za-z0-9_\-]', '_', seq_id)[:80] or 'seq'
        dest = plot_dir / f'{safe}_plot.png'
        seq = sequences.get(seq_id) or sequences.get(pred.seq_id) or ''
        length = len(seq) if seq else max(pred.sp_end or 0, pred.cs_before or 0, 40)
        # Truncate display to N-term window used by SignalP plots (~70 aa)
        window = min(max(length, 20), 70)
        series = _synthetic_region_series(pred, window)
        _render_probability_plot(
            dest,
            title=f'{pred.seq_id} · {pred.prediction}',
            series=series,
            sequence=(seq[:window] if seq else None),
            cs_before=pred.cs_before if pred.cs_before and pred.cs_before <= window else None,
            subtitle='MOCK — install SignalP 6.0 CLI for authentic model plots',
        )
        pred.plot_path = str(dest)


def run_signalp(
    sequences: dict[str, str],
    *,
    organism: str = 'eukarya',
    job_id: str | None = None,
    want_plots: bool = True,
) -> dict[str, SignalPPrediction]:
    """
    Run SignalP 6.0 fast on a batch of sequences.

    sequences: map of safe FASTA ids -> amino-acid sequence
    """
    if not sequences:
        return {}

    cleaned: dict[str, str] = {}
    id_map: dict[str, str] = {}
    for i, (orig, seq) in enumerate(sequences.items()):
        safe = re.sub(r'[^A-Za-z0-9_\-]', '_', orig)[:80] or f'seq_{i}'
        base = safe
        n = 1
        while safe in cleaned:
            safe = f'{base}_{n}'
            n += 1
        cleaned[safe] = ''.join(c for c in seq.upper() if c.isalpha())
        id_map[safe] = orig

    if not signalp_available():
        preds = _mock_predict(cleaned)
        remapped = {}
        seq_by_orig = {}
        for k, v in preds.items():
            orig = id_map.get(k, k)
            v.seq_id = orig
            remapped[orig] = v
            seq_by_orig[orig] = cleaned[k]
        if want_plots:
            _attach_mock_plots(remapped, job_id=job_id, sequences=seq_by_orig)
        return remapped

    job = job_id or uuid.uuid4().hex
    work = signalp_cache_dir() / job
    work.mkdir(parents=True, exist_ok=True)
    fasta = work / 'input.fasta'
    out_dir = work / 'out'
    out_dir.mkdir(exist_ok=True)

    with fasta.open('w', encoding='utf-8') as fh:
        for sid, seq in cleaned.items():
            fh.write(f'>{sid}\n{seq}\n')

    cmd = [
        signalp_bin(),
        '--fastafile', str(fasta),
        '--organism', 'eukarya' if organism == 'eukarya' else 'other',
        '--output_dir', str(out_dir),
        '--format', 'all' if want_plots else 'txt',
        '--mode', 'fast',
    ]
    model_dir = signalp_model_dir()
    if model_dir:
        cmd.extend(['--model_dir', model_dir])

    try:
        proc = subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True,
            timeout=int(getattr(settings, 'SIGNALP_TIMEOUT_SEC', 600)),
        )
    except subprocess.CalledProcessError as exc:
        detail = (exc.stderr or exc.stdout or '').strip()
        if 'signalp6-register' in detail or 'has not been installed yet' in detail:
            raise RuntimeError(
                'SignalP 6.0 conda stub is not registered with the licensed package. '
                'Download signalp-6.0h.fast.tar.gz from DTU Health Tech and run: '
                'signalp6-register /path/to/signalp-6.0h.fast.tar.gz '
                f'(details: {detail[:500]})'
            ) from exc
        raise RuntimeError(
            f'SignalP 6.0 failed (exit {exc.returncode}): {detail[:800] or exc}'
        ) from exc
    except (subprocess.TimeoutExpired, FileNotFoundError) as exc:
        raise RuntimeError(f'SignalP 6.0 failed: {exc}') from exc

    results_path = out_dir / 'prediction_results.txt'
    if not results_path.is_file():
        # Some versions nest outputs
        candidates = list(out_dir.rglob('prediction_results.txt'))
        if not candidates:
            raise RuntimeError('SignalP finished but prediction_results.txt was not found')
        results_path = candidates[0]
        out_dir = results_path.parent

    preds = parse_prediction_results(results_path)
    regions = parse_region_gff3(out_dir / 'region_output.gff3')
    if not regions:
        regions = parse_region_gff3(out_dir / 'output.gff3')
    _attach_regions(preds, regions)
    if want_plots:
        _collect_plots(out_dir, preds, work / 'plots')

    # Remap to original ids
    remapped: dict[str, SignalPPrediction] = {}
    for sid, pred in preds.items():
        orig = id_map.get(sid, sid)
        pred.seq_id = orig
        remapped[orig] = pred
    return remapped


def prediction_to_dict(pred: SignalPPrediction) -> dict[str, Any]:
    return {
        'seq_id': pred.seq_id,
        'prediction': pred.prediction,
        'probabilities': pred.probabilities,
        'sp_prob': pred.sp_prob,
        'cs_before': pred.cs_before,
        'cs_after': pred.cs_after,
        'cs_prob': pred.cs_prob,
        'sp_start': pred.sp_start,
        'sp_end': pred.sp_end,
        'regions': pred.regions,
        'plot_path': pred.plot_path,
        'plot_url': pred.plot_url,
    }


def compute_deltas(wt: SignalPPrediction, mut: SignalPPrediction) -> dict[str, Any]:
    """Delta metrics relative to WT-predicted class probability."""
    wt_class = wt.prediction
    wt_class_prob = wt.class_prob(wt_class) if wt_class else None
    mut_same_class_prob = mut.class_prob(wt_class) if wt_class else None
    delta_wt_class = None
    if wt_class_prob is not None and mut_same_class_prob is not None:
        delta_wt_class = mut_same_class_prob - wt_class_prob

    mut_pred_prob = mut.sp_prob
    delta_pred_prob = None
    if wt.sp_prob is not None and mut_pred_prob is not None:
        delta_pred_prob = mut_pred_prob - wt.sp_prob

    wt_had_sp = wt.prediction not in (None, '', 'OTHER')
    sp_lost = bool(wt_had_sp and mut.prediction == 'OTHER')

    cs_delta = None
    if wt.cs_before is not None and mut.cs_before is not None:
        cs_delta = mut.cs_before - wt.cs_before

    return {
        'delta_wt_class_prob': delta_wt_class,
        'delta_sp_prob': delta_pred_prob,
        'sp_lost': sp_lost,
        'cs_delta': cs_delta,
        'wt_class': wt_class,
        'mut_prob_of_wt_class': mut_same_class_prob,
    }
