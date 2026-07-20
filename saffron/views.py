"""SAFFRON views: home → loading → results, exports, SignalP plot serving."""

from __future__ import annotations

import csv
import io
from pathlib import Path

from django.http import FileResponse, Http404, HttpResponse
from django.shortcuts import redirect, render
from django.views.decorators.http import require_GET
from openpyxl import Workbook

from designer.services import UserInputError

from . import result_store
from .services import run_saffron_analysis
from .signalp_client import signalp_cache_dir, signalp_status


def home(request):
    return render(request, 'saffron/home.html', {'signalp': signalp_status()})


def loading(request):
    if request.method != 'POST':
        return redirect('saffron_home')
    return render(request, 'saffron/loading.html', {'request': request})


def _form_from_post(post) -> dict:
    return {
        'input_mode': post.get('input_mode', 'uniprot'),
        'uniprot_id': post.get('uniprot_id', ''),
        'aa_sequence': post.get('aa_sequence', ''),
        'mutation': post.get('mutation', ''),
        'organism': post.get('organism', 'eukarya'),
        'editor': post.get('editor', 'BOTH'),
        'pam_type': post.get('pam_type', 'NGG'),
        'window_min': post.get('window_min', '4'),
        'window_max': post.get('window_max', '8'),
        'top_sgrnas': post.get('top_sgrnas', '5'),
    }


def results(request):
    if request.method != 'POST':
        data = result_store.load_analysis_results(request)
        if not data:
            return redirect('saffron_home')
        return render(request, 'saffron/results.html', _results_context(data))

    form = _form_from_post(request.POST)
    try:
        payload = run_saffron_analysis(form)
    except UserInputError as exc:
        return render(request, 'saffron/error.html', {'message': str(exc)}, status=400)
    except Exception as exc:  # noqa: BLE001 — show friendly page for unexpected SignalP/runtime errors
        return render(
            request,
            'saffron/error.html',
            {'message': f'Analysis failed: {exc}'},
            status=500,
        )

    payload['form_data'] = form
    result_store.save_analysis_results(request, payload=payload)
    return render(request, 'saffron/results.html', _results_context(payload))


def _fmt(val, digits=3):
    if val is None:
        return '—'
    if isinstance(val, float):
        return f'{val:.{digits}f}'
    return val


def _basename(path: str | None) -> str | None:
    if not path:
        return None
    return Path(path).name


def _plot_ref(path: str | None) -> str | None:
    """Prefer cache-relative path so plot URLs stay unique across jobs."""
    if not path:
        return None
    p = Path(path).resolve()
    try:
        root = signalp_cache_dir().resolve()
        rel = p.relative_to(root)
        # Keep URL path portable (forward slashes)
        return rel.as_posix()
    except (ValueError, OSError):
        return p.name


def _results_context(data: dict) -> dict:
    wt = dict(data.get('wt_signalp') or {})
    wt['plot_file'] = _plot_ref(wt.get('plot_path'))
    rows = []
    patho_count = 0
    for raw in data.get('guide_rows') or []:
        r = dict(raw)
        r['sp_prob_display'] = _fmt(r.get('sp_prob'))
        r['delta_display'] = _fmt(r.get('delta_wt_class_prob'))
        r['cs_display'] = (
            f"{r['cs_before']}-{r.get('cs_after')}" if r.get('cs_before') else '—'
        )
        r['cs_delta_display'] = _fmt(r.get('cs_delta'), digits=0) if r.get('cs_delta') is not None else '—'
        r['plot_file'] = _plot_ref(r.get('plot_path'))
        clin = r.get('clin_sig')
        r['clin_sig_display'] = clin if clin else '—'
        if r.get('paper_pathogenic'):
            r['paper_pathogenic_display'] = 'Yes'
            patho_count += 1
        elif r.get('in_patho_catalogue'):
            r['paper_pathogenic_display'] = 'Other'
        else:
            r['paper_pathogenic_display'] = '—'
        rows.append(r)

    focus = data.get('focus_mutation')
    if focus:
        focus = dict(focus)
        focus['plot_file'] = _plot_ref(focus.get('plot_path'))
        clin = focus.get('clin_sig')
        focus['clin_sig_display'] = clin if clin else '—'
        if focus.get('paper_pathogenic'):
            focus['paper_pathogenic_display'] = 'Yes'
        elif focus.get('in_patho_catalogue'):
            focus['paper_pathogenic_display'] = 'Other'
        else:
            focus['paper_pathogenic_display'] = '—'

    return {
        'data': data,
        'form_data': data.get('form_data') or data.get('form_meta') or {},
        'wt': wt,
        'guide_rows': rows,
        'guide_count': len(rows),
        'patho_count': patho_count,
        'focus': focus,
        'sp_span': data.get('sp_span'),
        'mode': data.get('mode'),
        'message': data.get('message') or '',
        'no_sp': data.get('no_sp', False),
        'guides_available': data.get('guides_available', False),
        'signalp': signalp_status(),
        'signalp_backend': data.get('signalp_backend') or signalp_status()['backend'],
    }



@require_GET
def plot_file(request, job_file: str):
    """Serve a SignalP plot PNG from the SAFFRON cache directory."""
    if '..' in job_file or job_file.startswith(('/', '\\')):
        raise Http404('Not found')
    root = signalp_cache_dir().resolve()
    candidate = (root / job_file).resolve()
    path = None
    if candidate.is_file() and candidate.suffix.lower() == '.png':
        try:
            candidate.relative_to(root)
            path = candidate
        except ValueError:
            path = None
    if path is None:
        # Backward-compatible: basename search, prefer newest (avoids stale mock plots)
        safe = Path(job_file).name
        if not safe.endswith('.png'):
            raise Http404('Not found')
        matches = [p for p in root.rglob(safe) if p.is_file()]
        if not matches:
            raise Http404('Plot not found')
        path = max(matches, key=lambda p: p.stat().st_mtime).resolve()
        try:
            path.relative_to(root)
        except ValueError:
            raise Http404('Invalid path') from None
    return FileResponse(path.open('rb'), content_type='image/png')


def _export_rows(data: dict) -> list[dict]:
    return data.get('guide_rows') or []


@require_GET
def download_csv(request):
    data = result_store.load_analysis_results(request)
    if not data:
        return redirect('saffron_home')
    rows = _export_rows(data)
    buf = io.StringIO()
    fields = [
        'position', 'wt_aa', 'mut_aa', 'editor_used', 'sgrna_seq', 'protospacer',
        'pam', 'strand', 'sp_prediction', 'sp_prob', 'delta_wt_class_prob',
        'cs_before', 'cs_prob', 'cs_delta', 'sp_lost', 'highlighted',
        'clin_sig', 'paper_pathogenic', 'in_patho_catalogue',
    ]
    writer = csv.DictWriter(buf, fieldnames=fields, extrasaction='ignore')
    writer.writeheader()
    for r in rows:
        writer.writerow(r)
    resp = HttpResponse(buf.getvalue(), content_type='text/csv')
    resp['Content-Disposition'] = 'attachment; filename="saffron_results.csv"'
    return resp


@require_GET
def download_excel(request):
    data = result_store.load_analysis_results(request)
    if not data:
        return redirect('saffron_home')
    rows = _export_rows(data)
    wb = Workbook()
    ws = wb.active
    ws.title = 'SAFFRON'
    headers = [
        'Position', 'WT AA', 'Mut AA', 'Editor', 'sgRNA', 'Protospacer', 'PAM', 'Strand',
        'SP prediction', 'SP prob', 'Δ WT-class prob', 'CS before', 'CS prob', 'CS Δ', 'SP lost',
        'ClinVar/VEP CLIN_SIG', 'Paper potentially pathogenic',
    ]
    ws.append(headers)
    for r in rows:
        ws.append([
            r.get('position'), r.get('wt_aa'), r.get('mut_aa'), r.get('editor_used'),
            r.get('sgrna_seq'), r.get('protospacer'), r.get('pam'), r.get('strand'),
            r.get('sp_prediction'), r.get('sp_prob'), r.get('delta_wt_class_prob'),
            r.get('cs_before'), r.get('cs_prob'), r.get('cs_delta'), r.get('sp_lost'),
            r.get('clin_sig') or '',
            'Yes' if r.get('paper_pathogenic') else ('Other' if r.get('in_patho_catalogue') else ''),
        ])
    out = io.BytesIO()
    wb.save(out)
    out.seek(0)
    resp = HttpResponse(
        out.getvalue(),
        content_type='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet',
    )
    resp['Content-Disposition'] = 'attachment; filename="saffron_results.xlsx"'
    return resp


def about(request):
    return render(request, 'saffron/about.html')
