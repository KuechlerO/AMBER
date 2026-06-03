# for excel:
from openpyxl import Workbook
from openpyxl.styles import Font
from openpyxl.styles import Font, Alignment, PatternFill

from django.shortcuts import render, redirect
from django.http import HttpResponse
import csv
from .services import run_analysis, UserInputError # for errors
from .result_store import save_analysis_results, get_full_rows, get_form_data, load_analysis_results
from .pipeline import apply_duplicate_mode, normalize_duplicate_mode
import math



DEFAULT_COLUMNS = [
    'position', 'wt_codon', 'wt_aa', 'mut_aa',
    'outcomes', 'avg_alpha_score', 'alpha_score',
    'sgrna_seq', 'protospacer',
    'pam', 'strand', 'target_position',
]

COLUMN_LABELS = {
    'position': 'Position',
    'wt_codon': 'WT Codon',
    'wt_aa': 'WT AA',
    'mut_aa': 'Mutated AA',
    'outcomes': 'Guide RNA Outcomes',
    'avg_alpha_score': 'avg. AlphaMissense Score',
    'alpha_score': 'AlphaMissense Score',
    'editor_used': 'Editor',
    'sgrna_seq': 'Guide RNA',
    'protospacer': 'Protospacer',
    'pam': 'PAM',
    'strand': 'Strand',
    'target_position': 'Target Position',
}

DISPLAY_COLUMNS_MAP = [
    ('position', 'Position'),
    ('wt_codon', 'WT Codon'),
    ('wt_aa', 'WT AA'),
    ('mut_aa', 'Mutated AA'),
    ('outcomes', 'Outcomes'),
    ('avg_alpha_score', 'avg. AlphaMissense Score'),
    ('alpha_score', 'AlphaMissense Score'),
    ('editor_used', 'Editor'),
    ('sgrna_seq', 'Guide RNA'),
    ('protospacer', 'Protospacer'),
    ('pam', 'PAM'),
    ('strand', 'Strand'),
    ('target_position', 'Target Pos.'),
]

EXPORT_COLUMN_MAPPING = {
    'position': ('Position', 'position'),
    'wt_codon': ('WT Codon', 'wt_codon'),
    'wt_aa': ('WT AA', 'wt_aa'),
    'mut_aa': ('Mutated AA', 'mut_aa'),
    'alpha_score': ('AlphaMissense Score', 'alpha_score'),
    'editor_used': ('Editor', 'editor_used'),
    'sgrna_seq': ('sgRNA Sequence', 'sgrna_seq'),
    'protospacer': ('Protospacer', 'protospacer'),
    'pam': ('PAM', 'pam'),
    'strand': ('Strand', 'strand'),
    'target_position': ('Target Position', 'target_position'),
    'outcomes': ('Guide RNA Outcomes', 'outcomes'),
    'avg_alpha_score': ('avg. AlphaMissense Score', 'avg_alpha_score'),
}


def parse_editing_window(post):
    window_min_raw = post.get('window_min')
    window_max_raw = post.get('window_max')
    try:
        window_min = int(window_min_raw) if window_min_raw not in [None, '', 'None'] else None
        window_max = int(window_max_raw) if window_max_raw not in [None, '', 'None'] else None
    except (TypeError, ValueError):
        window_min = None
        window_max = None
    return window_min, window_max


def resolve_score_filter(post, default_cutoff='0.5'):
    """Return (score_filter_mode, display_score_cutoff) for results display."""
    mode = post.get('score_filter_mode')
    cutoff_raw = post.get('display_score_cutoff')

    if mode is None:
        legacy = post.get('filter_threshold')
        if legacy == 'true':
            mode = 'above'
        elif legacy == 'false':
            mode = 'all'
        else:
            output_mode = post.get('output')
            if output_mode == 'only_over_threshold':
                mode = 'above'
            elif output_mode == 'all':
                mode = 'all'
            else:
                mode = 'above'

    if mode not in ('all', 'above', 'below'):
        mode = 'all'

    if cutoff_raw in [None, '']:
        cutoff_raw = post.get('alpha_threshold', default_cutoff)
    try:
        cutoff = float(cutoff_raw)
    except (TypeError, ValueError):
        cutoff = float(default_cutoff)

    return mode, cutoff


def resolve_selected_columns(post, editor):
    selected_columns = post.getlist('columns')
    if not selected_columns:
        selected_columns = list(DEFAULT_COLUMNS)
    if editor == 'BOTH' and 'editor_used' not in selected_columns:
        if 'alpha_score' in selected_columns:
            idx = selected_columns.index('alpha_score') + 1
            selected_columns.insert(idx, 'editor_used')
        else:
            selected_columns.append('editor_used')
    return selected_columns


def _row_has_valid_score(row):
    score = row.get('alpha_score')
    if score is None:
        return False
    if isinstance(score, float) and math.isnan(score):
        return False
    return True


def filter_rows_by_score(rows, score_filter_mode, display_score_cutoff):
    if score_filter_mode == 'all':
        return rows
    try:
        cutoff = float(display_score_cutoff)
    except (TypeError, ValueError):
        return rows

    if score_filter_mode == 'below':
        return [
            row for row in rows
            if _row_has_valid_score(row) and float(row['alpha_score']) < cutoff
        ]

    # default: 'above'
    return [
        row for row in rows
        if _row_has_valid_score(row) and float(row['alpha_score']) >= cutoff
    ]


def build_form_data(**kwargs):
    return kwargs


def _build_results_context(result_rows, form_data):
    unique_positions = sorted({row['position'] for row in result_rows}) if result_rows else []
    return {
        'results': result_rows,
        'position_count': len(unique_positions),
        'guide_count': len(result_rows),
        'display_columns': [COLUMN_LABELS.get(col, col) for col in form_data.get('selected_columns', [])],
        'display_columns_map': DISPLAY_COLUMNS_MAP,
        'form_data': form_data,
    }


def _rows_for_display(full_rows, form_data, post=None):
    """Apply duplicate handling then score filter for the results table."""
    post_data = post if post is not None else form_data
    duplicate_mode = normalize_duplicate_mode(
        post_data.get('duplicate_mode', form_data.get('duplicate_mode', 'best'))
    )
    rows = apply_duplicate_mode(full_rows, duplicate_mode)
    score_filter_mode, display_score_cutoff = resolve_score_filter(
        post_data,
        default_cutoff=form_data.get('alpha_threshold', '0.5'),
    )
    return filter_rows_by_score(rows, score_filter_mode, display_score_cutoff)


def _render_results(request, result_rows, full_rows, form_data):
    save_analysis_results(
        request,
        full_rows=full_rows,
        filtered_rows=result_rows,
        form_data=form_data,
    )
    return render(request, 'designer/results.html', _build_results_context(result_rows, form_data))


def _refresh_results_from_cache(request):
    cached = load_analysis_results(request)
    full_rows = cached.get('full_rows') or cached.get('raw_rows')
    if not full_rows:
        return redirect('home')

    form_data = dict(cached.get('form_data') or {})
    editor = request.POST.get('editor', form_data.get('editor', ''))
    alpha_threshold = request.POST.get('alpha_threshold', form_data.get('alpha_threshold', '0.5'))
    score_filter_mode, display_score_cutoff = resolve_score_filter(
        request.POST, default_cutoff=form_data.get('alpha_threshold', alpha_threshold)
    )
    selected_columns = resolve_selected_columns(request.POST, editor)
    editor_display = 'ABE & CBE' if editor == 'BOTH' else editor

    result_rows = _rows_for_display(full_rows, form_data, post=request.POST)

    form_data.update({
        'editor': editor,
        'editor_display': editor_display,
        'alpha_threshold': alpha_threshold,
        'selected_columns': selected_columns,
        'score_filter_mode': score_filter_mode,
        'display_score_cutoff': display_score_cutoff,
        'filter_threshold': 'true' if score_filter_mode == 'above' else 'false',
    })

    return _render_results(request, result_rows, full_rows, form_data)


def _rows_for_export(request):
    form_data = get_form_data(request)
    full_rows = get_full_rows(request)
    if not full_rows:
        return [], form_data
    duplicate_mode = normalize_duplicate_mode(form_data.get('duplicate_mode', 'best'))
    rows = apply_duplicate_mode(full_rows, duplicate_mode)
    rows = apply_threshold_filter(rows, form_data)
    return rows, form_data


def home(request):
    return render(request, 'designer/home.html')

def loading(request):
    return render(request, 'designer/loading.html')

def results(request):
    if request.method != 'POST':
        return redirect('home')

    if request.POST.get('refresh_only') == '1':
        return _refresh_results_from_cache(request)

    uniprot_id = request.POST.get('uniprot_id', '')
    editor = request.POST.get('editor', '')
    alpha_threshold = request.POST.get('alpha_threshold', '0.5')
    top_sgrnas = request.POST.get('top_sgrnas', '')
    duplicate_mode = normalize_duplicate_mode(request.POST.get('duplicate_mode', 'best'))
    pam_type = request.POST.get('pam_type', 'NGG')

    score_filter_mode, display_score_cutoff = resolve_score_filter(
        request.POST, default_cutoff=alpha_threshold
    )
    window_min, window_max = parse_editing_window(request.POST)
    selected_columns = resolve_selected_columns(request.POST, editor)
    editor_display = 'ABE & CBE' if editor == 'BOTH' else editor

    try:
        full_rows = run_analysis(
            uniprot_id=uniprot_id,
            editor=editor,
            alpha_threshold=alpha_threshold,
            top_sgrnas=top_sgrnas,
            window_min=window_min,
            window_max=window_max,
            pam_type=pam_type,
        )
        form_data = build_form_data(
            uniprot_id=uniprot_id,
            editor=editor,
            editor_display=editor_display,
            alpha_threshold=alpha_threshold,
            top_sgrnas=top_sgrnas,
            gene_name=uniprot_id,
            selected_columns=selected_columns,
            score_filter_mode=score_filter_mode,
            display_score_cutoff=display_score_cutoff,
            filter_threshold='true' if score_filter_mode == 'above' else 'false',
            window_min=window_min,
            window_max=window_max,
            duplicate_mode=duplicate_mode,
            pam_type=pam_type,
        )
        result_rows = _rows_for_display(full_rows, form_data, post=request.POST)

        return _render_results(request, result_rows, full_rows, form_data)

    except UserInputError as e:
        return render(request, 'designer/error.html', {'error': str(e), 'error_type': 'user'})
    except Exception as e:
        import traceback
        print(traceback.format_exc())
        return render(request, 'designer/error.html', {'error': str(e), 'error_type': 'server'})


def apply_threshold_filter(rows, form_data):
    mode = form_data.get('score_filter_mode')
    if mode is None:
        mode = 'above' if form_data.get('filter_threshold') == 'true' else 'all'
    cutoff = form_data.get(
        'display_score_cutoff',
        form_data.get('alpha_threshold', '0.5'),
    )
    return filter_rows_by_score(rows, mode, cutoff)

def download_excel(request):
    rows, form_data = _rows_for_export(request)
    duplicate_mode = normalize_duplicate_mode(form_data.get('duplicate_mode', 'best'))

    gene_name = form_data.get('gene_name') or form_data.get('uniprot_id', 'Unknown target')
    editor_display = form_data.get('editor_display', form_data.get('editor', ''))

    output_mode = form_data.get('output_mode', 'all')
    alpha_threshold = form_data.get('alpha_threshold')

    # Excel Workbook
    wb = Workbook()
    ws = wb.active
    ws.title = 'Results'

    row_idx = 1


    # ===== titel =====
    title = f"CRISPR Base Editing Designer Results for {gene_name}"
    ws.cell(row=row_idx, column=1, value=title)
    ws.cell(row=row_idx, column=1).font = Font(bold=True)
    row_idx += 2

    # ===== Tabelle 1: Zusammenfassung =====
    headers_1 = ['UniProt ID', 'Editor', 'PAM', 'AlphaMissense Threshold', 'Top Guide RNAs']
    values_1 = [
        form_data.get('uniprot_id', ''),
        editor_display,
        form_data.get('pam_type', 'NGG'),
        form_data.get('alpha_threshold', ''),
        form_data.get('top_sgrnas', ''),
    ]

    for col, header in enumerate(headers_1, start=1):
        ws.cell(row=row_idx, column=col, value=header).font = Font(bold=True)

    row_idx += 1

    for col, value in enumerate(values_1, start=1):
        ws.cell(row=row_idx, column=col, value=value)

    row_idx += 3

    #tabel 2: results
    column_mapping = {
        'position': ('Position', 'position'),
        'wt_codon': ('WT Codon', 'wt_codon'),
        'wt_aa': ('WT AA', 'wt_aa'),
        'mut_aa': ('Mutated AA', 'mut_aa'),
        'alpha_score': ('AlphaMissense Score', 'alpha_score'),
        'editor_used': ('Editor', 'editor_used'),
        'sgrna_seq': ('sgRNA Sequence', 'sgrna_seq'),
        'protospacer': ('Protospacer', 'protospacer'),
        'pam': ('PAM', 'pam'),
        'strand': ('Strand', 'strand'),
        'target_position': ('Target Position', 'target_position'),
        'outcomes': ('Guide RNA Outcomes', 'outcomes'),
        'avg_alpha_score': ('avg. AlphaMissense Score','avg_alpha_score'),
    }

    selected_columns = form_data.get('selected_columns', [])

    # if nothing selected -> show all
    if not selected_columns:
        selected_columns = list(column_mapping.keys())

    # Header
    for col_idx, col_key in enumerate(selected_columns, start=1):
        header = column_mapping[col_key][0]
        ws.cell(row=row_idx, column=col_idx, value=header).font = Font(bold=True)

    row_idx += 1

    # colors for groups of occurrences of same sgRNA
    group_colors = [
        '86EFAC',  # soft green
        '93C5FD',  # soft blue
        'D8B4FE',  # soft purple
        'FDB474',  # soft orange
        'F9A8D4',  # soft pink
    ]

    # counts every sgRNA
    from collections import Counter
    seq_counts = Counter(r.get('sgrna_seq') for r in rows)

    # map color to same sgRNAs
    seq_color_map = {}
    color_index = 0
    for r in rows:
        seq = r.get('sgrna_seq')
        if seq and seq_counts[seq] > 1 and seq not in seq_color_map:
            seq_color_map[seq] = group_colors[color_index % len(group_colors)]
            color_index += 1

    for r in rows:
        seq = r.get('sgrna_seq')
        row_color = seq_color_map.get(seq) if duplicate_mode == 'group' else None

        for col_idx, col_key in enumerate(selected_columns, start=1):
            _, data_key = column_mapping[col_key]
            value = r.get(data_key, '')

            if data_key == 'target_position' and isinstance(value, list):
                value = ', '.join(str(x) for x in value)

            if data_key == 'outcomes' and isinstance(value, list):
                value = '\n'.join(value)

            cell = ws.cell(row=row_idx, column=col_idx, value=value)

            if data_key == 'outcomes':
                cell.alignment = Alignment(wrap_text=True)

            if row_color:
                cell.fill = PatternFill(
                    start_color=row_color,
                    end_color=row_color,
                    fill_type='solid'
                )

        row_idx += 1

    for i in range(1, len(selected_columns) + 1):
        ws.column_dimensions[chr(64 + i)].width = 20

    # File:
    filename = f"{gene_name}_crispr_results.xlsx"

    response = HttpResponse(
        content_type='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'
    )
    response['Content-Disposition'] = f'attachment; filename="{filename}"'

    wb.save(response)

    return response


def download_csv(request):
    rows, form_data = _rows_for_export(request)

    output_mode = form_data.get('output_mode', 'all')
    alpha_threshold = form_data.get('alpha_threshold')


    column_mapping = {
        'position': ('Position', 'position'),
        'wt_codon': ('WT Codon', 'wt_codon'),
        'wt_aa': ('WT AA', 'wt_aa'),
        'mut_aa': ('Mutated AA', 'mut_aa'),
        'alpha_score': ('AlphaMissense Score', 'alpha_score'),
        'editor_used': ('Editor', 'editor_used'),
        'sgrna_seq': ('sgRNA Sequence', 'sgrna_seq'),
        'protospacer': ('Protospacer', 'protospacer'),
        'pam': ('PAM', 'pam'),
        'strand': ('Strand', 'strand'),
        'target_position': ('Target Position', 'target_position'),
        'outcomes': ('Guide RNA Outcomes', 'outcomes'),
        'avg_alpha_score': ('avg. AlphaMissense Score','avg_alpha_score'),
    }

    selected_columns = form_data.get('selected_columns', [])
    if not selected_columns:
        selected_columns = list(column_mapping.keys())


    gene_name = form_data.get('gene_name') or form_data.get('uniprot_id', 'Unknown target')

    filename = f"{gene_name}_crispr_results.csv"
    response = HttpResponse(content_type='text/csv')
    response['Content-Disposition'] = f'attachment; filename="{filename}"'

    writer = csv.writer(response)

    # header
    headers = [column_mapping[col][0] for col in selected_columns]
    writer.writerow(headers)

    # data
    for r in rows:
        csv_row = []
        for col_key in selected_columns:
            _, data_key = column_mapping[col_key]
            value = r.get(data_key, '')

            if data_key == 'target_position':
                if isinstance(value, list):
                    value = 'pos ' + ' / '.join(str(x) for x in value)
                elif value is None:
                    value = ''
                else:
                    value = f'pos {value}'

            if data_key == 'outcomes' and isinstance(value, list):
                value = ' | '.join(value)

            csv_row.append(value)

        writer.writerow(csv_row)

    return response


def tutorial(request):
    return render(request, 'designer/tutorial.html')

def about(request):
    return render(request, 'designer/about.html')
