"""Lazy-loaded screen data (Figures.ipynb NGG pipeline)."""

from __future__ import annotations

import json
import threading
import time
from pathlib import Path

import numpy as np
import pandas as pd
from django.conf import settings

from .gene_lookup import SCREEN_MARKERS
from .helpers import combine_pvals, ext_loc_splice, is_KO, is_syn


def _screen_data_root() -> Path:
    return Path(getattr(settings, 'SCREEN_DATA_DIR', settings.BASE_DIR / 'files-archive-dir'))


class ScreenDataStore:
    """Thread-safe singleton; loads tables once on first plot request."""

    _instance = None
    _lock = threading.Lock()

    def __init__(self):
        self._loaded = False
        self._loading = False
        self._load_error: str | None = None
        self._background_started = False
        self._domains_only_loaded = False
        self._domains_lock = threading.Lock()
        self._load_lock = threading.Lock()
        self.load_stage: str = 'idle'
        self._load_started_at: float | None = None
        self.root = _screen_data_root()
        self.protein_length: dict[str, int] = {}
        self.data_cutoff: pd.DataFrame | None = None
        self.data_grouped: pd.DataFrame | None = None
        self.domains: dict = {}
        self.domains_avg: pd.DataFrame | None = None

    @classmethod
    def get(cls) -> 'ScreenDataStore':
        if cls._instance is None:
            with cls._lock:
                if cls._instance is None:
                    cls._instance = cls()
        return cls._instance

    @classmethod
    def reset_for_tests(cls):
        with cls._lock:
            cls._instance = None

    def is_ready(self) -> bool:
        return self._loaded

    def is_loading(self) -> bool:
        return self._loading and not self._loaded

    def load_elapsed_sec(self) -> int | None:
        if self._load_started_at is None:
            return None
        return int(time.time() - self._load_started_at)

    def get_load_error(self) -> str | None:
        return self._load_error

    def start_background_load(self) -> bool:
        """Start loading screen tables in a daemon thread. Returns True if a new load was started."""
        if self._loaded:
            return False
        with self._lock:
            if self._loaded or self._background_started:
                return False
            self._background_started = True
            self._loading = True
            self._load_started_at = time.time()
            self.load_stage = 'table1'
        thread = threading.Thread(target=self._load_in_background, daemon=True)
        thread.start()
        return True

    def _load_in_background(self):
        try:
            self.ensure_loaded()
        except Exception as exc:
            self._load_error = str(exc)
            self.load_stage = 'error'
            with self._lock:
                self._background_started = False
        finally:
            if not self._loaded:
                self._loading = False

    def ensure_domains_only(self):
        """Load protein lengths and UniProt domain JSON (fast path for overview track)."""
        if self._domains_only_loaded:
            return
        # Separate lock from ensure_loaded so overview is not blocked by background full load.
        with self._domains_lock:
            if self._domains_only_loaded:
                return
            root = self.root
            protein_data = pd.read_csv(root / 'Supplementary Table 1.txt', sep='\t', index_col=0)
            if 'gene' in protein_data.columns:
                self.protein_length = {
                    str(row['gene']): int(row['uniprot_len'])
                    for _, row in protein_data.iterrows()
                }
            else:
                self.protein_length = {
                    str(gene): int(row['uniprot_len'])
                    for gene, row in protein_data.iterrows()
                }
            domain_path = root / 'libraries' / 'domain_data.json'
            if domain_path.is_file():
                with domain_path.open(encoding='utf-8') as fh:
                    self.domains = json.load(fh)
            self._domains_only_loaded = True

    def ensure_loaded(self):
        if self._loaded:
            return
        with self._load_lock:
            if not self._loaded:
                self._loading = True
                try:
                    self._load_all()
                    self._loaded = True
                    self._domains_only_loaded = True
                    self.load_stage = 'ready'
                except Exception:
                    self.load_stage = 'error'
                    raise
                finally:
                    self._loading = False

    def _load_all(self):
        root = self.root
        self.load_stage = 'table1'

        protein_data = pd.read_csv(root / 'Supplementary Table 1.txt', sep='\t', index_col=0)
        if 'gene' in protein_data.columns:
            self.protein_length = {
                str(row['gene']): int(row['uniprot_len'])
                for _, row in protein_data.iterrows()
            }
        else:
            self.protein_length = {
                str(gene): int(row['uniprot_len'])
                for gene, row in protein_data.iterrows()
            }

        lib_ngg = pd.read_csv(root / 'Supplementary Table 2.txt', sep='\t', index_col=0)
        lib_ngg['sites'] = lib_ngg['sites'].str.split(',')
        lib_ngg['sites_mapped'] = lib_ngg['sites_mapped'].str.split(',')

        self.load_stage = 'table3'
        data_all = pd.read_csv(root / 'Supplementary Table 3.txt', sep='\t', low_memory=False)
        data_all['control_count'] = data_all['control_count'].str.split('/')
        data_all['treatment_count'] = data_all['treatment_count'].str.split('/')
        data = data_all[data_all['screen'] == 'NGG'].copy()
        data = data.merge(lib_ngg, on=['sgrna', 'be_type'], how='left')
        drop_cols = [c for c in ['Gene', 'control_var', 'adj_var', 'p.low', 'p.high', 'high_in_treatment', 'type'] if c in data.columns]
        data = data.drop(columns=drop_cols, errors='ignore')
        data = data.fillna('None')

        control_counts = data['control_mean']
        count_cutoff = control_counts.quantile(q=0.025)
        cutoff_mask = (data['control_mean'] > count_cutoff) | (data['treat_mean'] > count_cutoff)
        self.data_cutoff = data[cutoff_mask].copy()

        self.load_stage = 'stats'
        base_level_stats = pd.read_csv(
            root / 'Supplementary Table 6 base_level_stats.txt', sep='\t'
        )
        residue_level_stats = base_level_stats.groupby(
            ['gene', 'marker', 'be_type', 'loc_base'], as_index=False
        ).agg({
            'p-value_base': list,
            'effect_size_base': list,
            'SE_base': list,
            'base': list,
            'codon_pos': list,
        }).rename({'loc_base': 'loc'}, axis=1)
        residue_level_stats['p-value_loc'] = residue_level_stats['p-value_base'].apply(combine_pvals)
        residue_level_stats['avg_effect_size_loc'] = residue_level_stats['effect_size_base'].apply(np.mean)
        residue_level_stats['best_p'] = residue_level_stats['p-value_base'].apply(min)

        data_no_KO = self.data_cutoff[
            ~self.data_cutoff['sites_mapped'].apply(
                lambda sites: any(is_KO(s) for s in sites)
            )
        ]
        data_exp = data_no_KO.explode(
            column=['treatment_count', 'control_count']
        ).explode(column='sites_mapped')
        data_exp = data_exp.copy()
        data_exp['sites'] = data_exp['sites_mapped']
        data_exp = data_exp[~data_exp['sites'].apply(is_syn)]
        data_exp['loc'] = data_exp['sites_mapped'].apply(ext_loc_splice)
        data_exp = data_exp[data_exp['loc'] != 'None']
        p_min = data_exp['p.twosided'][data_exp['p.twosided'] > 0].min()
        data_exp['p.twosided'] = data_exp['p.twosided'].replace(0.0, p_min)

        self.load_stage = 'grouping'
        data_grouped = data_exp.groupby(['gene', 'marker', 'be_type', 'loc']).agg({
            'sgrna': set,
            'LFC': 'mean',
            'control_mean': 'mean',
            'treat_mean': 'mean',
            'p.twosided': set,
        }).reset_index()
        data_grouped['p_values'] = data_grouped['p.twosided']
        data_grouped['p.twosided'] = data_grouped['p_values'].apply(lambda lst: combine_pvals(list(lst)))
        data_grouped = data_grouped.merge(
            residue_level_stats[
                ['gene', 'marker', 'be_type', 'loc', 'p-value_base', 'SE_base', 'effect_size_base', 'codon_pos']
            ],
            on=['gene', 'marker', 'be_type', 'loc'],
            how='left',
        )
        self.data_grouped = data_grouped.set_index(['gene', 'marker', 'be_type', 'loc'])

        self.load_stage = 'domains'
        domain_path = root / 'libraries' / 'domain_data.json'
        if domain_path.is_file():
            with domain_path.open(encoding='utf-8') as fh:
                self.domains = json.load(fh)
        # Skip global domains_avg precompute (~minutes); plots use structural domains only.
        self.domains_avg = None

    def _build_domains_avg(self, data_exp: pd.DataFrame) -> pd.DataFrame:
        domains_avg = pd.DataFrame(
            columns=['gene', 'marker', 'be_type', 'domain', 'avg_LFC', 'domain_type', 'start', 'end']
        )
        for gene in data_exp['gene'].unique():
            domain_list = self.domains.get(gene)
            if not domain_list:
                continue
            for domain in domain_list:
                start = int(domain['location']['start']['value'])
                end = int(domain['location']['end']['value'])
                name = domain['description']
                typ = domain['type']
                try:
                    domain_data = (
                        self.data_grouped[['LFC']]
                        .loc[gene, slice(None), slice(None), start:end]
                        .groupby(['marker', 'be_type'])
                        .mean()
                        .reset_index()[['marker', 'be_type', 'LFC']]
                        .rename({'LFC': 'avg_LFC'}, axis=1)
                    )
                except (KeyError, IndexError):
                    continue
                domain_data['domain_type'] = typ
                domain_data['gene'] = gene
                domain_data['domain'] = name
                domain_data['start'] = start
                domain_data['end'] = end
                domains_avg = pd.concat([domains_avg, domain_data], ignore_index=True)
        return domains_avg

    def get_structural_domain_layout(self, gene: str) -> list[dict]:
        """UniProt domain boundaries without screen LFC coloring."""
        self.ensure_domains_only()
        layout = []
        for domain in self.domains.get(gene) or []:
            layout.append({
                'name': str(domain['description']),
                'start': int(domain['location']['start']['value']),
                'end': int(domain['location']['end']['value']),
                'avg_lfc': None,
                'domain_type': domain.get('type', ''),
            })
        return layout

    def get_domain_layout(self, gene: str, marker: str) -> list[dict]:
        """
        UniProt domain boundaries with optional screen avg_LFC coloring for marker.

        Each entry: name, start, end, avg_lfc (float|None), domain_type.
        """
        self.ensure_loaded()
        domain_list = self.domains.get(gene) or []
        lfc_by_name: dict[str, float] = {}
        if self.domains_avg is not None and not self.domains_avg.empty:
            subset = self.domains_avg[
                (self.domains_avg['gene'] == gene) & (self.domains_avg['marker'] == marker)
            ]
            for _, row in subset.iterrows():
                lfc_by_name[str(row['domain'])] = float(row['avg_LFC'])

        layout = []
        for domain in domain_list:
            start = int(domain['location']['start']['value'])
            end = int(domain['location']['end']['value'])
            name = str(domain['description'])
            layout.append({
                'name': name,
                'start': start,
                'end': end,
                'avg_lfc': lfc_by_name.get(name),
                'domain_type': domain.get('type', ''),
            })
        return layout

    def screen_residue_locs_by_marker(self, gene: str) -> dict[str, set[int]]:
        """Residue positions with screen data per marker (any editor type)."""
        self.ensure_loaded()
        result = {marker: set() for marker in SCREEN_MARKERS}
        if self.data_grouped is None:
            return result
        try:
            gene_slice = self.data_grouped.loc[gene]
        except KeyError:
            return result
        for marker in SCREEN_MARKERS:
            try:
                marker_slice = gene_slice.loc[marker]
            except KeyError:
                continue
            locs: set[int] = set()
            for be_type in ('ABE', 'CBE'):
                try:
                    be_df = marker_slice.loc[be_type]
                except KeyError:
                    continue
                for loc in be_df.index:
                    try:
                        locs.add(int(loc))
                    except (TypeError, ValueError):
                        continue
            result[marker] = locs
        return result


def get_screen_data_store() -> ScreenDataStore:
    return ScreenDataStore.get()
