"""AlphaMissense Postgres access (no Django)."""

from __future__ import annotations

import os
from typing import Iterable

import pandas as pd
import psycopg2

DB_CONFIG_BASE = {
    'dbname': os.environ.get('AMBER_DB_NAME', 'swpdb'),
    'user': os.environ.get('AMBER_DB_USER', 'swp2026'),
    'password': os.environ.get('AMBER_DB_PASSWORD', 'swp2026'),
    'port': os.environ.get('AMBER_DB_PORT', '5432'),
    'options': '-c search_path=alpha_missense_db',
}


def connect_am_db():
    hosts = []
    env_host = os.environ.get('AMBER_DB_HOST')
    if env_host:
        hosts.append(env_host)
    hosts.extend(['127.0.0.1', 'medgen_pgsql_db'])
    last_err = None
    seen = set()
    for host in hosts:
        if host in seen:
            continue
        seen.add(host)
        try:
            return psycopg2.connect(host=host, connect_timeout=15, **DB_CONFIG_BASE)
        except psycopg2.OperationalError as exc:
            last_err = exc
    raise last_err or RuntimeError('Could not connect to AlphaMissense database')


def list_accessions() -> list[str]:
    """All distinct UniProt IDs present in the AM table."""
    sql = """
        SELECT DISTINCT uniprot_id
        FROM alpha_missense_db.alpha_missense
        WHERE uniprot_id IS NOT NULL AND uniprot_id <> ''
        ORDER BY uniprot_id
    """
    conn = connect_am_db()
    try:
        with conn.cursor() as cur:
            cur.execute(sql)
            return [str(r[0]).strip().upper() for r in cur.fetchall() if r[0]]
    finally:
        conn.close()


def fetch_variants(uniprot_ids: str | Iterable[str]) -> pd.DataFrame:
    """Load AM variants (no score threshold). Accepts one ID or a list."""
    if isinstance(uniprot_ids, str):
        ids = [uniprot_ids.strip().upper()]
    else:
        ids = [str(u).strip().upper() for u in uniprot_ids if u]
    if not ids:
        return pd.DataFrame(
            columns=['uniprot_id', 'a.a.1', 'position', 'a.a.2', 'pathogenicity score']
        )
    sql = """
        SELECT uniprot_id, protein_variant, am_pathogenicity
        FROM alpha_missense_db.alpha_missense
        WHERE uniprot_id = ANY(%s)
    """
    conn = connect_am_db()
    try:
        with conn.cursor() as cur:
            cur.execute(sql, (ids,))
            rows = cur.fetchall()
            cols = [d[0] for d in cur.description]
        df = pd.DataFrame.from_records(rows, columns=cols)
    finally:
        conn.close()
    if df.empty:
        return pd.DataFrame(
            columns=['uniprot_id', 'a.a.1', 'position', 'a.a.2', 'pathogenicity score']
        )
    parsed = df['protein_variant'].str.extract(r'([A-Z])(\d+)([A-Z])')
    parsed.columns = ['a.a.1', 'position', 'a.a.2']
    df = pd.concat([df[['uniprot_id', 'am_pathogenicity']], parsed], axis=1)
    df = df.dropna(subset=['position', 'a.a.1'])
    df['position'] = df['position'].astype(int)
    df['uniprot_id'] = df['uniprot_id'].astype(str).str.strip().str.upper()
    return df.rename(columns={'am_pathogenicity': 'pathogenicity score'})
