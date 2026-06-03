# SWP2026-AMBER Designer

This is the repository of the **AMBER Designer** project (SWP 2026 @BTG) – the **A**lpha**M**issense **B**ase **E**diting **R**NA Designer.

It contains a Django-based web application for designing base-editing-compatible guide RNAs and ranking them with Alpha Missense scores for target proteins.

---

## Description

The application allows users to:

- Enter a **UniProt ID** or **Ensembl transcript ID**
- Choose a base editor: **ABE** (A→G), **CBE** (C→T), or **both**
- Set an **AlphaMissense pathogenicity threshold** to filter variants
- Define the **editing window** (default: positions 4–8)
- Select the **number of top guide RNA candidates** per position
- Choose how to handle **duplicate guide RNAs** (show once or group all occurrences)
- Decide which columns should be shown in the results table
- Decide if the results table should show all values or only the ones over the choosen threshold

The tool returns ranked guide RNA candidates with protospacer, PAM, strand, target position, guide RNA outcomes annotated with AlphaMissense scores and their average Alpha Missense scores.

Results can be **downloaded as CSV or Excel**. If the results table is sorted by a column in the browser, only the **CSV download** preserves that sort order – the Excel download always uses the default order, which is by position.

Guide RNAs that span **exon-intron boundaries** are automatically excluded.

AlphaMissense scores are retrieved from a **PostGreSQL database**.

---

## Repository Structure

* `manage.py` – Django entry point
* `crispr_webapp/` – Django project configuration (settings, urls, wsgi)
* `designer/` – main Django application
* `designer/static/designer/` – static CSS files
* `designer/templates/designer/` – HTML templates (base, home, results, loading, tutorial, error, about)
* `designer/views.py` – request handling and result rendering
* `designer/services.py` – input validation and pipeline entry point
* `designer/pipeline.py` – core analysis pipeline (guide RNA search, exon filtering, outcome annotation)
* `designer/models.py` – AlphaMissense database model

---

> **For Charité users:** The application is accessible via the Charité intranet – no local installation required.
> The setup instructions below are intended for developers working on the codebase.

---

## Requirements

Make sure you have **Conda** installed before setting up the project.

---

## Setup

Clone the repository and create the Conda environment:

```bash
git clone https://github.com/KuechlerO/SWP2026-CRISPR-Tool.git
cd SWP2026-CRISPR-Tool
conda env create -f environment.yml
conda activate crispr_tool_env
```

---

## Run the Application

Start the Django development server:

```bash
python manage.py runserver
```

Then open the application in your browser at:

```
http://127.0.0.1:8000/
```

---

## How It Works

1. The user enters a UniProt or Ensembl ID and selects analysis parameters
2. The CDS is taken from the **Ensembl REST API**
3. AlphaMissense pathogenicity scores are pulled from a **PostGreSQL database**
4. For each position guide RNAs are searched on both strands
5. Guide RNAs spanning exon-intron boundaries are excluded (via Ensembl mapping)
6. Possible codon outcomes are enumerated and annotated with AlphaMissense scores
7. Results are sorted by position (ascending) by default and can be re-sorted in the browser

---

## Notes

- The AlphaMissense data is served from a **PostGreSQL database** (`legacy_db` in Django settings) – not from local files
- The Django project configuration is in `crispr_webapp/`
- The main analysis logic is in `designer/pipeline.py`
