/**
 * Gene symbol → UniProt accession picker for AMBER / SAFFRON home forms.
 *
 * Expects elements:
 *   #uniprot_id (input)
 *   #gene-lookup-btn (button)
 *   #gene-picker-modal (dialog)
 *   form that contains #uniprot_id
 *
 * data-search-url on #gene-picker-modal supplies the JSON endpoint.
 */
(function () {
    'use strict';

    var UNIPROT_RE = /^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-\d+)?$/i;

    function normalizeId(raw) {
        return String(raw || '').replace(/\s+/g, '').toUpperCase();
    }

    function looksLikeUniprot(raw) {
        return UNIPROT_RE.test(normalizeId(raw));
    }

    function looksLikeEnsembl(raw) {
        return normalizeId(raw).indexOf('ENST') === 0;
    }

    function looksLikeGeneSymbol(raw) {
        var s = normalizeId(raw);
        if (!s || s.length < 2) return false;
        if (looksLikeUniprot(s) || looksLikeEnsembl(s)) return false;
        return /^[A-Z0-9][A-Z0-9_-]{0,31}$/i.test(s);
    }

    function escapeHtml(text) {
        var d = document.createElement('div');
        d.textContent = text == null ? '' : String(text);
        return d.innerHTML;
    }

    function isModalOpen(modal) {
        return !!(modal && (modal.open || modal.hasAttribute('open')));
    }

    function initGenePicker(options) {
        options = options || {};
        var input = document.getElementById(options.inputId || 'uniprot_id');
        var lookupBtn = document.getElementById(options.lookupBtnId || 'gene-lookup-btn');
        var modal = document.getElementById(options.modalId || 'gene-picker-modal');
        if (!input || !lookupBtn || !modal) return;

        var form = input.closest('form');
        var searchUrl = modal.getAttribute('data-search-url') || options.searchUrl || '';
        var statusEl = modal.querySelector('.gene-picker-status');
        var listEl = modal.querySelector('.gene-picker-list');
        var titleEl = modal.querySelector('.gene-picker-modal-title');
        var closeBtns = modal.querySelectorAll('[data-gene-picker-close]');
        var pendingChoice = false;

        function setStatus(msg, isError) {
            if (!statusEl) return;
            statusEl.textContent = msg || '';
            statusEl.classList.toggle('gene-picker-status-error', !!isError);
        }

        function closeModal() {
            pendingChoice = false;
            if (listEl) listEl.innerHTML = '';
            setStatus('');
            if (typeof modal.close === 'function' && isModalOpen(modal)) {
                modal.close();
            } else {
                modal.removeAttribute('open');
                modal.hidden = true;
            }
        }

        function openModal() {
            modal.hidden = false;
            if (typeof modal.showModal === 'function') {
                if (!isModalOpen(modal)) modal.showModal();
            } else {
                modal.setAttribute('open', '');
            }
        }

        function selectAccession(accession, meta) {
            input.value = accession;
            pendingChoice = false;
            closeModal();
            var note = document.getElementById('gene-picker-selected');
            if (note) {
                var bits = [accession];
                if (meta && meta.gene) bits.push(meta.gene);
                if (meta && meta.reviewed) bits.push('Swiss-Prot');
                note.textContent = 'Selected: ' + bits.join(' · ');
                note.hidden = false;
            }
            input.focus();
        }

        function renderResults(results, query) {
            if (!listEl) return;
            listEl.innerHTML = '';
            if (titleEl) {
                titleEl.textContent = 'UniProt matches for “' + query + '”';
            }
            if (!results || !results.length) {
                setStatus('No human UniProt entries found for “' + query + '”.', true);
                pendingChoice = false;
                return;
            }

            var reviewed = results.filter(function (r) { return r.reviewed; });
            if (results.length === 1 && results[0].reviewed) {
                selectAccession(results[0].accession, results[0]);
                var note = document.getElementById('gene-picker-selected');
                if (note) {
                    note.textContent = 'Selected ' + results[0].accession +
                        (results[0].gene ? ' (' + results[0].gene + ')' : '') +
                        ' — sole Swiss-Prot match.';
                    note.hidden = false;
                }
                return;
            }

            pendingChoice = true;
            setStatus(
                'Choose one accession to continue' +
                (reviewed.length ? ' (Swiss-Prot listed first).' : '.'),
                false
            );

            results.forEach(function (row, idx) {
                var id = 'gene-pick-' + idx;
                var label = document.createElement('label');
                label.className = 'gene-picker-item';
                label.setAttribute('for', id);

                var radio = document.createElement('input');
                radio.type = 'radio';
                radio.name = 'gene_picker_choice';
                radio.id = id;
                radio.value = row.accession;
                radio.className = 'gene-picker-radio';

                var body = document.createElement('span');
                body.className = 'gene-picker-item-body';

                var title = document.createElement('span');
                title.className = 'gene-picker-item-title';
                title.innerHTML =
                    '<strong class="gene-picker-accession">' + escapeHtml(row.accession) + '</strong>' +
                    (row.reviewed
                        ? ' <span class="gene-picker-badge gene-picker-badge-reviewed">Swiss-Prot</span>'
                        : ' <span class="gene-picker-badge">TrEMBL</span>');

                var detail = document.createElement('span');
                detail.className = 'gene-picker-item-detail';
                var parts = [];
                if (row.gene) parts.push(row.gene);
                if (row.protein_name) parts.push(row.protein_name);
                if (row.length) parts.push(row.length + ' aa');
                detail.textContent = parts.join(' · ');

                body.appendChild(title);
                body.appendChild(detail);
                label.appendChild(radio);
                label.appendChild(body);
                listEl.appendChild(label);

                radio.addEventListener('change', function () {
                    if (radio.checked) selectAccession(row.accession, row);
                });
            });
        }

        function runSearch() {
            var q = normalizeId(input.value);
            if (!q) {
                openModal();
                if (titleEl) titleEl.textContent = 'Look up gene';
                setStatus('Enter a gene symbol first (e.g. TRBC1).', true);
                return;
            }
            if (looksLikeUniprot(q) || looksLikeEnsembl(q)) {
                closeModal();
                var note = document.getElementById('gene-picker-selected');
                if (note) note.hidden = true;
                pendingChoice = false;
                return;
            }
            if (!looksLikeGeneSymbol(q)) {
                openModal();
                if (titleEl) titleEl.textContent = 'Look up gene';
                setStatus('That does not look like a gene symbol.', true);
                return;
            }

            openModal();
            if (titleEl) titleEl.textContent = 'Looking up “' + q + '”…';
            setStatus('Searching UniProt…');
            if (listEl) listEl.innerHTML = '';
            lookupBtn.disabled = true;

            var url = searchUrl + (searchUrl.indexOf('?') >= 0 ? '&' : '?') + 'q=' + encodeURIComponent(q);
            fetch(url, { headers: { Accept: 'application/json' } })
                .then(function (res) {
                    return res.json().then(function (data) {
                        return { ok: res.ok, status: res.status, data: data };
                    });
                })
                .then(function (pack) {
                    if (!pack.ok) {
                        setStatus((pack.data && pack.data.error) || 'Lookup failed.', true);
                        pendingChoice = false;
                        return;
                    }
                    renderResults(pack.data.results || [], q);
                })
                .catch(function () {
                    setStatus('Network error during gene lookup.', true);
                    pendingChoice = false;
                })
                .finally(function () {
                    lookupBtn.disabled = false;
                });
        }

        lookupBtn.addEventListener('click', function (e) {
            e.preventDefault();
            runSearch();
        });

        closeBtns.forEach(function (btn) {
            btn.addEventListener('click', function (e) {
                e.preventDefault();
                closeModal();
            });
        });

        modal.addEventListener('click', function (e) {
            if (e.target === modal) closeModal();
        });

        modal.addEventListener('cancel', function (e) {
            e.preventDefault();
            closeModal();
        });

        input.addEventListener('keydown', function (e) {
            if (e.key === 'Enter' && looksLikeGeneSymbol(input.value)) {
                e.preventDefault();
                runSearch();
            }
        });

        input.addEventListener('input', function () {
            pendingChoice = false;
            var note = document.getElementById('gene-picker-selected');
            if (note) note.hidden = true;
        });

        if (form) {
            form.addEventListener('submit', function (e) {
                var q = normalizeId(input.value);
                if (!q) return;
                if (looksLikeUniprot(q) || looksLikeEnsembl(q)) {
                    pendingChoice = false;
                    return;
                }
                if (looksLikeGeneSymbol(q)) {
                    e.preventDefault();
                    if (isModalOpen(modal) && pendingChoice) {
                        setStatus('Please select a UniProt accession from the list before running.', true);
                        return;
                    }
                    runSearch();
                }
            });
        }
    }

    window.initGenePicker = initGenePicker;
})();
