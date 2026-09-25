/* Mode 2: what motif interactions do Puffin and PuffinD learn? */

const Mode2 = (() => {

  const ORIENTS = [['+', '+'], ['+', '-'], ['-', '+'], ['-', '-']];
  let orientSel = new Set(['++']);

  function selectedModels() {
    return $$('#m2Models input[type=checkbox]:checked').map(c => c.value);
  }

  function baseBody() {
    return {
      motif_a: $('#m2A').value,
      motif_b: $('#m2B').value,
      strand_a: $('#m2StrandA').value,
      strand_b: $('#m2StrandB').value,
      anchor: parseInt($('#m2Anchor').value, 10),
      spacing: parseInt($('#m2Spacing').value, 10),
      spacing_mode: $('#m2SpacingMode').value,
      profile_models: selectedModels(),
      track: $('#m2Track').value,
      residual_scale: $('#m2Scale').value,
      n_backgrounds: parseInt($('#m2NBg').value || '1', 10),
      background_params: {
        length: parseInt($('#m2Length').value, 10),
        gc: parseFloat($('#m2GC').value),
        cpg_oe: parseFloat($('#m2CpG').value),
        seed: parseInt($('#m2Seed').value, 10),
      },
    };
  }

  function updateGeom() {
    const a = State.motifById[$('#m2A').value];
    const b = State.motifById[$('#m2B').value];
    if (!a || !b) return;
    const anchor = parseInt($('#m2Anchor').value, 10);
    const sp = parseInt($('#m2Spacing').value, 10);
    let pb;
    if ($('#m2SpacingMode').value === 'edge') pb = anchor + a.width + sp;
    else pb = Math.round(anchor + a.width / 2 + sp - b.width / 2);
    const overlap = Math.max(0, Math.min(anchor + a.width, pb + b.width) - Math.max(anchor, pb));
    const gap = pb >= anchor ? pb - (anchor + a.width) : anchor - (pb + b.width);
    const maxW = Math.max(...State.motifs.map(m => m.width));
    let flag = '';
    if (overlap > 0) {
      flag = `<br><span style="color:var(--danger)">⚠ overlap ${overlap} bp — the A+B construct will not contain two intact sites</span>`;
    } else if (gap < maxW) {
      flag = `<br><span style="color:var(--warn)">⚠ gap ${gap} bp &lt; widest motif (${maxW} bp) — scan windows span the junction, so part of any residual is a sequence artefact</span>`;
    }
    $('#m2Geom').innerHTML =
      `A <b>${a.name}</b> at ${anchor}…${anchor + a.width - 1} · ` +
      `B <b>${b.name}</b> at ${pb}…${pb + b.width - 1} · gap ${gap} bp` + flag;
  }

  function updateScaleNote() {
    const models = State.models ? State.models.profile : [];
    const sel = selectedModels();
    const spaces = Array.from(new Set(sel.map(n => (models.find(m => m.name === n) || {}).output_space)));
    const chosen = $('#m2Scale').value;
    const note = chosen === 'linear'
      ? 'I is computed after back-transforming to linear counts: it tests departure from <b>additivity</b> in count-like units.'
      : `I is computed on the model's own output (${spaces.join(', ') || '?'}). On a log-like output, I = 0 means the two effects combine <b>multiplicatively</b> in count units.`;
    $('#m2ScaleNote').innerHTML = note;
  }

  function progress(prefix) {
    const bar = $('#m2Progress'), txt = $('#m2ProgressText');
    bar.style.display = 'block';
    return {
      update(p) {
        const pct = p.total ? Math.round(100 * p.done / p.total) : 0;
        bar.firstElementChild.style.width = pct + '%';
        txt.textContent = `${prefix} ${p.done}/${p.total} ${p.detail || ''}`;
      },
      done() { bar.style.display = 'none'; txt.textContent = ''; },
    };
  }

  /* ===================== single pair ===================== */
  async function runPair() {
    const btn = $('#m2RunPair'); btn.disabled = true;
    const pg = progress('pair');
    try {
      const res = await API.post('/api/mode2/pair', baseBody());
      State.lastMode2.pair = res;
      renderPair(res);
      showSub('pair');
    } catch (e) { fail(e); }
    finally { btn.disabled = false; pg.done(); }
  }

  function renderPair(res) {
    const box = $('#m2Pair');
    clear(box);

    const A = State.elementMeta(res.motif_a).name;
    const B = State.elementMeta(res.motif_b).name;

    box.appendChild(h('div', { class: 'stat-strip' },
      statEl('Pair', `${A} ${res.strand_a}  ×  ${B} ${res.strand_b}`, ''),
      statEl('A at', String(res.position_a), ''),
      statEl('B at', String(res.position_b), ''),
      statEl('Spacing', `${res.spacing} bp (${res.spacing_mode})`, ''),
      statEl('Backgrounds', String(res.n_backgrounds || 1), ''),
      statEl('Gap', res.gap_bp + ' bp', res.overlap_bp ? 'bad' : res.junction_risk ? 'mid' : '')));

    if (res.overlap_warning) {
      box.appendChild(h('div', { class: 'note warn' }, h('strong', { text: '⚠ ' }), res.overlap_warning));
    }
    if (res.junction_warning) {
      const additive = Object.entries(res.models).find(
        ([, m]) => m.meta && m.meta.additive_by_construction);
      box.appendChild(h('div', { class: 'note warn' },
        h('strong', { text: '⚠ junction artefact · ' }), res.junction_warning,
        additive ? h('div', { style: { marginTop: '5px' } },
          `Here ${additive[1].label} reports Σ|I| = ${fmtNum(additive[1].summary.I_sum_abs, 4)}. `,
          'That model cannot be non-additive, so this is a direct measure of the artefact: ',
          'any interaction you want to claim has to be bigger than it.') : null));
    }

    /* summary across models */
    const modelNames = Object.keys(res.models);
    const rows = modelNames.map(n => {
      const m = res.models[n], s = m.summary, sd = m.summary_sd || {};
      return h('tr', {},
        h('td', {}, m.label,
          h('span', { class: 'badge ' + (m.is_mock ? 'mock' : 'real'),
                      style: { marginLeft: '6px', fontSize: '9px', padding: '1px 5px' },
                      text: m.is_mock ? 'mock' : 'real' })),
        h('td', { class: 'num', text: fmtNum(s.effect_A_max, 3) }),
        h('td', { class: 'num', text: fmtNum(s.effect_B_max, 3) }),
        h('td', { class: 'num', text: fmtNum(s.effect_sum_max, 3) }),
        h('td', { class: 'num', text: fmtNum(s.effect_AB_max, 3) }),
        h('td', { class: 'num', text: fmtNum(s.I_sum_abs, 4) + (sd.I_sum_abs ? ' ±' + fmtNum(sd.I_sum_abs, 3) : '') }),
        h('td', { class: 'num', text: fmtNum(s.I_at_max_abs, 4) }),
        h('td', { class: 'num', text: s.I_position_at_max_abs }),
        h('td', { class: 'num', text: fmtNum(s.relative_I, 4) }));
    });
    box.appendChild(h('div', { class: 'panel' },
      h('h3', {}, 'Interaction summary'),
      h('table', { class: 'data' },
        h('thead', {}, h('tr', {},
          h('th', { text: 'model' }),
          h('th', { class: 'num', text: 'max A' }), h('th', { class: 'num', text: 'max B' }),
          h('th', { class: 'num', text: 'max A+B expected' }), h('th', { class: 'num', text: 'max A+B observed' }),
          h('th', { class: 'num', text: 'Σ|I|' }), h('th', { class: 'num', text: 'I at max|I|' }),
          h('th', { class: 'num', text: 'at pos' }), h('th', { class: 'num', text: 'rel |I|' }))),
        h('tbody', {}, ...rows)),
      h('div', { class: 'hint', text:
        '"expected" is f(A) + f(B) − f(background), i.e. what a purely additive model would give. ' +
        'Σ|I| accumulates across the whole window, so it grows with window width — compare it only between runs with the same window.' })));

    /* per-model plots */
    for (const n of modelNames) {
      const m = res.models[n];
      const pos = m.positions;

      const panel = h('div', { class: 'panel' },
        h('div', { class: 'plot-title' }, m.label,
          h('span', { class: 'badge ' + (m.is_mock ? 'mock' : 'real'),
                      style: { fontSize: '9px', padding: '1px 6px' }, text: m.is_mock ? 'mock' : 'real' })),
        h('p', { class: 'plot-sub', text: `output: ${m.scale} · residual computed on the ${m.scale_used} scale` }));

      const w1 = h('div', { class: 'plotwrap' }); const c1 = h('canvas'); w1.appendChild(c1);
      const w2 = h('div', { class: 'plotwrap' }); const c2 = h('canvas'); w2.appendChild(c2);
      panel.appendChild(w1); panel.appendChild(w2);
      box.appendChild(panel);

      Plot.managed(c1, () => Plot.line(c1, {
        height: 180, legendRight: true,
        xlabel: 'position relative to TSS (bp)', ylabel: 'prediction',
        vlines: [{ x: 0, label: 'TSS' },
                 { x: res.position_a, color: '#4dd4ac', dash: [2, 2], label: 'A' },
                 { x: res.position_b, color: '#b06bd6', dash: [2, 2], label: 'B' }],
        series: [
          { x: pos, y: m.profiles.background, color: '#6b7784', label: 'background' },
          { x: pos, y: m.profiles.A, color: '#4dd4ac', label: 'A only' },
          { x: pos, y: m.profiles.B, color: '#b06bd6', label: 'B only' },
          { x: pos, y: m.profiles.AB, color: '#f0a13a', label: 'A + B', width: 2.2 },
          { x: pos, y: m.expected_additive, color: '#58a6ff', label: 'additive expectation', dash: [5, 3] },
        ],
      }));

      Plot.managed(c2, () => Plot.line(c2, {
        height: 150, legendRight: true, ySymmetric: true,
        xlabel: 'position relative to TSS (bp)',
        ylabel: 'I',
        vlines: [{ x: 0, label: 'TSS' },
                 { x: res.position_a, color: '#4dd4ac', dash: [2, 2] },
                 { x: res.position_b, color: '#b06bd6', dash: [2, 2] }],
        series: [{ x: pos, y: m.residual, sd: m.residual_sd, color: '#f2545b',
                   label: 'I = AB − A − B + bg', width: 2, fill: true }],
      }));

      panel.appendChild(h('div', { class: 'hint', text: m.residual_meaning }));
      if (m.meta && m.meta.additive_by_construction) {
        panel.appendChild(h('div', { class: 'note' },
          'This adapter is additive by construction, so I should be flat at zero. ' +
          'Any visible deviation means the two elements are close enough to share a ' +
          'PWM scan window — a useful check that the pipeline is doing what you think.'));
      }
    }

    box.appendChild(interpretationNote());
  }

  function statEl(lab, val, cls) {
    return h('div', { class: 'stat' },
      h('div', { class: 'lab', text: lab }),
      h('div', { class: 'num ' + (cls || ''), text: val }));
  }

  function interpretationNote() {
    return h('div', { class: 'note warn' },
      h('strong', { text: 'What a non-zero I does and does not show' }),
      h('ul', {},
        h('li', { text: 'It shows this model does not combine these two elements additively on the stated scale. That is a statement about the model.' }),
        h('li', { text: 'It is not evidence of biological cooperativity, and not evidence that the two factors touch.' }),
        h('li', { text: 'Profiles that merely look different between constructs establish even less — that is why all four constructs are differenced.' }),
        h('li', { text: 'The scale matters: zero residual on a log output means multiplicative combination in counts, not no interaction. Switch "Residual on" to compare.' }),
        h('li', { text: 'Saturation alone produces a non-zero I. Compare against "PuffinD (mock, pair terms off)" to separate a saturating response from a genuine pair term.' }),
        h('li', {}, h('strong', { text: 'Keep an additive model in the comparison. ' }),
          'It cannot be non-additive, so whatever I it reports is pure sequence-composition ' +
          'artefact from the junction between the two sites. That number is the floor any real ' +
          'interaction claim has to clear.'),
        h('li', { text: 'Repeat across backgrounds. The shaded band on the residual is the spread across replicates.' })));
  }

  /* ===================== grid ===================== */
  async function runGrid() {
    const btn = $('#m2RunGrid'); btn.disabled = true;
    const pg = progress('grid');
    try {
      const body = Object.assign(baseBody(), { statistic: $('#m2Stat').value });
      const res = await API.job('/api/mode2/grid', body, pg.update);
      State.lastMode2.grid = res;
      renderGrid(res);
      showSub('grid');
    } catch (e) { fail(e); }
    finally { btn.disabled = false; pg.done(); }
  }

  function renderGrid(res) {
    const box = $('#m2Grid');
    clear(box);
    const labels = res.motif_ids.map(i => State.elementMeta(i).name);

    box.appendChild(h('div', { class: 'note' },
      `Every pair from the 10-motif library, including same-motif pairs on the diagonal. `,
      `Statistic: ${$('#m2Stat').selectedOptions[0].textContent}. `,
      `Anchor ${res.anchor}, spacing ${res.spacing} bp (${res.spacing_mode}), `,
      `orientation ${res.strand_a}${res.strand_b}, ${res.n_backgrounds} background(s). `,
      h('strong', { text: 'Click a cell to open that pair in the detail view.' })));

    const models = Object.keys(res.grids);
    /* one shared colour scale so the panels are comparable */
    let vmax = 0;
    for (const m of models) for (const row of res.grids[m]) for (const v of row)
      if (isFinite(v)) vmax = Math.max(vmax, Math.abs(v));

    for (const m of models) {
      const info = (State.models.profile || []).find(x => x.name === m) || {};
      const wrap = h('div', { class: 'plotwrap' });
      const canvas = h('canvas');
      wrap.appendChild(canvas);
      const panel = h('div', { class: 'panel' },
        h('div', { class: 'plot-title' }, info.label || m,
          h('span', { class: 'badge ' + (info.is_mock ? 'mock' : 'real'),
                      style: { fontSize: '9px', padding: '1px 6px' }, text: info.is_mock ? 'mock' : 'real' })),
        h('p', { class: 'plot-sub', text: info.description || '' }),
        wrap);
      box.appendChild(panel);

      Plot.managed(canvas, () => Plot.heatmap(canvas, {
        matrix: res.grids[m], rowLabels: labels, colLabels: labels,
        scheme: 'diverging', vmin: -vmax, vmax: vmax,
        plotHeight: Math.max(220, labels.length * 26),
        valueLabel: res.statistic,
        onCell: (i, j) => {
          $('#m2A').value = res.motif_ids[i];
          $('#m2B').value = res.motif_ids[j];
          updateGeom();
          runPair();
        },
      }));
    }

    /* ranked table */
    const flat = [];
    res.cells.forEach(c => {
      models.forEach(m => flat.push({
        pair: `${State.elementMeta(c.a).name} × ${State.elementMeta(c.b).name}`,
        a: c.a, b: c.b, model: m, junction: c.junction_risk,
        value: c.values[m].value, sd: c.values[m].sd, overlap: c.overlap_bp,
      }));
    });
    const top = flat.filter(f => isFinite(f.value))
      .sort((a, b) => Math.abs(b.value) - Math.abs(a.value)).slice(0, 25);
    box.appendChild(h('div', { class: 'panel' },
      h('h3', {}, 'Strongest 25 cells'),
      h('table', { class: 'data' },
        h('thead', {}, h('tr', {},
          h('th', { text: 'pair' }), h('th', { text: 'model' }),
          h('th', { class: 'num', text: res.statistic }), h('th', { class: 'num', text: 'sd' }),
          h('th', { class: 'num', text: 'overlap' }), h('th', {}))),
        h('tbody', {}, ...top.map(f => h('tr', { class: f.overlap ? 'target' : f.junction ? 'hl' : '' },
          h('td', { text: f.pair }), h('td', { text: f.model }),
          h('td', { class: 'num', text: fmtNum(f.value, 4) }),
          h('td', { class: 'num', text: fmtNum(f.sd, 4) }),
          h('td', { class: 'num', text: f.overlap + ' bp' }),
          h('td', {}, h('button', {
            class: 'btn sm', text: 'open',
            onclick: () => { $('#m2A').value = f.a; $('#m2B').value = f.b; updateGeom(); runPair(); },
          })))))),
      h('div', { class: 'hint', text:
        'Blue rows: sites overlap at this spacing — not an interaction measurement. ' +
        'Green rows: gap smaller than the widest motif, so part of the residual is a junction artefact.' })));

    box.appendChild(h('div', { class: 'note warn' },
      h('strong', { text: 'One grid is one geometry. ' }),
      'These values are for a single anchor, spacing and orientation. A pair that looks inert here ' +
      'may interact strongly at another spacing — use the spacing curve before concluding anything.'));
  }

  /* ===================== spacing ===================== */
  async function runSpacing() {
    const btn = $('#m2RunSpacing'); btn.disabled = true;
    const pg = progress('spacing');
    try {
      const from = parseInt($('#m2SpFrom').value, 10);
      const to = parseInt($('#m2SpTo').value, 10);
      const step = Math.max(1, parseInt($('#m2SpStep').value, 10));
      const spacings = [];
      for (let v = Math.min(from, to); v <= Math.max(from, to); v += step) spacings.push(v);

      const orientations = ORIENTS.filter(o => orientSel.has(o[0] + o[1]));
      const body = Object.assign(baseBody(), {
        spacings, statistic: $('#m2Stat').value,
        orientations: orientations.length ? orientations : [['+', '+']],
      });
      const res = await API.job('/api/mode2/spacing', body, pg.update);
      State.lastMode2.spacing = res;
      renderSpacing(res);
      showSub('spacing');
    } catch (e) { fail(e); }
    finally { btn.disabled = false; pg.done(); }
  }

  function renderSpacing(res) {
    const box = $('#m2SpacingView');
    clear(box);
    const A = State.elementMeta(res.motif_a).name;
    const B = State.elementMeta(res.motif_b).name;

    box.appendChild(h('div', { class: 'note' },
      `${A} × ${B}, ${res.statistic} against ${res.spacing_mode === 'center' ? 'centre-to-centre' : 'edge-to-edge'} spacing, `,
      `anchor ${res.anchor}, ${res.n_backgrounds} background(s). `,
      h('strong', { text: 'A periodic ripple near 10–11 bp would indicate helical-face dependence — check it is not aliasing from your step size.' })));

    /* Which spacings put the two sites on top of each other? There the A+B
     * construct does not contain two intact sites, so the residual is not an
     * interaction measurement. Those points are excluded from the curve --
     * otherwise their (large, meaningless) values set the y-axis and flatten
     * everything worth looking at. They stay in the table below. */
    const overlapAt = res.spacings.map((_, k) =>
      Object.values(res.series).some(b => (b.overlap_bp || [])[k] > 0));
    const nExcluded = overlapAt.filter(Boolean).length;
    const junctionAt = res.spacings.map((_, k) =>
      !overlapAt[k] && Object.values(res.series).some(b => (b.junction_risk || [])[k]));

    const overlapBands = [];
    for (let k = 0; k < overlapAt.length; k++) {
      if (!overlapAt[k]) continue;
      let j = k;
      while (j + 1 < overlapAt.length && overlapAt[j + 1]) j++;
      const step = res.spacings.length > 1 ? (res.spacings[1] - res.spacings[0]) : 1;
      overlapBands.push({
        from: res.spacings[k] - step / 2, to: res.spacings[j] + step / 2,
        color: 'rgba(242,84,91,.13)',
      });
      k = j;
    }
    const step0 = res.spacings.length > 1 ? (res.spacings[1] - res.spacings[0]) : 1;
    for (let k = 0; k < junctionAt.length; k++) {
      if (!junctionAt[k]) continue;
      let j = k;
      while (j + 1 < junctionAt.length && junctionAt[j + 1]) j++;
      overlapBands.push({
        from: res.spacings[k] - step0 / 2, to: res.spacings[j] + step0 / 2,
        color: 'rgba(240,161,58,.12)',
      });
      k = j;
    }

    const series = [];
    let ci = 0;
    for (const [okey, byModel] of Object.entries(res.series)) {
      for (const m of res.models) {
        if (!byModel[m]) continue;
        series.push({
          x: res.spacings,
          y: byModel[m].value.map((v, k) => (overlapAt[k] ? null : v)),
          sd: byModel[m].sd,
          color: Plot.PALETTE[ci % Plot.PALETTE.length],
          label: `${m.replace('mock_', '')} ${okey}`,
          width: m.includes('puffind') ? 2 : 1.5,
          dash: okey === '++' ? [] : okey === '+-' ? [5, 3] : okey === '-+' ? [2, 3] : [7, 2, 2, 2],
          points: res.spacings.length < 45,
        });
        ci++;
      }
    }

    const wrap = h('div', { class: 'plotwrap' });
    const canvas = h('canvas');
    wrap.appendChild(canvas);
    box.appendChild(h('div', { class: 'panel' },
      h('div', { class: 'plot-title', text: `Interaction versus spacing — ${A} × ${B}` }),
      h('p', { class: 'plot-sub', text:
        'Shaded ribbon = spread across background replicates. Red band = spacings where the two ' +
        'sites overlap; those points are left out of the curve because the A+B construct does not ' +
        'contain two intact sites there' +
        (nExcluded ? ` (${nExcluded} spacing${nExcluded === 1 ? '' : 's'} excluded — still listed in the table below).` : '.') +
        ' Amber band = sites too close for scan windows to stay on one site: an additive model' +
        ' showing a non-zero value there is measuring the sequence artefact, not an interaction.' }),
      wrap));

    Plot.managed(canvas, () => Plot.line(canvas, {
      height: 260, legendRight: true, yZero: true,
      xlabel: `${res.spacing_mode === 'center' ? 'centre-to-centre' : 'edge-to-edge'} spacing (bp)`,
      ylabel: res.statistic,
      bands: overlapBands,
      series,
    }));

    /* table */
    const rows = res.spacings.map((s, k) => {
      const tds = [h('td', { class: 'num', text: s })];
      for (const [okey, byModel] of Object.entries(res.series)) {
        for (const m of res.models) {
          if (!byModel[m]) continue;
          tds.push(h('td', { class: 'num', text: fmtNum(byModel[m].value[k], 4) }));
        }
      }
      return h('tr', { class: overlapAt[k] ? 'target' : junctionAt[k] ? 'hl' : '' }, ...tds);
    });
    const heads = [h('th', { class: 'num', text: 'spacing' })];
    for (const [okey, byModel] of Object.entries(res.series))
      for (const m of res.models) if (byModel[m]) heads.push(h('th', { class: 'num', text: `${m.replace('mock_', '')} ${okey}` }));

    box.appendChild(h('div', { class: 'panel' },
      h('h3', {}, 'Values'),
      h('div', { style: { maxHeight: '300px', overflow: 'auto' } },
        h('table', { class: 'data' },
          h('thead', {}, h('tr', {}, ...heads)),
          h('tbody', {}, ...rows)))));

    box.appendChild(interpretationNote());
  }

  /* ===================== sub-tabs & exports ===================== */
  function showSub(name) {
    $$('#m2Subtabs button').forEach(b => b.classList.toggle('active', b.dataset.sub === name));
    $('#m2Pair').style.display = name === 'pair' ? '' : 'none';
    $('#m2Grid').style.display = name === 'grid' ? '' : 'none';
    $('#m2SpacingView').style.display = name === 'spacing' ? '' : 'none';
    Plot.redrawVisible();
  }

  function exportCsv() {
    const { pair, grid, spacing } = State.lastMode2;
    if (grid) {
      const rows = [];
      grid.cells.forEach(c => Object.entries(c.values).forEach(([m, v]) => rows.push({
        motif_a: c.a, motif_b: c.b, model: m, statistic: grid.statistic,
        value: v.value, sd: v.sd, overlap_bp: c.overlap_bp,
        anchor: grid.anchor, spacing: grid.spacing, spacing_mode: grid.spacing_mode,
        strand_a: grid.strand_a, strand_b: grid.strand_b,
        n_backgrounds: grid.n_backgrounds, residual_scale: grid.residual_scale,
      })));
      download(`mode2_grid_${stamp()}.csv`, toCSV(rows), 'text/csv');
      return;
    }
    if (spacing) {
      const rows = [];
      Object.entries(spacing.series).forEach(([okey, byModel]) => {
        spacing.models.forEach(m => {
          if (!byModel[m]) return;
          spacing.spacings.forEach((s, k) => rows.push({
            motif_a: spacing.motif_a, motif_b: spacing.motif_b, orientation: okey,
            model: m, spacing: s, statistic: spacing.statistic,
            value: byModel[m].value[k], sd: byModel[m].sd[k],
            overlap_bp: (byModel.overlap_bp || [])[k],
          }));
        });
      });
      download(`mode2_spacing_${stamp()}.csv`, toCSV(rows), 'text/csv');
      return;
    }
    if (pair) {
      const names = Object.keys(pair.models);
      const pos = pair.models[names[0]].positions;
      const rows = pos.map((p, i) => {
        const r = { position: p };
        for (const n of names) {
          const m = pair.models[n];
          r[`${n}_bg`] = m.profiles.background[i];
          r[`${n}_A`] = m.profiles.A[i];
          r[`${n}_B`] = m.profiles.B[i];
          r[`${n}_AB`] = m.profiles.AB[i];
          r[`${n}_expected`] = m.expected_additive[i];
          r[`${n}_I`] = m.residual[i];
        }
        return r;
      });
      download(`mode2_pair_${stamp()}.csv`, toCSV(rows), 'text/csv');
      return;
    }
    toast('run something first', true);
  }

  async function exportFasta() {
    const p = State.lastMode2.pair;
    if (!p) return toast('run a single pair first', true);
    const bp = baseBody().background_params;
    const sets = [
      ['background', []],
      ['A_only', [{ element_id: p.motif_a, position: p.position_a, strand: p.strand_a }]],
      ['B_only', [{ element_id: p.motif_b, position: p.position_b, strand: p.strand_b }]],
      ['A_plus_B', [{ element_id: p.motif_a, position: p.position_a, strand: p.strand_a },
                    { element_id: p.motif_b, position: p.position_b, strand: p.strand_b }]],
    ];
    try {
      let out = '';
      for (const [name, placements] of sets) {
        const r = await API.post('/api/export/fasta',
          Object.assign({}, bp, { placements, name: `${p.motif_a}_${p.motif_b}_${name}` }));
        out += r.fasta;
      }
      download(`mode2_constructs_${stamp()}.fa`, out, 'text/plain');
    } catch (e) { fail(e); }
  }

  function init() {
    /* motif selects */
    for (const sel of ['#m2A', '#m2B']) {
      const el = $(sel);
      clear(el);
      State.motifs.forEach(m => el.appendChild(h('option', { value: m.id, text: m.name })));
    }
    $('#m2A').value = 'tata';
    $('#m2B').value = 'inr';

    /* model checkboxes */
    const box = $('#m2Models');
    clear(box);
    (State.models.profile || []).forEach(m => {
      const id = 'm2m_' + m.name;
      const on = ['mock_puffin', 'mock_puffind'].includes(m.name);
      box.appendChild(h('div', { class: 'row', style: { marginBottom: '4px' } },
        h('input', { type: 'checkbox', id, value: m.name, checked: on && m.available,
                     disabled: !m.available, style: { flex: '0 0 auto', width: 'auto' },
                     onchange: updateScaleNote }),
        h('label', { for: id, class: 'wide', style: { flex: '1', cursor: 'pointer' },
                     title: m.available ? m.description : m.unavailable_reason,
                     text: m.label + (m.available ? '' : ' — unavailable') })));
    });

    /* orientation toggles */
    const ob = $('#m2Orients');
    clear(ob);
    ORIENTS.forEach(([a, b]) => {
      const key = a + b;
      const btn = h('button', {
        class: 'btn sm' + (orientSel.has(key) ? ' primary' : ''),
        text: `${a}${b}`,
        onclick: () => {
          orientSel.has(key) ? orientSel.delete(key) : orientSel.add(key);
          if (!orientSel.size) orientSel.add('++');
          init_orients();
        },
      });
      ob.appendChild(btn);
    });
    function init_orients() {
      $$('#m2Orients button').forEach(b =>
        b.classList.toggle('primary', orientSel.has(b.textContent)));
    }
    init_orients();

    ['#m2A', '#m2B', '#m2Anchor', '#m2Spacing', '#m2SpacingMode'].forEach(s =>
      $(s).addEventListener('input', updateGeom));
    $('#m2SpacingMode').addEventListener('change', updateGeom);
    $('#m2Scale').addEventListener('change', updateScaleNote);

    $('#m2RunPair').onclick = runPair;
    $('#m2RunGrid').onclick = runGrid;
    $('#m2RunSpacing').onclick = runSpacing;
    $$('#m2Subtabs button').forEach(b => b.onclick = () => showSub(b.dataset.sub));

    $('#m2ExCsv').onclick = exportCsv;
    $('#m2ExFasta').onclick = exportFasta;
    $('#m2ExJson').onclick = () => {
      const payload = State.lastMode2;
      if (!payload.pair && !payload.grid && !payload.spacing) return toast('run something first', true);
      download(`mode2_${stamp()}.json`, JSON.stringify(payload, null, 1), 'application/json');
    };
    $('#m2ExBundle').onclick = async () => {
      try {
        const r = await API.post('/api/export/bundle',
          { name: 'mode2', payload: State.lastMode2, settings: baseBody() });
        toast('saved ' + r.filename);
      } catch (e) { fail(e); }
    };

    updateGeom();
    updateScaleNote();
  }

  return { init, runPair, runGrid, runSpacing, showSub };
})();
