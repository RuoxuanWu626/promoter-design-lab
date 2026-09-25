/* Mode 1: where is cell-type specificity encoded? */

const Mode1 = (() => {

  let heatMode = 'delta';   // 'delta' | 'absolute'

  function elementSpec() {
    const id = $('#m1Element').value;
    const spec = { element_id: id, strand: $('#m1Strand').value };
    if (id === 'custom') spec.custom_sequence = $('#m1Custom').value.toUpperCase();
    if (id === 'cpg_segment') spec.seg_length = parseInt($('#m1SegLen').value || '60', 10);
    return spec;
  }

  function offsets() {
    const a = parseInt($('#m1From').value, 10);
    const b = parseInt($('#m1To').value, 10);
    const s = Math.max(1, parseInt($('#m1Step').value, 10));
    const out = [];
    for (let v = Math.min(a, b); v <= Math.max(a, b); v += s) out.push(v);
    return out;
  }

  function updateCount() {
    const n = offsets().length;
    const nb = parseInt($('#m1NBg').value || '1', 10);
    $('#m1Count').textContent = `${n} positions × ${nb} background${nb > 1 ? 's' : ''} = ${n * nb} predictions`;
  }

  async function run() {
    const btn = $('#m1Run');
    btn.disabled = true;
    const bar = $('#m1Progress'); bar.style.display = 'block';
    const txt = $('#m1ProgressText');
    try {
      const body = {
        element: elementSpec(),
        offsets: offsets(),
        celltype_model: $('#m1Model').value,
        target_cell_type: $('#m1Target').value || null,
        activity_method: $('#m1ActMethod').value,
        activity_window: [parseInt($('#m1ActFrom').value, 10), parseInt($('#m1ActTo').value, 10)],
        n_backgrounds: parseInt($('#m1NBg').value || '1', 10),
        background_params: {
          length: parseInt($('#m1Length').value, 10),
          gc: parseFloat($('#m1GC').value),
          cpg_oe: parseFloat($('#m1CpG').value),
          seed: parseInt($('#m1Seed').value, 10),
        },
      };
      const res = await API.job('/api/mode1/scan', body, (p) => {
        const pct = p.total ? Math.round(100 * p.done / p.total) : 0;
        bar.firstElementChild.style.width = pct + '%';
        txt.textContent = `${p.stage || ''} ${p.done}/${p.total} ${p.detail || ''}`;
      });
      State.lastMode1 = res;
      render(res);
      toast(`scan done — best position ${res.best_offset >= 0 ? '+' : ''}${res.best_offset} (tau ${fmtNum(res.best_tau, 3)})`);
    } catch (e) {
      fail(e);
      $('#m1Content').innerHTML = '';
      $('#m1Content').appendChild(h('div', { class: 'err', text: (e.message || '') + '\n' + (e.detail || '') }));
    } finally {
      btn.disabled = false;
      bar.style.display = 'none';
      txt.textContent = '';
    }
  }

  /* ---------------- rendering ---------------- */
  function render(res) {
    const box = $('#m1Content');
    clear(box);

    const cts = res.cell_types;
    const act = res.activity;                    // (cells × offsets)
    const del = res.delta_activity;
    const ti = res.target_index;
    const label = (c) => (res.cell_type_labels?.[c] || c).split(' (')[0];

    /* ---- headline ---- */
    const strongestSet = Array.from(new Set(res.strongest_cell_type));
    box.appendChild(h('div', { class: 'stat-strip' },
      stat('Best position', (res.best_offset >= 0 ? '+' : '') + res.best_offset + ' bp', 'good'),
      stat('Tau there', fmtNum(res.best_tau, 3), 'good'),
      stat('Tau range', `${fmtNum(Math.min(...res.tau), 3)} – ${fmtNum(Math.max(...res.tau), 3)}`, ''),
      stat('Element', State.elementMeta(res.element.element_id).name + ' (' + res.element.strand + ')', ''),
      stat('Strongest cell types seen', String(strongestSet.length), strongestSet.length > 1 ? 'mid' : ''),
      stat('Backgrounds', String(res.n_backgrounds || 1), '')));

    box.appendChild(h('div', { class: 'note mock' },
      h('strong', { text: 'Mock model. ' }),
      `y-values are ${res.scale}. `,
      'Tau here is computed on the activity scale shown, over this exact panel of ' +
      cts.length + ' cell types — it is not comparable to a tau from a different panel.'));

    /* ---- specificity vs position ---- */
    const targetSeries = [];
    if (ti !== null && ti !== undefined) {
      const offMean = res.offsets.map((_, j) => {
        let s = 0, n = 0;
        for (let i = 0; i < cts.length; i++) { if (i === ti) continue; s += act[i][j]; n++; }
        return n ? s / n : 0;
      });
      const offMax = res.offsets.map((_, j) => {
        let m = -Infinity;
        for (let i = 0; i < cts.length; i++) { if (i === ti) continue; m = Math.max(m, act[i][j]); }
        return m;
      });
      targetSeries.push(
        { x: res.offsets, y: act[ti], color: '#4dd4ac', label: label(cts[ti]) + ' (target)', width: 2.2, points: true },
        { x: res.offsets, y: offMean, color: '#58a6ff', label: 'off-target mean', dash: [5, 3] },
        { x: res.offsets, y: offMax, color: '#f2545b', label: 'off-target max', dash: [2, 3] });
    }

    box.appendChild(plotPanel(
      'Specificity versus position',
      'Tau rises when the panel becomes uneven. Read it together with the arms below: ' +
      'the same tau can come from raising the target or from losing everything else.',
      (canvas) => Plot.line(canvas, {
        height: 190, legendRight: true,
        x: res.offsets,
        xlabel: "element 5′ position relative to TSS (bp)",
        ylabel: 'tau',
        vlines: [{ x: 0, label: 'TSS' }],
        series: [{ x: res.offsets, y: res.tau, sd: res.tau_sd, color: '#f0a13a',
                   label: 'tau', width: 2.2, points: true }],
      })));

    if (targetSeries.length) {
      box.appendChild(plotPanel(
        'Target versus off-target activity',
        'The two arms of specificity, plotted on the same axis so you can see which one moves.',
        (canvas) => Plot.line(canvas, {
          height: 190, legendRight: true,
          xlabel: "element 5′ position relative to TSS (bp)",
          ylabel: `activity (${res.activity_method})`,
          vlines: [{ x: 0, label: 'TSS' }],
          series: targetSeries,
        })));
    }

    /* ---- mechanism strip ---- */
    if (res.decomposition && res.decomposition.some(Boolean)) {
      box.appendChild(mechanismPanel(res));
    }

    /* ---- heatmap ---- */
    const order = cts.map((c, i) => i)
      .sort((a, b) => Math.max(...act[b]) - Math.max(...act[a]));
    const matrix = order.map(i => (heatMode === 'delta' ? del[i] : act[i]));
    const rowLabels = order.map(i => label(cts[i]) + (i === ti ? ' ◀' : ''));

    const hmPanel = h('div', { class: 'panel' },
      h('h3', {}, 'Cell type × position'),
      h('div', { class: 'subtabs' },
        h('button', { class: heatMode === 'delta' ? 'active' : '', text: 'change vs background',
                      onclick: () => { heatMode = 'delta'; render(res); } }),
        h('button', { class: heatMode === 'absolute' ? 'active' : '', text: 'absolute activity',
                      onclick: () => { heatMode = 'absolute'; render(res); } })),
      h('p', { class: 'plot-sub', text: heatMode === 'delta'
        ? 'Activity minus the same background with nothing placed. Red = the element raises this cell type, blue = lowers it.'
        : 'Absolute predicted activity. Rows sorted by peak activity; the target is marked ◀.' }));
    const hmWrap = h('div', { class: 'plotwrap' });
    const hmCanvas = h('canvas');
    hmWrap.appendChild(hmCanvas);
    hmPanel.appendChild(hmWrap);
    box.appendChild(hmPanel);
    Plot.managed(hmCanvas, () => Plot.heatmap(hmCanvas, {
      matrix, rowLabels, xticksFrom: res.offsets,
      scheme: heatMode === 'delta' ? 'diverging' : 'viridis',
      annotate: false, plotHeight: Math.max(150, cts.length * 19),
      xlabel: "element 5′ position relative to TSS (bp)",
      valueLabel: heatMode === 'delta' ? 'Δ activity' : 'activity',
    }));

    /* ---- per-position table ---- */
    box.appendChild(tablePanel(res, label));

    /* ---- caveats ---- */
    box.appendChild(h('div', { class: 'note warn' },
      h('strong', { text: 'Reading this honestly' }),
      h('ul', {}, ...(res.caveats || []).map(c => h('li', { text: c })),
        h('li', { text: 'A position effect seen in one background may be a property of that background. ' +
                        'Raise "Replicates" and check the shaded spread on the tau curve.' }))));
  }

  function stat(lab, val, cls) {
    return h('div', { class: 'stat' },
      h('div', { class: 'lab', text: lab }),
      h('div', { class: 'num ' + (cls || ''), text: val }));
  }

  function plotPanel(title, sub, draw) {
    const wrap = h('div', { class: 'plotwrap' });
    const canvas = h('canvas');
    wrap.appendChild(canvas);
    const panel = h('div', { class: 'panel' },
      h('div', { class: 'plot-title', text: title }),
      h('p', { class: 'plot-sub', text: sub }),
      wrap);
    Plot.managed(canvas, () => draw(canvas));
    return panel;
  }

  const MECH_COLOR = {
    'target gain': '#4dd4ac',
    'target gain outpaces off-target': '#2e9e7f',
    'off-target loss': '#58a6ff',
    'off-target falls faster than target': '#3f7fd0',
    'off-target rises faster than target': '#f0a13a',
    'off-target gain': '#f0a13a',
    'target loss': '#f2545b',
    'target falls faster than off-target': '#c0392b',
    'no change': '#3a434e',
  };

  function mechanismPanel(res) {
    const wrap = h('div', { class: 'plotwrap' });
    const canvas = h('canvas');
    wrap.appendChild(canvas);

    const panel = h('div', { class: 'panel' },
      h('div', { class: 'plot-title', text: 'How is specificity being won?' }),
      h('p', { class: 'plot-sub', text:
        'Top: the two arms separately — target change (green) and off-target mean change (blue). ' +
        'Bottom: which arm dominates at each position. Increased target activity and decreased ' +
        'off-target activity are different designs, even when tau is the same.' }),
      wrap);

    const dt = res.decomposition.map(d => d ? d.delta_target : 0);
    const dof = res.decomposition.map(d => d ? d.delta_offtarget_mean : 0);
    const contrast = res.decomposition.map(d => d ? d.contrast_gain : 0);

    Plot.managed(canvas, () => {
      const geo = Plot.line(canvas, {
        height: 170, legendRight: true, ySymmetric: true,
        xlabel: "element 5′ position relative to TSS (bp)",
        ylabel: 'Δ vs background',
        vlines: [{ x: 0, label: 'TSS' }],
        pad: { l: 54, r: 128, t: 10, b: 52 },
        series: [
          { x: res.offsets, y: dt, color: '#4dd4ac', label: 'Δ target', width: 2 },
          { x: res.offsets, y: dof, color: '#58a6ff', label: 'Δ off-target mean', width: 2 },
          { x: res.offsets, y: contrast, color: '#f0a13a', label: 'Δ contrast', dash: [4, 3] },
        ],
      });
      /* mechanism strip under the axis */
      const ctx = canvas.getContext('2d');
      const n = res.offsets.length;
      const bw = geo.PW / n;
      const y = geo.pad.t + geo.PH + 20;
      for (let j = 0; j < n; j++) {
        const d = res.decomposition[j];
        ctx.fillStyle = d ? (MECH_COLOR[d.mechanism] || '#666') : '#333';
        ctx.fillRect(geo.pad.l + j * bw, y, Math.ceil(bw), 9);
      }
      ctx.strokeStyle = Plot.CSS('--line');
      ctx.strokeRect(geo.pad.l, y, geo.PW, 9);
      ctx.fillStyle = Plot.CSS('--fg-faint');
      ctx.font = '9px ' + Plot.CSS('--sans');
      ctx.textAlign = 'right'; ctx.textBaseline = 'middle';
      ctx.fillText('mechanism', geo.pad.l - 6, y + 4.5);
    });

    const used = Array.from(new Set(res.decomposition.filter(Boolean).map(d => d.mechanism)));
    panel.appendChild(h('div', { class: 'legend' },
      ...used.map(m => h('div', { class: 'item' },
        h('span', { class: 'sw', style: { background: MECH_COLOR[m] || '#666', height: '8px', width: '12px' } }),
        h('span', { text: m })))));
    return panel;
  }

  function tablePanel(res, label) {
    const cts = res.cell_types;
    const rows = res.offsets.map((o, j) => {
      const d = res.decomposition[j];
      return h('tr', { class: o === res.best_offset ? 'hl' : '' },
        h('td', { class: 'num', text: (o >= 0 ? '+' : '') + o }),
        h('td', { class: 'num', text: fmtNum(res.tau[j], 4) }),
        h('td', { text: label(res.strongest_cell_type[j]) }),
        h('td', { class: 'num', text: res.target_index !== null && res.target_index !== undefined
          ? fmtNum(res.activity[res.target_index][j], 3) : '–' }),
        h('td', { class: 'num', text: d ? fmtNum(d.delta_target, 3) : '–' }),
        h('td', { class: 'num', text: d ? fmtNum(d.delta_offtarget_mean, 3) : '–' }),
        h('td', { text: d ? d.mechanism : '–' }));
    });
    return h('div', { class: 'panel' },
      h('h3', {}, 'Per-position results'),
      h('div', { style: { maxHeight: '320px', overflow: 'auto' } },
        h('table', { class: 'data' },
          h('thead', {}, h('tr', {},
            h('th', { class: 'num', text: 'pos' }), h('th', { class: 'num', text: 'tau' }),
            h('th', { text: 'strongest' }), h('th', { class: 'num', text: 'target' }),
            h('th', { class: 'num', text: 'Δ target' }), h('th', { class: 'num', text: 'Δ off-mean' }),
            h('th', { text: 'mechanism' }))),
          h('tbody', {}, ...rows))));
  }

  /* ---------------- exports ---------------- */
  function exportCsv() {
    const r = State.lastMode1;
    if (!r) return toast('run a scan first', true);
    const rows = r.offsets.map((o, j) => {
      const row = { offset: o, tau: r.tau[j], tau_sd: r.tau_sd ? r.tau_sd[j] : '',
                    strongest_cell_type: r.strongest_cell_type[j] };
      const d = r.decomposition[j];
      if (d) Object.assign(row, {
        delta_target: d.delta_target, delta_offtarget_mean: d.delta_offtarget_mean,
        contrast_gain: d.contrast_gain, mechanism: d.mechanism,
        fraction_from_target_arm: d.fraction_from_target_arm,
      });
      r.cell_types.forEach((c, i) => { row['act_' + c] = r.activity[i][j]; });
      r.cell_types.forEach((c, i) => { row['delta_' + c] = r.delta_activity[i][j]; });
      return row;
    });
    download(`mode1_scan_${stamp()}.csv`, toCSV(rows), 'text/csv');
  }

  function init() {
    ['#m1From', '#m1To', '#m1Step', '#m1NBg'].forEach(s =>
      $(s).addEventListener('input', updateCount));
    $('#m1Element').addEventListener('change', () => {
      const v = $('#m1Element').value;
      $('#m1CustomWrap').style.display = v === 'custom' ? '' : 'none';
      $('#m1SegWrap').style.display = v === 'cpg_segment' ? '' : 'none';
    });
    $('#m1Run').onclick = run;
    $('#m1ExCsv').onclick = exportCsv;
    $('#m1ExJson').onclick = () => {
      if (!State.lastMode1) return toast('run a scan first', true);
      download(`mode1_scan_${stamp()}.json`, JSON.stringify(State.lastMode1, null, 1), 'application/json');
    };
    $('#m1ExBundle').onclick = async () => {
      if (!State.lastMode1) return toast('run a scan first', true);
      try {
        const r = await API.post('/api/export/bundle',
          { name: 'mode1_scan', payload: State.lastMode1, settings: { mode: 'mode1' } });
        toast('saved ' + r.filename);
      } catch (e) { fail(e); }
    };
    updateCount();
  }

  return { init, run, render, updateCount };
})();
