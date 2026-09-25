/* The Design tab: motif palette, construct track with dragging, and the
 * prediction / metric panels. */

const Designer = (() => {

  let evalTimer = null;
  let busy = false;
  let pending = false;

  /* ================= palette ================= */
  function renderPalette() {
    const box = $('#palette');
    clear(box);
    const items = State.motifs.concat(State.extraElements);
    for (const m of items) {
      const chip = h('div', {
        class: 'chip', draggable: 'true', title: m.notes || '',
        style: { borderLeftColor: m.color || '#888' },
        onclick: () => { State.addPlacement(m.id); changed(); },
        ondragstart: (ev) => {
          ev.dataTransfer.setData('text/plain', m.id);
          ev.dataTransfer.effectAllowed = 'copy';
        },
      },
        h('span', { text: m.name }),
        m.consensus ? h('span', { class: 'cons', text: m.consensus }) : null);
      box.appendChild(chip);
    }
  }

  /* ================= element list ================= */
  function renderList() {
    const box = $('#elist');
    clear(box);
    $('#elCount').textContent = State.design.placements.length
      ? `(${State.design.placements.length})` : '';

    if (!State.design.placements.length) {
      box.appendChild(h('div', { class: 'hint', text: 'Nothing placed yet — click a motif above.' }));
      return;
    }

    for (const p of State.design.placements) {
      const meta = State.elementMeta(p.element_id);
      const w = State.elementWidth(p);
      const row = h('div', {
        class: 'el' + (State.design.selected === p.uid ? ' sel' : ''),
        onclick: () => { State.design.selected = p.uid; renderList(); drawTrack(); },
      },
        h('div', { class: 'swatch', style: { background: meta.color || '#888' } }),
        h('div', {},
          h('div', { class: 'name', text: meta.name }),
          h('div', { class: 'meta', text: `${p.position} … ${p.position + w - 1}  (${w} bp, ${p.strand})` }),
          h('div', { class: 'row', style: { marginTop: '3px', marginBottom: '0' } },
            h('input', {
              type: 'number', value: p.position, title: '5′ position relative to TSS',
              oninput: (e) => { p.position = parseInt(e.target.value || '0', 10); changed(); },
              onclick: (e) => e.stopPropagation(),
            }),
            p.element_id === 'cpg_segment' ? h('input', {
              type: 'number', value: p.seg_length, title: 'segment length (bp)',
              style: { width: '54px' },
              oninput: (e) => { p.seg_length = Math.max(5, parseInt(e.target.value || '60', 10)); changed(); },
              onclick: (e) => e.stopPropagation(),
            }) : null,
          ),
          p.element_id === 'custom' ? h('input', {
            type: 'text', value: p.custom_sequence, placeholder: 'ACGT…',
            style: { marginTop: '3px' },
            oninput: (e) => { p.custom_sequence = e.target.value.toUpperCase(); changed(); },
            onclick: (e) => e.stopPropagation(),
          }) : null,
        ),
        h('div', { class: 'acts' },
          h('button', {
            class: 'btn icon sm', title: 'reverse strand', text: p.strand === '+' ? '→' : '←',
            onclick: (e) => { e.stopPropagation(); p.strand = p.strand === '+' ? '-' : '+'; changed(); },
          }),
          h('button', {
            class: 'btn icon sm', title: 'duplicate', text: '⧉',
            onclick: (e) => { e.stopPropagation(); State.duplicatePlacement(p.uid); changed(); },
          }),
          h('button', {
            class: 'btn icon sm', title: 'remove', text: '×',
            onclick: (e) => { e.stopPropagation(); State.removePlacement(p.uid); changed(); },
          }),
        ));
      box.appendChild(row);
    }
  }

  /* ================= weight sliders ================= */
  function renderWeights() {
    const box = $('#weightSliders');
    clear(box);
    const w = State.design.weights;
    const shown = State.metricCatalogue.filter(m => m.name in w);
    const rest = State.metricCatalogue.filter(m => !(m.name in w));

    for (const m of shown) {
      box.appendChild(h('div', { class: 'row' },
        h('label', { class: 'wide', text: m.label, title: `${m.higher_is_better ? 'higher' : 'lower'} is better · normalised over [${m.lo}, ${m.hi}]` }),
        h('input', {
          type: 'range', min: '0', max: '3', step: '0.05', value: w[m.name],
          oninput: (e) => {
            w[m.name] = parseFloat(e.target.value);
            e.target.nextSibling.textContent = w[m.name].toFixed(2);
            $('#weightPreset').value = 'custom';
            scheduleEvaluate();
          },
        }),
        h('span', { class: 'val', text: Number(w[m.name]).toFixed(2) }),
        h('button', {
          class: 'btn icon sm', text: '×', title: 'drop from score',
          onclick: () => { delete w[m.name]; $('#weightPreset').value = 'custom'; renderWeights(); scheduleEvaluate(); },
        })));
    }

    if (rest.length) {
      const sel = h('select', {
        onchange: (e) => {
          if (!e.target.value) return;
          w[e.target.value] = 1.0;
          $('#weightPreset').value = 'custom';
          renderWeights(); scheduleEvaluate();
        },
      }, h('option', { value: '', text: '+ add a metric to the score' }),
         ...rest.map(m => h('option', { value: m.name, text: m.label })));
      box.appendChild(h('div', { class: 'row' }, sel));
    }
  }

  /* ================= construct track ================= */
  const track = { drag: null, lanes: [] };

  function viewRange() {
    const v = State.design.view;
    return [v.center - v.halfWidth, v.center + v.halfWidth];
  }

  function assignLanes() {
    const sorted = State.design.placements.slice()
      .sort((a, b) => a.position - b.position);
    const laneEnd = [];
    const lanes = {};
    for (const p of sorted) {
      const w = State.elementWidth(p);
      let L = 0;
      while (L < laneEnd.length && laneEnd[L] > p.position - 2) L++;
      laneEnd[L] = p.position + w;
      lanes[p.uid] = L;
    }
    track.lanes = lanes;
    return Math.max(1, laneEnd.length);
  }

  function drawTrack() {
    const canvas = $('#trackCanvas');
    const nLanes = assignLanes();
    const rulerH = 26, seqH = 16, laneH = 19;
    const height = rulerH + seqH + nLanes * laneH + 16;
    const { ctx, w, h } = Plot.setup(canvas, height);

    const [x0, x1] = viewRange();
    const padL = 6, padR = 6;
    const PW = w - padL - padR;
    const X = v => padL + (v - x0) / (x1 - x0) * PW;
    const bpPerPx = (x1 - x0) / PW;

    ctx.fillStyle = Plot.CSS('--bg');
    ctx.fillRect(0, 0, w, h);

    /* ruler */
    const tk = Plot.ticks(x0, x1, Math.max(4, Math.floor(w / 95)));
    ctx.font = '9.5px ' + Plot.CSS('--mono');
    ctx.textAlign = 'center'; ctx.textBaseline = 'top';
    for (const t of tk) {
      if (t < x0 || t > x1) continue;
      const x = X(t);
      ctx.strokeStyle = Plot.CSS('--line-soft');
      ctx.beginPath(); ctx.moveTo(x, 14); ctx.lineTo(x, h - 4); ctx.stroke();
      ctx.fillStyle = Plot.CSS('--fg-faint');
      ctx.fillText(t > 0 ? '+' + t : String(t), x, 2);
    }

    /* sequence letters / composition strip */
    const seqY = rulerH;
    const seq = State.lastEvaluation && State.lastEvaluation.sequence;
    const tssIdx = State.lastEvaluation && State.lastEvaluation.background.tss_index;
    if (seq && bpPerPx < 0.14) {
      ctx.font = '10px ' + Plot.CSS('--mono');
      ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
      const COL = { A: '#5fb878', C: '#4c8dd8', G: '#e0a63c', T: '#d9534f', N: '#666' };
      for (let p = Math.floor(x0); p <= Math.ceil(x1); p++) {
        const i = tssIdx + p;
        if (i < 0 || i >= seq.length) continue;
        const ch = seq[i];
        ctx.fillStyle = COL[ch] || '#888';
        ctx.fillText(ch, X(p + 0.5), seqY + seqH / 2);
      }
    } else if (seq) {
      /* CpG ticks + GC shading when too zoomed out for letters */
      const i0 = Math.max(0, tssIdx + Math.floor(x0));
      const i1 = Math.min(seq.length - 1, tssIdx + Math.ceil(x1));
      ctx.fillStyle = 'rgba(88,166,255,.42)';
      for (let i = i0; i < i1; i++) {
        if (seq[i] === 'C' && seq[i + 1] === 'G') {
          ctx.fillRect(X(i - tssIdx), seqY + 3, Math.max(1, 1 / bpPerPx), seqH - 6);
        }
      }
      ctx.fillStyle = Plot.CSS('--fg-faint');
      ctx.font = '9px ' + Plot.CSS('--sans');
      ctx.textAlign = 'left'; ctx.textBaseline = 'middle';
      ctx.fillText('CpG', padL + 2, seqY + seqH / 2);
    }

    /* TSS marker */
    if (0 >= x0 && 0 <= x1) {
      const xt = X(0);
      ctx.strokeStyle = Plot.CSS('--warn'); ctx.lineWidth = 1.4;
      ctx.beginPath(); ctx.moveTo(xt, 12); ctx.lineTo(xt, h - 4); ctx.stroke();
      ctx.beginPath();
      ctx.moveTo(xt, 15); ctx.lineTo(xt + 16, 15);
      ctx.moveTo(xt + 16, 15); ctx.lineTo(xt + 11, 12);
      ctx.moveTo(xt + 16, 15); ctx.lineTo(xt + 11, 18);
      ctx.stroke();
      ctx.fillStyle = Plot.CSS('--warn');
      ctx.font = 'bold 9px ' + Plot.CSS('--sans');
      ctx.textAlign = 'left'; ctx.textBaseline = 'bottom';
      ctx.fillText('TSS', xt + 19, 19);
      ctx.lineWidth = 1;
    }

    /* element boxes */
    const boxTop = rulerH + seqH + 4;
    for (const p of State.design.placements) {
      const meta = State.elementMeta(p.element_id);
      const wbp = State.elementWidth(p);
      const bx = X(p.position), bw = Math.max(3, wbp / bpPerPx);
      const by = boxTop + (track.lanes[p.uid] || 0) * laneH;
      const bh = laneH - 4;
      if (bx + bw < -20 || bx > w + 20) continue;

      const sel = State.design.selected === p.uid;
      ctx.fillStyle = (meta.color || '#888') + (sel ? 'ee' : 'aa');
      ctx.fillRect(bx, by, bw, bh);
      if (sel) {
        ctx.strokeStyle = Plot.CSS('--fg'); ctx.lineWidth = 1.4;
        ctx.strokeRect(bx - 0.5, by - 0.5, bw + 1, bh + 1);
        ctx.lineWidth = 1;
      }

      /* strand arrow */
      ctx.fillStyle = 'rgba(0,0,0,.62)';
      ctx.font = 'bold 9px ' + Plot.CSS('--mono');
      ctx.textBaseline = 'middle';
      if (bw > 13) {
        ctx.textAlign = p.strand === '+' ? 'right' : 'left';
        ctx.fillText(p.strand === '+' ? '▸' : '◂',
                     p.strand === '+' ? bx + bw - 2 : bx + 2, by + bh / 2);
      }
      if (bw > 42) {
        ctx.fillStyle = 'rgba(0,0,0,.8)';
        ctx.font = '9.5px ' + Plot.CSS('--sans');
        ctx.textAlign = 'center';
        const label = meta.name.length * 6 < bw - 12 ? meta.name : (meta.id || '').slice(0, 5);
        ctx.fillText(label, bx + bw / 2, by + bh / 2);
      }
    }

    /* overlap warnings */
    const ev = State.lastEvaluation;
    if (ev && ev.construct && ev.construct.overlaps && ev.construct.overlaps.length) {
      ctx.fillStyle = Plot.CSS('--danger');
      ctx.font = '9.5px ' + Plot.CSS('--sans');
      ctx.textAlign = 'right'; ctx.textBaseline = 'bottom';
      ctx.fillText('⚠ overlapping elements overwrite each other', w - 6, h - 3);
    }

    /* interaction */
    canvas.style.cursor = track.drag ? 'grabbing' : 'default';
    canvas.onmousedown = (ev2) => {
      const r = canvas.getBoundingClientRect();
      const mx = ev2.clientX - r.left, my = ev2.clientY - r.top;
      for (const p of State.design.placements.slice().reverse()) {
        const wbp = State.elementWidth(p);
        const bx = X(p.position), bw = Math.max(3, wbp / bpPerPx);
        const by = boxTop + (track.lanes[p.uid] || 0) * laneH;
        if (mx >= bx - 2 && mx <= bx + bw + 2 && my >= by && my <= by + laneH - 4) {
          track.drag = { uid: p.uid, grabOffset: mx - bx, moved: false };
          State.design.selected = p.uid;
          renderList(); drawTrack();
          ev2.preventDefault();
          return;
        }
      }
      State.design.selected = null;
      renderList(); drawTrack();
    };

    canvas.onmousemove = (ev2) => {
      const r = canvas.getBoundingClientRect();
      const mx = ev2.clientX - r.left;
      if (!track.drag) {
        const my = ev2.clientY - r.top;
        let over = false;
        for (const p of State.design.placements) {
          const bx = X(p.position), bw = Math.max(3, State.elementWidth(p) / bpPerPx);
          const by = boxTop + (track.lanes[p.uid] || 0) * laneH;
          if (mx >= bx - 2 && mx <= bx + bw + 2 && my >= by && my <= by + laneH - 4) { over = true; break; }
        }
        canvas.style.cursor = over ? 'grab' : 'default';
        return;
      }
      const p = State.getPlacement(track.drag.uid);
      if (!p) return;
      const newPos = Math.round(x0 + (mx - track.drag.grabOffset - padL) * bpPerPx);
      if (newPos !== p.position) { p.position = newPos; track.drag.moved = true; drawTrack(); renderList(); }
    };

    const endDrag = () => {
      if (track.drag) {
        const moved = track.drag.moved;
        track.drag = null;
        canvas.style.cursor = 'default';
        if (moved) changed();
      }
    };
    canvas.onmouseup = endDrag;
    canvas.onmouseleave = endDrag;

    canvas.ondragover = (e) => { e.preventDefault(); e.dataTransfer.dropEffect = 'copy'; };
    canvas.ondrop = (e) => {
      e.preventDefault();
      const id = e.dataTransfer.getData('text/plain');
      if (!id) return;
      const r = canvas.getBoundingClientRect();
      const pos = Math.round(x0 + (e.clientX - r.left - padL) * bpPerPx);
      State.addPlacement(id, pos);
      changed();
    };

    canvas.onwheel = (e) => {
      e.preventDefault();
      const v = State.design.view;
      if (e.shiftKey) {
        v.center += e.deltaY * bpPerPx * 1.4;
      } else {
        const r = canvas.getBoundingClientRect();
        const anchor = x0 + (e.clientX - r.left - padL) * bpPerPx;
        const factor = Math.exp(e.deltaY * 0.0016);
        const maxHalf = Math.max(60, Math.floor(State.design.background.length / 2));
        v.halfWidth = Math.max(25, Math.min(maxHalf, v.halfWidth * factor));
        v.center = anchor - (anchor - v.center) * factor;
      }
      clampView();
      syncZoomUI();
      drawTrack();
    };
  }

  function clampView() {
    const v = State.design.view;
    const half = Math.floor(State.design.background.length / 2);
    v.halfWidth = Math.max(25, Math.min(half, v.halfWidth));
    v.center = Math.max(-half + v.halfWidth, Math.min(half - v.halfWidth, v.center));
  }

  function syncZoomUI() {
    const v = State.design.view;
    const half = Math.max(60, Math.floor(State.design.background.length / 2));
    const t = Math.log(v.halfWidth / 25) / Math.log(half / 25);
    $('#zoomRange').value = String(Math.round(t * 100));
    $('#zoomLabel').textContent = '±' + Math.round(v.halfWidth);
  }

  /* ================= evaluation ================= */
  function scheduleEvaluate() {
    clearTimeout(evalTimer);
    evalTimer = setTimeout(evaluate, 170);
  }

  async function evaluate() {
    if (busy) { pending = true; return; }
    busy = true;
    try {
      const half = Math.max(50, Math.floor(State.design.background.length / 2) - 1);
      const win = Math.min(1000, half);
      const res = await API.post('/api/evaluate',
        State.designPayload({ window: [-win, win] }));
      State.lastEvaluation = res;
      renderStats(res);
      renderProfiles(res);
      renderCellTypes(res);
      renderMetrics(res);
      renderSeqBox(res);
      drawTrack();
    } catch (e) {
      fail(e);
    } finally {
      busy = false;
      if (pending) { pending = false; scheduleEvaluate(); }
    }
  }

  function changed() {
    renderList();
    drawTrack();
    scheduleEvaluate();
  }

  /* ================= result panels ================= */
  function renderStats(res) {
    const box = $('#designStats');
    clear(box);
    const m = res.metrics || {};
    const s = res.score || {};
    const items = [
      ['Weighted score', fmtNum(s.score, 3), 'good'],
      ['Peak height', fmtNum(m.peak_height, 2), ''],
      ['Peak position', m.peak_position !== undefined ? (m.peak_position > 0 ? '+' : '') + m.peak_position + ' bp' : '–', ''],
      ['Focus ±50 bp', m.focus_fraction !== undefined ? (100 * m.focus_fraction).toFixed(0) + '%' : '–', ''],
      ['FWHM', m.fwhm !== undefined ? m.fwhm + ' bp' : '–', ''],
      ['Tau', fmtNum(m.tau, 3), m.tau > 0.4 ? 'good' : m.tau > 0.2 ? 'mid' : ''],
    ];
    for (const [lab, val, cls] of items) {
      box.appendChild(h('div', { class: 'stat' },
        h('div', { class: 'lab', text: lab }),
        h('div', { class: 'num ' + cls, text: val })));
    }
    $('#bgStats').textContent =
      `GC ${fmtNum(res.background.stats.gc, 3)} · CpG o/e ${fmtNum(res.background.stats.cpg_oe, 3)} · ` +
      `construct GC ${fmtNum(res.construct.stats.gc, 3)} · hash ${res.sequence_hash}`;
  }

  function renderProfiles(res) {
    const box = $('#profilePlots');
    clear(box);
    const [vx0, vx1] = viewRange();

    for (const [name, pred] of Object.entries(res.profiles || {})) {
      if (pred.error) {
        box.appendChild(h('div', { class: 'err', text: `${name}: ${pred.error}` }));
        continue;
      }
      const wrap = h('div', { class: 'plotwrap' });
      const canvas = h('canvas');
      const modelInfo = (State.models.profile || []).find(m => m.name === name) || {};
      box.appendChild(h('div', {},
        h('div', { class: 'plot-title' },
          modelInfo.label || name,
          h('span', { class: 'badge ' + (pred.is_mock ? 'mock' : 'real'),
                      style: { fontSize: '9px', padding: '1px 6px' },
                      text: pred.is_mock ? 'mock' : 'real' })),
        h('p', { class: 'plot-sub', text: `y: ${pred.scale}` }),
        wrap));
      wrap.appendChild(canvas);

      Plot.managed(canvas, () => Plot.line(canvas, {
        height: 150, legendRight: true,
        xrange: [vx0, vx1],
        xlabel: 'position relative to TSS (bp)',
        ylabel: 'predicted',
        vlines: [{ x: 0, label: 'TSS' }],
        series: [
          { x: pred.positions, y: pred.tracks.plus, color: '#4dd4ac', label: 'plus', fill: true },
          { x: pred.positions, y: pred.tracks.minus, color: '#58a6ff', label: 'minus', dash: [4, 3] },
        ],
      }));
    }
  }

  function renderCellTypes(res) {
    const ct = res.celltype;
    const canvas = $('#ctBars');
    const sum = $('#ctSummary');
    clear(sum);
    if (!ct || ct.error) {
      Plot.setup(canvas, 40);
      sum.appendChild(h('div', { class: 'err', text: (ct && ct.error) || 'no cell-type model' }));
      return;
    }
    const m = ct.metrics;
    const items = ct.cell_types.map(c => ({
      label: (ct.meta.cell_type_labels?.[c] || c).split(' (')[0],
      value: m.activities[c],
      highlight: c === m.strongest_cell_type,
      color: c === m.target ? Plot.CSS('--accent-2')
           : c === m.strongest_cell_type ? Plot.CSS('--accent') : '#3d4b5c',
    })).sort((a, b) => b.value - a.value);

    Plot.managed(canvas, () => Plot.bars(canvas, {
      items, xlabel: `activity (${ct.activity_method} over ${ct.activity_window[0]}…${ct.activity_window[1]} bp)`,
      decimals: 3, labelWidth: 108,
    }));

    sum.appendChild(h('dl', { class: 'kv', style: { marginTop: '8px' } },
      h('dt', { text: 'strongest' }), h('dd', { text: m.strongest_cell_type + '  (' + fmtNum(m.strongest_activity) + ')' }),
      h('dt', { text: 'tau' }), h('dd', { text: fmtNum(m.tau, 4) }),
      h('dt', { text: 'gini' }), h('dd', { text: fmtNum(m.gini, 4) }),
      ...(m.target ? [
        h('dt', { text: 'target' }), h('dd', { text: `${m.target} — rank ${m.target_rank}` }),
        h('dt', { text: 'target − best off' }), h('dd', { text: fmtNum(m.target_minus_offtarget_max, 4) }),
        h('dt', { text: 'log2 t/off-mean' }), h('dd', { text: fmtNum(m.log2_target_over_offtarget_mean, 3) }),
      ] : [])));
    sum.appendChild(h('div', { class: 'hint', text: ct.scale }));
  }

  function renderMetrics(res) {
    const box = $('#metricsTable');
    clear(box);
    const contrib = res.score.contributions || {};
    const rows = Object.entries(contrib).map(([k, c]) => h('tr', {},
      h('td', { text: c.label }),
      h('td', { class: 'num', text: fmtNum(c.raw, 3) }),
      h('td', { class: 'num', text: fmtNum(c.normalised, 3) }),
      h('td', { class: 'num', text: fmtNum(c.weight, 2) }),
      h('td', { class: 'num', text: fmtNum(c.contribution, 3) })));

    box.appendChild(h('table', { class: 'data' },
      h('thead', {}, h('tr', {},
        h('th', { text: 'metric' }), h('th', { class: 'num', text: 'raw' }),
        h('th', { class: 'num', text: 'norm' }), h('th', { class: 'num', text: 'weight' }),
        h('th', { class: 'num', text: 'contrib' }))),
      h('tbody', {}, rows.length ? rows :
        [h('tr', {}, h('td', { colspan: '5', text: 'no metrics weighted' }))]),
      h('tfoot', {}, h('tr', { class: 'hl' },
        h('td', { text: 'total' }), h('td', {}), h('td', {}),
        h('td', { class: 'num', text: fmtNum(res.score.weight_sum, 2) }),
        h('td', { class: 'num', text: fmtNum(res.score.score, 3) })))));

    box.appendChild(h('div', { class: 'hint', text:
      `Normalised value = metric mapped onto [0,1] over its declared range, ` +
      `flipped where lower is better. Score = Σ weight × normalised.` }));
  }

  function renderSeqBox(res) {
    const box = $('#seqBox');
    const seq = res.sequence;
    const t = res.background.tss_index;
    const lo = Math.max(0, t - 120), hi = Math.min(seq.length, t + 60);
    clear(box);
    box.appendChild(h('span', { text: seq.slice(lo, t) }));
    box.appendChild(h('span', { style: { color: 'var(--warn)', fontWeight: '700' }, text: seq[t] || '' }));
    box.appendChild(h('span', { text: seq.slice(t + 1, hi) }));
    box.appendChild(h('div', { class: 'hint', text:
      `showing ${lo - t}…${hi - t - 1} relative to the TSS (highlighted base = position 0); ` +
      `full construct is ${seq.length} bp` }));
  }

  /* ================= keyboard ================= */
  function onKey(e) {
    if (!$('#view-design').classList.contains('active')) return;
    const tag = (e.target.tagName || '').toLowerCase();
    if (tag === 'input' || tag === 'textarea' || tag === 'select') return;
    const p = State.design.selected && State.getPlacement(State.design.selected);
    if (!p) return;
    const step = e.shiftKey ? 10 : 1;
    if (e.key === 'ArrowLeft') { p.position -= step; changed(); e.preventDefault(); }
    else if (e.key === 'ArrowRight') { p.position += step; changed(); e.preventDefault(); }
    else if (e.key === 'r' || e.key === 'R') { p.strand = p.strand === '+' ? '-' : '+'; changed(); }
    else if (e.key === 'd' || e.key === 'D') { State.duplicatePlacement(p.uid); changed(); }
    else if (e.key === 'Delete' || e.key === 'Backspace') { State.removePlacement(p.uid); changed(); e.preventDefault(); }
  }

  /* ================= wiring ================= */
  function init() {
    renderPalette();
    renderList();
    renderWeights();
    clampView();
    syncZoomUI();

    const bg = State.design.background;
    const bind = (sel, key, cast, label) => {
      const el = $(sel);
      el.value = bg[key];
      el.addEventListener('input', () => {
        bg[key] = cast(el.value);
        if (label) $(label).textContent = Number(bg[key]).toFixed(2);
        if (key === 'length') { clampView(); syncZoomUI(); }
        changed();
      });
    };
    bind('#bgLength', 'length', v => Math.max(200, parseInt(v || '1001', 10)));
    bind('#bgGC', 'gc', parseFloat, '#bgGCv');
    bind('#bgCpG', 'cpg_oe', parseFloat, '#bgCpGv');
    bind('#bgSeed', 'seed', v => parseInt(v || '0', 10));

    $('#btnNewBg').onclick = () => {
      bg.seed = Math.floor(Math.random() * 1e6);
      $('#bgSeed').value = bg.seed; changed();
    };
    $('#btnRerollBg').onclick = $('#btnNewBg').onclick;
    $('#btnClearEls').onclick = () => { State.design.placements = []; State.design.selected = null; changed(); };
    $('#btnSnapTypical').onclick = () => {
      for (const p of State.design.placements) {
        const m = State.motifById[p.element_id];
        if (m && m.typical_offset !== undefined) p.position = m.typical_offset;
      }
      changed();
    };

    $('#zoomRange').addEventListener('input', (e) => {
      const half = Math.max(60, Math.floor(bg.length / 2));
      const t = parseInt(e.target.value, 10) / 100;
      State.design.view.halfWidth = Math.round(25 * Math.pow(half / 25, t));
      clampView(); syncZoomUI(); drawTrack(); renderProfiles(State.lastEvaluation || { profiles: {} });
    });
    $('#btnZoomFit').onclick = () => {
      State.design.view = { center: 0, halfWidth: Math.floor(bg.length / 2) };
      clampView(); syncZoomUI(); drawTrack();
      if (State.lastEvaluation) renderProfiles(State.lastEvaluation);
    };
    $('#btnZoomTSS').onclick = () => {
      State.design.view = { center: 0, halfWidth: 100 };
      syncZoomUI(); drawTrack();
      if (State.lastEvaluation) renderProfiles(State.lastEvaluation);
    };

    $('#targetCT').addEventListener('change', (e) => {
      State.design.targetCellType = e.target.value || null;
      scheduleEvaluate();
    });

    $('#weightPreset').addEventListener('change', (e) => {
      const v = e.target.value;
      if (v === 'custom') return;
      State.design.weights = Object.assign({}, State.defaultWeights[v] || {});
      renderWeights(); scheduleEvaluate();
    });

    /* exports */
    $('#btnExFasta').onclick = async () => {
      try {
        const r = await API.post('/api/export/fasta', State.designPayload({ name: 'design_' + stamp() }));
        download(`design_${stamp()}.fa`, r.fasta, 'text/plain');
      } catch (e) { fail(e); }
    };
    $('#btnExJson').onclick = () => {
      download(`design_${stamp()}.json`, JSON.stringify({
        settings: State.designPayload(), evaluation: State.lastEvaluation,
      }, null, 1), 'application/json');
    };
    $('#btnExProfiles').onclick = () => {
      const ev = State.lastEvaluation;
      if (!ev) return toast('nothing to export yet', true);
      const names = Object.keys(ev.profiles);
      const pos = ev.profiles[names[0]].positions;
      const rows = pos.map((p, i) => {
        const r = { position: p };
        for (const n of names) {
          r[`${n}_plus`] = ev.profiles[n].tracks.plus[i];
          r[`${n}_minus`] = ev.profiles[n].tracks.minus[i];
        }
        if (ev.celltype && ev.celltype.profiles) {
          ev.celltype.cell_types.forEach((c, ci) => { r[`ct_${c}`] = ev.celltype.profiles[ci][i]; });
        }
        return r;
      });
      download(`profiles_${stamp()}.csv`, toCSV(rows), 'text/csv');
    };
    $('#btnExBundle').onclick = async () => {
      try {
        const r = await API.post('/api/export/bundle', {
          name: 'design', settings: State.designPayload(), payload: State.lastEvaluation,
        });
        toast('saved ' + r.filename);
      } catch (e) { fail(e); }
    };

    document.addEventListener('keydown', onKey);
    Bus.on('design:changed', changed);
  }

  return { init, evaluate, changed, drawTrack, renderWeights, renderList,
           renderPalette, syncZoomUI, clampView, viewRange };
})();
