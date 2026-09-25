/* Minimal canvas plotting: line charts, heatmaps, bar charts.
 *
 * Written from scratch rather than pulled from a CDN so the app works on a
 * compute node with no outbound network. Everything is device-pixel-ratio
 * aware and re-renders on resize.
 */
const Plot = (() => {

  const CSS = (n) => getComputedStyle(document.documentElement)
    .getPropertyValue(n).trim();

  const PALETTE = ['#4dd4ac', '#58a6ff', '#f0a13a', '#f2545b', '#b06bd6',
                   '#2ec4b6', '#ff9f1c', '#8ecae6', '#e07a5f', '#81b29a',
                   '#c77dff', '#48cae4'];

  /* ---------- canvas setup ---------- */
  function setup(canvas, height) {
    const dpr = window.devicePixelRatio || 1;
    const w = canvas.clientWidth || canvas.parentElement.clientWidth || 600;
    const h = height || canvas.clientHeight || 220;
    canvas.width = Math.max(1, Math.round(w * dpr));
    canvas.height = Math.max(1, Math.round(h * dpr));
    canvas.style.height = h + 'px';
    const ctx = canvas.getContext('2d');
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, w, h);
    return { ctx, w, h };
  }

  /* ---------- nice axis ticks ---------- */
  function ticks(lo, hi, target = 6) {
    if (!isFinite(lo) || !isFinite(hi) || lo === hi) return [lo];
    const span = hi - lo;
    const raw = span / target;
    const mag = Math.pow(10, Math.floor(Math.log10(raw)));
    const norm = raw / mag;
    const step = (norm < 1.5 ? 1 : norm < 3 ? 2 : norm < 7 ? 5 : 10) * mag;
    const out = [];
    for (let v = Math.ceil(lo / step) * step; v <= hi + step * 1e-9; v += step) {
      out.push(Math.abs(v) < step * 1e-9 ? 0 : v);
    }
    return out;
  }

  function fmt(v, step) {
    if (v === 0) return '0';
    const a = Math.abs(v);
    if (a >= 1e5 || (a < 1e-3 && a > 0)) return v.toExponential(1);
    const dec = step !== undefined && step < 1
      ? Math.min(4, Math.max(0, -Math.floor(Math.log10(step))))
      : (a < 1 ? 2 : a < 10 ? 2 : a < 100 ? 1 : 0);
    return v.toFixed(dec);
  }

  /* isFinite(null) is true in JavaScript, so a null gap would be drawn as a
   * zero. Every value check goes through this instead. */
  function ok(v) { return v !== null && v !== undefined && isFinite(v); }

  function extent(arrs) {
    let lo = Infinity, hi = -Infinity;
    for (const a of arrs) for (const v of a) {
      if (v === null || v === undefined || !isFinite(v)) continue;
      if (v < lo) lo = v; if (v > hi) hi = v;
    }
    if (!isFinite(lo)) { lo = 0; hi = 1; }
    if (lo === hi) { lo -= 0.5; hi += 0.5; }
    return [lo, hi];
  }

  /* ---------- colour maps ---------- */
  function lerp(a, b, t) { return a + (b - a) * t; }
  function rgb(c) { return `rgb(${c[0]|0},${c[1]|0},${c[2]|0})`; }

  const VIRIDIS = [[68,1,84],[72,40,120],[62,74,137],[49,104,142],[38,130,142],
                   [31,158,137],[53,183,121],[109,205,89],[180,222,44],[253,231,37]];
  const MAGMA = [[0,0,4],[28,16,68],[79,18,123],[129,37,129],[181,54,122],
                 [229,80,100],[251,135,97],[254,194,135],[252,253,191]];

  function rampColor(ramp, t) {
    t = Math.max(0, Math.min(1, t));
    const x = t * (ramp.length - 1);
    const i = Math.min(ramp.length - 2, Math.floor(x));
    const f = x - i;
    return rgb([lerp(ramp[i][0], ramp[i+1][0], f),
                lerp(ramp[i][1], ramp[i+1][1], f),
                lerp(ramp[i][2], ramp[i+1][2], f)]);
  }

  /* Diverging blue-grey-red, centred on zero. Used for interaction residuals,
   * where the sign is the whole point. */
  function diverging(t) {
    t = Math.max(-1, Math.min(1, t));
    const neg = [56, 122, 223], mid = [30, 36, 44], pos = [235, 78, 86];
    const a = t < 0 ? neg : pos;
    const k = Math.abs(t);
    return rgb([lerp(mid[0], a[0], k), lerp(mid[1], a[1], k), lerp(mid[2], a[2], k)]);
  }

  function colorFor(scheme, t) {
    if (scheme === 'magma') return rampColor(MAGMA, t);
    if (scheme === 'diverging') return diverging(t);
    return rampColor(VIRIDIS, t);
  }

  /* ---------- tooltip ---------- */
  function tooltipFor(canvas) {
    const wrap = canvas.closest('.plotwrap') || canvas.parentElement;
    if (getComputedStyle(wrap).position === 'static') wrap.style.position = 'relative';
    let tip = wrap.querySelector(':scope > .tooltip');
    if (!tip) { tip = document.createElement('div'); tip.className = 'tooltip'; wrap.appendChild(tip); }
    return tip;
  }

  /* =====================================================================
   * Line chart
   * ===================================================================== */
  function line(canvas, cfg) {
    const height = cfg.height || 220;
    const { ctx, w, h } = setup(canvas, height);
    const series = (cfg.series || []).filter(s => s && s.y && s.y.length);

    const pad = Object.assign(
      { l: 54, r: cfg.legendRight ? 128 : 12, t: 10, b: 34 }, cfg.pad || {});
    const PW = Math.max(10, w - pad.l - pad.r);
    const PH = Math.max(10, h - pad.t - pad.b);

    const xs = series.map(s => s.x);
    const ys = series.map(s => s.y);
    let [x0, x1] = cfg.xrange || extent(xs);
    let [y0, y1] = cfg.yrange || extent(ys.concat(cfg.extraY ? [cfg.extraY] : []));
    if (cfg.yZero) { y0 = Math.min(0, y0); y1 = Math.max(0, y1); }
    if (cfg.ySymmetric) { const m = Math.max(Math.abs(y0), Math.abs(y1)); y0 = -m; y1 = m; }
    const yPad = (y1 - y0) * 0.06; y0 -= yPad; y1 += yPad;

    const X = v => pad.l + (v - x0) / (x1 - x0) * PW;
    const Y = v => pad.t + PH - (v - y0) / (y1 - y0) * PH;

    /* grid + axes */
    const xt = ticks(x0, x1, cfg.xticks || 7);
    const yt = ticks(y0, y1, 5);
    ctx.font = '10px ' + CSS('--mono');
    ctx.lineWidth = 1;

    ctx.strokeStyle = CSS('--line-soft');
    for (const t of yt) { const y = Y(t); ctx.beginPath(); ctx.moveTo(pad.l, y); ctx.lineTo(pad.l + PW, y); ctx.stroke(); }
    for (const t of xt) { const x = X(t); ctx.beginPath(); ctx.moveTo(x, pad.t); ctx.lineTo(x, pad.t + PH); ctx.stroke(); }

    /* shaded bands */
    for (const b of (cfg.bands || [])) {
      ctx.fillStyle = b.color || 'rgba(88,166,255,.07)';
      const bx0 = X(Math.max(b.from, x0)), bx1 = X(Math.min(b.to, x1));
      ctx.fillRect(bx0, pad.t, Math.max(0, bx1 - bx0), PH);
    }

    /* error ribbons (drawn under the lines) */
    for (const s of series) {
      if (!s.sd) continue;
      ctx.fillStyle = (s.color || PALETTE[0]) + '22';
      ctx.beginPath();
      let started = false;
      for (let i = 0; i < s.x.length; i++) {
        if (!ok(s.y[i])) continue;
        const px = X(s.x[i]), py = Y(s.y[i] + (s.sd[i] || 0));
        started ? ctx.lineTo(px, py) : (ctx.moveTo(px, py), started = true);
      }
      for (let i = s.x.length - 1; i >= 0; i--) {
        if (!ok(s.y[i])) continue;
        ctx.lineTo(X(s.x[i]), Y(s.y[i] - (s.sd[i] || 0)));
      }
      ctx.closePath(); ctx.fill();
    }

    /* zero line */
    if (y0 < 0 && y1 > 0) {
      ctx.strokeStyle = CSS('--line'); ctx.setLineDash([]);
      ctx.beginPath(); ctx.moveTo(pad.l, Y(0)); ctx.lineTo(pad.l + PW, Y(0)); ctx.stroke();
    }

    /* vertical markers (TSS etc.) */
    for (const m of (cfg.vlines || [])) {
      if (m.x < x0 || m.x > x1) continue;
      ctx.strokeStyle = m.color || CSS('--warn');
      ctx.setLineDash(m.dash || [4, 3]);
      ctx.beginPath(); ctx.moveTo(X(m.x), pad.t); ctx.lineTo(X(m.x), pad.t + PH); ctx.stroke();
      ctx.setLineDash([]);
      if (m.label) {
        ctx.fillStyle = m.color || CSS('--warn');
        ctx.textAlign = 'left'; ctx.textBaseline = 'top';
        ctx.fillText(m.label, X(m.x) + 3, pad.t + 2);
      }
    }

    /* series */
    for (const s of series) {
      ctx.strokeStyle = s.color || PALETTE[0];
      ctx.lineWidth = s.width || 1.6;
      ctx.setLineDash(s.dash || []);
      if (s.fill) {
        ctx.fillStyle = (s.color || PALETTE[0]) + '1f';
        ctx.beginPath(); ctx.moveTo(X(s.x[0]), Y(Math.max(y0, 0)));
        for (let i = 0; i < s.x.length; i++) if (ok(s.y[i])) ctx.lineTo(X(s.x[i]), Y(s.y[i]));
        ctx.lineTo(X(s.x[s.x.length - 1]), Y(Math.max(y0, 0)));
        ctx.closePath(); ctx.fill();
      }
      ctx.beginPath();
      let started = false;
      for (let i = 0; i < s.x.length; i++) {
        if (!ok(s.y[i])) { started = false; continue; }
        const px = X(s.x[i]), py = Y(s.y[i]);
        started ? ctx.lineTo(px, py) : (ctx.moveTo(px, py), started = true);
      }
      ctx.stroke();
      ctx.setLineDash([]);
      if (s.points) {
        ctx.fillStyle = s.color || PALETTE[0];
        for (let i = 0; i < s.x.length; i++) {
          if (!ok(s.y[i])) continue;
          ctx.beginPath(); ctx.arc(X(s.x[i]), Y(s.y[i]), s.pointSize || 2.2, 0, 7); ctx.fill();
        }
      }
    }

    /* axis frame + labels */
    ctx.setLineDash([]);
    ctx.strokeStyle = CSS('--line'); ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(pad.l, pad.t); ctx.lineTo(pad.l, pad.t + PH);
    ctx.lineTo(pad.l + PW, pad.t + PH); ctx.stroke();

    ctx.fillStyle = CSS('--fg-faint');
    ctx.textAlign = 'right'; ctx.textBaseline = 'middle';
    const ystep = yt.length > 1 ? Math.abs(yt[1] - yt[0]) : undefined;
    for (const t of yt) { if (t < y0 || t > y1) continue; ctx.fillText(fmt(t, ystep), pad.l - 6, Y(t)); }
    ctx.textAlign = 'center'; ctx.textBaseline = 'top';
    const xstep = xt.length > 1 ? Math.abs(xt[1] - xt[0]) : undefined;
    for (const t of xt) { if (t < x0 || t > x1) continue; ctx.fillText(fmt(t, xstep), X(t), pad.t + PH + 6); }

    ctx.fillStyle = CSS('--fg-dim'); ctx.font = '10.5px ' + CSS('--sans');
    if (cfg.xlabel) { ctx.textAlign = 'center'; ctx.textBaseline = 'bottom';
                      ctx.fillText(cfg.xlabel, pad.l + PW / 2, h - 1); }
    if (cfg.ylabel) { ctx.save(); ctx.translate(11, pad.t + PH / 2); ctx.rotate(-Math.PI / 2);
                      ctx.textAlign = 'center'; ctx.textBaseline = 'top';
                      ctx.fillText(cfg.ylabel, 0, 0); ctx.restore(); }

    /* right-hand legend */
    if (cfg.legendRight) {
      let ly = pad.t + 4;
      ctx.font = '10px ' + CSS('--sans'); ctx.textAlign = 'left'; ctx.textBaseline = 'middle';
      for (const s of series) {
        if (s.hideLegend) continue;
        ctx.strokeStyle = s.color || PALETTE[0]; ctx.lineWidth = 2.4;
        ctx.setLineDash(s.dash || []);
        ctx.beginPath(); ctx.moveTo(pad.l + PW + 8, ly); ctx.lineTo(pad.l + PW + 22, ly); ctx.stroke();
        ctx.setLineDash([]);
        ctx.fillStyle = CSS('--fg-dim');
        const lbl = (s.label || '').slice(0, 18);
        ctx.fillText(lbl, pad.l + PW + 27, ly);
        ly += 14;
      }
    }

    /* hover readout */
    if (cfg.hover !== false && series.length) {
      const tip = tooltipFor(canvas);
      canvas.onmousemove = (ev) => {
        const r = canvas.getBoundingClientRect();
        const mx = ev.clientX - r.left, my = ev.clientY - r.top;
        if (mx < pad.l || mx > pad.l + PW || my < pad.t || my > pad.t + PH) {
          tip.classList.remove('on'); return;
        }
        const xv = x0 + (mx - pad.l) / PW * (x1 - x0);
        const ref = series[0];
        let bi = 0, bd = Infinity;
        for (let i = 0; i < ref.x.length; i++) {
          const d = Math.abs(ref.x[i] - xv);
          if (d < bd) { bd = d; bi = i; }
        }
        const lines = [`${cfg.xlabel || 'x'} = ${ref.x[bi]}`];
        for (const s of series) {
          if (s.hideLegend) continue;
          const v = s.y[bi];
          if (!ok(v)) { lines.push(`${(s.label || '?')}: excluded`); continue; }
          lines.push(`${(s.label || '?')}: ${v.toFixed(4)}`);
        }
        tip.textContent = lines.join('\n');
        tip.classList.add('on');
        const tw = tip.offsetWidth, th = tip.offsetHeight;
        tip.style.left = Math.min(w - tw - 4, Math.max(2, mx + 12)) + 'px';
        tip.style.top = Math.max(2, my - th - 8) + 'px';
      };
      canvas.onmouseleave = () => tip.classList.remove('on');
    }

    return { X, Y, pad, PW, PH, x0, x1, y0, y1 };
  }

  /* =====================================================================
   * Heatmap
   * ===================================================================== */
  function heatmap(canvas, cfg) {
    const m = cfg.matrix || [];
    const rows = m.length, cols = rows ? m[0].length : 0;
    if (!rows || !cols) { setup(canvas, cfg.height || 200); return; }

    const cell = cfg.cell || 0;
    const pad = Object.assign({ l: cfg.rowLabels ? 84 : 44, r: 62, t: 10,
                                b: cfg.colLabels ? 68 : 34 }, cfg.pad || {});
    const wantW = canvas.clientWidth || 600;
    const PW = Math.max(40, wantW - pad.l - pad.r);
    const PH = cell ? cell * rows : (cfg.plotHeight || Math.min(420, Math.max(120, cols ? PW * rows / cols : 200)));
    const h = PH + pad.t + pad.b;
    const { ctx, w } = setup(canvas, h);

    let flat = [];
    for (const r of m) for (const v of r) if (ok(v)) flat.push(v);
    let lo = cfg.vmin, hi = cfg.vmax;
    if (lo === undefined || hi === undefined) {
      const e = extent([flat]);
      lo = lo === undefined ? e[0] : lo;
      hi = hi === undefined ? e[1] : hi;
    }
    const div = cfg.scheme === 'diverging';
    if (div) { const mx = Math.max(Math.abs(lo), Math.abs(hi)) || 1; lo = -mx; hi = mx; }

    const cw = PW / cols, ch = PH / rows;
    for (let i = 0; i < rows; i++) {
      for (let j = 0; j < cols; j++) {
        const v = m[i][j];
        if (!ok(v)) { ctx.fillStyle = CSS('--panel'); }
        else {
          const t = div ? (hi ? v / hi : 0) : (hi === lo ? 0.5 : (v - lo) / (hi - lo));
          ctx.fillStyle = colorFor(cfg.scheme || 'viridis', t);
        }
        ctx.fillRect(pad.l + j * cw, pad.t + i * ch, Math.ceil(cw) , Math.ceil(ch));
      }
    }

    /* cell values, when there is room */
    if (cfg.annotate !== false && cw > 34 && ch > 17) {
      ctx.font = '9px ' + CSS('--mono'); ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
      for (let i = 0; i < rows; i++) for (let j = 0; j < cols; j++) {
        const v = m[i][j];
        if (!ok(v)) continue;
        const t = div ? Math.abs(hi ? v / hi : 0) : (hi === lo ? 0.5 : (v - lo) / (hi - lo));
        ctx.fillStyle = t > 0.55 ? 'rgba(0,0,0,.72)' : 'rgba(255,255,255,.78)';
        ctx.fillText(fmt(v), pad.l + (j + 0.5) * cw, pad.t + (i + 0.5) * ch);
      }
    }

    ctx.strokeStyle = CSS('--line'); ctx.lineWidth = 1;
    ctx.strokeRect(pad.l, pad.t, PW, PH);

    /* labels */
    ctx.fillStyle = CSS('--fg-faint'); ctx.font = '10px ' + CSS('--sans');
    if (cfg.rowLabels) {
      ctx.textAlign = 'right'; ctx.textBaseline = 'middle';
      for (let i = 0; i < rows; i++) {
        if (ch < 9) break;
        ctx.fillText(String(cfg.rowLabels[i]).slice(0, 14), pad.l - 6, pad.t + (i + 0.5) * ch);
      }
    }
    if (cfg.colLabels) {
      ctx.textAlign = 'right'; ctx.textBaseline = 'middle';
      for (let j = 0; j < cols; j++) {
        if (cw < 9) break;
        ctx.save();
        ctx.translate(pad.l + (j + 0.5) * cw, pad.t + PH + 6);
        ctx.rotate(-Math.PI / 4);
        ctx.fillText(String(cfg.colLabels[j]).slice(0, 14), 0, 0);
        ctx.restore();
      }
    } else if (cfg.xticksFrom) {
      ctx.textAlign = 'center'; ctx.textBaseline = 'top';
      const vals = cfg.xticksFrom;
      const t = ticks(vals[0], vals[vals.length - 1], 7);
      for (const tv of t) {
        const frac = (tv - vals[0]) / (vals[vals.length - 1] - vals[0]);
        if (frac < 0 || frac > 1) continue;
        ctx.fillText(fmt(tv), pad.l + frac * PW, pad.t + PH + 6);
      }
    }

    ctx.fillStyle = CSS('--fg-dim'); ctx.font = '10.5px ' + CSS('--sans');
    if (cfg.xlabel) { ctx.textAlign = 'center'; ctx.textBaseline = 'bottom';
                      ctx.fillText(cfg.xlabel, pad.l + PW / 2, h - 1); }
    if (cfg.ylabel) { ctx.save(); ctx.translate(11, pad.t + PH / 2); ctx.rotate(-Math.PI / 2);
                      ctx.textAlign = 'center'; ctx.textBaseline = 'top';
                      ctx.fillText(cfg.ylabel, 0, 0); ctx.restore(); }

    /* colour bar */
    const cbx = pad.l + PW + 16, cbw = 12, cbh = Math.min(PH, 150);
    const cby = pad.t + (PH - cbh) / 2;
    for (let k = 0; k < cbh; k++) {
      const t = 1 - k / (cbh - 1);
      ctx.fillStyle = colorFor(cfg.scheme || 'viridis', div ? (t * 2 - 1) : t);
      ctx.fillRect(cbx, cby + k, cbw, 1);
    }
    ctx.strokeStyle = CSS('--line'); ctx.strokeRect(cbx, cby, cbw, cbh);
    ctx.fillStyle = CSS('--fg-faint'); ctx.font = '9px ' + CSS('--mono');
    ctx.textAlign = 'left'; ctx.textBaseline = 'middle';
    ctx.fillText(fmt(hi), cbx + cbw + 4, cby);
    ctx.fillText(fmt(div ? 0 : (lo + hi) / 2), cbx + cbw + 4, cby + cbh / 2);
    ctx.fillText(fmt(lo), cbx + cbw + 4, cby + cbh);

    /* hover */
    const tip = tooltipFor(canvas);
    canvas.onmousemove = (ev) => {
      const r = canvas.getBoundingClientRect();
      const mx = ev.clientX - r.left, my = ev.clientY - r.top;
      const j = Math.floor((mx - pad.l) / cw), i = Math.floor((my - pad.t) / ch);
      if (i < 0 || i >= rows || j < 0 || j >= cols) { tip.classList.remove('on'); return; }
      const rl = cfg.rowLabels ? cfg.rowLabels[i] : i;
      const cl = cfg.colLabels ? cfg.colLabels[j] : (cfg.xticksFrom ? cfg.xticksFrom[j] : j);
      const v = m[i][j];
      tip.textContent = `${rl}  ×  ${cl}\n${cfg.valueLabel || 'value'}: ` +
        (!ok(v) ? 'n/a' : v.toFixed(5));
      tip.classList.add('on');
      tip.style.left = Math.min(w - tip.offsetWidth - 4, Math.max(2, mx + 12)) + 'px';
      tip.style.top = Math.max(2, my - tip.offsetHeight - 8) + 'px';
    };
    canvas.onmouseleave = () => tip.classList.remove('on');
    canvas.onclick = cfg.onCell ? (ev) => {
      const r = canvas.getBoundingClientRect();
      const j = Math.floor((ev.clientX - r.left - pad.l) / cw);
      const i = Math.floor((ev.clientY - r.top - pad.t) / ch);
      if (i >= 0 && i < rows && j >= 0 && j < cols) cfg.onCell(i, j, m[i][j]);
    } : null;
  }

  /* =====================================================================
   * Horizontal bar chart
   * ===================================================================== */
  function bars(canvas, cfg) {
    const items = cfg.items || [];
    const rowH = cfg.rowH || 17;
    const pad = Object.assign({ l: cfg.labelWidth || 104, r: 54, t: 6, b: 22 }, cfg.pad || {});
    const h = items.length * rowH + pad.t + pad.b;
    const { ctx, w } = setup(canvas, h);
    const PW = Math.max(30, w - pad.l - pad.r);

    const vals = items.map(i => i.value);
    let hi = cfg.vmax !== undefined ? cfg.vmax : Math.max(0, ...vals.filter(ok));
    let lo = cfg.vmin !== undefined ? cfg.vmin : Math.min(0, ...vals.filter(ok));
    if (hi === lo) hi = lo + 1;
    const X = v => pad.l + (v - lo) / (hi - lo) * PW;

    ctx.font = '10px ' + CSS('--mono');
    const xt = ticks(lo, hi, 5);
    ctx.strokeStyle = CSS('--line-soft');
    for (const t of xt) { ctx.beginPath(); ctx.moveTo(X(t), pad.t); ctx.lineTo(X(t), pad.t + items.length * rowH); ctx.stroke(); }

    items.forEach((it, i) => {
      const y = pad.t + i * rowH;
      const x = X(Math.min(it.value, 0) === it.value && lo < 0 ? it.value : Math.max(lo, 0));
      const xv = X(it.value);
      ctx.fillStyle = it.color || (it.highlight ? CSS('--accent') : CSS('--accent-2'));
      const bx = Math.min(x, xv), bw = Math.abs(xv - x);
      ctx.globalAlpha = it.dim ? 0.42 : 1;
      ctx.fillRect(bx, y + 2.5, Math.max(1, bw), rowH - 5);
      ctx.globalAlpha = 1;

      ctx.fillStyle = it.highlight ? CSS('--fg') : CSS('--fg-dim');
      ctx.font = (it.highlight ? 'bold ' : '') + '10.5px ' + CSS('--sans');
      ctx.textAlign = 'right'; ctx.textBaseline = 'middle';
      ctx.fillText(String(it.label).slice(0, 20), pad.l - 6, y + rowH / 2);

      ctx.fillStyle = CSS('--fg-faint'); ctx.font = '10px ' + CSS('--mono');
      ctx.textAlign = 'left';
      ctx.fillText(ok(it.value) ? it.value.toFixed(cfg.decimals ?? 3) : 'n/a',
                   pad.l + PW + 6, y + rowH / 2);
    });

    if (cfg.marker !== undefined && ok(cfg.marker)) {
      ctx.strokeStyle = CSS('--warn'); ctx.setLineDash([3, 3]);
      ctx.beginPath(); ctx.moveTo(X(cfg.marker), pad.t);
      ctx.lineTo(X(cfg.marker), pad.t + items.length * rowH); ctx.stroke();
      ctx.setLineDash([]);
    }

    ctx.fillStyle = CSS('--fg-faint'); ctx.font = '9.5px ' + CSS('--mono');
    ctx.textAlign = 'center'; ctx.textBaseline = 'top';
    for (const t of xt) ctx.fillText(fmt(t), X(t), pad.t + items.length * rowH + 4);
    if (cfg.xlabel) {
      ctx.fillStyle = CSS('--fg-dim'); ctx.font = '10px ' + CSS('--sans');
      ctx.fillText(cfg.xlabel, pad.l + PW / 2, pad.t + items.length * rowH + 14);
    }
  }

  /* =====================================================================
   * Scatter, with a least-squares line. Used for model-vs-model agreement.
   * ===================================================================== */
  function scatter(canvas, cfg) {
    const height = cfg.height || 220;
    const { ctx, w, h } = setup(canvas, height);
    const pad = Object.assign({ l: 56, r: 14, t: 18, b: 40 }, cfg.pad || {});
    const PW = Math.max(10, w - pad.l - pad.r);
    const PH = Math.max(10, h - pad.t - pad.b);

    const xs = [], ys = [];
    const n = Math.min(cfg.x.length, cfg.y.length);
    for (let i = 0; i < n; i++) {
      if (ok(cfg.x[i]) && ok(cfg.y[i])) { xs.push(cfg.x[i]); ys.push(cfg.y[i]); }
    }
    if (!xs.length) return;

    let [x0, x1] = extent([xs]);
    let [y0, y1] = extent([ys]);
    const xp = (x1 - x0) * 0.05, yp = (y1 - y0) * 0.05;
    x0 -= xp; x1 += xp; y0 -= yp; y1 += yp;
    const X = v => pad.l + (v - x0) / (x1 - x0) * PW;
    const Y = v => pad.t + PH - (v - y0) / (y1 - y0) * PH;

    ctx.strokeStyle = CSS('--line-soft');
    const xt = ticks(x0, x1, 5), yt = ticks(y0, y1, 5);
    for (const t of yt) { const y = Y(t); ctx.beginPath(); ctx.moveTo(pad.l, y); ctx.lineTo(pad.l + PW, y); ctx.stroke(); }
    for (const t of xt) { const x = X(t); ctx.beginPath(); ctx.moveTo(x, pad.t); ctx.lineTo(x, pad.t + PH); ctx.stroke(); }

    ctx.fillStyle = (cfg.color || PALETTE[0]) + 'aa';
    const r = xs.length > 600 ? 1.2 : xs.length > 200 ? 1.8 : 2.4;
    for (let i = 0; i < xs.length; i++) {
      ctx.beginPath(); ctx.arc(X(xs[i]), Y(ys[i]), r, 0, 7); ctx.fill();
    }

    /* least-squares fit */
    const mx = xs.reduce((a, b) => a + b, 0) / xs.length;
    const my = ys.reduce((a, b) => a + b, 0) / ys.length;
    let sxy = 0, sxx = 0;
    for (let i = 0; i < xs.length; i++) { sxy += (xs[i] - mx) * (ys[i] - my); sxx += (xs[i] - mx) ** 2; }
    if (sxx > 0) {
      const b1 = sxy / sxx, b0 = my - b1 * mx;
      ctx.strokeStyle = CSS('--warn'); ctx.lineWidth = 1.4; ctx.setLineDash([5, 3]);
      ctx.beginPath(); ctx.moveTo(X(x0), Y(b0 + b1 * x0)); ctx.lineTo(X(x1), Y(b0 + b1 * x1)); ctx.stroke();
      ctx.setLineDash([]); ctx.lineWidth = 1;
    }

    ctx.strokeStyle = CSS('--line');
    ctx.beginPath(); ctx.moveTo(pad.l, pad.t); ctx.lineTo(pad.l, pad.t + PH);
    ctx.lineTo(pad.l + PW, pad.t + PH); ctx.stroke();

    ctx.fillStyle = CSS('--fg-faint'); ctx.font = '9.5px ' + CSS('--mono');
    ctx.textAlign = 'right'; ctx.textBaseline = 'middle';
    for (const t of yt) { if (t >= y0 && t <= y1) ctx.fillText(fmt(t), pad.l - 6, Y(t)); }
    ctx.textAlign = 'center'; ctx.textBaseline = 'top';
    for (const t of xt) { if (t >= x0 && t <= x1) ctx.fillText(fmt(t), X(t), pad.t + PH + 5); }

    ctx.fillStyle = CSS('--fg-dim'); ctx.font = '10px ' + CSS('--sans');
    if (cfg.xlabel) { ctx.textAlign = 'center'; ctx.textBaseline = 'bottom';
                      ctx.fillText(String(cfg.xlabel).slice(0, 46), pad.l + PW / 2, h - 1); }
    if (cfg.ylabel) { ctx.save(); ctx.translate(12, pad.t + PH / 2); ctx.rotate(-Math.PI / 2);
                      ctx.textAlign = 'center'; ctx.textBaseline = 'top';
                      ctx.fillText(String(cfg.ylabel).slice(0, 34), 0, 0); ctx.restore(); }
    if (cfg.title) { ctx.fillStyle = CSS('--accent'); ctx.font = 'bold 11px ' + CSS('--mono');
                     ctx.textAlign = 'right'; ctx.textBaseline = 'top';
                     ctx.fillText(cfg.title, pad.l + PW, 2); }
  }

  /* ---------- re-render on resize ---------- */
  const _redraw = new Map();
  function managed(canvas, fn) {
    _redraw.set(canvas, fn);
    fn();
  }
  let _rt = null;
  window.addEventListener('resize', () => {
    clearTimeout(_rt);
    _rt = setTimeout(() => {
      for (const [c, fn] of _redraw) {
        if (!c.isConnected) { _redraw.delete(c); continue; }
        if (c.offsetParent === null) continue;   // hidden tab: redraw on show
        try { fn(); } catch (e) { /* ignore */ }
      }
    }, 130);
  });
  function redrawVisible() {
    for (const [c, fn] of _redraw) {
      if (c.isConnected && c.offsetParent !== null) { try { fn(); } catch (e) {} }
    }
  }

  return { line, heatmap, bars, scatter, setup, ticks, fmt, extent, ok, colorFor, diverging,
           PALETTE, managed, redrawVisible, CSS };
})();
