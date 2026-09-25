/* Bootstrap: load metadata, wire the tabs, render the About page. */

const App = (() => {

  function showView(name) {
    $$('nav.tabs button').forEach(b => b.classList.toggle('active', b.dataset.view === name));
    $$('.view').forEach(v => v.classList.toggle('active', v.id === 'view-' + name));
    setTimeout(() => { Plot.redrawVisible(); if (name === 'design') Designer.drawTrack(); }, 20);
  }

  function fillSelect(sel, items, valueKey, labelKey, selected) {
    const el = typeof sel === 'string' ? $(sel) : sel;
    clear(el);
    for (const it of items) {
      el.appendChild(h('option', { value: it[valueKey], text: it[labelKey] }));
    }
    if (selected !== undefined && selected !== null) el.value = selected;
    return el;
  }

  async function boot() {
    try {
      const [motifs, models, metrics] = await Promise.all([
        API.get('/api/motifs'), API.get('/api/models'), API.get('/api/metrics'),
      ]);

      State.motifs = motifs.motifs;
      State.motifById = Object.fromEntries(motifs.motifs.map(m => [m.id, m]));
      State.extraElements = motifs.extra_elements;
      State.provenance = motifs.provenance;
      State.models = models;
      State.metricCatalogue = metrics.metrics;
      State.defaultWeights = metrics.default_weights;

      State.design.weights = Object.assign({}, metrics.default_weights.design);
      State.design.profileModels = models.defaults.profile.filter(
        n => (models.profile.find(m => m.name === n) || {}).available);
      State.design.celltypeModel = models.defaults.celltype;

      const ct = await API.get('/api/celltypes?model=' + encodeURIComponent(State.design.celltypeModel));
      State.cellTypes = ct.cell_types;

      /* selects that depend on metadata */
      const ctOptions = [{ id: '', label: '— none —' }].concat(State.cellTypes);
      fillSelect('#targetCT', ctOptions, 'id', 'label', '');
      fillSelect('#m1Target', State.cellTypes, 'id', 'label', State.cellTypes[0]?.id);
      fillSelect('#m1Element',
        State.motifs.map(m => ({ id: m.id, label: m.name }))
          .concat(State.extraElements.map(e => ({ id: e.id, label: e.name }))),
        'id', 'label', 'ets');
      fillSelect('#m1Model',
        models.celltype.filter(m => m.available).map(m => ({ id: m.name, label: m.label })),
        'id', 'label', State.design.celltypeModel);

      /* badges */
      const anyReal = models.any_real_available;
      const badge = $('#mockBadge');
      badge.className = 'badge ' + (anyReal ? 'warn dot' : 'mock dot');
      badge.textContent = anyReal ? 'mixed: mock + real' : 'mock predictions';
      badge.title = anyReal
        ? 'Some adapters are real; check each plot for its own badge.'
        : 'No trained model is in use. Every number on screen comes from a hand-written stand-in.';

      const health = await API.get('/api/health');
      $('#serverStatus').textContent = 'connected';
      $('#serverStatus').style.color = 'var(--accent)';

      /* player name */
      const saved = localStorage.getItem('pdg_player') || '';
      $('#playerName').value = saved;
      $('#playerName').addEventListener('change', async () => {
        const name = $('#playerName').value.trim();
        localStorage.setItem('pdg_player', name);
        if (!name) return;
        try { State.player = await API.post('/api/players', { name }); } catch (e) { fail(e); }
      });
      if (saved) { try { State.player = await API.post('/api/players', { name: saved }); } catch (e) {} }

      /* modules */
      Designer.init();
      Mode1.init();
      Mode2.init();
      Board.init();
      renderAbout();

      $$('nav.tabs button').forEach(b => b.onclick = () => showView(b.dataset.view));

      /* a starting design so the first screen is not empty */
      State.addPlacement('tata', -31);
      State.addPlacement('inr', -2);
      State.addPlacement('sp1', -75);
      Designer.changed();

    } catch (e) {
      $('#serverStatus').textContent = 'offline';
      $('#serverStatus').style.color = 'var(--danger)';
      fail(e);
    }
  }

  function renderAbout() {
    const box = $('#aboutContent');
    clear(box);

    box.appendChild(h('div', { class: 'note mock' },
      h('strong', { text: 'Everything on screen is a mock prediction. ' }),
      'No trained model is running. The adapters below generate plausible-looking ' +
      'profiles from hand-written rules so that the interface and the experiment ' +
      'designs can be built and checked first. Do not read any number here as a ' +
      'statement about real promoters.'));

    box.appendChild(h('div', { class: 'panel' },
      h('h3', {}, 'Motif library'),
      h('p', { style: { marginTop: 0 }, text: State.provenance }),
      h('table', { class: 'data' },
        h('thead', {}, h('tr', {},
          h('th', { text: 'motif' }), h('th', { text: 'consensus' }),
          h('th', { class: 'num', text: 'width' }), h('th', { class: 'num', text: 'typical 5′' }),
          h('th', { text: 'notes' }))),
        h('tbody', {}, ...State.motifs.map(m => h('tr', {},
          h('td', {}, h('span', { style: { color: m.color, fontWeight: '600' }, text: m.name })),
          h('td', { class: 'num', text: m.consensus }),
          h('td', { class: 'num', text: m.width }),
          h('td', { class: 'num', text: (m.typical_offset > 0 ? '+' : '') + m.typical_offset }),
          h('td', { style: { fontSize: '11px', color: 'var(--fg-dim)' }, text: m.notes })))))));

    const modelRows = (State.models.profile || []).concat(State.models.celltype || []);
    box.appendChild(h('div', { class: 'panel' },
      h('h3', {}, 'Model adapters'),
      h('table', { class: 'data' },
        h('thead', {}, h('tr', {},
          h('th', { text: 'adapter' }), h('th', { text: 'kind' }), h('th', { text: 'status' }),
          h('th', { text: 'output scale' }), h('th', { text: 'what it does' }))),
        h('tbody', {}, ...modelRows.map(m => h('tr', {},
          h('td', {}, m.label,
            h('span', { class: 'badge ' + (m.is_mock ? 'mock' : 'real'),
                        style: { marginLeft: '6px', fontSize: '9px', padding: '1px 5px' },
                        text: m.is_mock ? 'mock' : 'real' })),
          h('td', { text: m.kind }),
          h('td', { text: m.available ? 'available' : (m.unavailable_reason || 'unavailable'),
                    style: { color: m.available ? 'var(--accent)' : 'var(--fg-faint)', fontSize: '11px' } }),
          h('td', { style: { fontSize: '11px' }, text: m.scale }),
          h('td', { style: { fontSize: '11px', color: 'var(--fg-dim)' }, text: m.description })))))));

    box.appendChild(h('div', { class: 'panel' },
      h('h3', {}, 'How the two modes are set up'),
      h('div', { class: 'note' },
        h('strong', { text: 'Mode 1 — cell-type specificity. ' }),
        'One element is slid across positions relative to the TSS on a fixed background. ' +
        'At every position the whole cell-type panel is predicted and collapsed to one activity ' +
        'per cell type. Tau summarises how uneven the panel is: ',
        h('code', { text: 'τ = Σ(1 − xᵢ/x_max)/(n−1)' }),
        '. Tau alone cannot distinguish "the target went up" from "everything else went down", ' +
        'so each position is also decomposed into those two arms and labelled with which one moved.'),
      h('div', { class: 'note' },
        h('strong', { text: 'Mode 2 — motif interactions. ' }),
        'Four constructs are built on one background: background, A only, B only, A+B. They are ' +
        'the same length because placing an element overwrites bases. The residual ',
        h('code', { text: 'I = f(A+B) − f(A) − f(B) + f(bg)' }),
        ' is then a difference-in-differences that cancels the background and each element\'s own effect.'),
      h('div', { class: 'note warn' },
        h('strong', { text: 'The scale is part of the claim. ' }),
        'Both mock profile adapters emit a log-like quantity. A zero residual there means the two ' +
        'effects combine multiplicatively in count units — not that there is "no interaction". ' +
        'Switching "Residual on" to linear back-transforms first and tests additivity in counts instead. ' +
        'Any real adapter must declare its own output space for this to stay meaningful.'),
      h('div', { class: 'note warn' },
        h('strong', { text: 'Controls that are built in. ' }),
        h('ul', {},
          h('li', {}, h('code', { text: 'Puffin (mock)' }), ' is additive by construction, so its residual is identically zero. If it ever is not, the two sites are close enough to share a PWM scan window — or something is wrong.'),
          h('li', {}, h('code', { text: 'PuffinD (mock, pair terms off)' }), ' keeps only the saturating nonlinearity. Comparing against it separates "the response saturates" from "the model has a real pair term".'),
          h('li', { text: 'Overlapping sites are detected and flagged; the residual is not an interaction measurement there.' }),
          h('li', { text: 'Replicate backgrounds are averaged and the spread is drawn as a ribbon.' })))));

    box.appendChild(h('div', { class: 'panel' },
      h('h3', {}, 'Wiring up real inference'),
      h('p', { style: { fontSize: '12px', color: 'var(--fg-dim)' } },
        'Adapters are the only place that knows about a model. Implement ',
        h('code', { text: 'backend/adapters/real_puffin.py' }), ' or ',
        h('code', { text: 'backend/adapters/real_alphagenome.py' }),
        ', register the class in ', h('code', { text: 'adapters/__init__.py' }),
        ', and it appears in the pickers with a "real" badge. Each file carries a checklist of ' +
        'the things that have to be right — input length, output resolution, output scale, ' +
        'strandedness and determinism — because getting any of them wrong changes what the ' +
        'residual and tau mean without changing how the plots look.')));
  }

  return { boot, showView };
})();

window.addEventListener('DOMContentLoaded', App.boot);
