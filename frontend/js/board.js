/* Challenges, submissions and the leaderboard. */

const Board = (() => {

  async function refreshChallenges() {
    try {
      const { challenges } = await API.get('/api/challenges');
      const box = $('#challengeList');
      clear(box);
      if (!challenges.length) {
        box.appendChild(h('div', { class: 'hint', text: 'none yet' }));
        return;
      }
      for (const c of challenges) {
        const active = State.activeChallenge && State.activeChallenge.id === c.id;
        box.appendChild(h('div', {
          class: 'el' + (active ? ' sel' : ''),
          style: { cursor: 'pointer', gridTemplateColumns: '1fr auto' },
          onclick: () => open(c.id),
        },
          h('div', {},
            h('div', { class: 'name', text: c.name }),
            h('div', { class: 'meta', text:
              `${c.mode} · ${c.n_submissions} submission${c.n_submissions === 1 ? '' : 's'} · ` +
              `seed ${c.background_params.seed ?? '–'} · ${c.background_params.length ?? '–'} bp` })),
          h('div', { class: 'acts' },
            h('button', {
              class: 'btn icon sm', text: '×', title: 'delete challenge',
              onclick: async (e) => {
                e.stopPropagation();
                if (!confirm(`Delete challenge "${c.name}" and all of its submissions?`)) return;
                try {
                  await API.post('/api/challenges/delete', { challenge_id: c.id });
                  if (State.activeChallenge && State.activeChallenge.id === c.id) {
                    State.activeChallenge = null;
                    clear($('#boardContent'));
                    $('#boardContent').appendChild(h('div', { class: 'empty', text: 'Create or select a challenge.' }));
                  }
                  refreshChallenges();
                } catch (err) { fail(err); }
              },
            }))));
      }
    } catch (e) { fail(e); }
  }

  async function create() {
    const name = prompt('Name this challenge:',
      `${State.design.targetCellType || 'Focused promoter'} — ${State.design.background.length} bp`);
    if (name === null) return;
    try {
      const ch = await API.post('/api/challenges', {
        name,
        mode: State.design.targetCellType ? 'mode1' : 'design',
        background_params: {
          length: +State.design.background.length,
          gc: +State.design.background.gc,
          cpg_oe: +State.design.background.cpg_oe,
          seed: +State.design.background.seed,
        },
        target_cell_type: State.design.targetCellType,
        weights: State.design.weights,
        models: {
          profile_models: State.design.profileModels,
          celltype_model: State.design.celltypeModel,
        },
        player_id: State.player && State.player.id,
        notes: '',
      });
      toast('challenge created');
      await refreshChallenges();
      open(ch.id);
    } catch (e) { fail(e); }
  }

  async function open(cid) {
    try {
      const res = await API.get('/api/challenges/leaderboard?challenge_id=' + encodeURIComponent(cid));
      State.activeChallenge = res.challenge;
      render(res);
      refreshChallenges();
    } catch (e) { fail(e); }
  }

  function render(res) {
    const box = $('#boardContent');
    clear(box);
    const c = res.challenge;
    if (!c) { box.appendChild(h('div', { class: 'empty', text: 'challenge not found' })); return; }

    const bp = c.background_params;
    box.appendChild(h('div', { class: 'panel' },
      h('div', { class: 'plot-title', text: c.name }),
      h('dl', { class: 'kv' },
        h('dt', { text: 'shared background' }),
        h('dd', { text: `${bp.length} bp · GC ${bp.gc} · CpG o/e ${bp.cpg_oe} · seed ${bp.seed}` }),
        h('dt', { text: 'target cell type' }), h('dd', { text: c.target_cell_type || '—' }),
        h('dt', { text: 'weights' }),
        h('dd', { text: Object.entries(c.weights || {}).map(([k, v]) => `${k}=${v}`).join('  ') || '—' })),
      h('div', { class: 'hint', text:
        'The background is regenerated from this seed on the server, so every player starts from the ' +
        'identical sequence. Submissions are re-scored server-side from these settings — a client ' +
        'cannot submit its own numbers.' }),
      h('div', { class: 'btnrow', style: { marginTop: '9px' } },
        h('button', { class: 'btn', text: 'Load background into Design', onclick: () => loadInto(c) }),
        h('button', { class: 'btn primary', text: 'Submit current design', onclick: () => submit(c) }))));

    const rows = res.leaderboard.map(e => h('tr', { class: e.rank === 1 ? 'hl' : '' },
      h('td', { class: 'num', text: e.rank }),
      h('td', { text: e.player_name }),
      h('td', { text: e.label || '—' }),
      h('td', { class: 'num', text: fmtNum(e.score, 3) }),
      h('td', { class: 'num', text: fmtNum(e.metrics.peak_height, 2) }),
      h('td', { class: 'num', text: e.metrics.peak_position ?? '–' }),
      h('td', { class: 'num', text: fmtNum(e.metrics.tau, 3) }),
      h('td', { text: e.metrics.strongest_cell_type || '–',
                style: { color: e.metrics.target_cell_type &&
                  e.metrics.strongest_cell_type === e.metrics.target_cell_type
                  ? 'var(--accent)' : 'var(--fg-dim)' } }),
      h('td', { class: 'num', text: e.n_elements }),
      h('td', { text: new Date(e.created * 1000).toLocaleString() })));

    box.appendChild(h('div', { class: 'panel' },
      h('h3', {}, 'Leaderboard'),
      rows.length
        ? h('table', { class: 'data' },
            h('thead', {}, h('tr', {},
              h('th', { class: 'num', text: '#' }), h('th', { text: 'player' }),
              h('th', { text: 'label' }), h('th', { class: 'num', text: 'score' }),
              h('th', { class: 'num', text: 'peak' }), h('th', { class: 'num', text: 'pos' }),
              h('th', { class: 'num', text: 'tau' }), h('th', { text: 'strongest' }),
              h('th', { class: 'num', text: 'elements' }),
              h('th', { text: 'when' }))),
            h('tbody', {}, ...rows))
        : h('div', { class: 'empty', text: 'no submissions yet — be first' }),
      h('div', { class: 'hint', text: 'Best entry per player is shown.' })));
  }

  function loadInto(c) {
    const bp = c.background_params;
    Object.assign(State.design.background, {
      length: bp.length, gc: bp.gc, cpg_oe: bp.cpg_oe, seed: bp.seed,
    });
    $('#bgLength').value = bp.length;
    $('#bgGC').value = bp.gc; $('#bgGCv').textContent = Number(bp.gc).toFixed(2);
    $('#bgCpG').value = bp.cpg_oe; $('#bgCpGv').textContent = Number(bp.cpg_oe).toFixed(2);
    $('#bgSeed').value = bp.seed;
    if (c.weights && Object.keys(c.weights).length) {
      State.design.weights = Object.assign({}, c.weights);
      $('#weightPreset').value = 'custom';
      Designer.renderWeights();
    }
    if (c.target_cell_type) {
      State.design.targetCellType = c.target_cell_type;
      $('#targetCT').value = c.target_cell_type;
    }
    Designer.clampView(); Designer.syncZoomUI();
    App.showView('design');
    Designer.changed();
    toast('loaded the challenge background — design away');
  }

  async function submit(c) {
    if (!State.design.placements.length && !confirm('Submit an empty design?')) return;
    const name = ($('#playerName').value || '').trim();
    if (!name) { toast('enter a player name first', true); $('#playerName').focus(); return; }
    const label = prompt('Label for this submission (optional):', '');
    if (label === null) return;
    try {
      const res = await API.post('/api/challenges/submit', {
        challenge_id: c.id,
        player_name: name,
        player_id: State.player && State.player.id,
        label,
        placements: State.designPayload().placements,
      });
      State.player = State.player || { id: res.submission.player_id, name };
      toast(`submitted — score ${fmtNum(res.submission.score.score, 3)}`);
      open(c.id);
    } catch (e) { fail(e); }
  }

  function init() {
    $('#btnNewChallenge').onclick = create;
    refreshChallenges();
  }

  return { init, refreshChallenges, open };
})();
