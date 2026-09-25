/* Shared application state. */

const State = {
  motifs: [],            // Puffin's 10 core promoter motifs
  motifById: {},
  celltypeElements: [],  // lineage TF sites (celltype_elements.py)
  extraElements: [],     // cpg_segment, custom
  models: null,          // /api/models
  cellTypes: [],         // [{id,label,...}]
  metricCatalogue: [],
  defaultWeights: {},
  provenance: '',

  player: null,

  // --- current design ---
  design: {
    background: { length: 2001, gc: 0.45, cpg_oe: 0.25, seed: 42, scrub: true },
    placements: [],      // [{uid, element_id, position, strand, seg_length, custom_sequence}]
    selected: null,      // uid
    weights: {},
    targetCellType: null,
    profileModels: [],
    celltypeModel: null,
    view: { center: 0, halfWidth: 500 },
  },

  lastEvaluation: null,
  lastMode1: null,
  lastMode2: { pair: null, grid: null, spacing: null },
  activeChallenge: null,

  uidCounter: 0,
  newUid() { return 'e' + (++this.uidCounter); },

  elementMeta(id) {
    return this.motifById[id] ||
      this.celltypeElements.find(e => e.id === id) ||
      this.extraElements.find(e => e.id === id) ||
      { id, name: id, color: '#7f8c8d' };
  },

  elementWidth(p) {
    if (p.element_id === 'custom') {
      return Math.max(1, (p.custom_sequence || '').replace(/[^ACGTNacgtn]/g, '').length);
    }
    if (p.element_id === 'cpg_segment') return Math.max(1, p.seg_length || 60);
    const m = this.elementMeta(p.element_id);
    return m && m.width ? m.width : 10;
  },

  /* Payload shared by every backend call that needs the current design. */
  designPayload(extra) {
    const d = this.design;
    return Object.assign({
      length: +d.background.length,
      gc: +d.background.gc,
      cpg_oe: +d.background.cpg_oe,
      seed: +d.background.seed,
      scrub_celltype_elements: !!d.background.scrub,
      placements: d.placements.map(p => ({
        element_id: p.element_id, position: Math.round(p.position),
        strand: p.strand, uid: p.uid,
        seg_length: p.seg_length, custom_sequence: p.custom_sequence,
        instance: p.instance || 'consensus',
      })),
      weights: d.weights,
      target_cell_type: d.targetCellType,
      profile_models: d.profileModels,
      celltype_model: d.celltypeModel,
    }, extra || {});
  },

  addPlacement(elementId, position, strand) {
    const meta = this.elementMeta(elementId);
    const p = {
      uid: this.newUid(),
      element_id: elementId,
      position: Math.round(position ?? (meta.typical_offset ?? -50)),
      strand: strand || '+',
      seg_length: 60,
      custom_sequence: '',
      instance: 'consensus',
    };
    this.design.placements.push(p);
    this.design.selected = p.uid;
    return p;
  },

  removePlacement(uid) {
    this.design.placements = this.design.placements.filter(p => p.uid !== uid);
    if (this.design.selected === uid) this.design.selected = null;
  },

  getPlacement(uid) { return this.design.placements.find(p => p.uid === uid); },

  duplicatePlacement(uid) {
    const p = this.getPlacement(uid);
    if (!p) return null;
    const copy = Object.assign({}, p, {
      uid: this.newUid(),
      position: p.position + this.elementWidth(p) + 10,
    });
    this.design.placements.push(copy);
    this.design.selected = copy.uid;
    return copy;
  },
};

/* Subscribers re-render when the design changes. */
const Bus = {
  handlers: {},
  on(evt, fn) { (this.handlers[evt] = this.handlers[evt] || []).push(fn); },
  emit(evt, data) { (this.handlers[evt] || []).forEach(fn => { try { fn(data); } catch (e) { console.error(e); } }); },
};
