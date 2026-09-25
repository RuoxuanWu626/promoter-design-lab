/* API client, small DOM helpers, and export utilities. */

const API = (() => {
  async function req(path, body) {
    const opts = body === undefined
      ? { method: 'GET' }
      : { method: 'POST', headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify(body) };
    const res = await fetch(path, opts);
    let data;
    try { data = await res.json(); }
    catch (e) { throw new Error(`${res.status} ${res.statusText} (bad JSON)`); }
    if (!res.ok || data.error) {
      const err = new Error(data.error || `${res.status} ${res.statusText}`);
      err.detail = data.detail;
      throw err;
    }
    return data;
  }

  const get = (p) => req(p);
  const post = (p, b) => req(p, b || {});

  /* Start a background job and poll it. onProgress({done,total,stage,detail}) */
  async function job(path, body, onProgress) {
    const started = await post(path, Object.assign({}, body, { async: true }));
    if (!started.job_id) return started;           // server ran it inline
    const id = started.job_id;
    let delay = 150;
    for (;;) {
      await new Promise(r => setTimeout(r, delay));
      delay = Math.min(700, delay * 1.25);
      const st = await get(`/api/jobs/${id}`);
      if (onProgress && st.progress) onProgress(st.progress);
      if (st.status === 'done') return st.result;
      if (st.status === 'error') {
        const e = new Error(st.error); e.detail = st.traceback; throw e;
      }
    }
  }

  return { get, post, job };
})();


/* ---------------- DOM helpers ---------------- */
function h(tag, attrs, ...kids) {
  const e = document.createElement(tag);
  for (const [k, v] of Object.entries(attrs || {})) {
    if (v === null || v === undefined || v === false) continue;
    if (k === 'class') e.className = v;
    else if (k === 'html') e.innerHTML = v;
    else if (k === 'text') e.textContent = v;
    else if (k === 'style' && typeof v === 'object') Object.assign(e.style, v);
    else if (k.startsWith('on') && typeof v === 'function') e.addEventListener(k.slice(2), v);
    else e.setAttribute(k, v === true ? '' : v);
  }
  for (const k of kids.flat()) {
    if (k === null || k === undefined || k === false) continue;
    e.appendChild(typeof k === 'string' || typeof k === 'number'
      ? document.createTextNode(String(k)) : k);
  }
  return e;
}

const $ = (sel, root) => (root || document).querySelector(sel);
const $$ = (sel, root) => Array.from((root || document).querySelectorAll(sel));

function clear(node) { while (node && node.firstChild) node.removeChild(node.firstChild); }

let _toastTimer = null;
function toast(msg, isError) {
  let t = $('#toast');
  if (!t) { t = h('div', { id: 'toast', class: 'toast' }); document.body.appendChild(t); }
  t.className = 'toast on' + (isError ? ' err' : '');
  t.textContent = msg;
  clearTimeout(_toastTimer);
  _toastTimer = setTimeout(() => t.classList.remove('on'), isError ? 7000 : 3200);
}

function fail(err) {
  console.error(err);
  toast(err.message || String(err), true);
}

const fmtNum = (v, d = 3) =>
  (v === null || v === undefined || !isFinite(v)) ? '–' : Number(v).toFixed(d);


/* ---------------- export helpers ---------------- */
function download(filename, text, mime) {
  const blob = new Blob([text], { type: mime || 'text/plain;charset=utf-8' });
  const url = URL.createObjectURL(blob);
  const a = h('a', { href: url, download: filename });
  document.body.appendChild(a); a.click();
  setTimeout(() => { URL.revokeObjectURL(url); a.remove(); }, 0);
}

function toCSV(rows, columns) {
  const cols = columns || (rows.length ? Object.keys(rows[0]) : []);
  const esc = (v) => {
    if (v === null || v === undefined) return '';
    const s = String(v);
    return /[",\n]/.test(s) ? '"' + s.replace(/"/g, '""') + '"' : s;
  };
  return [cols.join(',')]
    .concat(rows.map(r => cols.map(c => esc(r[c])).join(',')))
    .join('\n') + '\n';
}

function stamp() {
  const d = new Date(), p = (n) => String(n).padStart(2, '0');
  return `${d.getFullYear()}${p(d.getMonth() + 1)}${p(d.getDate())}-${p(d.getHours())}${p(d.getMinutes())}${p(d.getSeconds())}`;
}
