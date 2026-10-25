// Stage 2 equality test: the browser generator (index.html, <script id="diffuser-core">) against the stage 1 Python canon.
//
//   npm test            (= node --test scripts/, run from the generator folder)
//
// Reference side: scripts/fixtures/stage1/ (STL + JSON of the 29 stage 1 tiles, plus layouts.json) and the other files in scripts/fixtures/ (made by scripts/reference.py).
// Both the JS mesh and the reference STL are measured by the SAME code in this file (STL parser, Moller-Trumbore ray
// cast, signed volume, shell finder), so the comparison cannot hide a bug in either side's own tooling. The first
// group of tests checks that tooling against the Python sidecars and needs no index.html.
// Contract v2 (derive margin 25 mm, wallOptions / wallBand / features / ESTIMATE / estimate) is covered in sections 3 and 5.
// Written from CONTRACT.md and the Python canon only; the generator's JS code is not read.

import { describe, it } from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import path from 'node:path';
import vm from 'node:vm';
import { fileURLToPath } from 'node:url';
import Module from 'manifold-3d';

const HERE = path.dirname(fileURLToPath(import.meta.url));
const GEN = path.resolve(HERE, '..');
const STAGE1 = path.join(HERE, 'fixtures', 'stage1');
const FIX = path.join(HERE, 'fixtures');

// ---------------------------------------------------------------- tolerances (CONTRACT.md, "Equality")
const TOL_H = 0.01;        // mm, heights by ray cast
const TOL_BB = 0.01;       // mm, bounding box
const TOL_VOL = 0.001;     // relative, volume (0.1 %)
const TOL_NUM = 1e-9;      // relative, pure number ports (bands, phases, depths)
const EPS_VOL = 1e-3;      // mm^3, a shell below this volume counts as zero volume
const SIDECAR_VOL_CM3 = 0.01;   // self-check: sidecar volume is rounded to 0.01 cm3

// ---------------------------------------------------------------- loading
const readJson = (p) => JSON.parse(fs.readFileSync(p, 'utf8'));
const plain = (x) => JSON.parse(JSON.stringify(x));        // strips the vm realm: arrays/objects of another realm break deepStrictEqual

const needFile = (p, hint) => {
  if (!fs.existsSync(p)) throw new Error(`missing ${p} (${hint})`);
  return p;
};
const EXP = readJson(needFile(path.join(FIX, 'expected.json'), 'run: python scripts/reference.py'));
const CASES = readJson(needFile(path.join(FIX, 'cases.json'), 'fixture cases'));
const LAYOUTS_JSON = readJson(needFile(path.join(STAGE1, 'layouts.json'), 'stage 1 output'));

const M = await Module();
M.setup();

let core = null;
let coreError = null;
try {
  const indexPath = process.env.DIFFUSER_INDEX || path.join(GEN, 'index.html');   // DIFFUSER_INDEX: harness validation only
  if (!fs.existsSync(indexPath)) throw new Error(`index.html not found: ${indexPath}`);
  const html = fs.readFileSync(indexPath, 'utf8');
  const tag = /<script\b([^>]*)\bid\s*=\s*["']diffuser-core["']([^>]*)>([\s\S]*?)<\/script>/i.exec(html);
  if (!tag) throw new Error('index.html has no <script id="diffuser-core">');
  if (/type\s*=\s*["']?module/i.test(tag[1] + tag[2])) throw new Error('<script id="diffuser-core"> must be a classic script, not type="module"');
  const ctx = vm.createContext({ Math, console });
  // manifold-3d's JS wrapper tests vector arguments with `instanceof Array` of the MAIN realm; array literals created inside
  // a vm context belong to the vm realm and would be rejected ("Missing field: x"). In a browser there is one realm, so this
  // is an artefact of the vm loading only: make the vm Array.prototype inherit from the main one.
  Object.setPrototypeOf(vm.runInContext('Array', ctx).prototype, Array.prototype);
  const factory = vm.runInContext(tag[3] + '\n;(typeof createDiffuser === "function") ? createDiffuser : undefined',
    ctx, { filename: 'index.html#diffuser-core' });
  if (typeof factory !== 'function') throw new Error('the core script does not define a global function createDiffuser(M)');
  core = factory(M);
} catch (e) {
  coreError = e;
}
function needCore() {
  if (coreError) throw new Error(`generator core not available: ${coreError.message}`);
  return core;
}

// ---------------------------------------------------------------- mesh tooling (used for BOTH sides)
function weld(pos, tris) {
  // pos: float32 xyz triples (flat), tris: vertex index triples; vertices with identical float32 coordinates are merged
  const f32 = Float32Array.from(pos);
  const map = new Map();
  const remap = new Int32Array(f32.length / 3);
  const out = [];
  for (let i = 0; i < remap.length; i++) {
    const k = `${f32[3 * i]}|${f32[3 * i + 1]}|${f32[3 * i + 2]}`;
    let id = map.get(k);
    if (id === undefined) { id = out.length / 3; map.set(k, id); out.push(f32[3 * i], f32[3 * i + 1], f32[3 * i + 2]); }
    remap[i] = id;
  }
  const t = new Int32Array(tris.length);
  for (let i = 0; i < tris.length; i++) t[i] = remap[tris[i]];
  return { v: Float64Array.from(out), t };
}

function parseStl(buf) {
  if (buf.length < 84) throw new Error('STL too short');
  const n = buf.readUInt32LE(80);
  if (84 + 50 * n !== buf.length) throw new Error(`not a binary STL: ${n} triangles announced, size ${buf.length}`);
  const pos = new Float32Array(9 * n);
  const tris = new Int32Array(3 * n);
  for (let i = 0; i < n; i++) {
    const o = 84 + 50 * i + 12;                      // skip the normal
    for (let j = 0; j < 9; j++) pos[9 * i + j] = buf.readFloatLE(o + 4 * j);
    tris[3 * i] = 3 * i; tris[3 * i + 1] = 3 * i + 1; tris[3 * i + 2] = 3 * i + 2;
  }
  return weld(pos, tris);
}

function meshFromManifold(man) {
  const g = man.getMesh();
  const np = g.numProp;
  const nv = g.vertProperties.length / np;
  const pos = new Float32Array(3 * nv);
  for (let i = 0; i < nv; i++) for (let j = 0; j < 3; j++) pos[3 * i + j] = g.vertProperties[np * i + j];
  return weld(pos, g.triVerts);
}

function bbox(mesh) {
  const min = [Infinity, Infinity, Infinity], max = [-Infinity, -Infinity, -Infinity];
  for (let i = 0; i < mesh.v.length; i += 3) for (let a = 0; a < 3; a++) {
    min[a] = Math.min(min[a], mesh.v[i + a]); max[a] = Math.max(max[a], mesh.v[i + a]);
  }
  return { min, max };
}

function triVolume(v, a, b, c) {
  const ax = v[3 * a], ay = v[3 * a + 1], az = v[3 * a + 2];
  const bx = v[3 * b], by = v[3 * b + 1], bz = v[3 * b + 2];
  const cx = v[3 * c], cy = v[3 * c + 1], cz = v[3 * c + 2];
  return (ax * (by * cz - bz * cy) - ay * (bx * cz - bz * cx) + az * (bx * cy - by * cx)) / 6;
}

/** Shells = faces connected through shared edges. Volume > 0: outer solid, < 0: inner cavity (chamber), ~0: sheet. */
function analyze(mesh) {
  const { v, t } = mesh;
  const nt = t.length / 3, nv = v.length / 3;
  const parent = Int32Array.from({ length: nt }, (_, i) => i);
  const find = (x) => { while (parent[x] !== x) { parent[x] = parent[parent[x]]; x = parent[x]; } return x; };
  const union = (a, b) => { a = find(a); b = find(b); if (a !== b) parent[a] = b; };
  const undirected = new Map();                      // key -> first face seen
  const directed = new Map();                        // a*nv+b -> count
  const edgeFaces = new Map();                       // undirected key -> face count
  let degenerate = 0;
  for (let f = 0; f < nt; f++) {
    const idx = [t[3 * f], t[3 * f + 1], t[3 * f + 2]];
    if (idx[0] === idx[1] || idx[1] === idx[2] || idx[0] === idx[2]) degenerate++;
    for (let e = 0; e < 3; e++) {
      const a = idx[e], b = idx[(e + 1) % 3];
      const dk = a * nv + b;
      directed.set(dk, (directed.get(dk) || 0) + 1);
      const uk = a < b ? a * nv + b : b * nv + a;
      edgeFaces.set(uk, (edgeFaces.get(uk) || 0) + 1);
      const first = undirected.get(uk);
      if (first === undefined) undirected.set(uk, f); else union(first, f);
    }
  }
  const shellOf = new Map();
  const shells = [];
  for (let f = 0; f < nt; f++) {
    const r = find(f);
    let s = shellOf.get(r);
    if (!s) { s = { faces: 0, volume: 0, badEdges: 0 }; shellOf.set(r, s); shells.push(s); }
    s.faces++;
    s.volume += triVolume(v, t[3 * f], t[3 * f + 1], t[3 * f + 2]);
  }
  // closed + consistently oriented: every directed edge a->b appears once and b->a appears once
  let nonManifold = 0;
  for (let f = 0; f < nt; f++) {
    const s = shellOf.get(find(f));
    for (let e = 0; e < 3; e++) {
      const a = t[3 * f + e], b = t[3 * f + (e + 1) % 3];
      if (directed.get(a * nv + b) !== 1 || directed.get(b * nv + a) !== 1) { s.badEdges++; }
      const uk = a < b ? a * nv + b : b * nv + a;
      if (edgeFaces.get(uk) !== 2) nonManifold++;
    }
  }
  const outer = shells.filter((s) => s.volume > EPS_VOL);
  const inner = shells.filter((s) => s.volume < -EPS_VOL);
  const zero = shells.filter((s) => Math.abs(s.volume) <= EPS_VOL);
  return {
    triangles: nt, shells: shells.length, outer: outer.length, inner: inner.length, zero: zero.length,
    innerVolumes: inner.map((s) => -s.volume).sort((x, y) => y - x),
    open: shells.filter((s) => s.badEdges > 0).length, nonManifoldEdgeUses: nonManifold, degenerate,
    volume: shells.reduce((a, s) => a + s.volume, 0),
  };
}

/** First hit of the ray o + t d (t > 0) with the mesh; Moller-Trumbore, edges inclusive. Returns the hit point or null. */
function firstHit(mesh, o, d) {
  const { v, t } = mesh;
  const EPS = 1e-7;
  let best = Infinity;
  for (let f = 0; f < t.length; f += 3) {
    const a = 3 * t[f], b = 3 * t[f + 1], c = 3 * t[f + 2];
    const e1x = v[b] - v[a], e1y = v[b + 1] - v[a + 1], e1z = v[b + 2] - v[a + 2];
    const e2x = v[c] - v[a], e2y = v[c + 1] - v[a + 1], e2z = v[c + 2] - v[a + 2];
    const px = d[1] * e2z - d[2] * e2y, py = d[2] * e2x - d[0] * e2z, pz = d[0] * e2y - d[1] * e2x;
    const det = e1x * px + e1y * py + e1z * pz;
    if (Math.abs(det) < 1e-12) continue;
    const inv = 1 / det;
    const sx = o[0] - v[a], sy = o[1] - v[a + 1], sz = o[2] - v[a + 2];
    const u = (sx * px + sy * py + sz * pz) * inv;
    if (u < -EPS || u > 1 + EPS) continue;
    const qx = sy * e1z - sz * e1y, qy = sz * e1x - sx * e1z, qz = sx * e1y - sy * e1x;
    const w = (d[0] * qx + d[1] * qy + d[2] * qz) * inv;
    if (w < -EPS || u + w > 1 + EPS) continue;
    const tt = (e2x * qx + e2y * qy + e2z * qz) * inv;
    if (tt > 1e-9 && tt < best) best = tt;
  }
  return best === Infinity ? null : [o[0] + best * d[0], o[1] + best * d[1], o[2] + best * d[2]];
}

// ---------------------------------------------------------------- measurement points (same points as the Python reference tests)
/** Height rays: deep - per well n at fractions 0.2/0.5/0.8 of the length, looking down -Y from above the front plane;
 *  blocks - per cell centre, looking down -Z. Returns tags + origins + direction + the axis of the height. */
function heightPoints(series, p) {
  const origins = [], tags = [];
  if (series === 'deep') {
    const T = p.dmax + p.back;
    for (let n = 0; n < p.N; n++) for (const f of [0.2, 0.5, 0.8]) {
      origins.push([(n + 0.5) * p.pitch, T + 10, f * p.length]); tags.push(`well ${n} @${f}`);
    }
    return { origins, tags, dir: [0, -1, 0], axis: 1 };
  }
  const top = p.base + p.hmax + 10;
  for (let n = 0; n < p.Nx; n++) for (let m = 0; m < p.Ny; m++) {
    origins.push([(n + 0.5) * p.W / p.Nx, (m + 0.5) * p.H / p.Ny, top]); tags.push(`cell n=${n} m=${m}`);
  }
  return { origins, tags, dir: [0, 0, -1], axis: 2, nm: tags.map((_, i) => [Math.floor(i / p.Ny), i % p.Ny]) };
}

function heights(mesh, series, p) {
  const pts = heightPoints(series, p);
  const hs = pts.origins.map((o) => { const h = firstHit(mesh, o, pts.dir); return h ? h[pts.axis] : NaN; });
  return { hs, tags: pts.tags, pts };
}

/** Label: grid of rays from outside the back (deep) / the bed face (blocks) through the sunk digit + arrow.
 *  The first hit is at 0 outside the label and at 0.4 inside it, so this measures the label shape (volume cannot: 0.03 %). */
function labelHits(mesh, series, p) {
  const out = [];
  const step = 0.7, off = 0.13;
  let cx, cy;
  if (series === 'deep') { cx = p.N * p.pitch / 2; cy = p.length / 2; } else { cx = p.W / 2; cy = p.H / 2; }
  for (let x = cx - 9 + off; x < cx + 9; x += step) for (let y = cy - 20 + off; y < cy + 20; y += step) {
    const o = series === 'deep' ? [x, -10, y] : [x, y, -10];
    const d = series === 'deep' ? [0, 1, 0] : [0, 0, 1];
    const h = firstHit(mesh, o, d);
    out.push(h ? h[series === 'deep' ? 1 : 2] : NaN);
  }
  return out;
}

const count = (a, pred) => a.reduce((n, x) => n + (pred(x) ? 1 : 0), 0);
const maxAbsDiff = (a, b) => {
  let worst = 0, at = -1;
  for (let i = 0; i < a.length; i++) {
    const d = (Number.isNaN(a[i]) || Number.isNaN(b[i])) ? (Number.isNaN(a[i]) && Number.isNaN(b[i]) ? 0 : Infinity) : Math.abs(a[i] - b[i]);
    if (d > worst) { worst = d; at = i; }
  }
  return { worst, at };
};
const near = (a, b, rel = TOL_NUM) => Math.abs(a - b) <= rel * Math.max(1, Math.abs(a), Math.abs(b));

// ---------------------------------------------------------------- the 29 stage 1 tiles and the fixture cases
const TILES = [];
for (const size of ['S', 'L']) for (let t = 1; t <= 7; t++) {
  TILES.push({ id: EXP.names.deep[size][t - 1], series: 'deep', size, tile: t, params: EXP.presets.deep[size] });
}
for (const size of ['S', 'L']) LAYOUTS_JSON.blocks[size].shifts.forEach((kl, i) => {
  TILES.push({ id: EXP.names.blocks[size][i], series: 'blocks', size, tile: i + 1, params: EXP.presets.blocks[size], shift: kl });
});
for (const t of TILES) { t.stl = path.join(STAGE1, `${t.id}.stl`); t.sidecar = path.join(STAGE1, `${t.id}.json`); }
const CASE_ITEMS = CASES.map((c) => ({
  id: c.id, series: c.series, size: null, tile: c.tile, params: c.params, shift: c.shift,
  stl: path.join(FIX, EXP.cases[c.id].file), sidecar: path.join(FIX, `${c.id}.json`), note: c.note,
}));

const refCache = new Map();
function reference(item) {
  if (!refCache.has(item.id)) {
    const mesh = parseStl(fs.readFileSync(needFile(item.stl, 'reference STL')));
    refCache.set(item.id, { mesh, an: analyze(mesh), bb: bbox(mesh) });
  }
  return refCache.get(item.id);
}

// ---------------------------------------------------------------- reporting helpers
class Report {
  constructor(label) { this.label = label; this.fails = []; }
  check(measure, ok, detail) { if (!ok) this.fails.push(`${measure}: ${typeof detail === 'function' ? detail() : detail}`); }
  done() { if (this.fails.length) assert.fail(`${this.label}\n  - ${this.fails.join('\n  - ')}`); }
}
const STATS = {};
function rec(measure, value, label) {
  if (!STATS[measure] || value > STATS[measure].value) STATS[measure] = { value, label };
}
const f4 = (x) => (Number.isFinite(x) ? x.toFixed(4) : String(x));

// ================================================================ 1. harness self-check (no index.html needed)
describe('harness self-check: STL parser, ray cast, volume, shells reproduce the Python sidecars', () => {
  for (const item of [...TILES, ...CASE_ITEMS]) {
    it(`${item.id}`, () => {
      const r = reference(item);
      const side = readJson(needFile(item.sidecar, 'Python sidecar'));
      const rep = new Report(item.id);

      // heights vs sidecar (depths / cells.top_z)
      const { hs, pts } = heights(r.mesh, item.series, item.params);
      if (item.series === 'deep') {
        const T = side.front_plane_y;
        let worst = 0;
        hs.forEach((h, i) => { const n = Math.floor(i / 3); worst = Math.max(worst, Math.abs(h - (T - side.depths[n]))); });
        rec('self-check heights vs sidecar (mm)', worst, item.id);
        rep.check('heights vs sidecar depths', worst <= TOL_H, `worst ${f4(worst)} mm`);
      } else {
        const top = new Map(side.cells.map((c) => [`${c.n},${c.m}`, c.top_z]));
        let worst = 0;
        pts.nm.forEach(([n, m], i) => { worst = Math.max(worst, Math.abs(hs[i] - top.get(`${n},${m}`))); });
        rec('self-check heights vs sidecar (mm)', worst, item.id);
        rep.check('heights vs sidecar top_z', worst <= TOL_H, `worst ${f4(worst)} mm`);
        rep.check('sidecar cells', side.cells.length === item.params.Nx * item.params.Ny, `${side.cells.length} cells`);
      }

      // volume
      const volCm3 = r.an.volume / 1000;
      const dv = Math.abs(volCm3 - side.volume_cm3);
      rec('self-check volume vs sidecar (cm3)', dv, item.id);
      rep.check('volume vs volume_cm3', dv <= SIDECAR_VOL_CM3, `STL ${volCm3.toFixed(4)} cm3, sidecar ${side.volume_cm3}`);

      // bounding box
      const ext = side.stl_extent_xyz;
      const bb = r.bb;
      let dbb = 0;
      for (let a = 0; a < 3; a++) dbb = Math.max(dbb, Math.abs(bb.min[a]), Math.abs(bb.max[a] - ext[a]));
      rec('self-check bbox vs sidecar (mm)', dbb, item.id);
      rep.check('bounding box', dbb <= TOL_BB, `min ${bb.min.map(f4)} max ${bb.max.map(f4)} want 0..${ext}`);

      // shells
      rep.check('triangle count vs sidecar', r.an.triangles === side.triangles, `${r.an.triangles} vs ${side.triangles}`);
      rep.check('one closed outer solid', r.an.outer === 1 && r.an.open === 0 && r.an.zero === 0 && r.an.degenerate === 0,
        `outer ${r.an.outer}, open ${r.an.open}, zero-volume ${r.an.zero}, degenerate ${r.an.degenerate}`);
      if (item.series === 'deep') {
        const p = item.params, T = side.front_plane_y;
        const chambers = count(side.depths, (d) => d < p.dmax && (T - d - p.floor) - p.back >= 4);
        rep.check('inner chamber shells', r.an.inner === chambers, `${r.an.inner} found, ${chambers} expected from sidecar depths`);
      } else {
        rep.check('blocks have no inner shells', r.an.inner === 0 && r.an.shells === 1, `${r.an.shells} shells, ${r.an.inner} inner`);
      }

      // label: both levels (0 outside, 0.4 inside) must be present, nothing else
      const lab = labelHits(r.mesh, item.series, item.params);
      const lo = count(lab, (h) => Math.abs(h) < 1e-4), hi = count(lab, (h) => Math.abs(h - 0.4) < 1e-4);
      rep.check('label rays', lo > 0 && hi > 0 && lo + hi === lab.length, `${lo} at 0, ${hi} at 0.4, ${lab.length} rays`);
      rep.done();
    });
  }

  it('synthetic: cube with a cavity (ray cast, volume, shell classification)', () => {
    const cube = (s, o, inward) => {
      const p = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]]
        .flatMap(([x, y, z]) => [o + s * x, o + s * y, o + s * z]);
      let f = [[0, 2, 1], [0, 3, 2], [4, 5, 6], [4, 6, 7], [0, 1, 5], [0, 5, 4], [1, 2, 6], [1, 6, 5], [2, 3, 7], [2, 7, 6], [3, 0, 4], [3, 4, 7]];
      if (inward) f = f.map(([a, b, c]) => [a, c, b]);
      return { p, f };
    };
    const outer = cube(10, 0, false), cav = cube(2, 4, true);
    const pos = [...outer.p, ...cav.p];
    const tris = [...outer.f.flat(), ...cav.f.flat().map((i) => i + 8)];
    const mesh = weld(pos, tris);
    const an = analyze(mesh);
    assert.equal(an.shells, 2); assert.equal(an.outer, 1); assert.equal(an.inner, 1); assert.equal(an.open, 0);
    assert.ok(Math.abs(an.volume - (1000 - 8)) < 1e-9, `volume ${an.volume}`);
    assert.ok(Math.abs(an.innerVolumes[0] - 8) < 1e-9);
    const hit = firstHit(mesh, [5, 5, 20], [0, 0, -1]);
    assert.ok(hit && Math.abs(hit[2] - 10) < 1e-9, `first hit ${hit}`);                 // top of the outer cube, not the cavity
    const edge = firstHit(mesh, [5, 5, -20], [0, 0, 1]);
    assert.ok(edge && Math.abs(edge[2]) < 1e-9);
    assert.equal(firstHit(mesh, [50, 50, 20], [0, 0, -1]), null);
    const open = analyze(weld(outer.p, outer.f.flat().slice(3)));                          // one face missing -> not closed
    assert.ok(open.open > 0);
  });

  it('synthetic: STL parser rejects a truncated file', () => {
    assert.throws(() => parseStl(Buffer.alloc(84 + 49)), /not a binary STL/);
  });
});

// ================================================================ 2. numeric ports and tables (need index.html)
describe('core: members, tables and pure number ports', () => {
  it('createDiffuser(M) returns every member listed in CONTRACT.md', () => {
    const c = needCore();
    const want = ['PRESETS', 'LAYOUTS', 'residues', 'shifted', 'deepShift', 'deepDepths', 'phaseFraction', 'phaseMax', 'blockHeights',
      'deepBand', 'blockBand', 'deepDims', 'blockDims', 'deepTile', 'blockTile', 'name', 'derive', 'limits'];
    const missing = want.filter((k) => c[k] === undefined);
    assert.deepEqual(missing, [], `missing members: ${missing.join(', ')}`);
  });

  it('PRESETS equal qrd.DEEP / qrd.BLOCKS field for field (bed string included)', () => {
    const c = needCore();
    assert.deepStrictEqual(plain(c.PRESETS), { deep: EXP.presets.deep, blocks: EXP.presets.blocks });
  });

  it('LAYOUTS equal the relevant parts of layouts.json (only room_left_to_right / shifts / grids)', () => {
    const c = needCore();
    const want = { deep: {}, blocks: {} };
    for (const s of ['S', 'L']) {
      want.deep[s] = { room_left_to_right: LAYOUTS_JSON.deep[s].room_left_to_right };
      want.blocks[s] = { shifts: LAYOUTS_JSON.blocks[s].shifts, grids: LAYOUTS_JSON.blocks[s].grids };
    }
    assert.deepStrictEqual(plain(c.LAYOUTS), want);
  });

  it('name() equals the stage 1 file stems (29 tiles, checked against the files in scripts/fixtures/stage1)', () => {
    const c = needCore();
    const rep = new Report('name()');
    for (const t of TILES) {
      const got = c.name(t.series, t.params, t.tile);
      rep.check(t.id, got === t.id, `got "${got}"`);
      rep.check(t.id, fs.existsSync(t.stl), 'no such STL in scripts/fixtures/stage1');
    }
    for (const t of CASE_ITEMS) {
      const got = c.name(t.series, t.params, t.tile);
      rep.check(t.id, got === EXP.cases[t.id].name, `got "${got}", Python qrd.name gives "${EXP.cases[t.id].name}"`);
    }
    rep.done();
  });

  it('residues, shifted, deepShift, phaseMax, phaseFraction equal qrd.py', () => {
    const c = needCore();
    const rep = new Report('sequences');
    for (const [n, want] of Object.entries(EXP.residues)) rep.check(`residues(${n})`, JSON.stringify(plain(c.residues(+n))) === JSON.stringify(want), `got ${JSON.stringify(c.residues(+n))}`);
    for (const [key, want] of Object.entries(EXP.shifted)) {
      const [n, k] = key.split(':').map(Number);
      rep.check(`shifted(${n},${k})`, JSON.stringify(plain(c.shifted(n, k))) === JSON.stringify(want), `got ${JSON.stringify(c.shifted(n, k))}`);
    }
    for (const [t, want] of Object.entries(EXP.deep_shift)) rep.check(`deepShift(${t})`, c.deepShift(+t) === want, `got ${c.deepShift(+t)}`);
    for (const [key, want] of Object.entries(EXP.phase_max)) {
      const [nx, ny] = key.split('x').map(Number);
      rep.check(`phaseMax(${nx},${ny})`, near(c.phaseMax(nx, ny), want), `got ${c.phaseMax(nx, ny)}, want ${want}`);
    }
    for (const q of EXP.phase_fraction) {
      const got = c.phaseFraction(q.Nx, q.Ny, q.n, q.m, q.k, q.l);
      rep.check(`phaseFraction(${q.Nx},${q.Ny},${q.n},${q.m},${q.k},${q.l})`, near(got, q.v), `got ${got}, want ${q.v}`);
    }
    rep.done();
  });

  it('deepDepths and blockHeights equal qrd.py (presets and fixture cases)', () => {
    const c = needCore();
    const rep = new Report('depths/heights');
    for (const size of ['S', 'L']) for (let t = 1; t <= 7; t++) {
      const got = plain(c.deepDepths(EXP.presets.deep[size], t)), want = EXP.deep_depths[size][t];
      rep.check(`deepDepths(${size},${t})`, got.length === want.length && got.every((x, i) => near(x, want[i])), `got ${JSON.stringify(got)}, want ${JSON.stringify(want)}`);
    }
    for (const size of ['S', 'L']) LAYOUTS_JSON.blocks[size].shifts.forEach(([k, l], i) => {
      const got = plain(c.blockHeights(EXP.presets.blocks[size], k, l)), want = EXP.block_heights[size][i];
      const ok = got.length === want.length && got.every((row, m) => row.length === want[m].length && row.every((x, n) => near(x, want[m][n])));
      rep.check(`blockHeights(${size}, ${k},${l})`, ok, 'differs ([m][n] layout, rows along H)');
    });
    for (const it2 of CASE_ITEMS) {
      const ex = EXP.cases[it2.id];
      if (it2.series === 'deep') {
        const got = plain(c.deepDepths(it2.params, it2.tile));
        rep.check(`case ${it2.id} deepDepths`, got.length === ex.depths.length && got.every((x, i) => near(x, ex.depths[i])), `got ${JSON.stringify(got)}`);
      } else {
        const got = plain(c.blockHeights(it2.params, it2.shift[0], it2.shift[1]));
        const ok = got.length === ex.heights.length && got.every((row, m) => row.every((x, n) => near(x, ex.heights[m][n])));
        rep.check(`case ${it2.id} blockHeights`, ok, 'differs');
      }
    }
    rep.done();
  });

  it('deepBand / blockBand / deepDims / blockDims equal the Python values (presets and fixture cases)', () => {
    const c = needCore();
    const rep = new Report('bands');
    const cmp = (label, got, want) => {
      got = plain(got);
      rep.check(`${label} keys`, JSON.stringify(Object.keys(got).sort()) === JSON.stringify(Object.keys(want).sort()),
        `got [${Object.keys(got).sort()}], want [${Object.keys(want).sort()}]`);
      for (const k of Object.keys(want)) rep.check(`${label}.${k}`, typeof got[k] === 'number' && near(got[k], want[k]), `got ${got[k]}, want ${want[k]}`);
    };
    const dims = (label, got, want) => rep.check(label, got.length === 3 && want.every((x, i) => near(got[i], x)), `got ${JSON.stringify(plain(got))}, want ${JSON.stringify(want)}`);
    for (const size of ['S', 'L']) {
      cmp(`deepBand(${size})`, c.deepBand(EXP.presets.deep[size]), EXP.presets_band.deep[size].band);
      dims(`deepDims(${size})`, c.deepDims(EXP.presets.deep[size]), EXP.presets_band.deep[size].dims);
      cmp(`blockBand(${size})`, c.blockBand(EXP.presets.blocks[size]), EXP.presets_band.blocks[size].band);
      dims(`blockDims(${size})`, c.blockDims(EXP.presets.blocks[size]), EXP.presets_band.blocks[size].dims);
    }
    for (const it2 of CASE_ITEMS) {
      const ex = EXP.cases[it2.id];
      if (it2.series === 'deep') { cmp(`case ${it2.id} deepBand`, c.deepBand(it2.params), ex.band); dims(`case ${it2.id} deepDims`, c.deepDims(it2.params), ex.dims); }
      else { cmp(`case ${it2.id} blockBand`, c.blockBand(it2.params), ex.band); dims(`case ${it2.id} blockDims`, c.blockDims(it2.params), ex.dims); }
    }
    rep.done();
  });
});

// ================================================================ 3. derive / limits (CONTRACT.md v2 rules)
const BEDS = {
  a1mini: { x: 180, y: 180, z: 180 },
  mk4: { x: 250, y: 210, z: 220 },
  p1s: { x: 256, y: 256, z: 256 },
  custom: { x: 200, y: 160, z: 150 },
  narrow: { x: 230, y: 155, z: 200 },     // H/W < 10/13: W is reduced (no .5 rounding case)
  slim: { x: 240, y: 190, z: 210 },       // 10/13 <= H/W < 12/13: Ny = 11
};
const MARGIN = 25;                         // v2: 25 mm on every bed axis
// independent implementation of the contract rules
function expectDerive(series, bed, depth) {
  if (series === 'deep') {
    const pitch = Math.floor((bed.x - MARGIN) / 7);
    const up01 = (v) => Math.ceil(v * 10 - 1e-9) / 10;      // rule v3: thickness grows with the depth
    return { N: 7, pitch, length: Math.min(bed.y - MARGIN, bed.z - MARGIN), dmax: depth, fin: Math.max(1.6, up01(depth / 82.5)), edge: Math.max(2.0, up01(depth / 66)), floor: 1.6, back: 2.0, cap: 1.2, chamfer: pitch < 30 ? 3 : 4 };
  }
  let W = bed.x - MARGIN;
  const H = Math.min(bed.y - MARGIN, W);
  const ratio = H / W;
  const Ny = ratio >= 12 / 13 ? 13 : 11;
  if (ratio < 10 / 13) W = Math.round(H * 13 / 10);
  return { Nx: 13, Ny, W, H, hmax: depth, base: 3.0 };
}
function expectLimits(series, bed) {
  if (series === 'deep') {
    const pitch = Math.floor((bed.x - MARGIN) / 7);
    const max = bed.y - MARGIN - 2.0;
    return { depthMin: 40, depthMax: max, depthDefault: Math.min(4 * pitch, max) };
  }
  const p = expectDerive('blocks', bed, 48);
  const max = Math.min(bed.z - MARGIN - 3.0, 4 * Math.min(p.W / p.Nx, p.H / p.Ny));
  return { depthMin: 20, depthMax: max, depthDefault: Math.min(p.Ny === 13 ? 48 : 60, Math.floor(max)) };
}
const withoutBed = (p) => { const q = { ...plain(p) }; delete q.bed; return q; };

describe('core: derive and limits (v2, margin 25 mm)', () => {
  it('"Results" lines of CONTRACT.md v2 on A1 mini, MK4 and P1S (hard-coded numbers)', () => {
    const c = needCore();
    const rep = new Report('v2 Results');
    const eq = (label, got, want) => rep.check(label, typeof got === 'number' && Math.abs(got - want) <= 0.005, `got ${got}, want ${want}`);
    const deep = [['A1 mini', BEDS.a1mini, 22, 154, 155, 153, 88], ['MK4', BEDS.mk4, 32, 224, 185, 183, 128], ['P1S', BEDS.p1s, 33, 231, 231, 229, 132]];
    for (const [bn, bed, pitch, w, len, max, def] of deep) {
      const lim = plain(c.limits('deep', bed));
      const p = plain(c.derive('deep', bed, def));
      eq(`deep ${bn} pitch`, p.pitch, pitch);
      eq(`deep ${bn} tile width`, c.deepDims(p)[0], w);
      eq(`deep ${bn} tile length`, c.deepDims(p)[1], len);
      eq(`deep ${bn} depthMax`, lim.depthMax, max);
      eq(`deep ${bn} depthDefault`, lim.depthDefault, def);
      eq(`deep ${bn} depthMin`, lim.depthMin, 40);
    }
    const blocks = [['A1 mini', BEDS.a1mini, 155, 155, 13, 47.69, 47], ['MK4', BEDS.mk4, 225, 185, 11, 67.27, 60], ['P1S', BEDS.p1s, 231, 231, 13, 71.08, 48]];
    for (const [bn, bed, w, h, ny, max, def] of blocks) {
      const lim = plain(c.limits('blocks', bed));
      const p = plain(c.derive('blocks', bed, def));
      eq(`blocks ${bn} W`, p.W, w); eq(`blocks ${bn} H`, p.H, h); eq(`blocks ${bn} Ny`, p.Ny, ny); eq(`blocks ${bn} Nx`, p.Nx, 13);
      eq(`blocks ${bn} depthMax`, lim.depthMax, max);
      eq(`blocks ${bn} depthDefault`, lim.depthDefault, def);
      eq(`blocks ${bn} depthMin`, lim.depthMin, 20);
    }
    rep.done();
  });

  it('rule v3: fin and edge of the deep series, hard-coded numbers from the Python canon (qrd.derive, rule 3)', () => {
    const c = needCore();
    const rep = new Report('v3 fin / edge');
    const eq = (label, got, want) => rep.check(label, typeof got === 'number' && Math.abs(got - want) <= 1e-9, `got ${got}, want ${want}`);
    // [depth, fin, edge]: defaults of the three beds, 132 = exact multiple of 66 (edge stays 2.0), the first depth above it,
    // the maximum depth of each bed, and 330 = exact multiple of both divisors (4.0 / 5.0)
    const table = [[40, 1.6, 2.0], [88, 1.6, 2.0], [128, 1.6, 2.0], [132, 1.6, 2.0], [132.5, 1.7, 2.1], [153, 1.9, 2.4], [165, 2.0, 2.5],
      [165.1, 2.1, 2.6], [183, 2.3, 2.8], [229, 2.8, 3.5], [330, 4.0, 5.0], [330.1, 4.1, 5.1]];
    for (const [d, fin, edge] of table) {
      const p = plain(c.derive('deep', BEDS.p1s, d));
      eq(`depth ${d} fin`, p.fin, fin); eq(`depth ${d} edge`, p.edge, edge);
    }
    for (const [bn, bed, max, fin, edge] of [['A1 mini', BEDS.a1mini, 153, 1.9, 2.4], ['MK4', BEDS.mk4, 183, 2.3, 2.8], ['P1S', BEDS.p1s, 229, 2.8, 3.5]]) {
      const lim = plain(c.limits('deep', bed)), p = plain(c.derive('deep', bed, lim.depthMax));
      eq(`${bn} depthMax`, lim.depthMax, max); eq(`${bn} fin at depthMax`, p.fin, fin); eq(`${bn} edge at depthMax`, p.edge, edge);
      const pd = plain(c.derive('deep', bed, lim.depthDefault));
      eq(`${bn} fin at default depth`, pd.fin, 1.6); eq(`${bn} edge at default depth`, pd.edge, 2.0);
    }
    for (let d = 40; d <= 573; d++) {                      // every integer depth the slider allows on a 600 mm custom bed
      const p = plain(c.derive('deep', BEDS.custom, d));
      if (p.fin < 1.6 || p.edge < 2.0 || p.fin * 82.5 < d - 1e-6 || p.edge * 66 < d - 1e-6 || p.fin > Math.max(1.6, d / 82.5 + 0.1) + 1e-9 || p.edge > Math.max(2.0, d / 66 + 0.1) + 1e-9) {
        rep.check(`depth ${d} bounds`, false, `fin ${p.fin}, edge ${p.edge}`);
      }
    }
    rep.done();
  });

  it('every standard bed is at least 20 mm smaller than the stage 1 tile in each tile dimension (owner request)', () => {
    const c = needCore();
    const rep = new Report('20 mm smaller');
    for (const [bn, bed] of Object.entries({ a1mini: BEDS.a1mini, mk4: BEDS.mk4, p1s: BEDS.p1s })) {
      for (const series of ['deep', 'blocks']) {
        const lim = plain(c.limits(series, bed));
        const p = c.derive(series, bed, lim.depthDefault);
        const d = series === 'deep' ? c.deepDims(p) : c.blockDims(p);
        rep.check(`${series} ${bn} fits the bed with margin`, d[0] <= bed.x - 20 && d[1] <= bed.y - 20,
          `tile ${plain(d)} on bed ${bed.x}x${bed.y}`);
      }
    }
    rep.done();
  });

  it('derive follows the v2 rules on all beds incl. custom, slim (Ny 11) and narrow (W reduced)', () => {
    const c = needCore();
    const rep = new Report('derive rules');
    for (const series of ['deep', 'blocks']) for (const [bn, bed] of Object.entries(BEDS)) {
      const depth = series === 'deep' ? 90 : 40;
      let got;
      try { got = withoutBed(c.derive(series, bed, depth)); } catch (e) { rep.check(`derive(${series}, ${bn})`, false, `threw ${e.message}`); continue; }
      const want = expectDerive(series, bed, depth);
      rep.check(`derive(${series}, ${bn}=${bed.x}x${bed.y}x${bed.z}, ${depth})`, JSON.stringify(Object.entries(got).sort()) === JSON.stringify(Object.entries(want).sort()),
        `got ${JSON.stringify(got)}, want ${JSON.stringify(want)}`);
    }
    rep.done();
  });

  it('limits(): min, max and default per the v2 rules on all beds', () => {
    const c = needCore();
    const rep = new Report('limits');
    for (const series of ['deep', 'blocks']) for (const [bn, bed] of Object.entries(BEDS)) {
      const got = plain(c.limits(series, bed)), want = expectLimits(series, bed);
      for (const k of ['depthMin', 'depthMax', 'depthDefault']) {
        rep.check(`${series} ${bn} ${k}`, typeof got[k] === 'number' && near(got[k], want[k], 1e-9), `got ${got[k]}, want ${want[k]}`);
      }
    }
    rep.done();
  });
});

// ================================================================ 4. geometry equality
function compareGeometry(item) {
  const c = needCore();
  const rep = new Report(`${item.id} (${item.series} tile ${item.tile}${item.shift ? ` shift ${item.shift}` : ''})`);
  const ref = reference(item);
  const params = plain(item.params);
  const t0 = performance.now();
  let man = null;
  try {
    try {
      man = item.series === 'deep' ? c.deepTile(params, item.tile) : c.blockTile(params, item.tile, item.shift[0], item.shift[1]);
    } catch (e) {
      rep.check('build', false, `${item.series === 'deep' ? 'deepTile' : 'blockTile'} threw: ${e && e.stack ? e.stack.split('\n').slice(0, 3).join(' | ') : e}`);
      rep.done();
    }
    const ms = performance.now() - t0;
    rec('slowest build (ms)', ms, item.id);

    // manifold health
    const status = man.status();
    rep.check('manifold status', status === 'NoError', `status ${status}`);
    rep.check('manifold not empty', !man.isEmpty(), 'empty');
    const js = meshFromManifold(man);
    const an = analyze(js);
    const bb = bbox(js);
    rep.check('manifold own volume() consistent with mesh', near(man.volume(), an.volume, 1e-4), `manifold ${man.volume()} vs mesh ${an.volume}`);

    // bounding box
    let dbb = 0, dbbWhere = '';
    for (let a = 0; a < 3; a++) for (const [nm, x, y] of [['min', bb.min[a], ref.bb.min[a]], ['max', bb.max[a], ref.bb.max[a]]]) {
      const d = Math.abs(x - y);
      if (d > dbb) { dbb = d; dbbWhere = `${nm} ${'xyz'[a]}`; }
    }
    rec('bounding box (mm)', dbb, item.id);
    rep.check('bounding box', dbb <= TOL_BB, () => `${dbbWhere} differs by ${f4(dbb)} mm; js min ${bb.min.map(f4)} max ${bb.max.map(f4)}, ref min ${ref.bb.min.map(f4)} max ${ref.bb.max.map(f4)}`);

    // volume
    const rel = Math.abs(an.volume - ref.an.volume) / ref.an.volume;
    rec('volume (relative)', rel, item.id);
    rep.check('volume', rel <= TOL_VOL, `js ${(an.volume / 1000).toFixed(3)} cm3, ref ${(ref.an.volume / 1000).toFixed(3)} cm3, ${(rel * 100).toFixed(4)} %`);

    // heights
    const hj = heights(js, item.series, params), hr = heights(ref.mesh, item.series, params);
    const { worst, at } = maxAbsDiff(hj.hs, hr.hs);
    rec('heights (mm)', Number.isFinite(worst) ? worst : 1e9, item.id);
    rep.check('heights (ray cast)', worst <= TOL_H, () => `${hj.tags[at]}: js ${f4(hj.hs[at])}, ref ${f4(hr.hs[at])}, diff ${f4(worst)} mm; ${count(hj.hs, Number.isNaN)} js rays without hit`);

    // shells
    rep.check('one closed outer solid', an.outer === 1 && an.open === 0 && an.nonManifoldEdgeUses === 0 && an.zero === 0 && an.degenerate === 0,
      `outer ${an.outer}, shells with open/inconsistent edges ${an.open}, non-2-manifold edge uses ${an.nonManifoldEdgeUses}, zero-volume shells ${an.zero}, degenerate triangles ${an.degenerate}`);
    rep.check('inner chamber shells', an.inner === ref.an.inner, `js ${an.inner}, reference ${ref.an.inner}${item.series === 'blocks' ? ' (blocks must have none)' : ''}`);
    if (an.inner === ref.an.inner) {
      let wi = 0, wiAt = -1;
      an.innerVolumes.forEach((x, i) => { const d = Math.abs(x - ref.an.innerVolumes[i]) / ref.an.innerVolumes[i]; if (d > wi) { wi = d; wiAt = i; } });
      rec('inner chamber volume (relative)', wi, item.id);
      rep.check('inner chamber volumes', wi <= TOL_VOL, () => `chamber #${wiAt} (sorted by size): js ${an.innerVolumes[wiAt].toFixed(2)} mm3, ref ${ref.an.innerVolumes[wiAt].toFixed(2)} mm3`);
    }
    rep.check('total shells', an.shells === ref.an.shells, `js ${an.shells}, reference ${ref.an.shells}`);

    // label (tile number + arrow)
    const lj = labelHits(js, item.series, params), lr = labelHits(ref.mesh, item.series, params);
    const lw = maxAbsDiff(lj, lr);
    rec('label depth (mm)', Number.isFinite(lw.worst) ? lw.worst : 1e9, item.id);
    rep.check('label (tile number)', lw.worst <= TOL_H, () => `${count(lj.map((x, i) => (Math.abs(x - lr[i]) > TOL_H ? 1 : 0)), (x) => x === 1)} of ${lj.length} label rays differ; js has ${count(lj, (x) => Math.abs(x - 0.4) < 1e-3)} rays inside the label, ref ${count(lr, (x) => Math.abs(x - 0.4) < 1e-3)}`);
  } finally {
    if (man) man.delete();
  }
  rep.done();
}

describe('equality: the 29 stage 1 tiles (JS generator vs stage 1 STL)', () => {
  for (const item of TILES) it(item.id, () => compareGeometry(item));
});

describe('equality: fixture cases (random params, Python reference STL from scripts/reference.py)', () => {
  for (const item of CASE_ITEMS) it(`${item.id} - ${item.note}`, () => compareGeometry(item));
});

// ================================================================ 5. v2 members: walls, wall band, features, estimates
const SPEED = 343;
const surfaceArea = (mesh) => {
  const { v, t } = mesh;
  let s = 0;
  for (let f = 0; f < t.length; f += 3) {
    const a = 3 * t[f], b = 3 * t[f + 1], c = 3 * t[f + 2];
    const ux = v[b] - v[a], uy = v[b + 1] - v[a + 1], uz = v[b + 2] - v[a + 2];
    const wx = v[c] - v[a], wy = v[c + 1] - v[a + 1], wz = v[c + 2] - v[a + 2];
    s += 0.5 * Math.hypot(uy * wz - uz * wy, uz * wx - ux * wz, ux * wy - uy * wx);
  }
  return s;
};

/** Params sets with a known Python band: stage 1 presets and fixture cases. */
const WALL_SETS = [];
for (const size of ['S', 'L']) {
  WALL_SETS.push({ id: `deep ${size} preset`, series: 'deep', params: EXP.presets.deep[size], band: EXP.presets_band.deep[size].band, dims: EXP.presets_band.deep[size].dims });
  WALL_SETS.push({ id: `blocks ${size} preset`, series: 'blocks', params: EXP.presets.blocks[size], band: EXP.presets_band.blocks[size].band, dims: EXP.presets_band.blocks[size].dims });
}
for (const it2 of CASE_ITEMS) WALL_SETS.push({ id: `case ${it2.id}`, series: it2.series, params: it2.params, band: EXP.cases[it2.id].band, dims: EXP.cases[it2.id].dims });

/** Wall options expected from layouts.json (independent of the generator's LAYOUTS copy). */
function expectedWalls(series, p) {
  const out = new Map();
  out.set('1x1', [[1]]);
  if (series === 'deep') {
    const room = LAYOUTS_JSON.deep[p.pitch < 30 ? 'S' : 'L'].room_left_to_right;
    for (const w of Object.keys(room).map(Number).sort((a, b) => a - b)) for (const rows of [1, 2, 3]) {
      out.set(`${w}x${rows}`, Array.from({ length: rows }, () => room[String(w)].slice()));
    }
  } else {
    const grids = LAYOUTS_JSON.blocks[p.Ny === 13 ? 'S' : 'L'].grids;
    for (const g of Object.values(grids)) out.set(`${g[0].length}x${g.length}`, g.map((r) => r.slice()));
  }
  return out;                                           // key "cols x rows"
}

describe('core v2: wallOptions', () => {
  const sets = [...WALL_SETS];
  for (const [bn, bed] of Object.entries({ a1mini: BEDS.a1mini, mk4: BEDS.mk4, p1s: BEDS.p1s })) {
    for (const series of ['deep', 'blocks']) sets.push({ id: `${series} derived on ${bn}`, series, bedDerive: [series, bed] });
  }
  for (const s of sets) {
    it(s.id, () => {
      const c = needCore();
      const params = s.params ? plain(s.params) : plain(c.derive(s.bedDerive[0], s.bedDerive[1], c.limits(...s.bedDerive).depthDefault));
      const got = plain(c.wallOptions(s.series, params));
      const want = expectedWalls(s.series, params);
      const rep = new Report(`wallOptions(${s.id})`);
      rep.check('is a list', Array.isArray(got) && got.length > 0, `got ${JSON.stringify(got).slice(0, 80)}`);
      if (!Array.isArray(got) || !got.length) return rep.done();
      rep.check('first entry is tile 1 alone', got[0].cols === 1 && got[0].rows === 1 && got[0].tiles === 1 && JSON.stringify(got[0].grid) === '[[1]]',
        `first entry ${JSON.stringify(got[0])}`);
      const seen = new Set();
      const ids = new Set();
      for (const w of got) {
        const key = `${w.cols}x${w.rows}`;
        rep.check(`${key} unique`, !seen.has(key), 'duplicate cols x rows');
        seen.add(key);
        rep.check(`${key} id`, (typeof w.id === 'string' || typeof w.id === 'number') && !ids.has(w.id), `id ${JSON.stringify(w.id)} missing or duplicated`);
        ids.add(w.id);
        const g = w.grid;
        rep.check(`${key} grid shape`, Array.isArray(g) && g.length === w.rows && g.every((r) => Array.isArray(r) && r.length === w.cols),
          `grid ${JSON.stringify(g)} is not ${w.rows} rows of ${w.cols}`);
        rep.check(`${key} tiles = cols x rows`, w.tiles === w.cols * w.rows, `tiles ${w.tiles}`);
        rep.check(`${key} grid content`, JSON.stringify(g) === JSON.stringify(want.get(key)), `got ${JSON.stringify(g)}, want ${JSON.stringify(want.get(key))}`);
      }
      const missing = [...want.keys()].filter((k) => !seen.has(k)), extra = [...seen].filter((k) => !want.has(k));
      rep.check('set of walls', !missing.length && !extra.length, `missing ${missing}, unexpected ${extra}`);
      rep.done();
    });
  }
});

describe('core v2: wallBand', () => {
  const picks = [['deep S preset', 4, 1], ['deep S preset', 7, 3], ['deep L preset', 3, 1], ['deep L preset', 7, 2], ['blocks S preset', 2, 2],
    ['blocks S preset', 3, 3], ['blocks L preset', 2, 3], ['case deep-a', 5, 1], ['case deep-b', 5, 2], ['case blocks-b', 3, 2], ['case blocks-c', 2, 2],
    ['deep S preset', 1, 1], ['blocks L preset', 1, 1]];
  for (const [id, cols, rows] of picks) {
    it(`${id}, ${cols} x ${rows}`, () => {
      const c = needCore();
      const s = WALL_SETS.find((x) => x.id === id);
      const got = plain(c.wallBand(s.series, plain(s.params), cols, rows));
      const width = cols * s.dims[0];
      const height = rows * s.dims[1];                    // deep: length along the wells; blocks: H (assumption, not in the contract)
      const rep = new Report(`wallBand(${id}, ${cols}, ${rows})`);
      const eq = (k, want) => rep.check(k, typeof got[k] === 'number' && near(got[k], want, 1e-9), `got ${got[k]}, want ${want}`);
      eq('high', s.band.f_high);
      // deep walls narrower than 3.5 design wavelengths start half an octave above f0 (2D BEM check of the old order, kept as a cautious start)
      const start = s.series === 'deep' && width < 3.5 * s.band.lambda0 - 1e-6 ? Math.SQRT2 * s.band.f0 : s.band.f0;
      eq('low', Math.max(start, SPEED / (width / 1000)));
      eq('wallWidthMm', width);
      eq('wallHeightMm', height);
      rep.check('low < high', got.low < got.high, `low ${got.low}, high ${got.high}`);
      rep.done();
    });
  }
  it('wider walls reach lower (low end falls with wall width until f0, never below it)', () => {
    const c = needCore();
    const s = WALL_SETS.find((x) => x.id === 'deep L preset');
    const lows = [1, 3, 5, 7].map((cols) => c.wallBand('deep', plain(s.params), cols, 1).low);
    const rep = new Report('low vs width');
    for (let i = 1; i < lows.length; i++) rep.check(`cols step ${i}`, lows[i] <= lows[i - 1] + 1e-9, `lows ${lows}`);
    rep.check('never below f0', lows.every((l) => l >= s.band.f0 - 1e-9), `lows ${lows}, f0 ${s.band.f0}`);
    rep.done();
  });
});

describe('core v2: features(manifold) vs the same quantities measured on the reference STL', () => {
  for (const item of [...TILES, ...CASE_ITEMS]) {
    it(item.id, () => {
      const c = needCore();
      const ref = reference(item);
      const rep = new Report(`features ${item.id}`);
      let man = null;
      try {
        man = item.series === 'deep' ? c.deepTile(plain(item.params), item.tile) : c.blockTile(plain(item.params), item.tile, item.shift[0], item.shift[1]);
        const f = plain(c.features(man));
        const want = {
          area_cm2: surfaceArea(ref.mesh) / 100,
          volume_cm3: ref.an.volume / 1000,
          height_mm: ref.bb.max[2] - ref.bb.min[2],
          footprint_cm2: (ref.bb.max[0] - ref.bb.min[0]) * (ref.bb.max[1] - ref.bb.min[1]) / 100,
        };
        for (const k of Object.keys(want)) {
          const ok = typeof f[k] === 'number' && (k === 'height_mm' ? Math.abs(f[k] - want[k]) <= TOL_H : Math.abs(f[k] - want[k]) <= TOL_VOL * want[k]);
          rec(`features ${k} (${k === 'height_mm' ? 'mm' : 'relative'})`, typeof f[k] === 'number' ? (k === 'height_mm' ? Math.abs(f[k] - want[k]) : Math.abs(f[k] - want[k]) / want[k]) : 1e9, item.id);
          rep.check(k, ok, `got ${f[k]}, reference STL gives ${want[k]}`);
        }
      } finally { if (man) man.delete(); }
      rep.done();
    });
  }
});

// ---------------------------------------------------------------- estimates (need scripts/calibration/{data,estimate,fit-report}.json from the calibration agent)
const CALIB = process.env.DIFFUSER_CALIB || path.join(HERE, 'calibration');
const ESTIMATE_PATH = path.join(CALIB, 'estimate.json');
const DATA_PATH = path.join(CALIB, 'data.json');
const skipNoEstimate = fs.existsSync(ESTIMATE_PATH) ? false : `scripts/calibration/estimate.json not found yet (calibration agent) - estimate tests skipped`;
const REPORT_PATH = path.join(CALIB, 'fit-report.json');
const skipNoData = [DATA_PATH, ESTIMATE_PATH, REPORT_PATH].every((f) => fs.existsSync(f)) ? false
  : `scripts/calibration/data.json, estimate.json and fit-report.json (the error it reports) not all there yet (calibration agent) - skipped`;
const FEATURE_KEYS = ['area_cm2', 'volume_cm3', 'height_mm', 'footprint_cm2'];
const applyCoef = (coef, feats) => Object.entries(coef).reduce((s, [k, v]) => s + (k === 'const' ? v : v * (feats[k] ?? 0)), 0);

/** data.json: array of rows {id, series, features, grams, minutes, ...} (calibration agent). */
function calibRows(data) {
  const rows = Array.isArray(data) ? data : (data.rows || data.data || data.samples || data.measurements || data.tiles || []);
  return rows.filter((r) => r.in_fit !== false).map((r) => ({       // in_fit false: deep rule v1/v2 rows kept as a record, not fitted
    id: r.id ?? r.name ?? r.file ?? '?', series: r.series,
    feats: r.features ?? Object.fromEntries(FEATURE_KEYS.filter((k) => k in r).map((k) => [k, r[k]])),
    grams: r.grams ?? r.grams_measured ?? r.measured?.grams ?? r.g,
    minutes: r.minutes ?? r.minutes_measured ?? r.measured?.minutes ?? r.min,
  }));
}
/** Reported fit error of the calibration agent: fit-report.json, [series][quantity].train_max_pct (worst relative error of
 *  the fit on its own data rows, rounded to 0.01 percentage points). Returns {tol (relative), note} or null. */
function reportedTolerance(report, series, q) {
  const r = report?.[series]?.[q];
  if (r && typeof r.train_max_pct === 'number') return { tol: r.train_max_pct / 100 + 5e-5, note: `fit-report.json ${series}.${q}.train_max_pct = ${r.train_max_pct} %` };
  return null;
}

describe('core v2: ESTIMATE and estimate()', () => {
  it('ESTIMATE equals scripts/calibration/estimate.json', { skip: skipNoEstimate }, () => {
    const c = needCore();
    const file = readJson(ESTIMATE_PATH);
    assert.deepStrictEqual(plain(c.ESTIMATE), file);
  });

  it('estimate(series, feats) = sum of coefficient x feature (+ const), both series, grams and minutes', { skip: skipNoEstimate }, () => {
    const c = needCore();
    const file = readJson(ESTIMATE_PATH);
    const rep = new Report('estimate formula');
    const featSets = [
      { area_cm2: 400, volume_cm3: 250, height_mm: 120, footprint_cm2: 300 },
      { area_cm2: 900, volume_cm3: 1100, height_mm: 60, footprint_cm2: 480 },
      { area_cm2: 1500, volume_cm3: 700, height_mm: 220, footprint_cm2: 520 },
    ];
    for (const series of ['deep', 'blocks']) for (const f of featSets) {
      const got = plain(c.estimate(series, f));
      for (const q of ['grams', 'minutes']) {
        const want = applyCoef(file[series][q], f);
        rep.check(`${series} ${q} @${f.area_cm2}`, typeof got[q] === 'number' && near(got[q], want, 1e-9), `got ${got[q]}, want ${want}`);
      }
    }
    rep.done();
  });

  it('estimate() reproduces every row of data.json within the error fit-report.json reports', { skip: skipNoData }, () => {
    const c = needCore();
    const data = readJson(DATA_PATH);
    const report = readJson(REPORT_PATH);
    const rows = calibRows(data);
    assert.ok(rows.length > 0, `data.json: no rows recognised (keys: ${Object.keys(data)})`);
    const rep = new Report('calibration rows');
    let worst = {};
    for (const r of rows) {
      for (const q of ['grams', 'minutes']) {
        if (typeof r[q] !== 'number') { rep.check(`${r.id} ${q}`, false, 'row has no measured value (data.json layout not understood)'); continue; }
        const rt = reportedTolerance(report, r.series, q);
        if (!rt) { rep.check(`${r.id} ${q}`, false, 'fit-report.json has no train_max_pct for this series/quantity'); continue; }
        const est = plain(c.estimate(r.series, r.feats))[q];
        const rel = Math.abs(est - r[q]) / r[q];
        const key = `${r.series} ${q}`;
        if (!worst[key] || rel > worst[key].rel) worst[key] = { rel, tol: rt.tol, id: r.id };
        rep.check(`${r.series} ${r.id} ${q}`, rel <= rt.tol + 1e-9, `estimate ${est.toFixed(2)}, measured ${r[q]}, off by ${(rel * 100).toFixed(2)} % (${rt.note})`);
      }
    }
    for (const [k, w] of Object.entries(worst)) rec(`estimate vs data.json ${k} (relative, tol ${w.tol.toFixed(3)})`, w.rel, w.id);
    rep.done();
  });

  it('estimate is positive and increasing from P1S size to 1.5x (default depth, tile 1)', { skip: skipNoEstimate }, () => {
    const c = needCore();
    const rep = new Report('estimate vs size');
    for (const series of ['deep', 'blocks']) {
      const grams = [], minutes = [];
      for (const s of [1, 1.25, 1.5]) {
        const bed = { x: 256 * s, y: 256 * s, z: 256 * s };
        const params = plain(c.derive(series, bed, plain(c.limits(series, bed)).depthDefault));
        const kl = series === 'blocks' ? LAYOUTS_JSON.blocks[params.Ny === 13 ? 'S' : 'L'].shifts[0] : null;
        const man = series === 'deep' ? c.deepTile(params, 1) : c.blockTile(params, 1, kl[0], kl[1]);
        try {
          const e = plain(c.estimate(series, plain(c.features(man))));
          grams.push(e.grams); minutes.push(e.minutes);
        } finally { man.delete(); }
      }
      rep.check(`${series} positive`, grams.every((x) => x > 0) && minutes.every((x) => x > 0), `grams ${grams.map((x) => x.toFixed(1))}, minutes ${minutes.map((x) => x.toFixed(1))}`);
      rep.check(`${series} grams increase`, grams[0] < grams[1] && grams[1] < grams[2], `grams ${grams.map((x) => x.toFixed(1))}`);
      rep.check(`${series} minutes increase`, minutes[0] < minutes[1] && minutes[1] < minutes[2], `minutes ${minutes.map((x) => x.toFixed(1))}`);
    }
    rep.done();
  });
});

describe('summary', () => {
  it('worst deviation per measurement', (t) => {
    const keys = Object.keys(STATS);
    const lines = keys.length ? keys.map((k) => `${k}: ${STATS[k].value.toPrecision(4)} (${STATS[k].label})`) : ['no geometry comparison ran'];
    t.diagnostic(lines.join('\n'));
    console.log('\nWorst deviation per measurement:\n  ' + lines.join('\n  '));
    console.log(`Tolerances: heights ${TOL_H} mm, bbox ${TOL_BB} mm, volume ${TOL_VOL * 100} %.`);
  });
});
