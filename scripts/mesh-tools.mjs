// Mesh tooling shared by the bass traps tests (not a test file: the npm test glob is scripts/*.test.mjs).
// The functions weld, parseStl, meshFromManifold, bbox, triVolume, analyze and firstHit are copied VERBATIM from
// scripts/equality.test.mjs (the measuring code of the diffuser tiles), so the JS side and the Python reference STL are measured
// by the same code as there; equality.test.mjs itself is untouched. Added at the end: the vm loader of the bass cores.

import fs from 'node:fs';
import path from 'node:path';
import vm from 'node:vm';
import { fileURLToPath } from 'node:url';
import Module from 'manifold-3d';

const EPS_VOL = 1e-3;      // mm^3, a shell below this volume counts as zero volume

export function weld(pos, tris) {
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

export function parseStl(buf) {
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

export function meshFromManifold(man) {
  const g = man.getMesh();
  const np = g.numProp;
  const nv = g.vertProperties.length / np;
  const pos = new Float32Array(3 * nv);
  for (let i = 0; i < nv; i++) for (let j = 0; j < 3; j++) pos[3 * i + j] = g.vertProperties[np * i + j];
  return weld(pos, g.triVerts);
}

export function bbox(mesh) {
  const min = [Infinity, Infinity, Infinity], max = [-Infinity, -Infinity, -Infinity];
  for (let i = 0; i < mesh.v.length; i += 3) for (let a = 0; a < 3; a++) {
    min[a] = Math.min(min[a], mesh.v[i + a]); max[a] = Math.max(max[a], mesh.v[i + a]);
  }
  return { min, max };
}

export function triVolume(v, a, b, c) {
  const ax = v[3 * a], ay = v[3 * a + 1], az = v[3 * a + 2];
  const bx = v[3 * b], by = v[3 * b + 1], bz = v[3 * b + 2];
  const cx = v[3 * c], cy = v[3 * c + 1], cz = v[3 * c + 2];
  return (ax * (by * cz - bz * cy) - ay * (bx * cz - bz * cx) + az * (bx * cy - by * cx)) / 6;
}

/** Shells = faces connected through shared edges. Volume > 0: outer solid, < 0: inner cavity (chamber), ~0: sheet. */
export function analyze(mesh) {
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
export function firstHit(mesh, o, d) {
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

// ---------------------------------------------------------------- loading the bass cores like the diffuser core
const HERE = path.dirname(fileURLToPath(import.meta.url));
export const GEN = path.resolve(HERE, '..');
export const FIX = path.join(HERE, 'fixtures', 'bass');
export const readJson = (p) => JSON.parse(fs.readFileSync(p, 'utf8'));
export const plain = (x) => JSON.parse(JSON.stringify(x));      // strips the vm realm

function scriptText(html, id) {
  const tag = new RegExp(`<script\\b([^>]*)\\bid\\s*=\\s*["']${id}["']([^>]*)>([\\s\\S]*?)<\\/script>`, 'i').exec(html);
  if (!tag) throw new Error(`index.html has no <script id="${id}">`);
  if (/type\s*=\s*["']?module/i.test(tag[1] + tag[2])) throw new Error(`<script id="${id}"> must be a classic script, not type="module"`);
  return tag[3];
}

/** Loads <script id="bass-physics"> (and, with withGeometry, <script id="bass-geometry"> and the manifold module) in a node:vm
 *  context holding only Math and console, the way equality.test.mjs loads diffuser-core. Returns {physics, geometry, M, html, sources}. */
export async function loadBass({ withGeometry = false } = {}) {
  const html = fs.readFileSync(process.env.DIFFUSER_INDEX || path.join(GEN, 'index.html'), 'utf8');
  const sources = { physics: scriptText(html, 'bass-physics'), geometry: scriptText(html, 'bass-geometry') };
  const ctx = vm.createContext({ Math, console });
  Object.setPrototypeOf(vm.runInContext('Array', ctx).prototype, Array.prototype);   // see equality.test.mjs: vm realm arrays vs manifold-3d
  const physicsFactory = vm.runInContext(sources.physics + '\n;(typeof createBassPhysics === "function") ? createBassPhysics : undefined', ctx, { filename: 'index.html#bass-physics' });
  if (typeof physicsFactory !== 'function') throw new Error('the physics script does not define a global function createBassPhysics()');
  const physics = physicsFactory();
  if (!withGeometry) return { physics, html, sources };
  const M = await Module();
  M.setup();
  const geometryFactory = vm.runInContext(sources.geometry + '\n;(typeof createBassGeometry === "function") ? createBassGeometry : undefined', ctx, { filename: 'index.html#bass-geometry' });
  if (typeof geometryFactory !== 'function') throw new Error('the geometry script does not define a global function createBassGeometry(M, physics)');
  return { physics, geometry: geometryFactory(M, physics), M, html, sources };
}
