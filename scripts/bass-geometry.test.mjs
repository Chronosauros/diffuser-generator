// Bass traps, geometry core: <script id="bass-geometry"> of index.html (createBassGeometry) against its Python twin
// (reference/bass_cell.py; fixtures in scripts/fixtures/bass/: geometry.json and the binary STL written by scripts/bass_reference.py).
//
//   npm test            (= node --test scripts/*.test.mjs)
//
// Tolerances (spec of the bass core): rectangle plans 1e-6 mm, STL volume 0.1 %, bounding box 0.01 mm, sections and heights by ray
// cast 0.01 mm. The JS mesh and the reference STL are measured by the SAME code (scripts/mesh-tools.mjs = the measuring code of
// equality.test.mjs). The volumes of the research package STL (w = 6 mm blocks for 42.9 and 95.2 Hz) are constants below.

import { describe, it } from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import path from 'node:path';
import vm from 'node:vm';
import { FIX, GEN, readJson, loadBass, plain, parseStl, meshFromManifold, bbox, analyze, firstHit } from './mesh-tools.mjs';

const GEO = readJson(path.join(FIX, 'geometry.json'));
const { physics: P, geometry: G, M, html, sources } = await loadBass({ withGeometry: true });

const TOL_PLAN = 1e-6;     // mm
const TOL_NUM = 1e-9;      // relative, numbers of the design (the fixtures hold 10 digits)
const TOL_VOL = 1e-3;      // relative, volume of a mesh
const TOL_BB = 0.01;       // mm, bounding box
const TOL_H = 0.01;        // mm, ray cast

// STL of the research package (cell/stl-v2, w = 6 mm, the current defaults), volume of the float32 mesh in mm3: the numbers the report quotes
const RESEARCH = {
  'bass-f042.9-1x1-k1-block.stl': { file: 'qwr2-m100-f042.9-w6-block.stl', volume: 354810.1005055583, size: [123.6, 120, 120] },
  'bass-f095.2-1x2-k2-block.stl': { file: 'qwr2-m111-f095.2-w6-block.stl', volume: 352235.42079070956, size: [58.8, 120, 240] },
  'bass-f095.2-inserts-plate.stl': { file: 'qwr2-m111-f095.2-w6-inserts-plate.stl', volume: 252491.09130664569, size: [197.2, 119, 53.469] },
};

const near = (a, b, rel = TOL_NUM, abs = 1e-12) => Math.abs(a - b) <= abs + rel * Math.max(Math.abs(a), Math.abs(b));
class Report {
  constructor(label) { this.label = label; this.fails = []; this.n = 0; }
  check(name, a, b, rel = TOL_NUM, abs = 1e-12) { this.n++; if (!near(a, b, rel, abs)) this.fails.push(`${name}: js ${a} python ${b}`); }
  equal(name, a, b) { this.n++; if (JSON.stringify(a) !== JSON.stringify(b)) this.fails.push(`${name}: js ${JSON.stringify(a)} python ${JSON.stringify(b)}`); }
  done() { if (this.fails.length) assert.fail(`${this.label}: ${this.fails.length} of ${this.n} differ\n  - ${this.fails.slice(0, 12).join('\n  - ')}`); }
}
const designs = new Map();
const designOf = (f) => { if (!designs.has(f)) designs.set(f, G.design(f)); return designs.get(f); };

/** First-hit z of rays from above on a grid in x-y (the section of the block seen from the top: solid roof or channel floor). */
function topGrid(mesh, size, step = 1.3, off = 0.37) {
  const out = [];
  for (let x = off; x < size[0]; x += step) for (let y = off; y < size[1]; y += step) {
    const h = firstHit(mesh, [x, y, size[2] + 10], [0, 0, -1]);
    out.push(h ? h[2] : NaN);
  }
  return out;
}
function maxDiff(a, b) {
  let worst = 0;
  for (let i = 0; i < a.length; i++) {
    const d = (Number.isNaN(a[i]) || Number.isNaN(b[i])) ? (Number.isNaN(a[i]) && Number.isNaN(b[i]) ? 0 : Infinity) : Math.abs(a[i] - b[i]);
    if (d > worst) worst = d;
  }
  return worst;
}
const stlOf = (file) => parseStl(fs.readFileSync(path.join(FIX, file)));
const freeAll = (...ms) => ms.forEach((m) => m && m.delete());

describe('bass-geometry: the script is a classic core that needs the physics core', () => {
  it('has no import, DOM or fetch, and defines createBassGeometry(M, physics)', () => {
    assert.match(sources.geometry, /function createBassGeometry\(M, physics\)/);
    assert.doesNotMatch(sources.geometry.replace(/\/\*[\s\S]*?\*\//g, '').replace(/\/\/.*$/gm, '').replace(/\.window\b|\bwindow:|"window"/g, ''), /\bimport\b|\brequire\s*\(|\bwindow\b|\bdocument\b|\bfetch\b/);
  });
  it('refuses to run without the physics object instead of failing later', () => {
    const ctx = vm.createContext({ Math, console });
    const factory = vm.runInContext(sources.geometry + '\n;createBassGeometry', ctx);
    assert.throws(() => factory(M), /physics/);
  });
  it('printer beds are the beds of the diffuser generator', () => {
    const m = html.match(/<script id="diffuser-core">([\s\S]*?)<\/script>/)[1];
    const ctx = vm.createContext({ Math, console });
    Object.setPrototypeOf(vm.runInContext('Array', ctx).prototype, Array.prototype);
    const core = vm.runInContext(m + '\n;createDiffuser', ctx)(M);
    assert.deepEqual(plain(G.BEDS), plain(core.BEDS));
    assert.equal(G.BED_MARGIN, 5);
  });
  it('the print estimate coefficients are the ones fitted by scripts/calibration/bass-calibrate.py', () => {
    assert.deepEqual(plain(G.ESTIMATE), readJson(path.join(GEN, 'scripts', 'calibration', 'bass-estimate.json')));
  });
});

describe('bass-geometry: design of the cell against the Python twin', () => {
  for (const [key, fx] of Object.entries(GEO.designs)) {
    it(`mode ${key} Hz: channel length, legs, inserts`, () => {
      const d = designOf(Number(key)), r = new Report(`design ${key}`);
      assert.equal(d.ok, true);
      r.check('L0', d.L0, fx.L0);
      r.check('f0', d.f0, fx.f0);
      r.equal('legs', d.plan.n, fx.legs);
      r.equal('bends', d.plan.bends, fx.bends);
      r.check('y_in', d.plan.y_in, fx.y_in);
      r.check('p_nom', d.plan.p_nom, fx.p_nom);
      r.check('L_geom', d.plan.L_geom, fx.L_geom);
      r.check('depth', d.depth, fx.depth);
      r.check('p_limit', d.p_limit, fx.p_limit);
      r.check('bend_dL', d.bend_dL, fx.bend_dL);
      r.equal('ns_down', d.ns_down, fx.ns_down);
      r.equal('ns_up', d.ns_up, fx.ns_up);
      r.equal('warnings', d.warn.length, fx.warn.length);
      r.equal('insert keys', Object.keys(d.pk).map(Number).sort((a, b) => a - b), Object.keys(fx.pk).map(Number).sort((a, b) => a - b));
      for (const k of Object.keys(fx.pk)) { r.check(`L_k ${k}`, d.Lk[k], fx.Lk[k]); r.check(`p_k ${k}`, d.pk[k], fx.pk[k]); }
      r.done();
    });
  }
  it('rules of the design: f cell = f mode x 1.005, bend correction -0.245 w, inserts -1 ... +4 % in 0.5 % steps, w = 6, wall 1.2', () => {
    const d = designOf(42.9);
    assert.ok(Math.abs(d.f0 - 42.9 * 1.005) < 1e-12);
    assert.equal(d.w, 6); assert.equal(d.wall, 1.2); assert.equal(d.width, 120); assert.equal(d.height, 120);
    assert.ok(Math.abs(d.bend_dL - (-0.245 * 6)) < 1e-12);
    assert.equal(d.ns_down, 2); assert.equal(d.ns_up, 8); assert.equal(d.step_pct, 0.5);
    // the acoustic length reproduces the cell frequency on the grid of first_res
    assert.ok(Math.abs(P.firstRes(6, d.L0) - d.f0) < 1e-6);
    // the physics cell of a mode uses the same L0 (one source of the channel length)
    assert.equal(P.cellForMode(42.9).L, d.L0);
  });
  it('the inserts get longer with k (shorter channel, higher f) and k = 0 is the plan\'s nominal insert', () => {
    const d = designOf(95.2);
    let prev = -Infinity;
    for (let k = -d.ns_down; k <= d.ns_up; k++) { assert.ok(d.pk[k] > prev); prev = d.pk[k]; assert.ok(d.Lk[k] <= d.L0 + 1e-9 || k < 0); }
    assert.equal(d.pk[0], d.plan.p_nom);
    assert.ok(d.pk[d.ns_up] <= d.p_limit);
  });
  it('out-of-range modes and silly input give {ok: false}, not an exception', () => {
    assert.equal(G.design(20).ok, false);
    assert.equal(G.design(500).ok, false);
    assert.equal(G.design(NaN).ok, false);
    assert.match(G.design(20).error, /35 - 120/);
  });
});

describe('bass-geometry: rectangle plans of the blocks against the Python twin', () => {
  const rectEq = (r, a, b, tag) => {
    r.equal(`${tag} kind`, a.kind, b.kind);
    for (const ax of ['x', 'y', 'z']) for (let i = 0; i < 2; i++) r.check(`${tag} ${ax}${i}`, a[ax][i], b[ax][i], 0, TOL_PLAN);
  };
  for (const fx of GEO.plans) {
    it(`mode ${fx.f_mode} Hz, block ${fx.ny}x${fx.nz}, ${fx.k} cell${fx.k > 1 ? 's' : ''}`, () => {
      const d = designOf(fx.f_mode), bp = G.planBlock(d, fx.ny, fx.nz, fx.k), r = new Report(`plan ${fx.f_mode} ${fx.ny}x${fx.nz} k${fx.k}`);
      const fb = fx.block;
      fb.size_mm.forEach((v, i) => r.check(`size ${i}`, bp.size_mm[i], v, 0, TOL_PLAN));
      r.check('L_axis', bp.L_axis, fb.L_axis, TOL_NUM);
      r.check('L_acoustic', bp.L_acoustic, fb.L_acoustic, TOL_NUM);
      r.check('volume', bp.volume_mm3, fb.volume_mm3, TOL_NUM);
      r.check('x_last', bp.x_last, fb.x_last, 0, TOL_PLAN);
      r.equal('cells', bp.cells.length, fb.cells.length);
      fb.cells.forEach((fc, i) => {
        const c = bp.cells[i], tag = `cell ${i}`;
        r.equal(`${tag} position`, [c.index_y, c.index_z, c.mirror_y], [fc.index_y, fc.index_z, fc.mirror_y]);
        r.equal(`${tag} air count`, c.air.length, fc.air.length);
        fc.air.forEach((a, j) => rectEq(r, c.air[j], a, `${tag} air ${j}`));
        r.equal(`${tag} void count`, c.voids.length, fc.voids.length);
        fc.voids.forEach((a, j) => rectEq(r, c.voids[j], a, `${tag} void ${j}`));
        if (fc.walls) { r.equal(`${tag} wall count`, c.walls.length, fc.walls.length); fc.walls.forEach((a, j) => rectEq(r, c.walls[j], a, `${tag} wall ${j}`)); }
        r.equal(`${tag} axis count`, c.axis_points_xy.length, fc.axis_points_xy.length);
        fc.axis_points_xy.forEach((q, j) => { r.check(`${tag} axis ${j} x`, c.axis_points_xy[j][0], q[0], 0, TOL_PLAN); r.check(`${tag} axis ${j} y`, c.axis_points_xy[j][1], q[1], 0, TOL_PLAN); });
        r.check(`${tag} window plane`, c.window.outer_wall_plane_y, fc.window.outer_wall_plane_y, 0, TOL_PLAN);
        r.equal(`${tag} insert axis`, c.window.insert_axis_into_cell, fc.window.insert_axis_into_cell);
      });
      r.done();
    });
  }
  it('the nominal insert (k = 0) closes the channel: its body spans the wall plane to the end of the axis', () => {
    for (const f of [35, 42.9, 68.7, 95.2, 120]) {
      const d = designOf(f), bp = G.planBlock(d, 2, 2, 4);
      for (const c of bp.cells) {
        const end = c.axis_points_xy[c.axis_points_xy.length - 1][1];
        assert.ok(Math.abs(Math.abs(c.window.outer_wall_plane_y - end) - d.pk[0]) < 2e-4, `f ${f} cell ${c.index}: ${Math.abs(c.window.outer_wall_plane_y - end)} vs p ${d.pk[0]}`);
      }
      assert.ok(Math.abs(bp.L_acoustic - d.L0) < 5e-4, `acoustic length ${bp.L_acoustic} vs ${d.L0}`);
    }
  });
});

describe('bass-geometry: solids of the blocks (manifold) against their rectangle plans', () => {
  for (const fx of GEO.plans) {
    it(`mode ${fx.f_mode} Hz, block ${fx.ny}x${fx.nz} k${fx.k}: one closed solid, pockets, volume`, () => {
      const d = designOf(fx.f_mode), bp = G.planBlock(d, fx.ny, fx.nz, fx.k), m = G.blockSolid(d, bp);
      try {
        const mesh = meshFromManifold(m), an = analyze(mesh), bb = bbox(mesh);
        const pockets = bp.cells.reduce((s, c) => s + c.voids.length, 0);
        assert.equal(an.outer, 1, 'one outer solid');
        assert.equal(an.inner, pockets, 'the closed pockets are the only cavities');
        assert.equal(an.zero, 0);
        assert.equal(an.open, 0);
        assert.equal(an.nonManifoldEdgeUses, 0);
        assert.equal(an.degenerate, 0);
        assert.ok(Math.abs(an.volume / bp.volume_mm3 - 1) <= TOL_VOL, `volume ${an.volume} vs plan ${bp.volume_mm3}`);
        for (let a = 0; a < 3; a++) {
          assert.ok(Math.abs(bb.min[a]) <= TOL_BB, `min ${a} ${bb.min[a]}`);
          assert.ok(Math.abs(bb.max[a] - bp.size_mm[a]) <= TOL_BB, `max ${a} ${bb.max[a]} vs ${bp.size_mm[a]}`);
        }
      } finally { freeAll(m); }
    });
  }
  it('planFeatures (no solid) agree with the features of the solid', () => {
    const d = designOf(95.2), bp = G.planBlock(d, 1, 2, 2), m = G.blockSolid(d, bp);
    try {
      const a = G.planFeatures(bp), b = G.features(m);
      assert.ok(Math.abs(a.volume_cm3 / b.volume_cm3 - 1) < 1e-6);
      assert.ok(Math.abs(a.height_mm - b.height_mm) < 1e-6);
      assert.ok(Math.abs(a.footprint_cm2 / b.footprint_cm2 - 1) < 1e-6);
      assert.ok(b.area_cm2 > 0);
    } finally { freeAll(m); }
  });
});

describe('bass-geometry: JS meshes against the reference STL (same measuring code)', () => {
  const entries = Object.entries(GEO.stl);
  it('the fixture set holds blocks 1x1, 2x1, 1x2, 2x2 k3 and k4, inserts and the plate', () => {
    assert.ok(entries.filter(([, e]) => e.kind === 'block').length >= 5);
    assert.ok(entries.some(([, e]) => e.kind === 'insert') && entries.some(([, e]) => e.kind === 'plate'));
  });
  for (const [file, e] of entries) {
    it(`${file}`, () => {
      const ref = stlOf(file), refAn = analyze(ref), refBB = bbox(ref), d = designOf(e.f_mode);
      let m, bp, size;
      if (e.kind === 'block') { bp = G.planBlock(d, e.ny, e.nz, e.k); m = G.blockSolid(d, bp); size = bp.size_mm; }
      else if (e.kind === 'insert') { bp = G.planBlock(d, 1, 2, 2); m = G.insertSolid(d, bp, e.k); }
      else {
        bp = G.planBlock(d, 1, 2, 2);
        m = G.platesSolid(d, bp, { rotated: false, items: G.platePlan(d, bp).map((it) => ({ k: it.k, p: it.p, x0: it.x0, y0: 0 })) });
      }
      try {
        const mesh = meshFromManifold(m), an = analyze(mesh), bb = bbox(mesh), r = new Report(file);
        r.check('volume (fixture sidecar vs STL)', e.volume_mm3, refAn.volume, 1e-6);
        r.check('volume js vs reference STL', an.volume, refAn.volume, TOL_VOL, 0);
        for (let a = 0; a < 3; a++) { r.check(`bbox min ${a}`, bb.min[a], refBB.min[a], 0, TOL_BB); r.check(`bbox max ${a}`, bb.max[a], refBB.max[a], 0, TOL_BB); }
        r.equal('outer solids', an.outer, refAn.outer);
        r.equal('inner cavities', an.inner, refAn.inner);
        r.equal('open shells', an.open, refAn.open);
        r.equal('non-manifold edges', an.nonManifoldEdgeUses, refAn.nonManifoldEdgeUses);
        refAn.innerVolumes.forEach((v, i) => r.check(`cavity ${i}`, an.innerVolumes[i], v, TOL_VOL, 0));
        if (e.kind === 'block') {
          r.equal('cavities = pockets of the plan', an.inner, e.pockets);
          const ext = [bb.max[0] - bb.min[0], bb.max[1] - bb.min[1], bb.max[2] - bb.min[2]];
          const h = maxDiff(topGrid(mesh, ext), topGrid(ref, ext));
          if (!(h <= TOL_H)) r.fails.push(`section / heights by ray cast: worst ${h} mm`);
        } else {
          // heights of the flange and body tops along a grid from above
          const ext = [bb.max[0], bb.max[1], bb.max[2]];
          const h = maxDiff(topGrid(mesh, ext, 1.1, 0.23), topGrid(ref, ext, 1.1, 0.23));
          if (!(h <= TOL_H)) r.fails.push(`heights by ray cast: worst ${h} mm`);
        }
        r.done();
      } finally { freeAll(m); }
    });
  }
  it('REPORT: w = 6 mm blocks for 42.9 and 95.2 Hz against the research STL (volume of the float32 mesh, mm3)', (t) => {
    const lines = [];
    for (const [file, rs] of Object.entries(RESEARCH)) {
      const e = GEO.stl[file], d = designOf(e.f_mode);
      let m;
      if (e.kind === 'block') m = G.blockSolid(d, G.planBlock(d, e.ny, e.nz, e.k));
      else { const bp = G.planBlock(d, 1, 2, 2); m = G.platesSolid(d, bp, { rotated: false, items: G.platePlan(d, bp).map((it) => ({ k: it.k, p: it.p, x0: it.x0, y0: 0 })) }); }
      try {
        const v = analyze(meshFromManifold(m)).volume;
        lines.push(`${rs.file}: js ${v.toFixed(2)}, research STL ${rs.volume.toFixed(2)}, regenerated STL ${e.volume_mm3.toFixed(2)} (js/research - 1 = ${((v / rs.volume - 1) * 100).toExponential(2)} %)`);
        assert.ok(Math.abs(v / rs.volume - 1) <= TOL_VOL, rs.file);
        assert.ok(Math.abs(e.volume_mm3 / rs.volume - 1) < 1e-9, 'the regenerated reference STL equals the research STL');
      } finally { freeAll(m); }
    }
    t.diagnostic('\n' + lines.join('\n'));
  });
});

describe('bass-geometry: inserts and the plate', () => {
  it('plate plan against the Python twin (positions, flange and body dimensions)', () => {
    const r = new Report('plate plans');
    for (const [key, items] of Object.entries(GEO.plates)) {
      const d = designOf(Number(key)), bp = G.planBlock(d, 1, 1, 1), pp = G.platePlan(d, bp);
      r.equal(`${key} count`, pp.length, items.length);
      items.forEach((fx, i) => {
        r.equal(`${key}[${i}] k`, pp[i].k, fx.k);
        r.check(`${key}[${i}] x0`, pp[i].x0, fx.x0, 1e-9);
        r.check(`${key}[${i}] p`, pp[i].p, fx.p, 1e-9);
        fx.flange.forEach((v, a) => r.check(`${key}[${i}] flange ${a}`, pp[i].dims.flange[a], v, 1e-9));
        fx.body.forEach((v, a) => r.check(`${key}[${i}] body ${a}`, pp[i].dims.body[a], v, 1e-9));
        fx.body_origin.forEach((v, a) => r.check(`${key}[${i}] origin ${a}`, pp[i].dims.body_origin[a], v, 1e-9));
      });
    }
    r.done();
  });
  it('file names in the style of the diffuser generator', () => {
    const d = designOf(95.2), bp = G.planBlock(d, 1, 2, 2);
    assert.equal(G.stem(d), 'soundblooms-bass-trap-qwr-f095.2-w6');
    assert.equal(G.blockName(d, bp), 'soundblooms-bass-trap-qwr-f095.2-w6-block-1x2-k2-58.8x120x240');
    assert.equal(G.insertName(d, -2), 'soundblooms-bass-trap-qwr-f095.2-w6-insert-m2');
    assert.equal(G.insertName(d, 0), 'soundblooms-bass-trap-qwr-f095.2-w6-insert-0');
    assert.equal(G.insertName(d, 3), 'soundblooms-bass-trap-qwr-f095.2-w6-insert-p3');
    assert.equal(G.platesName(d, 0, 1), 'soundblooms-bass-trap-qwr-f095.2-w6-inserts-plate');
    assert.equal(G.platesName(d, 1, 3), 'soundblooms-bass-trap-qwr-f095.2-w6-inserts-plate-2of3');
    assert.equal(G.stem(designOf(42.9)), 'soundblooms-bass-trap-qwr-f042.9-w6');
  });
  it('the plates of a P1S: one plate with all 11 inserts, 197.2 mm long; the nominal insert is the end cap', () => {
    const d = designOf(95.2), bp = G.planBlock(d, 1, 1, 1), res = G.packPlates(d, bp, G.BEDS[2], { caps: 1, tuning: true });
    assert.equal(res.ok, true);
    assert.equal(res.plates.length, 1);
    assert.equal(res.plates[0].items.length, 11);
    assert.ok(Math.abs(res.plates[0].size_mm[0] - 197.2) < 1e-6 && Math.abs(res.plates[0].size_mm[1] - 119) < 1e-6);
    assert.ok(res.plates[0].minutes > 120 && res.plates[0].minutes < 260, `${res.plates[0].minutes} min (research: 3 h 02)`);
    const caps = G.packPlates(d, bp, G.BEDS[2], { caps: 1, tuning: false });
    assert.equal(caps.plates[0].items.length, 1);
    assert.equal(caps.plates[0].items[0].k, 0);
  });
  it('plates fit every bed, split where needed, keep every insert, and turn 90 degrees when that holds more', () => {
    for (const bed of [...G.BEDS, { name: 'Custom', x: 150, y: 220, z: 130 }]) for (const f of [42.9, 95.2]) {
      const d = designOf(f), bp = G.planBlock(d, 1, 1, 1);
      for (const caps of [1, 8, 16]) {
        const res = G.packPlates(d, bp, bed, { caps, tuning: true });
        assert.equal(res.ok, true, `${bed.name} ${f} ${res.reasons}`);
        assert.equal(res.plates.reduce((s, p) => s + p.items.length, 0), 10 + caps, `${bed.name} caps ${caps}`);
        res.plates.forEach((pl, i) => {
          assert.ok(G.fitsBed(pl.size_mm, bed), `${bed.name} ${f} caps ${caps} plate ${i}: ${pl.size_mm}`);
          const m = G.platesSolid(d, bp, pl);
          try {
            const bb = bbox(meshFromManifold(m));
            const ext = [bb.max[0] - bb.min[0], bb.max[1] - bb.min[1], bb.max[2] - bb.min[2]];
            assert.ok(Math.abs(bb.min[0]) < 1e-3 && Math.abs(bb.min[1]) < 1e-3, 'plate starts at the bed corner');
            pl.size_mm.forEach((v, a) => assert.ok(Math.abs(ext[a] - v) < 0.01, `${bed.name} plate ${i} axis ${a}: ${ext[a]} vs ${v}`));
          } finally { freeAll(m); }
        });
      }
    }
    const mk4 = G.packPlates(designOf(95.2), G.planBlock(designOf(95.2), 1, 1, 1), G.BEDS[1], { caps: 16, tuning: true });
    assert.equal(mk4.plates[0].rotated, true);
  });
  it('an insert that does not fit the bed gives ok: false with a reason', () => {
    const d = designOf(42.9), bp = G.planBlock(d, 1, 1, 1);
    const low = G.packPlates(d, bp, { name: 'Custom', x: 250, y: 250, z: 80 }, { caps: 1, tuning: true });
    assert.equal(low.ok, false);
    assert.ok(low.reasons.length >= 1);
    const narrow = G.packPlates(d, bp, { name: 'Custom', x: 100, y: 100, z: 250 }, { caps: 1, tuning: true });
    assert.equal(narrow.ok, false);
  });
});

describe('bass-geometry: packing blocks on a printer bed', () => {
  const beds = G.BEDS;
  const freqs = [35, 42.9, 53.6, 68.7, 77.3, 85.8, 95.2, 107.2, 120];
  it('every block fits the bed (5 mm margin per axis), the cells add up, estimates are positive', () => {
    for (const bed of beds) for (const f of freqs) for (const N of [1, 2, 3, 4, 5, 8, 16]) {
      const d = designOf(f), res = G.packBlocks(d, N, bed);
      if (!res.ok) { assert.ok(res.reasons.length > 0); continue; }
      assert.equal(res.blocks.reduce((s, b) => s + b.k * b.count, 0), N, `${bed.key} ${f} Hz N=${N}`);
      for (const b of res.blocks) {
        assert.ok(G.fitsBed(b.size_mm, bed), `${bed.key} ${f} ${b.ny}x${b.nz}: ${b.size_mm}`);
        assert.ok(b.size_mm[2] <= bed.z - 5 && Math.min(...b.size_mm.slice(0, 2)) <= Math.min(bed.x, bed.y) - 5);
        assert.ok(b.minutes > 0 && b.grams > 0 && b.count >= 1);
        if (!b.overLimit) assert.ok(b.hours <= 12);
      }
      assert.equal(res.prints, res.blocks.reduce((s, b) => s + b.count, 0));
    }
  });
  it('prefers a block of up to 12 h; a layout over 12 h only when nothing quicker fits (and then says so)', () => {
    for (const bed of beds) for (const f of freqs) {
      const d = designOf(f), res = G.packBlocks(d, 16, bed);
      if (!res.ok) continue;
      const first = res.blocks[0];
      if (first.overLimit) {
        // every layout that fits the bed and has at most 16 cells needs more than 12 h
        for (const c of G.BLOCK_ORDER) {
          const size = [d.depth, c.ny * 120, c.nz * 120];
          if (!G.fitsBed(size, bed)) continue;
          const est = G.estimate('block', { volume_cm3: G.planBlock(d, 1, 1, 1).volume_mm3 * c.k / 1000, footprint_cm2: size[0] * size[1] / 100 });
          assert.ok(est.minutes / 60 > 12 || c.k === first.k, `${bed.key} ${f} Hz ${c.ny}x${c.nz} would have been quicker`);
        }
        assert.ok(res.notes.length >= 1);
      }
    }
  });
  it('known cases of the research set on a P1S: 42.9 Hz one cell per block, 95.2 Hz two cells per block, 107.2 Hz an L of three', () => {
    const p1s = beds[2];
    const a = G.packBlocks(designOf(42.9), 8, p1s);
    assert.deepEqual(plain(a.blocks.map((b) => [b.ny, b.nz, b.k, b.count])), [[1, 1, 1, 8]]);
    assert.ok(Math.abs(a.blocks[0].hours - 8.4) < 0.2, `${a.blocks[0].hours} h (research slice: 8 h 24)`);
    const b = G.packBlocks(designOf(95.2), 8, p1s);
    assert.deepEqual(plain(b.blocks.map((x) => [x.ny, x.nz, x.k, x.count])), [[1, 2, 2, 4]]);
    assert.ok(Math.abs(b.blocks[0].hours - 8.28) < 0.2, `${b.blocks[0].hours} h (research slice: 8 h 17)`);
    const c = G.packBlocks(designOf(107.2), 6, p1s);
    assert.deepEqual(plain(c.blocks.map((x) => [x.ny, x.nz, x.k, x.count])), [[2, 2, 3, 2]]);
    const odd = G.packBlocks(designOf(95.2), 5, p1s);
    assert.deepEqual(plain(odd.blocks.map((x) => [x.ny, x.nz, x.k, x.count])), [[1, 2, 2, 2], [1, 1, 1, 1]]);
    assert.equal(a.blocks[0].name, 'soundblooms-bass-trap-qwr-f042.9-w6-block-1x1-k1-123.6x120x120');
  });
  it('beds without room for a stack or a second column fall back to what fits', () => {
    const d = designOf(95.2);
    assert.deepEqual(plain(G.packBlocks(d, 8, beds[0]).blocks.map((x) => [x.ny, x.nz, x.k, x.count])), [[1, 1, 1, 8]]);     // A1 mini: 175 mm
    assert.deepEqual(plain(G.packBlocks(d, 8, beds[1]).blocks.map((x) => [x.ny, x.nz, x.k, x.count])), [[2, 1, 2, 4]]);     // MK4: 215 mm high, 240 wide fits across 245
  });
  it('"does not fit" is a state with reasons, never an exception', () => {
    const deep = G.packBlocks(designOf(35), 1, beds[2]);               // the 35 Hz cell is deeper than 150 mm
    assert.equal(deep.ok, false);
    assert.match(deep.reasons.join(' '), /150/);
    const low = G.packBlocks(designOf(95.2), 1, { name: 'Custom', x: 300, y: 300, z: 110 });
    assert.equal(low.ok, false);
    assert.match(low.reasons.join(' '), /too low/);
    const small = G.packBlocks(designOf(95.2), 1, { name: 'Custom', x: 60, y: 300, z: 300 });
    assert.equal(small.ok, false);
    assert.match(small.reasons.join(' '), /does not fit the bed/);
    for (const bad of [0, -1, 1.5, NaN, undefined]) assert.equal(G.packBlocks(designOf(95.2), bad, beds[2]).ok, false);
    assert.equal(G.packBlocks(G.design(500), 1, beds[2]).ok, false);
    assert.equal(G.packBlocks(null, 1, beds[2]).ok, false);
    // a custom bed that does hold the deep cell turns the state around
    assert.equal(G.packBlocks(designOf(35), 1, { name: 'Custom', x: 300, y: 300, z: 300 }).ok, false);   // still deeper than max_depth
    assert.equal(G.packBlocks(G.design(35, { max_depth: 400 }), 1, { name: 'Custom', x: 300, y: 300, z: 300 }).ok, true);
  });
});

describe('bass-geometry: print time and filament estimate', () => {
  const data = readJson(path.join(GEN, 'scripts', 'calibration', 'bass-data.json')), report = readJson(path.join(GEN, 'scripts', 'calibration', 'bass-fit-report.json'));
  it('reproduces the research slices within the training errors of the fit (blocks from plan features, plates from solid features)', () => {
    for (const kind of ['block', 'plate']) for (const q of ['minutes', 'grams']) {
      const rows = data.filter((r) => r.kind === kind), lim = report[kind][q].train_max_pct;
      let worst = 0;
      for (const r of rows) worst = Math.max(worst, Math.abs(G.estimate(kind, r.features)[q] / r[q] - 1) * 100);
      assert.ok(worst <= lim + 0.01,   // the report rounds to 0.01 %
         `${kind} ${q}: worst ${worst.toFixed(2)} % above the fit report ${lim} %`);
      assert.ok(lim < 3, `${kind} ${q}: training error ${lim} %`);
    }
    assert.equal(data.filter((r) => r.kind === 'block').length, 10);
    assert.equal(data.filter((r) => r.kind === 'plate').length, 10);
  });
  it('the rows of the fit are the research slices (names, v1 and v2)', () => {
    const ids = data.map((r) => r.id);
    assert.ok(ids.includes('qwr2-m100-f042.9-w6-block') && ids.includes('qwr-f107.2-inserts-plate'));
    assert.equal(data.find((r) => r.id === 'qwr2-m100-f042.9-w6-block').minutes, 504.517);
  });
  it('a block of a mode designed by the core is estimated near the slice of the same cell (42.9 Hz: 8 h 24, 425 g; 95.2 Hz 1x2: 8 h 17, 418 g)', () => {
    const a = G.packBlocks(designOf(42.9), 1, G.BEDS[2]).blocks[0], b = G.packBlocks(designOf(95.2), 2, G.BEDS[2]).blocks[0];
    assert.ok(Math.abs(a.minutes / 504.517 - 1) < 0.02 && Math.abs(a.grams / 424.67 - 1) < 0.02, `${a.minutes} min ${a.grams} g`);
    assert.ok(Math.abs(b.minutes / 497.5 - 1) < 0.02 && Math.abs(b.grams / 417.54 - 1) < 0.02, `${b.minutes} min ${b.grams} g`);
  });
  it('estimate() handles a missing feature as zero and never returns NaN for a block', () => {
    const e = G.estimate('block', { volume_cm3: 100 });
    assert.ok(Number.isFinite(e.minutes) && Number.isFinite(e.grams) && e.minutes > 0);
    for (const bed of G.BEDS) for (const f of [42.9, 95.2]) {
      const r = G.packBlocks(designOf(f), 4, bed);
      assert.ok(r.ok && Number.isFinite(r.totalMinutes) && Number.isFinite(r.totalGrams));
    }
  });
});
