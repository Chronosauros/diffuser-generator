// Bass traps, physics core: <script id="bass-physics"> of index.html (createBassPhysics) against its Python twin
// (reference/bass_room.py, fixtures in scripts/fixtures/bass/physics.json written by scripts/bass_reference.py).
//
//   npm test            (= node --test scripts/*.test.mjs)
//
// Tolerances (spec of the bass core): physics numbers 1e-6 relative, roots of the Newton polish 1e-9 relative. The fixtures hold 10
// significant digits, so 1e-9 is the honest floor. Everything is calculated, not measured: the FEM figures of the research
// package appear below as REPORTS (diagnostics), not as thresholds.

import { describe, it } from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import path from 'node:path';
import { FIX, readJson, loadBass, plain } from './mesh-tools.mjs';

const PHYS = readJson(path.join(FIX, 'physics.json'));
const { physics: P, sources } = await loadBass();

const TOL = 1e-6;          // relative, physics numbers
const TOL_ROOT = 1e-9;     // relative, roots of the single-mode equation and L from the bisection
const ABS = 1e-9;          // absolute floor for numbers that are zero in exact arithmetic (dB of a nodal anchor)
const close = (a, b, rel = TOL, abs = ABS) => {
  if (a === null || b === null || a === undefined || b === undefined) return a === b;
  return Math.abs(a - b) <= abs + rel * Math.max(Math.abs(a), Math.abs(b));
};
class Report {
  constructor(label) { this.label = label; this.fails = []; this.n = 0; }
  check(name, a, b, rel = TOL, abs = ABS) { this.n++; if (!close(a, b, rel, abs)) this.fails.push(`${name}: js ${a} python ${b}`); }
  equal(name, a, b) { this.n++; if (JSON.stringify(a) !== JSON.stringify(b)) this.fails.push(`${name}: js ${JSON.stringify(a)} python ${JSON.stringify(b)}`); }
  done() { if (this.fails.length) assert.fail(`${this.label}: ${this.fails.length} of ${this.n} differ\n  - ${this.fails.slice(0, 12).join('\n  - ')}`); }
}

describe('bass-physics: the script is a classic, self-contained core', () => {
  it('has no import, DOM, fetch or manifold, and defines createBassPhysics()', () => {
    assert.equal(typeof P.analyzeRoom, 'function');
    assert.doesNotMatch(sources.physics.replace(/\/\*[\s\S]*?\*\//g, '').replace(/\/\/.*$/gm, ''), /\bimport\b|\brequire\s*\(|\bwindow\b|\bdocument\b|\bfetch\b|\bManifold\b|\bCrossSection\b/);
    assert.match(sources.physics, /function createBassPhysics\(\)/);
  });
  it('exposes the documented members', () => {
    for (const k of ['analyzeRoom', 'evaluateMode', 'roomModes', 'baseDecay', 'sabineDecay', 't60', 'cluster', 'cellsPsi', 'loadingMass',
      'solveL', 'firstRes', 'zsQwr', 'makeCell', 'cellForMode', 'solveSingle', 'singleModeEffect', 'peakDropDb', 'overlapFlag', 'orient',
      'approximateReasons', 'uncertainReasons', 'anchorText', 'absorption', 'dlSlitArray']) assert.equal(typeof P[k], 'function', k);
    assert.deepEqual(plain(P.ANCHORS.map((a) => a.id)), PHYS.anchors);
    assert.deepEqual(plain(P.N_LIST), PHYS.ns);
    assert.equal(P.PEAK_POINTS, PHYS.peak.points);
    assert.deepEqual(plain(P.PEAK_SPAN), PHYS.peak.span);
  });
});

describe('bass-physics: building blocks against the Python twin', () => {
  const pr = PHYS.primitives;
  it('constants of the research package', () => {
    assert.equal(P.C0, PHYS.constants.C0);
    assert.equal(P.RHO0, PHYS.constants.RHO0);
    assert.ok(close(P.Z0, PHYS.constants.RHO0 * PHYS.constants.C0, 1e-12));
  });
  it('end correction of the slit array', () => {
    const r = new Report('dl_slit_array');
    r.check('6 mm', P.dlSlitArray(6e-3, 0.12), pr.dl_slit_array['6'], 1e-9);
    r.check('8 mm', P.dlSlitArray(8e-3, 0.12), pr.dl_slit_array['8'], 1e-9);
    r.done();
  });
  it('cell impedance at complex angular frequency', () => {
    const r = new Report('zs_qwr');
    for (const c of pr.zs_qwr) {
      const z = P.zsQwr(2 * Math.PI * c.f, c.om_im, c.w, c.L, 120, 1);
      r.check(`f ${c.f} + j${c.om_im} Re`, z[0], c.zs[0], TOL);
      r.check(`f ${c.f} + j${c.om_im} Im`, z[1], c.zs[1], TOL);
    }
    r.done();
  });
  it('first_res and solve_L on the 0.25 Hz grid (the same numbers as Python)', () => {
    const r = new Report('first_res / solve_L');
    for (const [l, f] of Object.entries(pr.first_res)) r.check(`first_res(6, ${l})`, P.firstRes(6, Number(l)), f, TOL_ROOT);
    for (const [f, L] of Object.entries(pr.solve_L)) r.check(`solve_L(${f}, 6)`, P.solveL(Number(f), 6), L, TOL_ROOT);
    r.done();
  });
  it('first_res is NaN outside the grid and L reproduces the frequency', () => {
    assert.ok(Number.isNaN(P.firstRes(6, 20)));            // a 20 mm channel resonates above 520 Hz
    const L = P.solveL(60, 6), f = P.firstRes(6, L);
    assert.ok(Math.abs(f - 60) < 0.05, `first_res(solve_L(60)) = ${f}`);
  });
  it('absorption of the cell, radiation mass, Sabine decay', () => {
    const r = new Report('misc');
    const cell = P.makeCell(43.1145, 6, 1908.7908590295256, 0);
    for (const [f, a] of Object.entries(pr.absorption)) r.check(`alpha(${f})`, cell.alphaAt(Number(f)), a, TOL);
    for (const [an, byN] of Object.entries(pr.loading_mass)) for (const [n, m] of Object.entries(byN)) r.check(`mass ${an} N=${n}`, P.loadingMass(Number(n), an), m, 1e-12);
    for (const [a, d] of Object.entries(pr.sabine_decay)) r.check(`sabine ${a}`, P.sabineDecay([4, 3.2, 2.6], Number(a)), d, 1e-12);
    r.done();
  });
  it('cell clusters (order of the cells) and mode shapes on the cells', () => {
    const L = [4, 3.2, 2.6], r = new Report('cluster');
    for (const [an, cells] of Object.entries(pr.cluster_example)) {
      const js = P.cluster(L, an, 20);
      r.equal(`${an} count`, js.length, cells.length);
      cells.forEach((c, i) => {
        r.equal(`${an}[${i}] wall`, js[i].wall, c.wall);
        c.center.forEach((v, a) => r.check(`${an}[${i}] centre ${a}`, js[i].center[a], v, 1e-12));
      });
    }
    const modes = P.roomModes(L, 120).slice(0, 4);
    for (const [an, psi] of Object.entries(pr.cells_psi_example)) {
      const js = P.cellsPsi(modes, L, P.cluster(L, an, 3));
      psi.forEach((row, s) => row.forEach((v, m) => r.check(`psi ${an} cell ${s} mode ${m}`, js[s][m], v, 1e-9, 1e-12)));
    }
    r.done();
  });
  it('anchor descriptions in English', () => {
    for (const [an, text] of Object.entries(pr.anchor_text_example)) assert.equal(P.anchorText(an, [4, 3.2, 2.6]), text);
    assert.equal(P.anchorText('corner', [4, 3.2, 2.6]), 'floor corner (the three walls meet)');
    assert.equal(P.anchorText('edge-z', [4, 3.2, 2.6]), 'vertical corner edge, half height (2.6 m)');
    assert.equal(P.anchorText('wall-x', [4, 3.2, 2.6]), 'middle of the 3.2 x 2.6 m wall');
  });
});

describe('bass-physics: the example room 4.0 x 3.2 x 2.6 m, alpha 0.10 (numbers of the research package)', () => {
  const room = P.analyzeRoom([4.0, 3.2, 2.6], 0.10, 120);
  it('mode (1,0,0) 42.9 Hz: base decay 5.136 1/s, T60 1.345 s', () => {
    const m = room.rows[0];
    assert.deepEqual(plain(m.n), [1, 0, 0]);
    assert.equal(m.kind, 'axial');
    assert.ok(Math.abs(m.f - 42.9) < 1e-9);
    assert.ok(Math.abs(m.deltaBase - 5.1356) < 1e-4, `delta ${m.deltaBase}`);          // REVIEW.md: 5.1356 (MES 5.1364)
    assert.ok(Math.abs(m.t60Before - 1.345) < 1e-3, `T60 ${m.t60Before}`);
    assert.equal(room.rows.length, 12);
  });
  it('REPORT: T60 before / after, decay factor and peak drop for N = 1 and 8 at the best wall position (not thresholds)', (t) => {
    const ev = P.evaluateMode(room, 0);
    const lines = [];
    for (const N of [1, 8]) {
      const an = ev.best[N], r = ev.anchors[an].byN[N];
      lines.push(`N = ${N}: best "${ev.anchors[an].text}" (${an}): T60 ${ev.t60Before.toFixed(3)} -> ${r.t60After.toFixed(3)} s, decay x${r.decayFactor.toFixed(3)}, peak ${r.peakDropDb.toFixed(2)} dB, approximate ${r.approximate}`);
      for (const id of Object.keys(ev.anchors)) {
        const q = ev.anchors[id].byN[N];
        lines.push(`    ${id.padEnd(7)} x${q.decayFactor.toFixed(3)}  ${q.peakDropDb.toFixed(2)} dB  nodal ${ev.anchors[id].nodal}`);
      }
    }
    lines.push('FEM (research package, alpha 0.10): 8 cells x1.37 - 1.50 and -1.3 to -1.7 dB; one cell T60 1.35 -> 1.30 s, -0.2 dB');
    t.diagnostic('\n' + lines.join('\n'));
    // plausibility only: a trap helps at the best place and never raises the peak there
    const b1 = ev.anchors[ev.best[1]].byN[1], b8 = ev.anchors[ev.best[8]].byN[8];
    assert.ok(b1.decayFactor > 1 && b1.decayFactor < 1.1);
    assert.ok(b8.decayFactor > 1.2 && b8.decayFactor < 2.0);
    assert.ok(b1.peakDropDb < 0 && b8.peakDropDb < b1.peakDropDb);
  });
  it('8 cells at 42.9 Hz stay in the neighbourhood of the FEM results (a loose sanity band, the report above has the numbers)', () => {
    const ev = P.evaluateMode(room, 0), corner = ev.anchors.corner.byN[8];
    assert.ok(corner.peakDropDb < -1.0 && corner.peakDropDb > -2.0, `${corner.peakDropDb}`);
    assert.ok(corner.decayFactor > 1.3 && corner.decayFactor < 1.6, `${corner.decayFactor}`);
  });
});

describe('bass-physics: every fixture case against the Python twin', () => {
  for (const room of PHYS.rooms) for (const cs of room.cases) {
    it(`${room.id} ${room.L.join(' x ')} m, alpha ${cs.alpha}: ${cs.modes.length} modes x 7 anchors x N ${PHYS.ns.join('/')}`, () => {
      const rm = P.analyzeRoom(room.L, cs.alpha, PHYS.fmax), r = new Report(`${room.id} alpha ${cs.alpha}`);
      r.equal('mode count', rm.rows.length, cs.modes.length);
      r.check('sabine T60', P.t60(P.sabineDecay(room.L, cs.alpha)), cs.sabine_t60, 1e-12);
      cs.modes.forEach((fx, i) => {
        const tag = `mode ${fx.n.join('')}`, row = rm.rows[i];
        r.equal(`${tag} n`, row.n, fx.n);
        r.equal(`${tag} kind`, row.kind, fx.kind);
        r.check(`${tag} f`, row.f, fx.f, 1e-9);
        r.check(`${tag} Lambda`, row.Lam, fx.Lam, 1e-9);
        r.check(`${tag} delta_base`, row.deltaBase, fx.delta_base, 1e-9);
        r.check(`${tag} T60 before`, row.t60Before, fx.t60_before, 1e-9);
        r.check(`${tag} spacing`, row.spacingHz, fx.spacing_Hz, 1e-9, 1e-9);
        const ev = P.evaluateMode(rm, i);
        r.check(`${tag} L0`, ev.L0, fx.L0, TOL_ROOT);
        r.check(`${tag} f cell`, ev.fCell, fx.f_cell, 1e-9);
        r.check(`${tag} gamma`, ev.gammaCell, fx.gamma_cell, TOL);
        r.equal(`${tag} no effect`, ev.noEffect, fx.no_effect);
        r.equal(`${tag} overlap`, ev.overlap, fx.overlap);
        r.equal(`${tag} uncertain`, ev.uncertain, fx.uncertain);
        r.equal(`${tag} best`, PHYS.ns.map((n) => ev.best[n]), fx.best);
        P.ANCHORS.forEach((a, ai) => {
          r.check(`${tag} ${a.id} psi2`, ev.anchors[a.id].psi2, fx.psi2[ai], 1e-9, 1e-12);
          r.equal(`${tag} ${a.id} nodal`, ev.anchors[a.id].nodal, fx.nodal[ai]);
          PHYS.ns.forEach((n, ni) => {
            const q = ev.anchors[a.id].byN[n];
            r.check(`${tag} ${a.id} N=${n} delta_eff`, q.deltaEff, fx.delta_eff[ai][ni], TOL_ROOT);
            r.check(`${tag} ${a.id} N=${n} peak dB`, q.peakDropDb, fx.peak_drop_db[ai][ni], TOL);
            if (q.deltaEff !== null) r.check(`${tag} ${a.id} N=${n} T60 after`, q.t60After, 6.907755278982137 / fx.delta_eff[ai][ni], 1e-8);
          });
        });
      });
      r.done();
    });
  }
});

describe('bass-physics: flags', () => {
  const find = (id, alpha) => PHYS.rooms.find((r) => r.id === id).cases.find((c) => c.alpha === alpha);
  it('"no effect" = base decay above the cell width gamma; the trap then does not shorten the decay', () => {
    let seen = 0;
    for (const room of PHYS.rooms) for (const cs of room.cases) for (const m of cs.modes) {
      assert.equal(m.no_effect, m.delta_base > m.gamma_cell, `${room.id} ${m.n}`);
      if (m.no_effect) {
        seen++;
        const best = Math.max(...m.delta_eff.flat());
        assert.ok(best <= m.delta_base, `${room.id} alpha ${cs.alpha} ${m.n}: delta_eff ${best} above the base ${m.delta_base}`);
      }
    }
    assert.ok(seen >= 2, 'the fixtures hold the known cases (oblique 114 Hz of the small room and 104 Hz of the near-cube room at alpha 0.20)');
    const small = find('small', 0.2).modes.find((m) => m.n.join('') === '111');
    assert.equal(small.no_effect, true);
    assert.ok(Math.min(...small.peak_drop_db.flat()) < -1, 'the peak still drops (REVIEW.md 3): the decay does not get shorter, the peak does');
  });
  it('alpha 0.20 at 42.9 Hz: 8 cells in the corner give no real gain (the fixture decides, not a guess)', () => {
    const m = find('example', 0.2).modes[0], corner = 0, n8 = PHYS.ns.indexOf(8);
    const factor = m.delta_eff[corner][n8] / m.delta_base;
    assert.ok(factor < 1.1, `decay factor ${factor}`);
    assert.ok(m.uncertain.includes('alpha >= 0.15'));
    assert.equal(m.no_effect, false);                       // delta 10.27 is below gamma 11.66, but above 0.8 gamma
    assert.ok(m.uncertain.includes('base decay above 0.8 of the cell width'));
  });
  it('nodal: the trap sits in a node of the mode', () => {
    const m = find('example', 0.1).modes[0];               // (1,0,0): x = Lx/2 is a node, so the edge-x anchor is nodal
    assert.equal(m.nodal[P.ANCHORS.findIndex((a) => a.id === 'edge-x')], true);
    assert.equal(m.nodal[P.ANCHORS.findIndex((a) => a.id === 'corner')], false);
    assert.ok(Math.max(...m.delta_eff[P.ANCHORS.findIndex((a) => a.id === 'edge-x')]) / m.delta_base < 1.01);
  });
  it('approximate: overlapping modes, more than 8 cells, tangential / oblique with 4 or more cells', () => {
    assert.equal(P.orient('axial', false, 8), false);
    assert.equal(P.orient('axial', false, 16), true);
    assert.equal(P.orient('axial', true, 1), true);
    assert.equal(P.orient('tangential', false, 2), false);
    assert.equal(P.orient('tangential', false, 4), true);
    assert.equal(P.orient('oblique', false, 4), true);
    assert.deepEqual(plain(P.approximateReasons('axial', false, 16, 'wall-x', false)), ['more than 8 cells']);
    assert.deepEqual(plain(P.approximateReasons('axial', false, 1, 'corner', true)), ['corner with a neighbouring mode within 10 % of the cell frequency']);
    const room = P.analyzeRoom([4.0, 3.2, 2.6], 0.10, 120), ev = P.evaluateMode(room, 0);
    assert.equal(ev.anchors.corner.byN[8].approximate, false);
    assert.equal(ev.anchors.corner.byN[16].approximate, true);
    const t = P.evaluateMode(room, room.rows.findIndex((r) => r.kind === 'tangential' && !r.overlap));
    assert.equal(t.anchors['wall-x'].byN[2].approximate, false);
    assert.equal(t.anchors['wall-x'].byN[4].approximate, true);
  });
  it('overlap: the near-cube room has overlapping modes and they are flagged approximate at every N', () => {
    const rm = P.analyzeRoom([3.0, 2.9, 2.7], 0.2, 120), idx = rm.rows.findIndex((r) => r.overlap);
    assert.ok(idx >= 0);
    const ev = P.evaluateMode(rm, idx);
    for (const id of Object.keys(ev.anchors)) for (const n of PHYS.ns) assert.equal(ev.anchors[id].byN[n].approximate, true);
    assert.equal(PHYS.rooms.flatMap((r) => r.cases.flatMap((c) => c.modes.filter((m) => m.overlap))).length, 145);
  });
});

describe('bass-physics: peak drop and model behaviour', () => {
  const room = P.analyzeRoom([4.0, 3.2, 2.6], 0.10, 120), mode = room.modes[0];
  const cell = P.cellForMode(mode.f);
  it('no coupling (G = 0, same decay) is 0 dB', () => {
    const d = P.baseDecay(mode, room.L, 0.10);
    assert.ok(Math.abs(P.peakDropDb(mode, d, d, 0, cell)) < 1e-9);
  });
  it('the peak drop is the ratio of the two peaks of the formula in the spec (independent evaluation)', () => {
    const d = P.baseDecay(mode, room.L, 0.10), cl = P.cluster(room.L, 'corner', 8), psi = P.cellsPsi([mode], room.L, cl).map((r) => r[0]);
    const G = 0.0144 * psi.reduce((s, p) => s + p * p, 0) / mode.Lam;
    const c8 = P.makeCell(cell.fDesign, cell.w, cell.L, P.loadingMass(8, 'corner'));
    const dEff = d - 0.5 * P.C0 * G * (0.10 / 8);
    const wn = 2 * Math.PI * mode.f;
    let mb = 0, mt = 0;
    for (let i = 0; i < 3001; i++) {
      const w = wn * (0.85 + 0.3 * i / 3000), zs = c8.zs(w, 0);
      const hb = 1 / Math.hypot(w * w - wn * wn, 2 * w * d);
      const pr = w * w - wn * wn, pi = -2 * w * dEff;                // (w^2 - wn^2 - 2 j w d')
      const dr = pr * zs[0] - pi * zs[1], di = pr * zs[1] + pi * zs[0] - w * P.C0 * P.Z0 * G;
      const ht = Math.hypot(zs[0], zs[1]) / Math.hypot(dr, di);
      mb = Math.max(mb, hb); mt = Math.max(mt, ht);
    }
    assert.ok(Math.abs(20 * Math.log10(mt / mb) - P.peakDropDb(mode, d, dEff, G, c8)) < 1e-9);
  });
  it('the decay factor grows with N at the best position, and the peak drop deepens', () => {
    const ev = P.evaluateMode(room, 0);
    let prevF = 1, prevDb = 0;
    for (const N of [1, 2, 4, 8]) {
      const r = ev.anchors[ev.best[N]].byN[N];
      assert.ok(r.decayFactor >= prevF - 1e-12 && r.peakDropDb <= prevDb + 1e-12, `N=${N}`);
      prevF = r.decayFactor; prevDb = r.peakDropDb;
    }
  });
  it('the best anchor holds the largest deltaEff of the seven, and has an English description', () => {
    const ev = P.evaluateMode(room, 0);
    for (const N of P.N_LIST) {
      const best = ev.anchors[ev.best[N]].byN[N].deltaEff;
      for (const id of Object.keys(ev.anchors)) assert.ok(ev.anchors[id].byN[N].deltaEff <= best + 1e-12);
    }
    assert.match(ev.anchors[ev.best[8]].text, /^middle of the 3\.2 x 2\.6 m wall$/);
  });
  it('never throws: too many cells for a wall, an empty mode list, a mode above the design range', () => {
    const ev = P.evaluateMode(room, 0, { ns: [1, 5000], anchors: ['wall-x'] });
    assert.equal(ev.anchors['wall-x'].byN[5000].available, false);
    assert.equal(ev.best[5000], null);
    const none = P.analyzeRoom([1.0, 1.0, 1.0], 0.1, 50);
    assert.equal(none.rows.length, 0);
    assert.doesNotThrow(() => P.evaluateMode(P.analyzeRoom([2.0, 1.9, 1.8], 0.1, 120), 0));
  });
  it('a custom wall absorption per wall (six values) is accepted and equals the scalar for equal values', () => {
    const a = P.baseDecay(mode, room.L, 0.1), b = P.baseDecay(mode, room.L, [0.1, 0.1, 0.1, 0.1, 0.1, 0.1]);
    assert.equal(a, b);
    assert.ok(P.baseDecay(mode, room.L, [0.3, 0.3, 0.1, 0.1, 0.1, 0.1]) > a);
  });
});

// ------------------------------------------------------------------ verdict(): the policy of the mode table (interface, not physics)
// No Python twin: the thresholds are the owner's policy of 25.10.2026 (x >= 1.25 worth trying, 1.08 <= x < 1.25 small effect, below that or the
// "no effect" flag: no effect). The cases below are built by hand, so the edges are exact, plus the example room as a plausibility check.
describe('bass-physics: verdict() sorts a mode into worth trying / small effect / no effect', () => {
  const rec = (byN, anchorOver = {}) => ({ text: 'x', psi2: 1, nodal: false, cellsAvailable: 16, byN, ...anchorOver });
  const one = (x, extra = {}) => ({ available: true, solved: true, decayFactor: x, approximate: false, approximateReasons: [], ...extra });
  const evOf = (anchors, over = {}) => ({ mode: {}, noEffect: false, overlap: false, uncertain: [], anchors, best: {}, ...over });
  const levelAt = (x, over) => P.verdict(evOf({ corner: rec({ 8: one(x) }) }, over)).level;

  it('thresholds are exactly 1.25 and 1.08, and are exported', () => {
    assert.equal(P.VERDICT_WORTH, 1.25);
    assert.equal(P.VERDICT_SMALL, 1.08);
    assert.equal(P.VERDICT_MAX_N, 8);
    assert.equal(levelAt(1.25), 'worth');
    assert.equal(levelAt(1.2499999), 'small');
    assert.equal(levelAt(1.08), 'small');
    assert.equal(levelAt(1.0799999), 'none');
    assert.equal(levelAt(1.0), 'none');
    assert.equal(levelAt(0.95), 'none');
    assert.equal(levelAt(2.4), 'worth');
  });
  it('the "no effect" flag of the core wins over any factor', () => {
    const v = P.verdict(evOf({ corner: rec({ 8: one(1.6) }) }, { noEffect: true }));
    assert.equal(v.level, 'none');
    assert.equal(v.label, 'No effect');
    assert.equal(v.noEffect, true);
    assert.equal(v.x, 1.6);
  });
  it('takes the best place and N up to 8 cells; more cells, unavailable and unsolved cases are ignored', () => {
    const v = P.verdict(evOf({
      corner: rec({ 1: one(1.04), 4: one(1.18), 8: one(1.30), 16: one(1.9) }),
      'wall-x': rec({ 8: one(1.41), 16: one(5.0) }),
      'edge-z': rec({ 8: { available: true, solved: false, decayFactor: null, approximate: false, approximateReasons: [] }, 4: { available: false } })
    }));
    assert.equal(v.level, 'worth');
    assert.equal(v.x, 1.41);
    assert.equal(v.N, 8);
    assert.equal(v.anchor, 'wall-x');
    assert.equal(P.verdict(evOf({ corner: rec({ 8: one(1.41), 16: one(5.0) }) }), { maxN: 4 }).level, 'none');
    assert.equal(P.verdict(evOf({ corner: rec({ 16: one(5.0) }) })).level, 'none');
  });
  it('a mode with nothing solved is "no effect" with no x, N or place; ties go to the smaller N', () => {
    const v = P.verdict(evOf({ corner: rec({ 8: { available: false } }) }));
    assert.deepEqual([v.level, v.x, v.N, v.anchor], ['none', null, null, null]);
    const t = P.verdict(evOf({ corner: rec({ 2: one(1.3), 8: one(1.3) }) }));
    assert.equal(t.N, 2);
  });
  it('markers: approximate follows the best case (with its reasons), uncertain follows the reasons of the core', () => {
    const a = P.verdict(evOf({ corner: rec({ 4: one(1.3, { approximate: true, approximateReasons: ['overlapping modes'] }), 8: one(1.2) }) }, { uncertain: ['alpha >= 0.15'] }));
    assert.equal(a.approximate, true);
    assert.deepEqual(a.approximateReasons, ['overlapping modes']);
    assert.equal(a.uncertain, true);
    assert.deepEqual(a.uncertainReasons, ['alpha >= 0.15']);
    const b = P.verdict(evOf({ corner: rec({ 8: one(1.3) }) }));
    assert.equal(b.approximate, false);
    assert.equal(b.uncertain, false);
    assert.deepEqual(b.uncertainReasons, []);
    assert.equal(P.verdict(evOf({ corner: rec({ 8: { available: false } }) }, { overlap: true })).approximate, true);
  });
  it('example room 4.0 x 3.2 x 2.6 m: 42.9 Hz is worth trying at alpha 0.05 and 0.10, not at 0.20 (flagged uncertain)', () => {
    const at = (alpha, i) => {
      const room = P.analyzeRoom([4.0, 3.2, 2.6], alpha, 120);
      return P.verdict(P.evaluateMode(room, i, { ns: [1, 2, 4, 8] }));
    };
    const v10 = at(0.10, 0);
    assert.equal(v10.level, 'worth');
    assert.equal(v10.N, 8);
    assert.equal(v10.anchor, 'wall-x');
    assert.ok(Math.abs(v10.x - 1.55) < 0.01, `x ${v10.x}`);
    assert.equal(v10.approximate, false);
    assert.equal(v10.uncertain, false);
    assert.equal(at(0.05, 0).level, 'worth');
    const v20 = at(0.20, 0);
    assert.equal(v20.level, 'none');
    assert.equal(v20.uncertain, true);
    assert.ok(v20.x < P.VERDICT_SMALL);
  });
  it('is a pure function of the evaluation: the same input gives the same answer and the input is not changed', () => {
    const room = P.analyzeRoom([4.0, 3.2, 2.6], 0.10, 120), ev = P.evaluateMode(room, 3, { ns: [1, 2, 4, 8] });
    const before = JSON.stringify(ev);
    assert.deepEqual(plain(P.verdict(ev)), plain(P.verdict(ev)));
    assert.equal(JSON.stringify(ev), before);
  });
});
