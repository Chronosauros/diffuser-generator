// Recreate published QRD recipes with the live generator core (index.html <script id="diffuser-core">).
//   node checks/recipes.mjs [path/to/index.html]
// Sources of the "known answers":
//   COX  - standard N = 7 QRD of T. J. Cox, J. A. S. Angus, P. D'Antonio, "The performance of Schroeder diffusers:
//          a wide band BEM investigation", Proc. IOA 25(5), 2003, and J. A. Hargreaves, T. J. Cox (Salford):
//          f0 490 Hz, unit depth 0.05 m, well width 0.075 m, period 0.525 m, first critical frequency about 3.4 kHz.
//   QRDUDE - subwoofer-builder.com/qrd.htm (QRDude calculator reference, after Cox & D'Antonio): residue tables
//          N = 7, 11, 13, build depth fractions 4/7, 9/11, 12/13, "N7 at 500 Hz -> plate frequency 3500 Hz",
//          2D depth (X^2 + Y^2) mod N, rules: well >= 25 mm up to 400 mm deep, period >= 2 x design depth.
//   GEARSPACE - forum worked example: 28 in period, 7 wells, design depth 14 in, depths 0 2 8 4 4 8 2 in.
//   SOUNDASSURED - web calculator: deepest well 171.5 mm for a 1000 Hz target.
import fs from 'node:fs';
import vm from 'node:vm';

const INDEX = process.argv[2] || new URL('../index.html', import.meta.url);
const html = fs.readFileSync(INDEX, 'utf8');
const tag = /<script\b[^>]*\bid\s*=\s*["']diffuser-core["'][^>]*>([\s\S]*?)<\/script>/i.exec(html);
const ctx = vm.createContext({ Math, console });
const core = vm.runInContext(tag[1] + '\n;createDiffuser', ctx)(null);   // no manifold: geometry builders unused here

const C = 343, IN = 25.4;
let fails = 0;
const ok = (cond) => { if (!cond) fails++; return cond ? 'PASS' : 'FAIL'; };
const rule = (cond) => (cond ? 'OK' : 'BELOW RULE');
const near = (a, b, tol) => Math.abs(a - b) <= tol;
const eqArr = (a, b, tol = 1e-9) => a.length === b.length && a.every((v, i) => near(v, b[i], tol));
const hz = (f) => f >= 1000 ? (f / 1000).toFixed(2) + ' kHz' : f.toFixed(1) + ' Hz';
const P = (o) => Object.assign({ fin: 0, edge: 0 }, o);
const TILE_SHIFT0 = 4;                                     // deepShift(4) = 0: the unshifted sequence 0 1 4 2 2 4 1

console.log('A. Published recipes through the generator core\n');

// COX standard N = 7
{
  const p = P({ N: 7, pitch: 75, dmax: 4 * 50 });
  const b = core.deepBand(p), d = core.deepDepths(p, TILE_SHIFT0);
  console.log('COX standard N=7 (unit 50 mm, well 75 mm):');
  console.log(`  depths ${d.join(' ')} mm  vs published 0 50 200 100 100 200 50  ${ok(eqArr(d, [0, 50, 200, 100, 100, 200, 50]))}`);
  console.log(`  f0 ${hz(b.f0)}  vs 490 Hz  ${ok(near(b.f0, 490, 0.5))}`);
  console.log(`  flat plate N*f0 ${hz(b.f_flat)}  vs about 3.4 kHz  ${ok(near(b.f_flat, 3430, 50))}`);
  console.log(`  period ${b.period} mm  vs 525 mm  ${ok(b.period === 525)}`);
  console.log(`  well-width limit c/2w ${hz(b.f_width)},  c/period ${hz(C / 0.525)}  (period / design depth = ${(525 / (7 * 50)).toFixed(2)})\n`);
}

// QRDude tables
{
  const pub = { 7: [0, 1, 4, 2, 2, 4, 1], 11: [0, 1, 4, 9, 5, 3, 3, 5, 9, 4, 1], 13: [0, 1, 4, 9, 3, 12, 10, 10, 12, 3, 9, 4, 1] };
  const frac = { 7: 4 / 7, 11: 9 / 11, 13: 12 / 13 };
  console.log('QRDUDE residue tables and build depth:');
  for (const N of [7, 11, 13]) {
    const r = core.residues(N);
    console.log(`  N=${N}: ${r.join(' ')}  ${ok(eqArr(r, pub[N]))};  s_max/N = ${(Math.max(...r) / N).toFixed(4)} vs ${frac[N].toFixed(4)}  ${ok(near(Math.max(...r) / N, frac[N], 1e-12))}`);
  }
  const dmax500 = 4 * C / (2 * 7 * 500) * 1000;            // deepest well of an N7 panel designed for 500 Hz
  const b = core.deepBand(P({ N: 7, pitch: 30, dmax: dmax500 }));
  console.log(`  N7 at 500 Hz (deepest ${dmax500.toFixed(1)} mm): f0 ${hz(b.f0)}, plate ${hz(b.f_flat)} vs 3500 Hz  ${ok(near(b.f_flat, 3500, 0.5))}`);
  let bad = 0;
  for (let n = 0; n < 13; n++) for (let m = 0; m < 13; m++)
    if (Math.round(core.phaseFraction(13, 13, n, m, 0, 0) * 13) !== (n * n + m * m) % 13) bad++;
  console.log(`  2D 13x13: (n^2 + m^2) mod 13 in all 169 cells, mismatches ${bad}  ${ok(bad === 0)}\n`);
}

// GEARSPACE forum example
{
  const p = P({ N: 7, pitch: 4 * IN, dmax: 8 * IN });
  const b = core.deepBand(p), d = core.deepDepths(p, TILE_SHIFT0).map((v) => +(v / IN).toFixed(6));
  console.log('GEARSPACE 28 in, 7 wells, design depth 14 in:');
  console.log(`  depths ${d.join(' ')} in  vs 0 2 8 4 4 8 2  ${ok(eqArr(d, [0, 2, 8, 4, 4, 8, 2], 1e-6))}`);
  console.log(`  f0 ${hz(b.f0)}, design wavelength ${(b.lambda0 / IN).toFixed(2)} in, design depth (N units) ${(7 * b.unit_depth / IN).toFixed(2)} in vs 14  ${ok(near(7 * b.unit_depth / IN, 14, 1e-6))}\n`);
}

// SOUNDASSURED calculator: a different convention
{
  const b = core.deepBand(P({ N: 7, pitch: 30, dmax: 171.5 }));
  const textbook = 4 * C / (2 * 7 * 1000) * 1000;
  console.log('SOUNDASSURED deepest well 171.5 mm "for 1000 Hz":');
  console.log(`  Cox & D'Antonio / QRDude give f0 ${hz(b.f0)} for that depth; their 1000 Hz panel needs ${textbook.toFixed(1)} mm.`);
  console.log('  That calculator scales the deepest well to lambda/2 (s_max instead of N in the denominator): not the textbook QRD.\n');
}

console.log('B. Our presets (default depth) against the QRDude rules\n');
const beds = core.BEDS || [
  { key: 'a1mini', short: 'A1 mini', x: 180, y: 180, z: 180 },
  { key: 'mk4', short: 'MK4 / Core One', x: 250, y: 210, z: 220 },
  { key: 'p1s', short: 'P1S / X1', x: 256, y: 256, z: 256 }];
const rows = [];
for (const bed of beds) {
  const lim = core.limits('deep', bed), p = core.derive('deep', bed, lim.depthDefault), b = core.deepBand(p);
  const wIn = p.pitch - p.fin, wEdge = p.pitch - p.fin / 2 - p.edge;
  const designDepth = 7 * p.dmax / 4, fPeriod = C / (b.period / 1000);
  rows.push({ bed: bed.short || bed.name, pitch: p.pitch, well: wIn.toFixed(1), edgeWell: wEdge.toFixed(1), dmax: p.dmax,
    f0: hz(b.f0), fWell: hz(b.f_width), flat: hz(b.f_flat), period: b.period, periodOverDesignDepth: (b.period / designDepth).toFixed(2),
    fPeriod: hz(fPeriod), depthMax: lim.depthMax, f0AtMax: hz(4 * C / (14 * lim.depthMax / 1000)),
    well25: rule(Math.min(wIn, wEdge) >= 25), depth16: rule(p.dmax <= 16 * Math.min(wIn, wEdge)), period2x: rule(b.period >= 2 * designDepth) });
}
console.table(rows);
console.log('  (rule columns are design findings against the QRDude rules of thumb, not errors in the math)');

console.log('\nC. What sets the low end (P1S): depth slider vs bed\n');
{
  const p1s = beds.find((b) => /P1S/.test(b.short || b.name));
  const lim = core.limits('deep', p1s);
  for (const d of [60, 100, 132, 200, lim.depthMax]) {
    const b = core.deepBand(core.derive('deep', p1s, d));
    console.log(`  depth ${String(d).padStart(3)} mm -> f0 ${hz(b.f0).padStart(9)}, c/2w ${hz(b.f_width)}, one tile c/period ${hz(C / (b.period / 1000))}`);
  }
  for (const bed of beds) {
    const b = core.deepBand(core.derive('deep', bed, 132));
    console.log(`  ${(bed.short || bed.name).padEnd(15)} at 132 mm -> f0 ${hz(b.f0)}, c/2w ${hz(b.f_width)}, one tile c/period ${hz(C / (b.period / 1000))}`);
  }
  console.log(`  deepest well for f0 = 200 Hz with N = 7: ${(4 * C / (14 * 200) * 1000).toFixed(0)} mm`);
}

console.log('\nD. Can the live page rebuild the COX standard diffuser?');
for (const bed of beds) {
  const lim = core.limits('deep', bed), p = core.derive('deep', bed, Math.min(200, lim.depthMax));
  console.log(`  ${(bed.short || bed.name).padEnd(15)} depth 200 mm ${lim.depthMax >= 200 ? 'reachable' : 'NOT reachable (max ' + lim.depthMax + ')'}, well pitch fixed at ${p.pitch} mm (recipe 75 mm), period ${7 * p.pitch} mm (recipe 525 mm)`);
}
console.log(`\n${fails === 0 ? 'ALL RECIPE CHECKS PASS' : fails + ' RECIPE CHECK(S) FAILED'}`);
process.exitCode = fails === 0 ? 0 : 1;
