// cPRA-MEX calculation engine (browser port of code/lib_cpramex.R).
//
// Each unacceptable specificity is resolved to the reference alleles of one
// locus, and cPRA = 1 - (1 - p)^2, where p is the frequency of haplotypes
// carrying at least one of them. Class and locus values restrict the target
// set; they are never added.
//
// Matching follows the primary analysis: an allele target matches every
// reference allele of its IPD-IMGT/HLA P group (identical antigen-recognition
// domain); antigens follow the WHO assignments within the P group; DP
// antigens are DPB1 first field; null alleles carry no antigen.
// tests/verify_engine.js checks this file against the analysis results.

(function (global) {
  'use strict';

  const GENE_LOCUS = {
    A: 'A', B: 'B', C: 'C', DRB1: 'DRB1', DRB3: 'DRB345', DRB4: 'DRB345',
    DRB5: 'DRB345', DQA1: 'DQA1', DQB1: 'DQB1', DPA1: 'DPA1', DPB1: 'DPB1'
  };
  const SCOPES = {
    Global: ['A', 'B', 'C', 'DRB345', 'DRB1', 'DQA1', 'DQB1', 'DPA1', 'DPB1'],
    ClaseI: ['A', 'B', 'C'],
    ClaseII: ['DRB345', 'DRB1', 'DQA1', 'DQB1', 'DPA1', 'DPB1'],
    A: ['A'], B: ['B'], C: ['C'], DRB1: ['DRB1'], DRB345: ['DRB345'],
    DQA1: ['DQA1'], DQB1: ['DQB1'], DPA1: ['DPA1'], DPB1: ['DPB1']
  };
  const SCOPE_NAMES = Object.keys(SCOPES);

  async function fetchBinary(baseUrl) {
    // The compressed copy is preferred; hosts that already decode it, and
    // browsers without DecompressionStream, are both handled.
    if (typeof DecompressionStream === 'function') {
      try {
        const response = await fetch(baseUrl + 'haplotypes.bin.gz');
        if (response.ok) {
          const buffer = await response.arrayBuffer();
          const head = new Uint8Array(buffer, 0, 2);
          if (head[0] !== 0x1f || head[1] !== 0x8b) return buffer;
          const stream = new Blob([buffer]).stream().pipeThrough(new DecompressionStream('gzip'));
          return await new Response(stream).arrayBuffer();
        }
      } catch (error) { /* fall back to the uncompressed file */ }
    }
    const response = await fetch(baseUrl + 'haplotypes.bin');
    if (!response.ok) throw new Error('No se pudo descargar la referencia (' + response.status + ').');
    return response.arrayBuffer();
  }

  async function loadReference(baseUrl) {
    const [metaResponse, buffer] = await Promise.all([fetch(baseUrl + 'reference.json'), fetchBinary(baseUrl)]);
    if (!metaResponse.ok) throw new Error('No se pudo descargar la referencia (' + metaResponse.status + ').');
    return buildReference(await metaResponse.json(), buffer);
  }

  function buildReference(meta, buffer) {
    const n = meta.nHaplotypes;
    const loci = meta.loci;
    let offset = 0;
    const haplotypes = loci.map(function () {
      const column = new Uint16Array(buffer.slice(offset, offset + 2 * n));
      offset += 2 * n;
      return column;
    });
    const pools = ['NACIONAL'].concat(meta.regions.map(function (region) { return region.id; }));
    const weights = pools.map(function () { return new Float64Array(n); });
    meta.regions.forEach(function (region, r) {
      const index = new Uint32Array(buffer.slice(offset, offset + 4 * region.nonzero));
      offset += 4 * region.nonzero;
      const value = new Float64Array(buffer.slice(offset, offset + 8 * region.nonzero));
      offset += 8 * region.nonzero;
      const regional = weights[r + 1];
      const national = weights[0];
      const factor = region.frequencySum * region.poolingWeight;
      for (let k = 0; k < index.length; k++) {
        regional[index[k]] = value[k];
        national[index[k]] += value[k] * factor;
      }
    });
    if (offset !== buffer.byteLength) throw new Error('El archivo de referencia está incompleto.');
    let total = 0;
    for (let i = 0; i < n; i++) total += weights[0][i];
    for (let i = 0; i < n; i++) weights[0][i] /= total;

    // Marginal allele frequency per pool; alleles of one locus are mutually
    // exclusive on a haplotype, so these add up exactly for a single target.
    const alleleFrequency = pools.map(function (pool, p) {
      return loci.map(function (locus, l) {
        const frequency = new Float64Array(meta.alleles[locus].length);
        const column = haplotypes[l];
        const weight = weights[p];
        for (let i = 0; i < n; i++) frequency[column[i]] += weight[i];
        return frequency;
      });
    });

    const alleleIndex = {};
    loci.forEach(function (locus) {
      const map = new Map();
      meta.alleles[locus].forEach(function (allele, i) { map.set(allele, i); });
      alleleIndex[locus] = map;
    });
    const knownAlleles = new Set();
    Object.keys(meta.knownAlleles).forEach(function (gene) {
      meta.knownAlleles[gene].split(',').forEach(function (fields) { knownAlleles.add(gene + '*' + fields); });
    });

    return {
      meta: meta, loci: loci, pools: pools, nHaplotypes: n, haplotypes: haplotypes, weights: weights,
      alleles: meta.alleles, alleleIndex: alleleIndex, alleleFrequency: alleleFrequency,
      antigenIndex: meta.antigenIndex, knownAntigens: new Set(meta.knownAntigens),
      knownAlleles: knownAlleles, groupOf: meta.groupOf, alleleGroup: meta.alleleGroup
    };
  }

  function splitTokens(text) {
    return String(text == null ? '' : text).split(/[,;|\s]+/).filter(function (token) { return token.length > 0; });
  }

  function normalizeToken(token) {
    let value = token.replace(/\s+/g, '').toUpperCase();
    value = value.replace(/^HLA-/, '');
    value = value.replace(/^QA1\*/, 'DQA1*');
    value = value.replace(/^QB1\*/, 'DQB1*');
    value = value.replace(/\.+$/, '');
    value = value.replace(/(\d)\.(\d)/g, '$1:$2');
    value = value.replace(/^C(?=[0-9])/, 'CW');
    // A*0201 is accepted as A*02:01.
    value = value.replace(/^([A-Z0-9]+)\*([0-9]{2})([0-9]{2})$/, '$1*$2:$3');
    return value;
  }

  function antigenLocus(key) {
    if (/^A[0-9]+$/.test(key)) return 'A';
    if (/^B[0-9]+$/.test(key) || key === 'BW4' || key === 'BW6') return 'B';
    if (/^CW[0-9]+$/.test(key)) return 'C';
    if (key === 'DR51' || key === 'DR52' || key === 'DR53') return 'DRB345';
    if (/^DR[0-9]+$/.test(key)) return 'DRB1';
    // Serologic DQ and DP are assigned to the beta chain, as in the Toolbox.
    if (/^DQ[0-9]+$/.test(key)) return 'DQB1';
    if (/^DP[0-9]+$/.test(key)) return 'DPB1';
    return null;
  }

  function twoFieldAllele(value) {
    const pieces = /^([A-Z0-9]+)\*([0-9]{2,4}):([0-9]{2,3})(.*)$/.exec(value);
    if (!pieces) return null;
    const rest = pieces[4];
    if (!/^(:[0-9]{2,3})*[NLSCAQ]?$/.test(rest)) return null;
    const suffix = /[NLSCAQ]$/.test(rest) ? rest.charAt(rest.length - 1) : '';
    return pieces[1] + '*' + pieces[2] + ':' + pieces[3] + suffix;
  }

  function displayToken(value) {
    return value.replace(/^CW(?=[0-9])/, 'Cw').replace(/^BW(?=[46]$)/, 'Bw');
  }

  function resolveToken(token, ref) {
    const value = normalizeToken(token);
    const out = { entrada: token, normalizado: value, tipo: null, locus: null, estado: 'no_reconocido', idx: [], nota: '' };
    if (value.indexOf('*') !== -1) {
      const gene = value.replace(/\*.*$/, '');
      if (!Object.prototype.hasOwnProperty.call(GENE_LOCUS, gene)) {
        out.nota = 'Locus no admitido.';
        return out;
      }
      const locus = GENE_LOCUS[gene];
      const pool = ref.alleles[locus];
      if (/^[A-Z0-9]+\*[0-9]{2,4}$/.test(value)) {
        const prefix = value + ':';
        out.tipo = 'Grupo (primer campo)';
        out.locus = locus;
        pool.forEach(function (allele, i) { if (allele.startsWith(prefix) && !/N$/.test(allele)) out.idx.push(i); });
        let known = false;
        if (!out.idx.length) ref.knownAlleles.forEach(function (allele) { if (!known && allele.startsWith(prefix)) known = true; });
        out.estado = out.idx.length ? 'ok' : (known ? 'sin_frecuencia' : 'no_reconocido');
        out.nota = out.idx.length ? 'Incluye todos los alelos expresados de la referencia con ese primer campo.'
          : (known ? 'Ningún alelo de este grupo tiene frecuencia en la referencia mexicana.'
            : 'Grupo alélico inexistente en IPD-IMGT/HLA 3.65.0.');
        return out;
      }
      const explicitGroup = /[0-9]P$/.test(value);
      if (/[0-9]G$/.test(value)) {
        out.nota = 'Los grupos G no se admiten; escriba el alelo a dos campos o el grupo P.';
        return out;
      }
      const allele = twoFieldAllele(explicitGroup ? value.slice(0, -1) : value);
      if (allele === null) {
        out.nota = 'Formato de alelo no válido.';
        return out;
      }
      out.normalizado = explicitGroup ? allele + 'P' : allele;
      out.tipo = 'Alelo (grupo P)';
      out.locus = locus;
      const groups = ref.alleleGroup[locus];
      let group = explicitGroup ? (gene + '*' + allele.split('*')[1] + 'P') : ref.groupOf[allele];
      if (group === undefined) {
        const own = ref.alleleIndex[locus].get(allele);
        if (own !== undefined) group = groups[own];
      }
      if (group !== undefined) groups.forEach(function (g, i) { if (g === group) out.idx.push(i); });
      if (out.idx.length) {
        out.estado = 'ok';
        const labels = out.idx.map(function (i) { return pool[i]; });
        const exact = labels.indexOf(allele) !== -1;
        out.nota = (group.endsWith('P') ? 'Grupo ' + group + ': ' : '') +
          (exact && labels.length === 1 ? '' : 'se cuentan las etiquetas de la referencia con el mismo dominio de reconocimiento (' + labels.join(', ') + ').');
      } else if (ref.knownAlleles.has(allele)) {
        out.estado = 'sin_frecuencia';
        out.nota = /N$/.test(allele) ? 'Alelo nulo: no se expresa.' : 'Ni el alelo ni su grupo P tienen frecuencia en la referencia mexicana.';
      } else {
        out.nota = 'Alelo inexistente en IPD-IMGT/HLA 3.65.0.';
      }
      return out;
    }
    let key = value.replace(/^DPB0*([0-9]+)$/, 'DP$1');
    const locus = antigenLocus(key);
    if (locus === null) {
      out.nota = 'No corresponde a un antígeno o alelo HLA admitido.';
      return out;
    }
    if (!ref.knownAntigens.has(key)) {
      // A02 or B07 are accepted as A2 or B7 when the padded form is not itself an antigen.
      const unpadded = key.replace(/^([A-Z]+)0+(?=[1-9])/, '$1');
      if (ref.knownAntigens.has(unpadded)) key = unpadded;
    }
    out.normalizado = key;
    out.tipo = 'Antígeno';
    out.locus = locus;
    const index = ref.antigenIndex[locus][key];
    if (index && index.length) {
      out.idx = index.slice();
      out.estado = 'ok';
    } else if (ref.knownAntigens.has(key)) {
      out.estado = 'sin_frecuencia';
      out.nota = 'Antígeno válido, pero ningún alelo equivalente tiene frecuencia en la referencia mexicana.';
    } else {
      out.nota = 'Antígeno no incluido en la tabla de equivalencias WHO.';
    }
    return out;
  }

  function parseUnacceptables(text, ref) {
    const seen = new Set();
    const tokens = [];
    splitTokens(Array.isArray(text) ? text.join(' ') : text).forEach(function (token) {
      const resolved = resolveToken(token, ref);
      const id = resolved.tipo + ' ' + resolved.normalizado;
      if (seen.has(id)) return;
      seen.add(id);
      tokens.push(resolved);
    });
    return tokens;
  }

  // Returns cPRA as a proportion for every scope and pool, plus the cPRA of
  // each individual target: result.cpra[scope][pool], result.target[i][pool].
  function computeCpra(tokens, ref) {
    const nPools = ref.pools.length;
    const nLoci = ref.loci.length;
    const masks = ref.loci.map(function () { return null; });
    tokens.forEach(function (token) {
      if (!token.idx.length) return;
      const l = ref.loci.indexOf(token.locus);
      if (!masks[l]) masks[l] = new Uint8Array(ref.alleles[token.locus].length);
      token.idx.forEach(function (i) { masks[l][i] = 1; });
    });
    const active = [];
    for (let l = 0; l < nLoci; l++) if (masks[l]) active.push(l);

    // Haplotype mass by pattern of loci carrying a target.
    const patterns = 1 << nLoci;
    const mass = new Float64Array(nPools * patterns);
    if (active.length) {
      for (let i = 0; i < ref.nHaplotypes; i++) {
        let pattern = 0;
        for (let a = 0; a < active.length; a++) {
          const l = active[a];
          if (masks[l][ref.haplotypes[l][i]]) pattern |= 1 << l;
        }
        if (pattern === 0) continue;
        for (let p = 0; p < nPools; p++) mass[p * patterns + pattern] += ref.weights[p][i];
      }
    }
    const cpra = {};
    SCOPE_NAMES.forEach(function (scope) {
      let scopeMask = 0;
      SCOPES[scope].forEach(function (locus) { scopeMask |= 1 << ref.loci.indexOf(locus); });
      cpra[scope] = {};
      ref.pools.forEach(function (pool, p) {
        let probability = 0;
        for (let pattern = 1; pattern < patterns; pattern++) {
          if (pattern & scopeMask) probability += mass[p * patterns + pattern];
        }
        probability = Math.min(Math.max(probability, 0), 1);
        cpra[scope][pool] = 1 - (1 - probability) * (1 - probability);
      });
    });
    const target = tokens.map(function (token) {
      const row = {};
      const l = token.locus === null ? -1 : ref.loci.indexOf(token.locus);
      ref.pools.forEach(function (pool, p) {
        let probability = 0;
        token.idx.forEach(function (i) { probability += ref.alleleFrequency[p][l][i]; });
        probability = Math.min(probability, 1);
        row[pool] = 1 - (1 - probability) * (1 - probability);
      });
      return row;
    });
    return { cpra: cpra, target: target };
  }

  global.CPRA = {
    SCOPE_NAMES: SCOPE_NAMES,
    loadReference: loadReference,
    buildReference: buildReference,
    splitTokens: splitTokens,
    resolveToken: resolveToken,
    parseUnacceptables: parseUnacceptables,
    computeCpra: computeCpra,
    displayToken: displayToken
  };
})(typeof window !== 'undefined' ? window : globalThis);
