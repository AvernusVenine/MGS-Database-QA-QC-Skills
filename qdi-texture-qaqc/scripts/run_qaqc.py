"""QA/QC of the QDI texture view (qdi.vw_qdtx). Read-only: findings go to CSVs of proposed changes.

Run with ArcGIS Pro's Python (arcpy is needed only when the cache is missing or --refresh is given):
    "C:\\Program Files\\ArcGIS\\Pro\\bin\\Python\\envs\\arcgispro-py3\\python.exe" run_qaqc.py --out <dir>
"""
import argparse
import os
import re
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qdi_load  # noqa: E402

SKILL_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CODES_PATH = os.path.join(SKILL_DIR, 'references', 'stratunit_codes.tsv')
FORMATIONS_PATH = os.path.join(os.path.dirname(SKILL_DIR), 'review-till-formations', 'references',
                               'stratunit_codes.csv')

OUT_COLUMNS = ['check', 'category', 'severity', 'objectid', 'relateid', 'sample_num', 'sampledepth', 'table',
               'field', 'current_value', 'proposed_value', 'reason']

# A valid code is "rare" when it is used on at most this many rows
RARE_ROWS = 10
# ... and it is called a likely typo when a code one keystroke away is used at least this many times more often
COMMON_RATIO = 10

# Minnesota land surface elevation range, ft (lowest: Lake Superior 602 ft; highest: Eagle Mountain 2,301 ft)
ELEV_MIN_FT, ELEV_MAX_FT = 600, 2310
# Deepest plausible sample, ft (Quaternary sediment in MN is < ~1,000 ft thick)
MAX_DEPTH_FT = 1000
# Typical lab sample is ~50 g; much heavier totals are almost always keying errors
MAX_SAMPLE_G = 500

PCT_COLUMNS = [
    'sand_pct', 'silt_pct', 'clay_pct', 'gravel_pct', 'precambrian_blk', 'paleozoic_blk', 'cretaceous_blk',
    'crystalline_blk', 'carbonate_blk', 'shale_blk', 'sand_12_pct', 'fine_pct', 'coarse_pct',
    'felic_cls', 'qtz_cls', 'clear_quartz_cls', 'mafic_cls', 'meta_cls', 'iron_cls', 'rhyolite_cls',
    'pc_arkosic_cls', 'pc_arenite_cls', 'pcother_cls', 'light_cls', 'dark_cls', 'red_cls', 'carbchert_cls',
    'ss_cls', 'shale_cls', 'pother_cls', 'chert_grp', 'sst_grp', 'psh_grp', 'gray_cls', 'speck_cls', 'lms_cls',
    'shell_cls', 'pyrite_cls', 'lignite_cls', 'ostrander_cls', 'cothers_cls', 'gsh_grp', 'misc_k_grp',
    'others_pct', 'ccs_others_pct']

# Groups of percentages that should add up to 100, with the slack that rounding each part to a whole number allows
SUM_GROUPS = {
    'sand/silt/clay': ['sand_pct', 'silt_pct', 'clay_pct'],
    '3 ages (blk)': ['precambrian_blk', 'paleozoic_blk', 'cretaceous_blk'],
    'xlline/carb/shale (blk)': ['crystalline_blk', 'carbonate_blk', 'shale_blk'],
    'PC crystalline class': ['felic_cls', 'qtz_cls', 'clear_quartz_cls', 'mafic_cls', 'meta_cls', 'iron_cls',
                             'rhyolite_cls', 'pc_arkosic_cls', 'pc_arenite_cls', 'pcother_cls'],
    'light/dark/red': ['light_cls', 'dark_cls', 'red_cls'],
    'Carb-Paleozoic class': ['carbchert_cls', 'ss_cls', 'shale_cls', 'pother_cls'],
    'Carb-Paleozoic group': ['chert_grp', 'sst_grp', 'psh_grp'],
    'Shale-Cretaceous class': ['gray_cls', 'speck_cls', 'lms_cls', 'shell_cls', 'pyrite_cls', 'lignite_cls',
                               'ostrander_cls', 'cothers_cls'],
    'Shale-Cretaceous group': ['gsh_grp', 'misc_k_grp'],
}
# Each part may be rounded to a whole number, so a group of k parts can be off by up to k/2

# Count totals and the counts they are the sum of
ROLLUPS = {
    'precambrian': ['felsic', 'quartzite', 'clear_quartz', 'mafic_igneous', 'metasedvol', 'iron', 'red_volcanic',
                    'pc_arkosic', 'pc_qtz_arenite', 'pc_other'],
    'paleozoic': ['carbonate', 'paleo_ss', 'paleo_shale', 'paleo_other'],
    'shale': ['gray_shale', 'speckled_shale'],
    'cret_total': ['shale', 'limestone', 'inoceramus', 'pyrite', 'lignite', 'ostrandersand', 'cret_other'],
    'other_total': ['chert', 'unknown', 'secondary', 'other_other'],
    'grand_total': ['precambrian', 'paleozoic', 'cret_total', 'other_total'],
}

WEIGHTS = ['totalweight', 'gravel', 'adjustweight', 'adjdryweight', 'sand_12', 'sand_a', 'sand_b', 'other_weight',
           'sand_total', 'clay', 'sand_clay', 'silt']


def derived_formulas(n):
    """Stored percentage -> the same value recomputed from weights or counts.

    A formula is only used to propose values if it reproduces the stored value on nearly every row (see
    VERIFY_SHARE); the agreement is written to summary.csv so a reviewer can see why each one was trusted.
    """
    ssc = n['sand_total'] + n['silt'] + n['clay']
    P, Z, K = n['precambrian'], n['paleozoic'], n['cret_total']
    ages = P + Z + K
    ccs = P + n['carbonate_total'] + n['shale']
    f = {
        'sand_pct': n['sand_total'] / ssc, 'silt_pct': n['silt'] / ssc, 'clay_pct': n['clay'] / ssc,
        'gravel_pct': n['gravel'] / n['totalweight'],
        'sand_12_pct': n['sand_12'] / n['sand_total'],
        'precambrian_blk': P / ages, 'paleozoic_blk': Z / ages, 'cretaceous_blk': K / ages,
        'crystalline_blk': P / ccs, 'carbonate_blk': n['carbonate_total'] / ccs, 'shale_blk': n['shale'] / ccs,
        'others_pct': n['other_total'] / n['grand_total'],
        'ccs_others_pct': n['ccs_other'] / n['grand_total'],
    }
    for cls, cnt in [('felic_cls', 'felsic'), ('qtz_cls', 'quartzite'), ('clear_quartz_cls', 'clear_quartz'),
                     ('mafic_cls', 'mafic_igneous'), ('meta_cls', 'metasedvol'), ('iron_cls', 'iron'),
                     ('rhyolite_cls', 'red_volcanic'), ('pc_arkosic_cls', 'pc_arkosic'),
                     ('pc_arenite_cls', 'pc_qtz_arenite'), ('pcother_cls', 'pc_other')]:
        f[cls] = n[cnt] / P
    for cls, cnt in [('carbchert_cls', 'carbonate'), ('ss_cls', 'paleo_ss'), ('shale_cls', 'paleo_shale'),
                     ('pother_cls', 'paleo_other')]:
        f[cls] = n[cnt] / Z
    for cls, cnt in [('gray_cls', 'gray_shale'), ('speck_cls', 'speckled_shale'), ('lms_cls', 'limestone'),
                     ('shell_cls', 'inoceramus'), ('pyrite_cls', 'pyrite'), ('lignite_cls', 'lignite'),
                     ('ostrander_cls', 'ostrandersand'), ('cothers_cls', 'cret_other')]:
        f[cls] = n[cnt] / K
    f['gsh_grp'] = n['shale'] / K
    return {k: 100 * v.replace([np.inf, -np.inf], np.nan) for k, v in f.items()}


VERIFY_SHARE = 0.97     # share of rows a formula must reproduce (within DERIVED_TOL) to be trusted
DERIVED_TOL = 1.5       # stored values are rounded, sometimes to whole numbers

MUNSELL = re.compile(r'^(?:(?:10|2\.5|5|7\.5)\s?(?:R|YR|Y|GY|G|BG|B|PB|P|RP)\s?\d+(?:\.\d+)?/\d+(?:\.\d+)?'
                     r'|N\s?\d+(?:\.\d+)?/0?|GLEY\s?[12]\s.*)$')


class Findings:
    def __init__(self, df):
        self.df = df
        self.parts = []

    def add(self, check, category, severity, mask, field, reason, proposed=None, table='qdtx', current=None):
        """One finding per row of df where mask is True. proposed/reason/current may be scalars or Series."""
        rows = self.df[mask]
        if rows.empty:
            return
        pick = lambda v: v[mask].to_numpy() if isinstance(v, pd.Series) else v  # noqa: E731
        out = pd.DataFrame({
            'check': check, 'category': category, 'severity': severity,
            'objectid': rows['objectid'].to_numpy(), 'relateid': rows['relateid'].to_numpy(),
            'sample_num': rows['sample_num'].to_numpy(), 'sampledepth': rows['sampledepth'].to_numpy(),
            'table': table, 'field': field,
            'current_value': pick(current) if current is not None else (rows[field].to_numpy()
                                                                         if field in rows else ''),
            'proposed_value': pick(proposed) if proposed is not None else '',
            'reason': pick(reason)})
        self.parts.append(out[OUT_COLUMNS])

    def add_rows(self, records):
        if records:
            self.parts.append(pd.DataFrame(records)[OUT_COLUMNS])

    def frame(self):
        if not self.parts:
            return pd.DataFrame(columns=OUT_COLUMNS)
        return pd.concat(self.parts, ignore_index=True)


def numeric(df):
    cols = set(PCT_COLUMNS) | set(WEIGHTS) | {c for k, v in ROLLUPS.items() for c in [k, *v]} | {
        'sampledepth', 'depth_top', 'depth_bottom', 'elevation', 'moist_correct', 'carbonate_total', 'ccs_other',
        'qdixid'}
    return pd.DataFrame({c: pd.to_numeric(df[c], errors='coerce') for c in cols if c in df.columns})


def fmt(x):
    return '' if pd.isna(x) else f'{x:.2f}'.rstrip('0').rstrip('.')


# ---------------------------------------------------------------- A. stratunit codes

def osa_distance(a, b):
    """Edit distance counting a swap of two neighbouring letters as one edit (QSHU -> QHSU = 1)."""
    d = [[0] * (len(b) + 1) for _ in range(len(a) + 1)]
    for i in range(len(a) + 1):
        d[i][0] = i
    for j in range(len(b) + 1):
        d[0][j] = j
    for i in range(1, len(a) + 1):
        for j in range(1, len(b) + 1):
            cost = a[i - 1] != b[j - 1]
            d[i][j] = min(d[i - 1][j] + 1, d[i][j - 1] + 1, d[i - 1][j - 1] + cost)
            if i > 1 and j > 1 and a[i - 1] == b[j - 2] and a[i - 2] == b[j - 1]:
                d[i][j] = min(d[i][j], d[i - 2][j - 2] + 1)
    return d[len(a)][len(b)]


def is_transposition(a, b):
    diff = [i for i in range(min(len(a), len(b))) if a[i] != b[i]]
    return len(a) == len(b) and len(diff) == 2 and diff[1] == diff[0] + 1 and sorted(a) == sorted(b)


def neighbours(code, valid, counts):
    """Valid codes one keystroke from code: swapped letters first, then the most used."""
    near = [v for v in valid if v != code and osa_distance(code, v) == 1]
    return sorted(near, key=lambda v: (not is_transposition(code, v), -counts.get(v, 0), v))


def same_formation(a, b):
    """Codes that differ only in the last letter name the same formation (member/provenance), e.g. QNUH/QNUD."""
    return a[:3] == b[:3]


def check_strat(df, F, codes, formation_of):
    raw = df['stratunit']
    s = raw.fillna('').astype(str)
    norm = s.str.strip().str.upper()
    valid = set(codes.index)
    coded = norm != ''
    counts = norm[coded].value_counts().to_dict()
    desc = codes['description'].to_dict()
    evidence = lambda c: f"{c} ({desc.get(c, 'not in table')}; {counts.get(c, 0)} rows)"  # noqa: E731

    # Formatting: whitespace-only, stray spaces, lower case
    blank_ws = (s != '') & (norm == '')
    F.add('STRAT_WHITESPACE', 'stratunit', 'medium', blank_ws, 'stratunit',
          'stratunit holds only spaces; blank means unclassified, so it should be NULL', proposed='NULL',
          current=s.map(repr))
    messy = coded & (s != norm)
    F.add('STRAT_FORMAT', 'stratunit', 'medium', messy, 'stratunit',
          'stratunit has stray spaces or lower case', proposed=norm, current=s.map(repr))

    # Codes not in the table: propose the most-used valid code one keystroke away
    unknown_codes = sorted(set(norm[coded]) - valid)
    proposal, why = {}, {}
    for c in unknown_codes:
        near = neighbours(c, valid, counts)
        proposal[c] = '|'.join(near[:4]) if near else ''
        why[c] = (f'{c} ({counts[c]} rows) is not in the strat code table'
                  + (f"; one keystroke from {', '.join(evidence(v) for v in near[:4])}" if near
                     else '; no valid code within one keystroke - needs a geologist'))
    unk = coded & norm.isin(unknown_codes)
    F.add('STRAT_UNKNOWN_CODE', 'stratunit', 'high', unk, 'stratunit', norm.map(why),
          proposed=norm.map(proposal))

    # Rare valid codes one keystroke from a much more common code of a different formation. A member that differs
    # from its own formation only in the last letter (QFRF vs QFRU) is a deliberate choice, not a typo.
    rare_prop, rare_why = {}, {}
    for c, k in counts.items():
        if k > RARE_ROWS or c not in valid:
            continue
        common = [v for v in neighbours(c, valid, counts) if counts.get(v, 0) >= COMMON_RATIO * k
                  and not same_formation(c, v)]
        if common:
            rare_prop[c] = common[0]
            rare_why[c] = (f'rare code {evidence(c)} is one keystroke from common code '
                           f'{", ".join(evidence(v) for v in common[:3])}')
    rare = coded & norm.isin(list(rare_prop))
    # Same-borehole support: does the borehole use the suggested common code elsewhere?
    hole_codes = df.assign(_n=norm).groupby('relateid')['_n'].agg(set)
    in_hole = pd.Series([rare_prop.get(c) in hole_codes.get(r, set()) for c, r in zip(norm, df['relateid'])],
                        index=df.index)
    # A rare code is still a real unit, so it is only a firm suspect when its borehole uses the common code too
    F.add('STRAT_RARE_TYPO', 'stratunit', 'medium', rare & in_hole, 'stratunit',
          norm.map(rare_why) + '; the suggested code is also used in this borehole', proposed=norm.map(rare_prop))
    F.add('STRAT_RARE_TYPO', 'stratunit', 'low', rare & ~in_hole, 'stratunit', norm.map(rare_why),
          proposed=norm.map(rare_prop))

    # In-borehole context: a sample whose code differs by one keystroke from the code shared by the samples
    # directly above and below it
    order = df.assign(_n=norm, _d=pd.to_numeric(df['sampledepth'], errors='coerce'))[coded & df.relateid.notna()]
    order = order.sort_values(['relateid', '_d'])
    above = order.groupby('relateid')['_n'].shift(1)
    below = order.groupby('relateid')['_n'].shift(-1)
    sandwich = (above == below) & above.notna() & (order['_n'] != above)
    sandwich &= [osa_distance(a, b) == 1 if s_ else False for a, b, s_ in zip(order['_n'], above.fillna(''),
                                                                              sandwich)]
    sandwich &= ~order['_n'].isin(list(rare_prop))           # already reported above
    idx = order.index[sandwich]
    member_only = pd.Series(False, index=df.index)
    member_only.loc[idx] = [same_formation(c, a) for c, a in zip(order['_n'][sandwich], above[sandwich])]
    prop = pd.Series('', index=df.index)
    prop.loc[idx] = above[sandwich].to_numpy()
    reason = pd.Series('', index=df.index)
    reason.loc[idx] = [f'{c} sits between two samples coded {a} in the same borehole and is one keystroke from it'
                       for c, a in zip(order['_n'][sandwich], above[sandwich])]
    in_idx = df.index.isin(idx)
    F.add('STRAT_CONTEXT_TYPO', 'stratunit', 'medium', in_idx & ~member_only, 'stratunit', reason, proposed=prop)
    F.add('STRAT_CONTEXT_MEMBER', 'stratunit', 'low', in_idx & member_only, 'stratunit',
          reason + ' (same formation, different member/provenance: confirm)', proposed=prop)

    # stratformation must be a valid code and the formation the stratunit belongs to
    sf = df['stratformation'].fillna('').astype(str).str.strip().str.upper()
    has_sf = sf != ''
    F.add('STRATFORMATION_UNKNOWN', 'stratunit', 'high', has_sf & ~sf.isin(valid), 'stratformation',
          'stratformation is not in the strat code table')
    expected = norm.map(formation_of)
    bad = has_sf & sf.isin(valid) & coded & expected.notna() & (sf != expected)
    F.add('STRATFORMATION_MISMATCH', 'stratunit', 'medium', bad, 'stratformation',
          'stratformation is not the formation of stratunit ' + norm + ' (' + expected.fillna('') + ')',
          proposed=expected)
    sm = df['stratmember'].fillna('').astype(str).str.strip().str.upper()
    F.add('STRATMEMBER_UNKNOWN', 'stratunit', 'high', (sm != '') & ~sm.isin(valid), 'stratmember',
          'stratmember is not in the strat code table')

    return pd.DataFrame({'code': list(counts), 'rows': list(counts.values())}).assign(
        in_table=lambda d: d.code.isin(valid),
        description=lambda d: d.code.map(desc).fillna('NOT IN TABLE')).sort_values('rows', ascending=False)


# ---------------------------------------------------------------- B. ranges, sums, weights, counts

def check_values(df, n, F, stats):
    for c in PCT_COLUMNS:
        if c not in n:
            continue
        F.add('PCT_OUT_OF_RANGE', 'ranges_sums', 'high', (n[c] < 0) | (n[c] > 100.5), c,
              f'{c} must be between 0 and 100')

    derived = derived_formulas(n)
    trusted = {}
    for c, est in derived.items():
        ok = n[c].notna() & est.notna()
        share = float(((n[c] - est).abs() <= DERIVED_TOL)[ok].mean()) if ok.any() else np.nan
        stats.append({'item': f'formula agreement: {c}', 'value': round(share, 4), 'rows': int(ok.sum())})
        if share >= VERIFY_SHARE:
            trusted[c] = est
    # Recomputed values can only be proposed when their inputs are sane
    inputs_ok = (n[WEIGHTS].fillna(0) >= 0).all(axis=1) & (n['totalweight'] <= MAX_SAMPLE_G)

    for name, cols in SUM_GROUPS.items():
        tol = 0.5 * len(cols)
        s = n[cols].sum(axis=1, min_count=1)
        off = s.notna() & (s != 0) & ((s - 100).abs() > tol)
        proposal = pd.Series('', index=n.index)
        if all(c in trusted for c in cols):
            rec = pd.concat([trusted[c] for c in cols], axis=1)
            good = inputs_ok & ((rec.sum(axis=1) - 100).abs() <= tol)
            proposal[good] = rec[good].apply(lambda r: '; '.join(f'{c}={fmt(v)}' for c, v in zip(cols, r)), axis=1)
        F.add('GROUP_SUM_NOT_100', 'ranges_sums', 'high' if name == 'sand/silt/clay' else 'medium', off,
              '+'.join(cols), f'{name} percentages add up to ' + s.map(fmt) + ' instead of 100',
              proposed=proposal, current=n[cols].apply(lambda r: '; '.join(fmt(v) for v in r), axis=1))

    # Texture stored as all zeros although weights exist -> "no data" coded as 0, or a lost calculation
    ssc = n[['sand_pct', 'silt_pct', 'clay_pct']].sum(axis=1, min_count=1)
    have_w = (n['sand_total'].fillna(0) + n['silt'].fillna(0) + n['clay'].fillna(0)) > 0
    zero_txt = (ssc == 0)
    F.add('TEXTURE_ZERO_WITH_WEIGHTS', 'ranges_sums', 'medium', zero_txt & have_w & inputs_ok,
          'sand_pct+silt_pct+clay_pct', 'texture % are all 0 although sand/silt/clay weights are present',
          proposed=pd.Series(['; '.join(f'{c}={fmt(trusted[c][i])}' for c in ['sand_pct', 'silt_pct', 'clay_pct']
                                        if c in trusted) for i in n.index], index=n.index),
          current='0; 0; 0')
    F.add('TEXTURE_ZERO_NO_DATA', 'ranges_sums', 'low', zero_txt & ~have_w, 'sand_pct+silt_pct+clay_pct',
          'texture % and weights are all 0: "not measured" is stored as 0, which reads as real data',
          proposed='NULL', current='0; 0; 0')

    # Stored value disagrees with the value recomputed from the weights / counts
    agreement = {d['item'].split(': ')[1]: d['value'] for d in stats}
    for c, est in trusted.items():
        diff = (n[c] - est).abs()
        bad = n[c].notna() & est.notna() & (diff > DERIVED_TOL)
        prop = est.map(fmt).where(inputs_ok, '')
        # Formulas exact on >= 99.5% of rows give firm proposals; the others are worth a look only
        sev = 'medium' if agreement[c] >= 0.995 else 'low'
        F.add('DERIVED_PCT_MISMATCH', 'ranges_sums', sev, bad, c,
              f'{c} does not match the value recomputed from the weights/counts (' + est.map(fmt)
              + f'; formula reproduces {agreement[c]:.1%} of stored values)', proposed=prop)

    # Weights
    for w in WEIGHTS:
        F.add('NEGATIVE_WEIGHT', 'ranges_sums', 'high', n[w] < 0, w, f'{w} (g) is negative')
    F.add('HEAVY_SAMPLE', 'ranges_sums', 'medium', n['totalweight'] > MAX_SAMPLE_G, 'totalweight',
          f'total weight > {MAX_SAMPLE_G} g (typical lab sample ~50 g); likely a keying/decimal error')
    F.add('GRAVEL_GT_TOTAL', 'ranges_sums', 'high', n['gravel'] > n['totalweight'], 'gravel',
          'gravel weight is larger than the total weight')
    expect = n['totalweight'] - n['gravel']
    F.add('ADJUSTWEIGHT_MISMATCH', 'ranges_sums', 'medium',
          n['adjustweight'].notna() & expect.notna() & ((n['adjustweight'] - expect).abs() > 0.5), 'adjustweight',
          'adjustweight should be totalweight - gravel (' + expect.map(fmt) + ')', proposed=expect.map(fmt))
    parts = n['sand_total'] + n['silt'] + n['clay']
    F.add('SSC_WEIGHTS_NE_ADJDRY', 'ranges_sums', 'medium',
          parts.notna() & n['adjdryweight'].notna() & ((parts - n['adjdryweight']).abs() > 0.5), 'adjdryweight',
          'sand_total + silt + clay weights (' + parts.map(fmt) + ') do not add up to adjdryweight')
    sc = n['sand_total'] + n['clay']
    F.add('SAND_CLAY_WEIGHT_MISMATCH', 'ranges_sums', 'low',
          sc.notna() & n['sand_clay'].notna() & ((sc - n['sand_clay']).abs() > 0.5), 'sand_clay',
          'sand_clay weight should be sand_total + clay (' + sc.map(fmt) + ')', proposed=sc.map(fmt))

    # Grain-count rollups
    for total, comps in ROLLUPS.items():
        s = n[comps].fillna(0).sum(axis=1)
        bad = n[total].notna() & ((n[total] - s).abs() > 0.5)
        lost = bad & (s == 0)
        F.add('COUNT_ROLLUP_MISMATCH', 'ranges_sums', 'medium', bad & ~lost, total,
              f'{total} should equal the sum of {"+".join(comps)} (' + s.map(fmt) + ')', proposed=s.map(fmt))
        F.add('COUNT_COMPONENTS_MISSING', 'ranges_sums', 'medium', lost, total,
              f'{total} is filled in but every one of {"+".join(comps)} is 0/NULL: component counts are missing')
    for c in sorted({c for k, v in ROLLUPS.items() for c in [k, *v]}):
        F.add('NEGATIVE_COUNT', 'ranges_sums', 'high', n[c] < 0, c, f'grain count {c} is negative')


def check_domains(df, F, stats_domains):
    for field, allowed, sev in [('deposittype', None, 'medium'), ('leached', {'Y', 'N'}, 'medium'),
                                ('type', None, 'low'), ('depth_unit', {'f', 'm'}, 'medium')]:
        s = df[field]
        txt = s.fillna('').astype(str)
        stats_domains.append(txt.replace('', '<blank/null>').value_counts().rename_axis('value')
                             .reset_index(name='rows').assign(field=field))
        F.add('VALUE_FORMAT', 'domains', sev, (txt != '') & (txt != txt.str.strip().str.upper()) &
              (txt.str.strip().str.upper() != 'T*') & ~((field == 'depth_unit') & txt.isin(['f', 'm'])),
              field, f'{field} has stray spaces or wrong case', proposed=txt.str.strip().str.upper(),
              current=txt.map(repr))
        if allowed is not None:
            F.add('VALUE_NOT_ALLOWED', 'domains', sev, (txt.str.strip() != '') & ~txt.str.strip().isin(allowed),
                  field, f'{field} must be one of {sorted(allowed)}')
    dt = df['deposittype'].fillna('').astype(str).str.strip()
    F.add('VALUE_NOT_ALLOWED', 'domains', 'medium', dt.str.fullmatch(r'\d+'), 'deposittype',
          'deposittype is a number; expected a deposit letter such as T (till) or L (lacustrine)')
    F.add('VALUE_NOT_ALLOWED', 'domains', 'low', (df['deposittype'].notna()) & (dt == ''), 'deposittype',
          'deposittype is an empty string; should be NULL or a deposit letter', proposed='NULL', current="''")
    for c in ['dry_color', 'wet_color']:
        v = df[c].fillna('').astype(str).str.strip()
        bad = (v != '') & ~v.str.upper().str.match(MUNSELL)
        F.add('MUNSELL_UNPARSEABLE', 'domains', 'low', bad, c,
              'not a Munsell notation (hue value/chroma, e.g. 10YR 5/3): free text, truncated or missing chroma')


# ---------------------------------------------------------------- C. depth and elevation

def check_depth_elev(df, n, F):
    sd, top, bot = n['sampledepth'], n['depth_top'], n['depth_bottom']
    F.add('DEPTH_SENTINEL', 'depth_elev', 'high', sd == -99, 'sampledepth',
          'sampledepth -99 is a "no value" placeholder stored as a number', proposed='NULL')
    F.add('DEPTH_NEGATIVE', 'depth_elev', 'high', (sd < 0) & (sd != -99), 'sampledepth',
          'sampledepth is negative', proposed=sd.abs().map(fmt))
    F.add('DEPTH_TOO_DEEP', 'depth_elev', 'medium', sd > MAX_DEPTH_FT, 'sampledepth',
          f'sampledepth > {MAX_DEPTH_FT} ft is deeper than Quaternary sediment in Minnesota')
    for c in ['depth_top', 'depth_bottom']:
        F.add('DEPTH_TOO_DEEP', 'depth_elev', 'medium', n[c] > MAX_DEPTH_FT, c,
              f'{c} > {MAX_DEPTH_FT} ft is deeper than Quaternary sediment in Minnesota')
    inv = top.notna() & bot.notna() & (top > bot)
    F.add('DEPTH_TOP_GT_BOTTOM', 'depth_elev', 'high', inv, 'depth_top+depth_bottom',
          'depth_top is below depth_bottom (swapped?)', current=top.map(fmt) + '; ' + bot.map(fmt),
          proposed=bot.map(fmt) + '; ' + top.map(fmt))
    # depth_top = 0 is a placeholder on most rows, so only real intervals are compared
    outside = (top > 0) & bot.notna() & ~inv & sd.notna() & (sd >= 0) & ((sd < top - 0.5) | (sd > bot + 0.5))
    F.add('SAMPLEDEPTH_OUTSIDE_INTERVAL', 'depth_elev', 'medium', outside, 'sampledepth',
          'sampledepth is outside depth_top-depth_bottom (' + top.map(fmt) + '-' + bot.map(fmt)
          + '); sampledepth may hold a sample sequence number', proposed=((top + bot) / 2).map(fmt) + ' (midpoint?)')
    F.add('DEPTH_UNIT_METRES', 'depth_elev', 'low', (df['depth_unit'] == 'm') & sd.notna(), 'depth_unit',
          'depth_unit is "m" but QDI depths are in feet: confirm whether sampledepth was converted '
          '(value in ft would be ' + (sd * 3.28084).map(fmt) + ')')
    F.add('DEPTH_UNIT_MISSING', 'depth_elev', 'low', df['depth_unit'].isna() & sd.notna(), 'depth_unit',
          'depth_unit is NULL', proposed='f')

    # Elevation comes from the qdix (borehole) table: report each borehole once
    holes = df.assign(_e=n['elevation']).dropna(subset=['relateid']).groupby('relateid').agg(
        elevation=('_e', 'first'), objectid=('objectid', 'first'), n=('objectid', 'size'))
    bad = holes.elevation.notna() & ((holes.elevation < ELEV_MIN_FT) | (holes.elevation > ELEV_MAX_FT))
    recs = []
    for rid, h in holes[bad].iterrows():
        e = h.elevation
        guess = ''
        if e <= 0:
            guess = 'NULL (missing elevation stored as a number)'
        elif e > ELEV_MAX_FT and ELEV_MIN_FT <= e / 10 <= ELEV_MAX_FT:
            guess = f'{e / 10:g} (extra digit?)'
        elif ELEV_MIN_FT <= e * 3.28084 <= ELEV_MAX_FT:
            guess = f'{e * 3.28084:.0f} (value looks like metres)'
        recs.append({'check': 'ELEVATION_OUT_OF_RANGE', 'category': 'depth_elev',
                     'severity': 'high' if e > 5000 else 'medium', 'objectid': h.objectid,
                     'relateid': rid, 'sample_num': f'{h.n} samples', 'sampledepth': '', 'table': 'qdix',
                     'field': 'elevation', 'current_value': e, 'proposed_value': guess,
                     'reason': f'borehole elevation outside Minnesota range {ELEV_MIN_FT}-{ELEV_MAX_FT} ft'})
    F.add_rows(recs)


# ---------------------------------------------------------------- D. integrity

def check_integrity(df, n, F):
    F.add('NULL_RELATEID', 'integrity', 'high', df['relateid'].isna(), 'relateid',
          'sample has no relateid (Q number): it is not tied to a borehole and is hidden by the published layer, '
          'whose object-ID field is relateid')
    # relateid should be a 10-character CWI number or Q number (e.g. 0000242265, 00Q0016567)
    rid = df['relateid'].fillna('').astype(str)
    padded = rid.str.zfill(10)
    placeholder = rid.str.fullmatch(r'\d{2}Q\d{2}')
    odd = (rid != '') & ~rid.str.fullmatch(r'\d{10}|\d{2}Q\d{7}')
    short = rid.str.fullmatch(r'\d{1,9}')
    pad_ok = odd & short & padded.isin(set(rid))
    F.add('RELATEID_PLACEHOLDER', 'integrity', 'high', placeholder, 'relateid',
          'relateid is a placeholder (e.g. 00Q00), not a real borehole id: samples from different boreholes '
          'are lumped together')
    F.add('RELATEID_FORMAT', 'integrity', 'medium', pad_ok, 'relateid',
          'relateid has lost its leading zeros; the padded id exists in the database', proposed=padded)
    F.add('RELATEID_FORMAT', 'integrity', 'medium', odd & ~pad_ok & ~placeholder, 'relateid',
          'relateid is not a 10-character CWI number or Q number', proposed=padded.where(short, ''))

    holes = df[df.relateid.notna() & n['qdixid'].isna()].groupby('relateid').agg(
        objectid=('objectid', 'first'), n=('objectid', 'size'))
    F.add_rows([{'check': 'NO_LOCATION', 'category': 'integrity', 'severity': 'medium', 'objectid': h.objectid,
                 'relateid': rid, 'sample_num': f'{h.n} samples', 'sampledepth': '', 'table': 'qdix',
                 'field': 'relateid', 'current_value': rid, 'proposed_value': '',
                 'reason': 'relateid has no qdix record, so its samples have no location'}
                for rid, h in holes.iterrows()])

    # Duplicates: same borehole, depth and sample name
    key = df.assign(_s=df['sample_num'].fillna('').astype(str).str.strip().str.upper(),
                    _d=n['sampledepth'])
    key = key[key.relateid.notna()]
    skip = {'objectid', 'qdixid', 'date_modify', 'importseq', 'groupid', 'x', 'y', 'importerrors1',
            'importerrors2', 'sourcefile', '_s', '_d'}
    value_cols = [c for c in key.columns if c not in skip]
    dup = key[key.duplicated(['relateid', '_d', '_s'], keep=False)]
    recs = []
    for (rid, d, s), g in dup.groupby(['relateid', '_d', '_s'], dropna=False):
        vals = g[value_cols].astype(str)
        differs = [c for c in value_cols if vals[c].nunique() > 1]
        same = not differs
        keep = g.objectid.min()
        for _, r in g.iterrows():
            if r.objectid == keep:
                continue
            recs.append({'check': 'DUPLICATE_EXACT' if same else 'DUPLICATE_CONFLICT', 'category': 'integrity',
                         'severity': 'medium' if same else 'high', 'objectid': r.objectid, 'relateid': rid,
                         'sample_num': r.sample_num, 'sampledepth': r.sampledepth, 'table': 'qdtx',
                         'field': '(row)', 'current_value': f'duplicate of objectid {keep}',
                         'proposed_value': 'delete row' if same else 'reconcile with objectid ' + str(keep),
                         'reason': ('identical copy of another sample (same relateid, depth, sample_num)' if same
                                    else 'same relateid, depth and sample_num as another row; differs in '
                                    + ', '.join(differs[:8]) + (' ...' if len(differs) > 8 else ''))})
    F.add_rows(recs)


# ---------------------------------------------------------------- main

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--lyr', default=qdi_load.DEFAULT_LYR)
    ap.add_argument('--out', default=r'C:\Users\sturm113\PycharmProjects\DatabaseQAQC\errors\qdi_texture')
    ap.add_argument('--cache', default=r'C:\Users\sturm113\PycharmProjects\DatabaseQAQC\cache\qdtx_full.pkl')
    ap.add_argument('--codes', default=CODES_PATH)
    ap.add_argument('--refresh', action='store_true', help='re-read the database instead of the cache')
    a = ap.parse_args(argv)

    df, meta = qdi_load.load(a.lyr, a.cache, a.refresh)
    df = df.reset_index(drop=True)
    n = numeric(df)
    codes = pd.read_csv(a.codes, sep='\t', dtype=str, keep_default_na=False).set_index('code')
    formation_of = {}
    if os.path.exists(FORMATIONS_PATH):
        formation_of = pd.read_csv(FORMATIONS_PATH, dtype=str).set_index('code')['formation_code'].to_dict()

    F = Findings(df)
    stats, domains = [], []
    code_counts = check_strat(df, F, codes, formation_of)
    check_values(df, n, F, stats)
    check_domains(df, F, domains)
    check_depth_elev(df, n, F)
    check_integrity(df, n, F)

    out = F.frame()
    out['severity'] = pd.Categorical(out['severity'], ['high', 'medium', 'low'], ordered=True)
    out = out.sort_values(['category', 'severity', 'check', 'relateid', 'sampledepth'], kind='stable')
    os.makedirs(a.out, exist_ok=True)
    files = {'stratunit': 'strat_codes.csv', 'ranges_sums': 'ranges_sums.csv', 'domains': 'domains.csv',
             'depth_elev': 'depth_elev.csv', 'integrity': 'integrity.csv'}
    for cat, name in files.items():
        out[out.category == cat].to_csv(os.path.join(a.out, name), index=False)
    out.to_csv(os.path.join(a.out, 'all_proposed_changes.csv'), index=False)
    code_counts.to_csv(os.path.join(a.out, 'stratunit_code_counts.csv'), index=False)
    pd.concat(domains).to_csv(os.path.join(a.out, 'field_value_counts.csv'), index=False)

    summary = (out.groupby(['category', 'check', 'severity'], observed=True)
               .agg(findings=('objectid', 'size'), boreholes=('relateid', 'nunique')).reset_index())
    summary.to_csv(os.path.join(a.out, 'summary.csv'), index=False)
    info = pd.concat([meta.rename('value').rename_axis('item').reset_index(), pd.DataFrame(stats)],
                     ignore_index=True)
    info.to_csv(os.path.join(a.out, 'run_info.csv'), index=False)
    print(meta.to_string())
    print(summary.to_string(index=False))
    print(f'{len(out)} findings on {out.objectid.nunique()} rows -> {a.out}')


if __name__ == '__main__':
    main()
