"""DRAFT. Round 30: the third-level prediction held against gravitational-wave data.
Run after third_level_prediction.json was committed on its own (82b48d2); the bounds below were looked up after.

Measured: GWTC-4.0 Tests of General Relativity II, Parameterized Tests (LVK, arXiv:2603.19020), Table 2,
fractional deviation of the non-spinning 3PN non-log phase coefficient, 66 events combined, 90% credible.
    python3 third_level_compare.py      # writes third_level_compare.json
"""
import json, os
HERE = os.path.dirname(os.path.abspath(__file__))
GWTC4 = {'source': 'GWTC-4.0 TGR II (arXiv:2603.19020), Table 2, d_phi6 (3PN non-log), 90% credible',
         'FTI_hierarchical': (0.00, -0.18, 0.18), 'FTI_restricted': (-0.01, -0.06, 0.06),
         'TIGER_hierarchical': (0.00, -0.29, 0.28), 'TIGER_restricted': (0.00, -0.14, 0.13)}


def compare():
    pr = json.load(open(os.path.join(HERE, 'third_level_prediction.json')))
    rows = {r['nu']: r for r in pr['floor_prediction']}
    out = {'prediction_written_utc': pr['written_utc'], 'prediction_sha256': pr['script_sha256'], 'measured': GWTC4,
           'by_mass_ratio': {}}
    for n in (0.2, 0.2222, 0.24, 0.25):
        f = rows[n]['fractional_d_phi6']
        out['by_mass_ratio'][str(n)] = {
            'q': rows[n]['mass_ratio_q'], 'floor_fraction': round(f, 3),
            'times_the_upper_bound': {k: round(f/v[2], 1) for k, v in GWTC4.items() if k != 'source'}}
    out['test_body'] = rows[0]['fractional_d_phi6']
    out['verdict'] = ('If the rules got every rational part of the third level right, their missing pi^2 alone would '
                      'move phi_6 by 0.65-0.69 of itself for the comparable-mass binaries that dominate the catalogue; '
                      'the 66-event bound allows 0.18 (hierarchical) or 0.06 (all events sharing one value). '
                      'The structural miss is certain; its size is a floor only if the rational parts agree.')
    json.dump(out, open(os.path.join(HERE, 'third_level_compare.json'), 'w'), indent=1)
    return out


if __name__ == '__main__':
    print(json.dumps(compare(), indent=1))
