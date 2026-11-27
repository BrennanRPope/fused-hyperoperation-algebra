"""DRAFT, round 29 (Nov 27, 2026): the echo prediction (echo_prediction.json) against Cassini.

Run after echo_predict.py. The prediction was written first (see its written_utc and the commit that holds it); this
script only reads it and compares with the measurement:
  Bertotti, Iess & Tortora 2003, Nature 425, 374: gamma = 1 + (2.1 +- 2.3) x 10^-5, from the Doppler tracking of the
  Cassini spacecraft during its June 2002 solar conjunction.
In the river gamma = lambda^2: the push's share of the delay over the time part's (whose size is fixed by Newton's mass).
Run: python3 echo_compare.py  (writes echo_compare.json, instant)
"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
CASSINI = {"gamma_minus_1": 2.1e-5, "sigma": 2.3e-5, "source": "Bertotti, Iess & Tortora 2003, Nature 425, 374"}


def compare():
    pred = json.load(open(os.path.join(HERE, "echo_prediction.json")))
    g, s = 1 + CASSINI["gamma_minus_1"], CASSINI["sigma"]
    p = pred["prediction"]
    lam_lo, lam_hi = pred["shadow_dial"]["lambda_range"]
    rt = p["round_trip_delay_s"]
    total = rt["total"]
    out = {"prediction_written_utc": pred["written_utc"], "prediction_sha256": pred["script_sha256"],
           "prediction_used_measured_delay": pred["uses_any_measured_delay_or_gamma"],
           "cassini": CASSINI,
           "derived_dials": {"gamma_predicted": p["gamma_predicted"],
                             "pull_in_sigma": (g - p["gamma_predicted"])/s},
           "push_over_time_part_measured": [g, s],
           "lambda_measured": [g**0.5, s/2],
           "shadow_only": {"gamma_range": pred["shadow_dial"]["gamma_range_lambda2"],
                           "cassini_inside": lam_lo**2 <= g <= lam_hi**2,
                           "narrowing_factor": (lam_hi**2 - lam_lo**2)/(2*s)},
           "time_part_alone": {"gamma": 0.0, "pull_in_sigma": (g - 0.0)/s},
           "newton_corpuscle": {"gamma": -2.0, "pull_in_sigma": (g + 2.0)/s},
           "twist": {"round_trip_s": rt["twist_round_trip"],
                     "cassini_delay_resolution_s": s/(1 + g)*total,
                     "twist_over_resolution": rt["twist_round_trip"]/(s/(1 + g)*total)}}
    return out


def main():
    res = compare()
    json.dump(res, open(os.path.join(HERE, "echo_compare.json"), "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
