"""Round 8: completeness from the injections (runs/round8_crowd/inject/*.json) -> completeness.csv + summary."""
import glob, json
import numpy as np, pandas as pd
from common import OUT
rows = []
for f in glob.glob(str(OUT / "inject" / "*.json")):
    r = json.load(open(f))
    rr = r.get("real")
    inj = r.get("inject", {})
    if not rr:
        rows.append({"id": r.get("id"), "error": r.get("error")}); continue
    H = rr["hold"]
    err = H.get("err_ppm") or np.nan
    rows.append({"id": r["id"], "stream": r["id"].split("_")[0], "mode": r["mode"],
                 "depth_inj": inj.get("depth_src_ppm"), "depth_nominal": inj.get("depth_ppm_nominal"),
                 "z_hold": H.get("z"), "z_exp": (inj.get("depth_src_ppm") or np.nan) / err if err else np.nan,
                 "his_P_ok": abs(rr["his"]["P"] / inj["period"] - 1) < 0.001, "pass": rr["pass"]["all"],
                 "flags": ";".join(k for k, v in rr["pass"].items() if not v and k != "all")})
d = pd.DataFrame(rows)
d.to_csv(OUT / "completeness.csv", index=False)
ok = d.dropna(subset=["z_exp"]).copy(); ok["pass"] = ok["pass"].astype(bool)
print("injections with a hold-out measurement:", len(ok), "of", len(d))
print("overall pass:", round(ok["pass"].mean(), 3))
bins = [-np.inf, 2, 4, 6, 10, np.inf]
ok["zbin"] = pd.cut(ok.z_exp, bins)
print(ok.groupby("zbin", observed=True)["pass"].agg(["size", "mean"]).round(3).to_string())
print(ok.groupby("mode")["pass"].agg(["size", "mean"]).round(3).to_string())
print(ok.groupby("stream")["pass"].agg(["size", "mean"]).round(3).to_string())
print("their-era BLS locked on injected period:", round(ok.his_P_ok.mean(), 3))
import collections
print(collections.Counter(x for s in ok.loc[~ok["pass"], "flags"].fillna("") for x in s.split(";") if x).most_common())
