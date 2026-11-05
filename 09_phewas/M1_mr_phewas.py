#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figure5 Panel C 数据采集: FinnGen R12 全表型 MR-PheWAS (候选基因 API)
=========================================================================
暴露 (Layer2 溯源): 胎盘 cis-eQTL 先导工具变量 (fig4_instruments.csv)
   · ERAP2  rs4869313   chr5:96888176  (Tier1 PPH4=0.925)   varid 5:96888176:T:G
   · TK2    rs12448207  chr16:66598482 (Tier1 PPH4=0.873)   varid 16:66598482:T:C
结局: FinnGen R12 全部 endpoint (manifest, 2469 endpoints)

数据源: FinnGen R12 公开 Results-API
   GET https://r12.finngen.fi/api/variant/{chr}:{pos}:{ref}:{alt}/{phenocode}
   返回该变异在该 endpoint 的 beta / pval / mlogp / maf_case / maf_control。
Layer4 对账: API 返回的 O15_PREECLAMPS(beta=-0.0578485, mlogp=3.82389,
   pval=1.50006e-4) 与本地 FinnGen R12 O15 区域汇总统计逐位一致,
   且 beta 均按 ALT 等位基因 (与暴露 ea 相同 -> 方向正号)。

方法: 单工具变量 Wald ratio MR
   b_mr = beta_out / b_exp ;  se_mr = se_out / |b_exp|
   se_out 由 beta 与 pval 精确反解: |z| = Phi^-1(1 - pval/2), se_out = |beta|/|z|
输出 -> /workspace/results/figure5/fig5C_phewas.csv
"""
import os, time, json, urllib.request, urllib.error
import threading
import numpy as np
import pandas as pd
from scipy import stats
from concurrent.futures import ThreadPoolExecutor

OUT = "/workspace/results/figure5"
MANIFEST = "/workspace/data/finngen_phewas/finngen_R12_manifest.tsv"
INST = "/workspace/results/figure4/fig4_instruments.csv"
CKPT = f"{OUT}/fig5C_phewas.partial.csv"
DONE = f"{OUT}/fig5C_done_pheno.txt"
API = "https://r12.finngen.fi/api/variant"
WORKERS = 20
STALL = 300

_OPENER = urllib.request.build_opener(
    urllib.request.ProxyHandler(urllib.request.getproxies()))

# ---- targets: varid (API) + lead IV rsid + exposure effect allele (== API alt)
TARGETS = {
    "ERAP2": dict(varid="5:96888176:T:G", rsid="rs4869313"),
    "TK2":   dict(varid="16:66598482:T:C", rsid="rs12448207"),
}
inst = pd.read_csv(INST)
for g, t in TARGETS.items():
    r = inst[(inst["symbol"] == g) & (inst["rsid"] == t["rsid"])].iloc[0]
    t["b_exp"] = float(r["b_exp"]); t["se_exp"] = float(r["se_exp"])
    t["ea"], t["oa"] = str(r["ea"]), str(r["oa"])
    chrom, pos, ref, alt = t["varid"].split(":")
    assert int(chrom) == int(r["chrom"]) and int(pos) == int(r["pos"])
    assert alt == t["ea"], f"{g}: API alt {alt} != exposure ea {t['ea']}"
    print(f"[target] {g} {t['rsid']} {t['varid']} ea={t['ea']} b_exp={t['b_exp']:.4f}",
          flush=True)

man = pd.read_csv(MANIFEST, sep="\t")
E = {r["phenocode"]: r for _, r in man.iterrows()}
print(f"[manifest] endpoints = {len(man)}", flush=True)


def _get(url):
    for a in range(4):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "curl/8"})
            with _OPENER.open(req, timeout=45) as r:
                return json.loads(r.read())
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
            time.sleep(0.6 * (a + 1))
        except Exception:
            time.sleep(0.6 * (a + 1))
    return "retry"


def _se_out(beta, pval, mlogp):
    """exact outcome SE from effect size + p-value."""
    if beta is None:
        return None
    if pval is not None and 0 < pval < 1:
        z = stats.norm.isf(pval / 2)
    elif mlogp is not None and mlogp > 0:
        z = stats.norm.isf((10 ** (-mlogp)) / 2)
    else:
        return None
    return abs(beta) / abs(z) if z > 0 else None


def _worker(pc):
    """fetch both target variants for one endpoint; return list of row dicts or 'nomatch'."""
    out = []
    for g, t in TARGETS.items():
        d = _get(f"{API}/{t['varid']}/{pc}")
        if d == "retry":
            return "error"
        res = (d or {}).get("results") if isinstance(d, dict) else None
        if not isinstance(res, dict) or res.get("beta") is None:
            continue
        beta = float(res["beta"])
        pval = res.get("pval"); mlogp = res.get("mlogp")
        se_out = _se_out(beta, pval, mlogp)
        if se_out is None or se_out <= 0:
            continue
        b_mr = beta / t["b_exp"]
        se_mr = se_out / abs(t["b_exp"])
        z = b_mr / se_mr
        p_mr = 2 * stats.norm.sf(abs(z))
        e = E[pc]
        out.append(dict(gene=g, snp=t["rsid"], phenocode=pc,
                        phenotype=e["phenotype"], category=e["category"],
                        num_cases=int(e["num_cases"]), num_controls=int(e["num_controls"]),
                        b_out=beta, se_out=se_out, p_out=(pval if pval is not None else np.nan),
                        b_mr=b_mr, se_mr=se_mr, p_mr=p_mr))
    return out if out else "nomatch"


# ---- resume
rows, done_all = [], set()
if os.path.exists(DONE):
    done_all |= {l.strip() for l in open(DONE) if l.strip()}
if os.path.exists(CKPT):
    _prev = pd.read_csv(CKPT)
    rows = _prev.to_dict("records")
    done_all |= set(_prev["phenocode"].unique())
    print(f"[resume] checkpoint rows={len(rows)}", flush=True)
print(f"[resume] endpoints already resolved = {len(done_all)}", flush=True)

man_run = man.loc[~man["phenocode"].isin(done_all)].copy()
print(f"[plan] endpoints to fetch = {len(man_run)}", flush=True)

_lock = threading.Lock()
_res = {}
_rows_prev = list(rows)
_done_prev = set(done_all)
_e_by_pc = {r["phenocode"]: r for _, r in man_run.iterrows()}
t0 = time.time()


def _materialize():
    all_rows = list(_rows_prev)
    dset = set(_done_prev)
    errs = nomatch = 0
    with _lock:
        items = list(_res.items())
    for pc, recs in items:
        if recs == "error":
            errs += 1
            continue
        dset.add(pc)
        if recs == "nomatch" or not recs:
            nomatch += 1
            continue
        all_rows.extend(recs)
    pd.DataFrame(all_rows).to_csv(CKPT, index=False)
    with open(DONE, "w") as fh:
        fh.write("\n".join(sorted(dset)))
    return all_rows, dset, errs, nomatch


def _run(pc):
    r = _worker(pc)
    with _lock:
        _res[pc] = r


print(f"[fetch] launching ThreadPool W={WORKERS} ...", flush=True)
with ThreadPoolExecutor(max_workers=WORKERS) as ex:
    futs = [ex.submit(_run, pc) for pc in man_run["phenocode"]]
    last_flush = time.time()
    while True:
        with _lock:
            k = len(_res)
        if k >= len(futs):
            break
        if time.time() - last_flush > 30:
            _, _dset, _, _ = _materialize()
            print(f"  {k}/{len(futs)} resolved={len(_dset)} t={time.time()-t0:.0f}s", flush=True)
            last_flush = time.time()
        time.sleep(2)

rows, done_all, errs, nomatch = _materialize()
df = pd.DataFrame(rows)
if len(df) == 0:
    raise SystemExit("no phewas rows collected — aborting")
print(f"[fetch] rows={len(df)} errors={errs} nomatch={nomatch} "
      f"resolved={len(done_all)}/{len(man)} wall={time.time()-t0:.0f}s", flush=True)


def bh(p):
    p = np.asarray(p, float); ok = ~np.isnan(p)
    q = np.full_like(p, np.nan); pv = p[ok]; n = len(pv); order = np.argsort(pv)
    qv = pv[order] * n / (np.arange(1, n + 1))
    qv = np.minimum.accumulate(qv[::-1])[::-1]
    out = np.empty(n); out[order] = np.clip(qv, 0, 1); q[ok] = out
    return q


df["fdr_mr"] = np.nan
for g, sub in df.groupby("gene"):
    df.loc[sub.index, "fdr_mr"] = bh(sub["p_mr"].values)
df = df.sort_values(["gene", "p_mr"])
df.to_csv(f"{OUT}/fig5C_phewas.csv", index=False)

for g, sub in df.groupby("gene"):
    sig = (sub["fdr_mr"] < 0.05).sum()
    print(f"[{g}] endpoints with MR={len(sub)}; FDR<0.05 = {sig}; "
          f"top: {sub.iloc[0]['phenotype']} (p={sub.iloc[0]['p_mr']:.2e})", flush=True)
