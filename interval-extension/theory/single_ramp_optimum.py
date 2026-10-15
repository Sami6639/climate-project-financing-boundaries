#!/usr/bin/env python3
"""Analytical global optimum for a uniform-noise, interval-wide single ramp.

The scalar candidate theorem is proved in single_ramp_theorem.md. Standard
library only. Negative algebraic capacity is retained and labelled infeasible.
All money is in the same currency; all emissions use the same physical unit.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
from pathlib import Path


def validate_inputs(*, h, delta, p, q, B, U, kappa=0.0, V=0.0, mu=0.0):
    vals = locals().copy()
    if not all(math.isfinite(float(v)) for v in vals.values()):
        raise ValueError("All inputs must be finite")
    if h <= 0 or p <= 0 or delta < 0 or q < 0 or kappa < 0 or V < 0:
        raise ValueError("Require h>0, p>0, delta,q,kappa,V>=0")


def ramp_reserves(h, delta, q, r):
    """Pointwise dominating ramp at fixed r, in physical units."""
    if q <= 0 or r < q:
        raise ValueError("Positive q and r>=q required")
    W, L = 2.0*h, 2.0*delta
    s = W*q/r
    a = W-s
    width = max(s, L+2.0*s-W)
    K = r*width
    if L <= W-s:
        R = q*L + W*q*q/(2.0*r)
        branch = "short"
    else:
        R = r*(L-W) + 2.0*W*q - W*q*q/(2.0*r)
        branch = "long"
    return {"r": r, "s": s, "relative_start": a, "width": width,
            "K": K, "R": R, "branch": branch}


def binding(B, U, R, K):
    x, y = B-R, U-K
    tol = 2e-12*max(1.0, abs(B), abs(U), abs(R), abs(K))
    if abs(x-y) <= tol:
        return "both"
    return "expected_budget" if x < y else "payment_span"


def class_frontier(h, delta, p, q, B, U, kappa=0.0, V=0.0):
    if q > p:
        return {"feasible": False, "reason": "q exceeds p"}
    W, L = 2*h, 2*delta
    n = math.floor(L/W)
    ell = L-n*W
    D = (n+1)*W*q + max(0.0, W*q-p*(W-ell))
    R = h*q*q/p + L*q
    floor = min(B-R, U-D)
    H = (floor-V)/(1+kappa)
    return {"D": D, "R": R, "n": n, "ell": ell,
            "H_algebraic": H, "feasible": H >= 0,
            "binding": binding(B,U,R,D)}


def conservative_frontier(h, delta, p, q, B, U, kappa=0.0, V=0.0):
    if q > p:
        return {"feasible": False, "reason": "q exceeds p"}
    K, R = 2*q*(h+delta), q*(h+2*delta)
    H = (min(B-R,U-K)-V)/(1+kappa)
    return {"K": K, "R": R, "r": q, "H_algebraic": H,
            "feasible": H >= 0, "binding": binding(B,U,R,K)}


def solve_single_ramp(*, h, delta, p, q, B, U, kappa=0.0, V=0.0, mu=0.0):
    """Return a finite-candidate certificate for the entire single-ramp class."""
    validate_inputs(h=h,delta=delta,p=p,q=q,B=B,U=U,kappa=kappa,V=V,mu=mu)
    if q > p:
        return {"feasible": False, "reason": "q exceeds p", "candidates": []}
    if q == 0:
        floor = min(B,U)
        H = (floor-V)/(1+kappa)
        return {"feasible": H>=0, "H_algebraic": H, "c": floor,
                "r": 0.0, "K": 0.0, "R": 0.0, "t": None,
                "binding": binding(B,U,0,0), "candidates": [],
                "certificate": "constant contract; zero target"}
    W, L = 2*h, 2*delta
    slopes = [("q_endpoint",q)]
    if p != q:
        slopes.append(("p_endpoint",p))
    if U > B:
        rc = W*q*q/(2*(U-B))
        if q <= rc <= p and L >= W-W*q/rc:
            if not any(rc == r for _,r in slopes):
                slopes.append(("interior_crossing",rc))
    candidates=[]
    for label,r in slopes:
        c = ramp_reserves(h,delta,q,r)
        c.update(label=label,t=mu-delta-h+c["relative_start"])
        c["floor_capacity"] = min(B-c["R"],U-c["K"])
        c["H_algebraic"] = (c["floor_capacity"]-V)/(1+kappa)
        c["binding"] = binding(B,U,c["R"],c["K"])
        candidates.append(c)
    best=max(candidates,key=lambda c:c["floor_capacity"])
    out={k:best[k] for k in ["r","K","R","t","width","relative_start","branch","binding","label","H_algebraic"]}
    out.update(feasible=best["H_algebraic"]>=0,c=best["floor_capacity"],
               candidates=candidates,certificate="finite analytical candidate theorem")
    return out


def zero_gap_region(h, delta, p, q, B, U):
    """Mathematical equality conditions with roundoff-only boundary tolerance.

    Monetary comparisons use 2e-13 times the relevant input scale; width
    comparisons use 2e-13. Nearby substantive scenarios are tested separately.
    """
    if q > p:
        return False
    if q == 0 or q == p:
        return True
    lam, alpha, W = delta/h, q/p, 2*h
    eps=2e-13
    money_eps=eps*max(1.0,abs(B),abs(U),W*q)
    if lam <= 1-alpha+eps:
        return True
    if lam < 1-eps:
        return U-B <= W*q*q/(2*p)+money_eps
    if abs(lam-round(lam))<=eps:
        return U-B <= W*q/2+money_eps
    return False


def comparison_row(label, **x):
    ramp=solve_single_ramp(**x)
    other={k:v for k,v in x.items() if k!="mu"}
    glob=class_frontier(**other)
    cons=conservative_frontier(**other)
    row={"scenario":label,"delta_over_h":x["delta"]/x["h"],"q_over_p":x["q"]/x["p"],
         **x,"global_H":glob.get("H_algebraic"),"ramp_H":ramp.get("H_algebraic"),
         "conservative_H":cons.get("H_algebraic"),"global_binding":glob.get("binding"),
         "ramp_binding":ramp.get("binding"),"optimal_r":ramp.get("r"),
         "ramp_K":ramp.get("K"),"ramp_R":ramp.get("R"),"ramp_t":ramp.get("t"),
         "ramp_feasible":ramp["feasible"]}
    if row["ramp_H"] is not None:
        row["global_minus_ramp"]=row["global_H"]-row["ramp_H"]
        row["global_minus_conservative"]=row["global_H"]-row["conservative_H"]
        row["predicted_zero_gap"]=zero_gap_region(x["h"],x["delta"],x["p"],x["q"],x["B"],x["U"])
    return row


def generate_results(out_dir):
    mu=316912.0
    base={"h":math.sqrt(3)*0.1*mu,"p":25.0,"q":18.75,"B":7925000.0,
          "U":1.1*7925000.0,"kappa":0.0,"V":0.0,"mu":mu}
    rows=[]
    for lam in [0,0.1,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25,3]:
        rows.append(comparison_row(f"central_lambda_{lam:g}",delta=lam*base["h"],**base))
    # Normalised money scale is B=20, W=2, p=1; B is large enough to
    # isolate differences without accidentally crossing zero-advance feasibility.
    for alpha in [0,0.1,0.25,0.5,0.75,0.9,1.0]:
        for lam in [0,0.1,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25,3]:
            for gap in [-0.25,0,0.1,0.25,0.5,0.75]:
                x={"h":1.0,"delta":lam,"p":1.0,"q":alpha,"B":20.0,
                   "U":20.0+gap*2*alpha,"kappa":0.0,"V":0.0,"mu":10.0}
                rows.append(comparison_row(f"normalised_a{alpha:g}_l{lam:g}_g{gap:g}",**x))
    # A strictly interior optimal slope: normalized slope r/q=1.2,
    # lambda=.75, alpha=.75, (U-B)/(Wq)=5/12.
    interior={"h":1.,"delta":.75,"p":1.,"q":.75,"B":20.,"U":20.625,
              "kappa":0.,"V":0.,"mu":10.}
    rows.append(comparison_row("interior_crossing_both_bind",**interior))
    payload={"model":"fixed-width uniform reporting; same ramp at every centre",
             "analytical_certificate":"single_ramp_theorem.md",
             "central_inputs":base,"central_delta_ratios":[0,.1,.25,.5,.75,1,1.25,1.5,1.75,2,2.25,3],
             "interior_crossing_inputs":interior,
             "interior_crossing_solution":solve_single_ramp(**interior),"rows":rows}
    out_dir=Path(out_dir);out_dir.mkdir(parents=True,exist_ok=True)
    (out_dir/"single_ramp_results.json").write_text(json.dumps(payload,indent=2)+"\n")
    with (out_dir/"single_ramp_comparison.csv").open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    central=rows[:12]
    with (out_dir/"central_three_way_comparison.csv").open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(central[0]));writer.writeheader();writer.writerows(central)
    return rows,payload


def self_check(rows,payload):
    # Closed-form algebraic checks, not a substitute for the independent validator.
    failed=[]
    tol=2e-10
    max_order_error=0.0
    for row in rows:
        scale=max(1,abs(row["B"]),abs(row["U"]))
        g,r,c=row["global_H"],row["ramp_H"],row["conservative_H"]
        err=max(0.,r-g,c-r)/scale
        max_order_error=max(max_order_error,err)
        if err>tol: failed.append((row["scenario"],"frontier ordering",err))
        actual_zero=abs(g-r)<=tol*scale
        if actual_zero != row["predicted_zero_gap"]:
            failed.append((row["scenario"],"equality region",g-r))
        if row["delta_over_h"]>=1 and abs(r-c)>tol*scale:
            failed.append((row["scenario"],"wide-interval ramp optimum",r-c))
    interior=payload["interior_crossing_solution"]
    if abs(interior["r"]-.9)>1e-12 or interior["binding"]!="both":
        failed.append(("interior_crossing","wrong slope or binding",interior))
    # q=0 constant case and nonnegative-advance infeasibility retain their meaning.
    x=dict(h=1.,delta=2.,p=1.,q=0.,B=2.,U=3.,V=.5,kappa=.5,mu=10.)
    z=solve_single_ramp(**x)
    if z["H_algebraic"]!=1 or z["K"]!=0:failed.append(("q0","constant failure",z))
    x.update(q=1.,B=.1,U=.1,V=0.)
    z=solve_single_ramp(**x)
    if z["feasible"] or z["H_algebraic"]>=0:failed.append(("infeasible","negative truncated",z))
    report={"status":"PASS" if not failed else "FAIL","row_count":len(rows),
            "tolerance_scaled":tol,"maximum_scaled_order_error":max_order_error,
            "random_seed":None,"failures":failed,
            "scope":"Algebraic implementation smoke checks; independent continuum validation is separate."}
    return report


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--output-dir",default=str(Path(__file__).parent))
    args=parser.parse_args();rows,payload=generate_results(args.output_dir)
    report=self_check(rows,payload)
    (Path(args.output_dir)/"single_ramp_self_check.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report,indent=2))
    raise SystemExit(0 if report["status"]=="PASS" else 1)
