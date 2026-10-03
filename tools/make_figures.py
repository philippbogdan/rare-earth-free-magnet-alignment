#!/usr/bin/env python3
"""Plot committed results only. Run from any directory with Python + Matplotlib.

No fitting, smoothing, resampling or outcome exclusions. The MnAl schematic is
derived from the stored alignment ceiling, not an image of a simulated specimen.
All five validation batches remain separate. See figures/provenance.json.
"""
import csv
import hashlib
import json
import math
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Arc, Patch

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "figures"
BLUE, ORANGE, GREY, INK = "#0072B2", "#B34E00", "#686868", "#202B33"
STYLE = {
    "font.family": "DejaVu Sans", "font.size": 11,
    "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK, "ytick.color": INK,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.edgecolor": "#808890", "axes.titleweight": "bold",
    "svg.fonttype": "none", "svg.hashsalt": "magnet-alignment",
    "savefig.facecolor": "white", "figure.facecolor": "white",
}


def save(fig, name):
    fig.savefig(OUT / f"{name}.png", dpi=180, facecolor="white")
    fig.savefig(OUT / f"{name}.svg", facecolor="white",
                metadata={"Date": None})
    plt.close(fig)


def mnal(rows):
    row = next(r for r in rows if (r["phase"], r["route"]) == ("tau-MnAl-C", "AM"))
    alternatives = json.loads(row["A_alt"])
    cap, random = alternatives["OR kept"], alternatives["random"]
    angle = math.degrees(math.acos(cap))
    fig = plt.figure(figsize=(10.8, 4.8))
    fig.suptitle("MnAl(C): inherited texture still limits alignment", x=.05,
                 y=.96, ha="left", fontsize=17, weight="bold")
    fig.text(.05, .885, "Ideal parent c-axis texture; no field or stress bias during transformation", fontsize=11)
    left = fig.add_axes([.04, .21, .36, .59])
    left.set(xlim=(-1.15, 1.2), ylim=(-.10, 1.30), aspect="equal")
    left.axis("off")
    left.annotate("", (0, 1.16), (0, 0), arrowprops={"arrowstyle": "->", "color": INK, "lw": 2})
    left.text(.06, 1.17, "Build direction\n(parent c-axis)", fontsize=10, va="top")
    x, z = math.sin(math.radians(angle)), cap
    for sign in [-1, 1]:
        left.annotate("", (sign*x, z), (0, 0),
                      arrowprops={"arrowstyle": "->", "color": BLUE, "lw": 2.5})
    left.plot([-x, x], [z, z], color=GREY, ls=":", lw=1.2)
    left.add_patch(Arc((0, 0), .86, .86, theta1=90-angle, theta2=90, color=INK))
    left.text(.20, .47, f"{angle:.1f}°", fontsize=11)
    left.text(-1.10, .81, "Product easy axes", color=BLUE, fontsize=11)
    left.text(0, -.09, "Section through the easy-axis cone (schematic)",
              ha="center", va="top", fontsize=9)
    ax = fig.add_axes([.61, .29, .31, .47])
    values = [random, cap, 1.0]
    colours = [GREY, BLUE, INK]
    for y, (v, col) in enumerate(zip(values, colours)):
        ax.barh(y, v, height=.44, color=col)
        ax.text(v+.025, y, f"{v:.2f}", va="center", fontsize=12, weight="bold")
    ax.set(yticks=range(3), yticklabels=["Random grains", "Inherited texture", "Perfect alignment"],
           xlim=(0, 1.15), xticks=[0, .5, 1], xlabel="Ideal alignment factor, A")
    ax.invert_yaxis()
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0, pad=8)
    fig.text(.05, .075, "A = mean |cos θ|. For independent uniaxial grains, A is the ideal remanence / saturation ratio.", fontsize=10)
    fig.text(.05, .032, "Equal variant weights; perfect orientation inheritance. Calculated limits, not measured magnet performance.", fontsize=10)
    save(fig, "mnal-texture")
    return {"A_random": random, "A_inherited": cap, "cone_angle_deg": angle}


def route_map(rows):
    phases = list(dict.fromkeys(r["phase"] for r in rows))
    routes = list(dict.fromkeys(r["route"] for r in rows))
    counts = Counter(r["verdict"] for r in rows)
    symbols = {"full": ("s", BLUE), "capped": ("o", ORANGE),
               "none": ("x", GREY), "open": ("$?$", INK)}
    names = {"tau-MnAl": "τ-MnAl", "tau-MnAl-C": "τ-MnAl-C",
             "AlNiCo-5": "AlNiCo 5 (Ti-free)", "AlNiCo-8-9": "AlNiCo 8/9",
             "Nd2Fe14B": "Nd2Fe14B (reference)",
             "HEA-FeCoNiAlCuTi": "FeCoNiAlCuTi alloy", "L10-HEA-Pt": "L10 Pt-based alloy"}
    fig = plt.figure(figsize=(11.6, 12.6))
    fig.suptitle("Which routes allow alignment in the model?", x=.04, y=.976,
                 ha="left", fontsize=18, weight="bold")
    fig.text(.04, .944, f"{len(phases)} phase entries × {len(routes)} routes. Each cell retains its assumptions and proposed test in account.csv.", fontsize=11)
    ax = fig.add_axes([.235, .175, .735, .72])
    for i in range(len(phases)):
        if i % 2 == 0:
            ax.axhspan(i-.5, i+.5, color="#F1F4F6", zorder=0)
    for verdict, (marker, colour) in symbols.items():
        selected = [r for r in rows if r["verdict"] == verdict]
        ax.scatter([routes.index(r["route"]) for r in selected],
                   [phases.index(r["phase"]) for r in selected],
                   marker=marker, c=colour, s=55, linewidths=1.3)
    ax.set(xlim=(-.6, len(routes)-.4), ylim=(len(phases)-.5, -.5),
           xticks=range(len(routes)), xticklabels=routes,
           yticks=range(len(phases)), yticklabels=[names.get(p, p) for p in phases])
    ax.xaxis.tick_top()
    ax.tick_params(length=0, axis="both", pad=8, labelsize=10)
    plt.setp(ax.get_xticklabels(), rotation=55, ha="left", rotation_mode="anchor")
    for s in ax.spines.values():
        s.set_visible(False)
    labels = {"full": f"A ceiling ≥ 0.95 ({counts['full']})",
              "capped": f"Geometry caps alignment ({counts['capped']})",
              "none": f"No gain from route ({counts['none']})",
              "open": f"Unresolved input ({counts['open']})"}
    handles = [Line2D([], [], color=c, marker=m, ls="", markersize=8, label=labels[v])
               for v, (m, c) in symbols.items()]
    fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(.03, .113),
               frameon=False, ncol=2, fontsize=11, columnspacing=2)
    fig.text(.04, .086, "FP  field-aligned particles   MO  shape-aligned particles   FS  field solidification   DS  directional solidification", fontsize=9.5)
    fig.text(.04, .065, "AM  additive manufacturing   FA  field ageing   SA  stress ageing   RX  rolled/recrystallised parent   SC  single crystal", fontsize=9.5)
    fig.text(.04, .044, "HD  hot deformation   HDDR  hydrogen disproportionation/recombination   EP  epitaxy (films)   TGG  templated grain growth", fontsize=9.5)
    fig.text(.04, .021, "These are conditional model limits. A high ceiling does not establish a viable bulk process or sufficient coercivity.", fontsize=10)
    save(fig, "route-map")


def validation():
    versions = [f"v{i}" for i in range(1, 6)]
    summaries = [json.loads((ROOT / f"results/scores_{v}.json").read_text())["summary"]["heldout"] for v in versions]
    keys = ["ceiling_ok", "ceiling_possible", "ceiling_violated"]
    colours, hatches = [BLUE, "#F0D5BC", ORANGE], ["", "///", "xx"]
    labels = ["Within ceiling", "Conversion range straddles", "Clear violation"]
    fig = plt.figure(figsize=(10.8, 5.0))
    fig.suptitle("Published measurements tested against frozen ceilings", x=.05,
                 y=.965, ha="left", fontsize=17, weight="bold")
    fig.text(.05, .885, "Each version uses its own held-out batch. Counts are measurements, not independent studies.", fontsize=11)
    ax = fig.add_axes([.16, .30, .51, .43])
    for y, s in enumerate(summaries):
        start = 0
        for key, colour, hatch in zip(keys, colours, hatches):
            width = 100*s[key]/s["n_scored"]
            ax.barh(y, width, left=start, height=.56, color=colour, hatch=hatch,
                    edgecolor=INK, linewidth=.5)
            start += width
        assert sum(s[k] for k in keys) == s["n_scored"]
        ax.text(107, y, " / ".join(str(s[k]) for k in keys), va="center", fontsize=11)
    ax.text(107, -.85, "Within / straddles / violates", fontsize=10, weight="bold")
    ax.set(yticks=range(5), yticklabels=[f"{v}   n = {s['n_scored']}" for v,s in zip(versions,summaries)],
           xlim=(0, 100), xticks=[0, 25, 50, 75, 100], xlabel="Share of scored measurements (%)")
    ax.invert_yaxis()
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    fig.legend(handles=[Patch(facecolor=c, hatch=h, edgecolor=INK, label=l)
                        for c,h,l in zip(colours,hatches,labels)],
               loc="lower left", bbox_to_anchor=(.04,.115), ncol=3, frameon=False, fontsize=10)
    fig.text(.05, .07, "Rule revisions are tested on later batches; the batches differ in size and composition. This is not a learning curve.", fontsize=10)
    fig.text(.05, .025, "Scoring tolerance: max(0.03 in A, twice stated uncertainty). Conversion ranges are not confidence intervals.", fontsize=10)
    save(fig, "validation")


def main():
    OUT.mkdir(exist_ok=True)
    inputs = [ROOT / "results/account.csv", ROOT / "results/stamp.json"]
    inputs += [ROOT / f"results/scores_v{i}.json" for i in range(1, 6)]
    with inputs[0].open(newline="") as f:
        rows = list(csv.DictReader(f))
    assert len({(r["phase"], r["route"]) for r in rows}) == len(rows)
    with plt.rc_context(STYLE):
        geometry = mnal(rows)
        route_map(rows)
        validation()
    manifest = {
        "result_stamp": json.loads(inputs[1].read_text()),
        "inputs_sha256": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "matplotlib": matplotlib.__version__,
        "destination": "GitHub README, web; MnAl figure also suitable as an email attachment",
        "exports": "PNG at 180 dpi plus SVG with editable text; opaque white background",
        "mnal_geometry": geometry,
        "transformations": [
            "MnAl-C AM: select A_alt random and OR kept; perfect alignment is the definitional reference A=1",
            "Schematic angle = acos(A_alt[OR kept]); representative section of the cone, not raw orientation samples",
            "Route map includes every account.csv row in file order; symbols encode stored verdicts",
            "Validation uses heldout summaries only, one bar per version; counts divided by n_scored",
        ],
        "missing_and_uncertainty": "Open route cells shown as question marks. No invented error bars. Conversion-straddling scores shown separately. No averaging across batches.",
        "scope": "Alignment ceilings under idealised assumptions; no new magnet experiments or process-feasibility validation",
    }
    (OUT / "provenance.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Wrote three PNG/SVG pairs and provenance to {OUT}")


if __name__ == "__main__":
    main()
