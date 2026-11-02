"""Render four static alignment figures from the included account and figure CSVs.

Run python tools/render_figures.py from the repository root.
Write PNGs to figures/. See REPRODUCE.md for commands and fonts.
"""

import argparse
import csv
import io
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Polygon, Rectangle
from matplotlib.patheffects import withStroke
from matplotlib.ticker import FormatStrFormatter
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "figures"
DATA = FIG / "data"
RESULTS = ROOT / "results"

BG = "#FAF8F3"
INK = "#2B2A28"
MUTED = "#77726A"
RULE = "#BDB7AC"
TEAL = "#00786E"
OCHRE = "#A77B1F"
VERMILION = "#D2361C"
NONE_GREY = "#D9D5CC"
HATCH = "#8F897F"
WHITE = "#FFFFFF"
GUIDE = "#F0EDE5"

WIDTH = 10.0
DPI = 200

PHASES = [
    "tau-MnAl", "tau-MnAl-C", "LTP-MnBi", "L10-MnGa", "D022-Mn3Ga", "D022-Mn3Ge",
    "Mn4N", "AlNiCo-5", "AlNiCo-8-9", "Fe-Cr-Co", "HEA-FeCoNiAlCuTi",
    "SrM-hexaferrite", "eps-Fe2O3", "Fe16N2", "L10-FeNi", "bct-FeCo", "(Fe,Co)2B",
    "Fe2P-based", "Fe5PB2", "Fe5SiB2", "Fe2Ta-C14", "Fe3Se4", "Zr2Co11", "HfCo7",
    "Co-nanorods", "Co3C", "L10-FePt", "L10-CoPt", "L10-FePd", "L10-HEA-Pt",
    "Nd2Fe14B",
]
ROUTES = [
    "FP", "MO", "FS", "DS", "AM", "FA", "DS+FA", "AM+FA", "RX+FA", "SA", "DS+SA",
    "RX+SA", "HD", "HDDR", "TGG", "SC+FA", "SC+SA", "EP",
]
REFERENCE = "Nd2Fe14B"
FILL = {"full": TEAL, "capped": OCHRE, "none": NONE_GREY, "open": WHITE}
EDGE = {"none": dict(edgecolor=HATCH, linewidth=1.0)}
CATEGORY = {"ok": TEAL, "possible violation": OCHRE, "violated": VERMILION}

MERIDIANS = 12
AZ0, EL0 = 265.0, 26.0
SELECTED_VARIANT = 0
MARK_K = 34000.0
CENTRES = (-1.4, 1.0)
XLIM = (-2.48, 3.12)
YLIM = (-0.8, -0.8 + 5.6 * 0.386)


def style():
    plt.rcParams.update({
        "font.family": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"],
        "font.size": 11,
        "text.color": INK,
        "axes.labelcolor": INK,
        "axes.edgecolor": MUTED,
        "axes.linewidth": 0.7,
        "xtick.color": INK,
        "ytick.color": INK,
        "xtick.major.width": 0.7,
        "ytick.major.width": 0.7,
        "xtick.major.size": 3.5,
        "ytick.major.size": 3.5,
        "figure.facecolor": BG,
        "axes.facecolor": BG,
        "savefig.facecolor": BG,
        "hatch.linewidth": 1.0,
    })


def save(fig, name):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=DPI)
    buf.seek(0)
    img = Image.open(buf).convert("RGBA")
    flat = Image.alpha_composite(Image.new("RGBA", img.size, BG), img).convert("RGB")
    flat.save(FIG / name, optimize=True, dpi=(DPI, DPI))


def table(path):
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def population(name):
    """Read axial endpoints and normalise each population to equal parent weight."""
    d = np.genfromtxt(DATA / name, delimiter=",", names=True)
    xyz = np.column_stack([d["x"], d["y"], d["z"]])
    parent = d["parent"].astype(int)
    variant = d["variant"].astype(int)
    assert (xyz[:, 2] >= 0).all()
    assert np.allclose(np.linalg.norm(xyz, axis=1), 1, atol=1e-5)
    w = d["weight"] / np.unique(parent).size
    return {"xyz": xyz, "w": w, "parent": parent, "variant": variant}


def hero_data():
    """Check the selected parent's variants and geometric cap against the summary."""
    random = population("mnal_c_axes_random_parent.csv")
    fibre = population("mnal_c_axes_fibre_parent.csv")
    pick = fibre["parent"] == 0
    ends = fibre["xyz"][pick][np.argsort(fibre["variant"][pick])]
    assert len(ends) == 6 and np.ptp(ends[:, 2]) < 1e-5
    z0 = float(ends[:, 2].mean())
    angle = float(np.degrees(np.arccos(z0)))
    spread = {r["parent_fibre_spread_deg"]: r["A_along_fibre"] for r in table(DATA / "mnal_c_A_vs_parent_spread.csv")}
    summary = json.loads((RESULTS / "summary.json").read_text())["mnal_c_textured_parent"]
    assert f"{angle:.1f}" == f"{summary['easy_axis_angle_to_fibre_deg']:.1f}"
    assert abs(float(spread["0"]) - summary["A_max_AM"]) < 1e-4
    assert abs(float(spread["random"]) - summary["A_random"]) < 1e-3
    return {"random": random, "fibre": fibre, "ends": ends, "z0": z0, "angle": angle}


def camera(az, el):
    """Return an orthographic camera basis for angles in degrees."""
    a, e = np.radians(az), np.radians(el)
    right = np.array([-np.sin(a), np.cos(a), 0.0])
    up = np.array([-np.sin(e) * np.cos(a), -np.sin(e) * np.sin(a), np.cos(e)])
    view = np.array([np.cos(e) * np.cos(a), np.cos(e) * np.sin(a), np.sin(e)])
    return right, up, view


def guide(ax, cx, cam):
    right, up, view = cam
    t = np.linspace(0, np.pi, 181)
    rim = np.column_stack([np.cos(t), np.sin(t)])
    a = np.arctan2(view[1], view[0])
    psi = np.linspace(-np.pi / 2, np.pi / 2, 181)
    eq = np.column_stack([np.cos(a + psi), np.sin(a + psi), np.zeros_like(psi)])
    front = np.column_stack([eq @ right, eq @ up])
    ax.add_patch(Polygon(np.vstack([rim, front]) + [cx, 0], closed=True, facecolor=GUIDE,
                         edgecolor="none", zorder=1))


def graticule():
    t = np.linspace(0, np.pi / 2, 61)
    lines = []
    for k in range(MERIDIANS):
        ph = 2 * np.pi * k / MERIDIANS
        lines.append(np.column_stack([np.sin(t) * np.cos(ph), np.sin(t) * np.sin(ph), np.cos(t)]))
    return lines


def wire(ax, cx, cam, scale):
    right, up, view = cam
    front, back = [], []
    for line in graticule():
        xy = np.column_stack([line @ right + cx, line @ up])
        depth = line @ view
        for i in range(len(line) - 1):
            (front if depth[i] + depth[i + 1] >= 0 else back).append(xy[i:i + 2])
    ax.add_collection(LineCollection(back, colors="#E2DDD3", linewidths=0.5 * scale, zorder=2))
    ax.add_collection(LineCollection(front, colors="#CFC9BE", linewidths=0.55 * scale, zorder=2.2))
    ph = np.linspace(0, 2 * np.pi, 361)
    eq = np.column_stack([np.cos(ph), np.sin(ph), np.zeros_like(ph)])
    xy = np.column_stack([eq @ right + cx, eq @ up])
    depth = eq @ view
    segs = [xy[i:i + 2] for i in range(len(ph) - 1)]
    seen = [depth[i] + depth[i + 1] >= 0 for i in range(len(ph) - 1)]
    ax.add_collection(LineCollection([s for s, f in zip(segs, seen) if not f],
                                     colors="#B3AC9F", linewidths=0.6 * scale, zorder=2.1))
    ax.add_collection(LineCollection([s for s, f in zip(segs, seen) if f],
                                     colors=MUTED, linewidths=0.85 * scale, zorder=2.4))
    t = np.linspace(0, np.pi, 181)
    ax.plot(cx + np.cos(t), np.sin(t), color=MUTED, lw=0.85 * scale, zorder=2.4,
            solid_capstyle="round")


def dots(ax, cx, pop, cam, scale):
    right, up, view = cam
    xyz = pop["xyz"]
    ax.scatter(xyz @ right + cx, xyz @ up, s=MARK_K * pop["w"] * scale ** 2, color=INK,
               linewidths=0, zorder=7)


def cone(ax, cx, hd, cam, scale):
    right, up, view = cam
    ends, z0 = hd["ends"], hd["z0"]

    def flat(p):
        p = np.atleast_2d(p)
        return p @ right + cx, p @ up

    top = 1.18
    x, y = flat(np.array([[0, 0, 0], [0, 0, top]]))
    ax.plot(x, y, color=MUTED, lw=1.0 * scale, zorder=4, solid_capstyle="butt")
    foot = np.array([0, 0, z0])
    for e in ends:
        x, y = flat(np.array([foot, e]))
        ax.plot(x, y, color=INK, lw=0.9 * scale, ls=(0, (1.0, 2.2)), zorder=4,
                dash_capstyle="round")
        x, y = flat(np.array([[0, 0, 0], e]))
        ax.plot(x, y, color=OCHRE, lw=1.5 * scale, zorder=4.2, solid_capstyle="round")
    fx, fy = flat(foot)
    ax.scatter(fx, fy, s=9 * scale ** 2, color=INK, linewidths=0, zorder=4.3)
    sel = ends[SELECTED_VARIANT]
    phi = np.arctan2(sel[1], sel[0])
    theta = np.arccos(z0)
    t = np.linspace(0, theta, 60)
    rho = 0.3
    arc = rho * np.column_stack([np.sin(t) * np.cos(phi), np.sin(t) * np.sin(phi), np.cos(t)])
    x, y = flat(arc)
    ax.plot(x, y, color=INK, lw=1.0 * scale, zorder=4.4)
    mid = 0.5 * np.array([np.sin(theta / 2) * np.cos(phi), np.sin(theta / 2) * np.sin(phi), np.cos(theta / 2)])
    x, y = flat(mid)
    ax.text(x[0], y[0], f"{hd['angle']:.1f}°", ha="center", va="center", fontsize=12.5,
            zorder=9, path_effects=[withStroke(linewidth=3.2, foreground=GUIDE)])
    _, ty = flat(np.array([0, 0, top]))
    ax.text(cx, ty[0] + 0.04, "Parent fibre", ha="center", va="bottom", fontsize=12.5, zorder=9)
    r0 = np.sqrt(1 - z0 ** 2)
    ly = z0 * up[2]
    ax.plot([cx + r0 + 0.015, cx + 1.1], [ly, ly], color=INK, lw=0.7 * scale, zorder=9)
    ax.text(cx + 1.14, ly, "Product easy axes", ha="left", va="center", fontsize=12.5, zorder=9)


def hero_frame(fig, hd, az, el, dpi):
    fig.clf()
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(*XLIM)
    ax.set_ylim(*YLIM)
    ax.axis("off")
    cam = camera(az, el)
    scale = 1.0
    for cx, key, label in zip(CENTRES, ("random", "fibre"), ("Random parent", "Perfect parent fibre")):
        pop = hd[key]
        guide(ax, cx, cam)
        wire(ax, cx, cam, scale)
        dots(ax, cx, pop, cam, scale)
        ax.text(cx, -0.62, label, ha="center", va="top", fontsize=14, zorder=9)
    cone(ax, CENTRES[1], hd, cam, scale)
    return ax


def hero_figure():
    return plt.figure(figsize=(WIDTH, WIDTH * 0.386))


def hero_still(hd):
    fig = hero_figure()
    hero_frame(fig, hd, AZ0, EL0, DPI)
    save(fig, "mnal-inheritance.png")
    plt.close(fig)


def atlas():
    """Draw every account cell, separating the rare-earth reference and specialist routes."""
    rows = table(RESULTS / "account.csv")
    verdict = {}
    for r in rows:
        key = (r["phase"], r["route"])
        assert key not in verdict
        verdict[key] = r["verdict"]
        a = float(r["A_max"])
        if r["verdict"] == "full":
            assert a >= 0.95
        if r["verdict"] == "capped":
            assert a < 0.95
        if r["verdict"] == "open":
            assert r["decides"] not in ("", "[]")
    assert len(rows) == len(PHASES) * len(ROUTES) == 558
    assert {p for p, _ in verdict} == set(PHASES) and {q for _, q in verdict} == set(ROUTES)
    assert len(set(PHASES)) == len(PHASES) and len(set(ROUTES)) == len(ROUTES)

    col_gap, row_gap = 0.75, 0.9
    split = ROUTES.index("SC+FA")
    xs = {q: j + (col_gap if j >= split else 0) for j, q in enumerate(ROUTES)}
    ys = {p: i + (row_gap if p == REFERENCE else 0) for i, p in enumerate(PHASES)}
    ncols = len(ROUTES) + col_gap
    nrows = len(PHASES) + row_gap

    left, right, top, bottom = 1.9, 0.12, 0.62, 0.62
    col_in = (WIDTH - left - right) / ncols
    row_in = 0.25
    height = top + bottom + nrows * row_in
    fig = plt.figure(figsize=(WIDTH, height))
    ax = fig.add_axes([left / WIDTH, bottom / height, ncols * col_in / WIDTH, nrows * row_in / height])
    ax.set_xlim(0, ncols)
    ax.set_ylim(nrows, 0)
    ax.axis("off")
    gx, gy = 0.035 / col_in, 0.035 / row_in
    for (p, q), v in verdict.items():
        x, y = xs[q], ys[p]
        kw = dict(facecolor=FILL[v], edgecolor="none", linewidth=0)
        kw.update(EDGE.get(v, {}))
        if v == "open":
            kw.update(hatch="/" * 5, hatchcolor=HATCH, edgecolor="#C4BEB3", linewidth=0.5)
        ax.add_patch(Rectangle((x + gx, y + gy), 1 - 2 * gx, 1 - 2 * gy, **kw))
    for q in ("AM", "DS"):
        x, y = xs[q], ys["tau-MnAl-C"]
        ax.add_patch(Rectangle((x + gx * 0.4, y + gy * 0.4), 1 - gx * 0.8, 1 - gy * 0.8,
                               fill=False, edgecolor=INK, linewidth=1.1, zorder=3))
    for p in PHASES:
        ax.text(-0.18, ys[p] + 0.5, f"{p} (reference)" if p == REFERENCE else p,
                ha="right", va="center", fontsize=11)
    for q in ROUTES:
        ax.text(xs[q] + 0.5, -0.22, q, ha="left", va="bottom", rotation=45,
                rotation_mode="anchor", fontsize=11)
    handles = [
        Patch(facecolor=TEAL, edgecolor="none", label="Full"),
        Patch(facecolor=OCHRE, edgecolor="none", label="Capped"),
        Patch(facecolor=NONE_GREY, label="None", **EDGE["none"]),
        Patch(facecolor=WHITE, hatch="/" * 5, hatchcolor=HATCH, edgecolor="#C4BEB3",
              linewidth=0.5, label="Open"),
    ]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.012), ncol=4,
              frameon=False, fontsize=11.5, handlelength=2.6, handleheight=1.3,
              columnspacing=2.2, handletextpad=0.6)
    save(fig, "alignment-atlas.png")
    plt.close(fig)


def heldout_rows():
    """Return every scored v1 held-out case in plotting order, retaining its case ID.

    Interval endpoints are measured alignment minus ceiling minus tolerance.
    Sort by ceiling, then relative interval midpoint, then case ID. Check the
    exported values and categories against results/scores_v1.json before drawing.
    """
    rows = [r for r in table(DATA / "measurements_vs_ceiling.csv")
            if r["version"] == "v1" and r["status"] == "heldout"]
    for r in rows:
        a, lo, hi, tol = (float(r[k]) for k in ("A_ceiling", "A_measured_lo", "A_measured_hi", "tol"))
        r["A"], r["lo"], r["hi"] = a, lo - a - tol, hi - a - tol
        r["mid"] = (r["lo"] + r["hi"]) / 2
    rows.sort(key=lambda r: (r["A"], r["mid"], r["case"]))
    scores = json.loads((RESULTS / "scores_v1.json").read_text())
    scored = {c["case"]: c for c in scores["cases"]
              if c["status"] == "heldout" and c["result"] == "scored"}
    assert len(scored) == len(rows) == len({r["case"] for r in rows})
    for r in rows:
        c = scored[r["case"]]
        assert c["phase"] == r["phase"] and c["ceiling"] == r["ceiling"]
        assert abs(c["A_hi"] - r["A"]) < 1e-9 and abs(c["tol"] - float(r["tol"])) < 1e-9
        assert abs(c["m_lo"] - float(r["A_measured_lo"])) < 1e-9
        assert abs(c["m_hi"] - float(r["A_measured_hi"])) < 1e-9
    return rows


def heldout():
    """Plot complete conversion intervals with the corresponding ceilings at right."""
    rows = heldout_rows()
    gap = 12
    low = sum(r["A"] < 0.97 for r in rows)
    y = np.array([i + (gap if i >= low else 0) for i in range(len(rows))], dtype=float)
    span = len(rows) + gap
    pitch = 0.034
    top, bottom = 0.55, 0.62
    height = top + bottom + span * pitch
    fig = plt.figure(figsize=(WIDTH, height))
    main = fig.add_axes([0.075, bottom / height, 0.745, span * pitch / height])
    strip = fig.add_axes([0.845, bottom / height, 0.13, span * pitch / height], sharey=main)
    pad = 2.5
    main.set_ylim(span - 1 + pad, -pad)
    main.set_xlim(-1.0, 0.4)
    main.axvline(0, color=INK, lw=0.8, zorder=1)
    for k, cat in enumerate(("ok", "possible violation", "violated")):
        idx = [i for i, r in enumerate(rows) if r["ceiling"] == cat]
        segs = [[(rows[i]["lo"], y[i]), (rows[i]["hi"], y[i])] for i in idx]
        main.add_collection(LineCollection(segs, colors=CATEGORY[cat], linewidths=1.2,
                                           capstyle="butt", zorder=2 + k))
        spans = [i for i in idx if float(rows[i]["A_measured_lo"]) != float(rows[i]["A_measured_hi"])]
        points = [i for i in idx if float(rows[i]["A_measured_lo"]) == float(rows[i]["A_measured_hi"])]
        main.scatter([rows[i]["mid"] for i in spans], y[spans], s=3.2, color=CATEGORY[cat],
                     linewidths=0, zorder=2.5 + k)
        big = cat == "violated"
        main.scatter([rows[i]["mid"] for i in points], y[points], s=24 if big else 8,
                     color=CATEGORY[cat], edgecolors=BG, linewidths=0.6 if big else 0.4,
                     zorder=6 if big else 2.6 + k)
    main.set_xticks([-1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4])
    main.set_xlabel("Measured A - ceiling - tolerance", fontsize=12)
    main.set_ylabel("Measurement", fontsize=12)
    main.tick_params(axis="y", left=False, labelleft=False)
    for side in ("top", "right", "left"):
        main.spines[side].set_visible(False)
    main.text(-0.985, -pad + 1, "Ceiling < 0.97", ha="left", va="top", fontsize=12)
    main.text(-0.985, low + gap / 2 - 0.5, "Ceiling ≥ 0.97", ha="left", va="center", fontsize=12)
    for part in (slice(0, low), slice(low, len(rows))):
        a = [r["A"] for r in rows[part]]
        edges = np.append(y[part] - 0.5, y[part][-1] + 0.5)
        strip.stairs(a, edges, orientation="horizontal", baseline=0, fill=True, color="#E2DED5", lw=0)
        strip.stairs(a, edges, orientation="horizontal", baseline=None, color=MUTED, lw=1.2)
    strip.set_xlim(0, 1.03)
    strip.set_xticks([0, 0.5, 1])
    strip.xaxis.set_major_formatter(FormatStrFormatter("%g"))
    strip.set_xlabel("Ceiling A", fontsize=12)
    strip.tick_params(axis="y", left=False, labelleft=False)
    for side in ("top", "right", "left"):
        strip.spines[side].set_visible(False)
    handles = [Line2D([], [], color=CATEGORY[c], lw=2.2, label=lab)
               for c, lab in (("ok", "Within"), ("possible violation", "Straddling"), ("violated", "Violated"))]
    main.legend(handles=handles, loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=3,
                frameon=False, fontsize=11.5, handlelength=1.8, columnspacing=1.6,
                borderaxespad=0.2)
    save(fig, "heldout-ceilings.png")
    plt.close(fig)
    return rows


def spread():
    """Connect sampled parent spreads without fitting or smoothing their alignment values."""
    rows = table(DATA / "mnal_c_A_vs_parent_spread.csv")
    pts = [(float(r["parent_fibre_spread_deg"]), float(r["A_along_fibre"]))
           for r in rows if r["parent_fibre_spread_deg"] != "random"]
    base = [float(r["A_along_fibre"]) for r in rows if r["parent_fibre_spread_deg"] == "random"]
    assert len(base) == 1 and len(pts) == len(rows) - 1
    x, a = np.array(pts).T
    fig = plt.figure(figsize=(WIDTH, 4.2))
    ax = fig.add_axes([0.085, 0.16, 0.89, 0.8])
    ax.axhline(1.0, color=MUTED, lw=1.0, zorder=1)
    ax.axhline(base[0], color=MUTED, lw=1.0, ls=(0, (4, 3)), zorder=1)
    ax.plot(x, a, color=OCHRE, lw=2.0, marker="o", ms=6.5, mfc=OCHRE, mec=BG, mew=1.2, zorder=3)
    ax.text(1.0, 0.988, "Perfect alignment", ha="left", va="top", fontsize=12, color=INK)
    ax.text(1.0, base[0] - 0.012, "Random parent", ha="left", va="top", fontsize=12, color=INK)
    ax.set_xlim(-2, 92)
    ax.set_ylim(0.44, 1.04)
    ax.set_xticks([0, 15, 30, 45, 60, 75, 90])
    ax.set_yticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    ax.set_xlabel("Parent spread σ (°)", fontsize=12.5)
    ax.set_ylabel("A along parent fibre", fontsize=12.5)
    ax.tick_params(labelsize=11)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    save(fig, "mnal-texture-spread.png")
    plt.close(fig)


def main():
    """Write four static figures from the saved numerical inputs."""
    argparse.ArgumentParser(description="Write four static alignment figures to figures/.").parse_args()
    style()
    atlas()
    heldout()
    spread()
    hd = hero_data()
    hero_still(hd)
    print("Wrote four static alignment figures to figures/.")


if __name__ == "__main__":
    main()
