"""Conditional predictions for published cases and their comparison with measurements.

A case spec (data/curated/cases.yaml) describes only the process conditions of a
published sample (never its outcome):

  - case: <outcome record id>            # value is read from the opened outcome record
    phase: <phase id>
    path: 0                              # formation path index (if a transformation is involved)
    texture: random | {fibre: <dir>} | single     # parent texture as made
    spread: [lo, hi]                     # deg, Gaussian sigma of the parent/particle texture
                                         #   (reported value -> [v, v]; unreported -> bounds)
    inherit: OR | random | either        # for incoherent (massive) transformations
    bias: none | field | stress          # applied during formation
    bias_acts: yes | no | phys           # 'phys': decided by the account's rules
    bias_axis_deg: <angle to fibre axis> # or 'best'
    free_body: none | single | internal  # particles/crystals rotated by a field
    axis: bias | fibre | perp | either   # measuring direction
    js: {ref: id} | {est:..}             # to convert Br to A when quantity is Br
    phi: [lo, hi]                        # hard-phase volume fraction when quantity is Br
    exchange_enhanced: false             # isotropic nanocrystalline Mr/Ms > 1/2 is not alignment

Prediction: [A_lo, A_hi] from the engine. Measurement: A_meas = Mr/Ms, or Br/(Js*phi)
(an interval when phi is an interval), or a texture metric converted to A.
"""
import numpy as np

from . import engine as E
from .routes import Path, _fibre

ISO = 0.5


def _parents(path, tex, spread):
    if tex == 'random':
        return E.texture_random(20000)
    if tex == 'single':
        return E.texture_single(n=2000, spread_deg=spread)
    d = tex['fibre']
    return E.texture_fibre(path.LP.direction(d), n_az=360 if spread else 720, spread_deg=spread)


def _A(path, parents, inherit, bias, bias_on, bias_axis, free_body, axis, sign=+1):
    ax, fr, per = path.axes(parents)
    if inherit == 'random':
        return ISO
    w = fr
    h = np.array([0.0, 0.0, 1.0])
    if bias != 'none' and bias_on:
        al = np.radians(bias_axis)
        h = np.array([np.sin(al), 0.0, np.cos(al)])
        w = E.bias_weights(ax, fr, len(path.V), per, 'field' if bias == 'field' else 'stress', h,
                           sign=sign)
    if free_body == 'internal':
        return E.body_alignment(ax, w)
    if axis == 'bias':
        return E.alignment(ax, w, h)
    if axis == 'fibre':
        return E.alignment(ax, w, (0, 0, 1))
    if axis == 'perp':
        return E.alignment(ax, w, (1, 0, 0))
    raise ValueError(axis)


def predict(ph, spec):
    """[A_lo, A_hi] for the case conditions; the extremes of every unreported parameter."""
    if spec.get('free_body') == 'single':
        lo_s, hi_s = spec.get('spread', [0, 0])
        vals = []
        for s in {lo_s, hi_s}:
            g = E.texture_single(n=4000, spread_deg=s)
            ax, fr, per = E.product_axes(g, np.eye(3)[None], np.array([0, 0, 1.0]))
            vals.append(E.alignment(ax, fr))
        return min(vals), max(vals)
    p = ph['paths'][spec.get('path', 0)] if spec.get('path') is not None else {'mode': 'direct'}
    path = Path(ph, p)
    inh = spec.get('inherit', 'OR')
    inherits = ['OR', 'random'] if inh == 'either' else [inh]
    axes = ['fibre', 'perp'] if spec.get('axis') == 'either' else [spec.get('axis', 'bias')]
    on = spec.get('bias_acts', 'no') == 'yes'
    vals = []
    for s in set(spec.get('spread', [0, 0])):
        parents = _parents(path, spec.get('texture', 'random'), s)
        for i in inherits:
            for a in axes:
                vals.append(_A(path, parents, i, spec.get('bias', 'none'), on,
                               spec.get('bias_axis_deg', 0), spec.get('free_body', 'none'), a,
                               p.get('stress_sign', +1)))
    return min(vals), max(vals)


def measured_A(rec, spec, js):
    """A_meas interval from an outcome record."""
    q, v = rec['quantity'], float(rec['value'])
    u = float(rec.get('uncertainty') or 0.0)
    if q == 'Mr_over_Ms':
        return v - u, v + u
    if q == 'Br':
        plo, phi = spec.get('phi', [1.0, 1.0])
        return (v - u) / (js * phi), (v + u) / (js * plo)
    if q == 'texture_metric' and rec.get('unit') == 'A_reported':
        return v - u, v + u
    if q == 'texture_metric' and rec.get('unit') == 'variant_fraction':
        return None
    return None
