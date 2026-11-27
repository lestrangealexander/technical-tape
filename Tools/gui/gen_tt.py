#!/usr/bin/env python3
"""Technical Tape (Special 1500) faceplate in UVP's faceplate-v4 style, from the user's sketch.

  gen_tt.py <root> canvas|atr|maxell
"""
import json, math, sys
import numpy as np

W, H = 1000, 660
SALMON, SALMON_DK, YELLOW, INK, WHITE = '#FEA993', '#EE917B', '#FDE6B0', '#3A302C', '#FFFFFF'
SKY = '#D1EEFE'
FONT = "'DM Sans'"
PEGASUS_URL = '/_blob/6febd93d1c650002fba59be3297df40f'

R = 110.0
REEL_L, REEL_R, VU_IN, VU_OUT = (138.0, 180.0), (379.0, 180.0), (621.0, 180.0), (862.0, 180.0)
RP_L, RP_R = 0.80 * R, 0.55 * R
BOX_W, BOX_H, BOX_Y = 150.0, 150.0, 310.0
# middle row, three boxes of one size, each centred between two of the top row: tape select between the
# reels, the badge between the take-up reel and the input meter (the faceplate's centre), speed select
# between the meters
SEL_X = (REEL_L[0] + REEL_R[0]) / 2 - BOX_W / 2      # 183.5
LOGO_X = (REEL_R[0] + VU_IN[0]) / 2 - BOX_W / 2      # 425
SPD_X = (VU_IN[0] + VU_OUT[0]) / 2 - BOX_W / 2       # 666.5
BOX_CY = BOX_Y + BOX_H / 2
SEL_HEAD = 38.0                                      # header row of a select box
# the wordmark under the badge, set justified: the first two lines exactly as wide as the box above them,
# "Special 1500" larger, 25 px past the box on each side (the user's mockup, 27 Nov)
SCRIPT_W = BOX_W + 50.0                              # 200
SCRIPT_X = LOGO_X + (BOX_W - SCRIPT_W) / 2           # 400
TEXT_BASELINES = (484.0, 504.0, 546.0)
ROW_Y = 548.0
SMALL_R, SMALL_Y = 28.0, 572.0
AG_C = (REEL_L[0] + 100.0, SMALL_Y)                  # right of INPUT
RESOURCES = '/home/claude/technical-tape/Resources'


def fit_size(font_file, s, width, spacing):
    """Font size at which `s` with this letter-spacing is `width` wide."""
    lo, hi = 4.0, 120.0
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        if text_width(font_file, mid, s, spacing) > width:
            hi = mid
        else:
            lo = mid
    return lo


def ink_edges(font_file, size, s):
    """(left bearing of the first glyph, right bearing of the last) at this size, from the outlines: what
    the advances add beyond the visible ink at each end (negative where a swash overhangs)."""
    from fontTools.ttLib import TTFont
    from fontTools.pens.boundsPen import BoundsPen
    font = TTFont(f'{RESOURCES}/{font_file}')
    upm, cmap, glyphs, hmtx = font['head'].unitsPerEm, font.getBestCmap(), font.getGlyphSet(), font['hmtx']

    def bounds(ch):
        name = cmap[ord(ch)]
        pen = BoundsPen(glyphs)
        glyphs[name].draw(pen)
        return hmtx[name][0], pen.bounds

    _, (x0, _, _, _) = bounds(s[0])
    adv, (_, _, x1, _) = bounds(s[-1])
    return x0 * size / upm, (adv - x1) * size / upm


def wordmark():
    """Three lines whose visible ends sit on their column's edges (a magazine column). The DM Sans lines are
    letter-spaced to the badge's width; the script is scaled as a whole until its ink spans its wider
    column (as the plug-in draws it). textLength makes each exact whatever the browser's font metrics,
    allowing for each line's end bearings so the ink, not the advances, lines up."""
    from PIL import ImageFont
    lines = [('DMSans-Medium.ttf', 'Sound Aquarium Recording', 0.3, 500, "'DM Sans'", 'spacing', LOGO_X, BOX_W),
             ('DMSans-Bold.ttf', 'TECHNICAL TAPE', 1.2, 700, "'DM Sans'", 'spacing', LOGO_X, BOX_W),
             ('Yellowtail-Regular.ttf', 'Special 1500', 0.0, 400, "'Yellowtail', cursive", 'spacingAndGlyphs', SCRIPT_X, SCRIPT_W)]
    o = []
    for (ff, txt, sp, wt, fam, adj, x0, width), y in zip(lines, TEXT_BASELINES):
        if adj == 'spacing':
            size = fit_size(ff, txt, width, sp)
            lsb, rsb = ink_edges(ff, size, txt)
            x, length = x0 - lsb, width + lsb + rsb
        else:
            # uniform scale: the size at which the ink (advance less the end bearings) spans the column
            def ink(sz):
                ft = ImageFont.truetype(f'{RESOURCES}/{ff}', int(round(sz * 64)))
                l, r = ink_edges(ff, sz, txt)
                return ft.getlength(txt) / 64.0 - l - r
            size = 40.0 * width / ink(40.0)
            ft = ImageFont.truetype(f'{RESOURCES}/{ff}', int(round(size * 64)))
            lsb, _ = ink_edges(ff, size, txt)
            x, length = x0 - lsb, ft.getlength(txt) / 64.0
        o.append(f'<text x="{x:.2f}" y="{f(y)}" font-family="{fam}" font-size="{size:.2f}" font-weight="{wt}" '
                 f'text-anchor="start" textLength="{length:.2f}" lengthAdjust="{adj}" fill="{WHITE}">{txt}</text>')
    return '\n'.join(o)


def text_width(font_file, size, s, spacing=0.0):
    from PIL import ImageFont
    return ImageFont.truetype(f'{RESOURCES}/{font_file}', int(round(size * 64))).getlength(s) / 64.0 + spacing * (len(s) - 1)
PACK, PACK_LINE = '#4A3C35', '#57483F'
REELS = {   # flange, rear flange seen through the windows, rim line, screws
    'atr': (SKY, '#A6CBE8', '#B5D8F1', '#9DC3E2'),
    'maxell': (WHITE, '#E6D8D1', '#EDE3DD', '#DCCBC3'),
    'black': (INK, '#231C19', '#2D2522', '#5C4C45'),
    'white': (WHITE, '#E6D8D1', '#EDE3DD', '#DCCBC3'),
}


def f(v):
    s = f'{v:.2f}'.rstrip('0').rstrip('.')
    return '0' if s in ('-0', '') else s


def polar(cx, cy, r, deg):
    a = math.radians(deg)
    return cx + r * math.sin(a), cy - r * math.cos(a)


def rot(p, deg):
    a = math.radians(deg)
    return p[0] * math.cos(a) - p[1] * math.sin(a), p[0] * math.sin(a) + p[1] * math.cos(a)


def arc(cx, cy, r, a0, a1):
    x0, y0 = polar(cx, cy, r, a0); x1, y1 = polar(cx, cy, r, a1)
    return f'M{f(x0)},{f(y0)} A{f(r)},{f(r)} 0 {1 if a1 - a0 > 180 else 0} 1 {f(x1)},{f(y1)}'


def text(x, y, s, size, weight, spacing=0.0, fill=INK, anchor='middle', opacity=None, italic=False):
    op = f' fill-opacity="{opacity}"' if opacity is not None else ''
    it = ' font-style="italic"' if italic else ''
    return (f'<text x="{f(x)}" y="{f(y)}" font-family="{FONT}" font-size="{f(size)}" font-weight="{weight}"{it} '
            f'text-anchor="{anchor}" letter-spacing="{f(spacing)}" fill="{fill}"{op}>{s}</text>')


def circle_path(cx, cy, r):
    return f'M{f(cx + r)},{f(cy)} A{f(r)},{f(r)} 0 1 0 {f(cx - r)},{f(cy)} A{f(r)},{f(r)} 0 1 0 {f(cx + r)},{f(cy)} Z'


def rounded_poly(pts, radii):
    P = [np.array(p, float) for p in pts]; n = len(P)

    def walk(i, dist, d):
        j, rem, passed = i, dist, []
        while True:
            k = (j + d) % n
            seg = float(np.linalg.norm(P[k] - P[j]))
            if seg >= rem - 1e-9:
                return P[j] + (P[k] - P[j]) * (rem / seg), passed
            rem -= seg; passed.append(k); j = k

    skip, special = set(), {}
    for c, r in radii.items():
        a, pa = walk(c, r, -1); b, pb = walk(c, r, 1)
        skip.update(pa); skip.update(pb); special[c] = (a, b)
    cmds = []
    for i in range(n):
        if i in skip:
            continue
        lead = 'M' if not cmds else 'L'
        if i in special:
            a, b = special[i]
            cmds.append(f'{lead}{f(a[0])},{f(a[1])} Q{f(P[i][0])},{f(P[i][1])} {f(b[0])},{f(b[1])}')
        else:
            cmds.append(f'{lead}{f(P[i][0])},{f(P[i][1])}')
    return ' '.join(cmds) + ' Z'


# ---------------------------------------------------------------- reels
def windows(style, cx, cy, turn):
    if style == 'three':      # rounded trapezoids, curved outer edge
        ri, ro, hi = 0.38 * R, 0.86 * R, 0.13 * R
        tho = math.degrees(math.asin(0.33 * R / ro))
        base = [(-hi, -ri), (hi, -ri)] + [polar(0, 0, ro, tho - 2 * tho * k / 23) for k in range(24)]
        out = []
        for a in (0, 120, 240):
            pts = [(x + cx, y + cy) for x, y in (rot(p, turn + a) for p in base)]
            out.append(rounded_poly(pts, {0: 5, 1: 5, 2: 7, len(pts) - 1: 7}))
        return out
    x0, x1, h, r = 0.38 * R, 0.84 * R, 0.12 * R, 6.0     # two straight slots
    out = []
    for side in (1, -1):
        def P(x, y):
            px, py = rot((side * x, side * y), turn)
            return f'{f(px + cx)},{f(py + cy)}'
        out.append(f'M{P(x0 + r, -h)} L{P(x1 - r, -h)} A{f(r)},{f(r)} 0 0 1 {P(x1, -h + r)} L{P(x1, h - r)} '
                   f'A{f(r)},{f(r)} 0 0 1 {P(x1 - r, h)} L{P(x0 + r, h)} A{f(r)},{f(r)} 0 0 1 {P(x0, h - r)} '
                   f'L{P(x0, -h + r)} A{f(r)},{f(r)} 0 0 1 {P(x0 + r, -h)} Z')
    return out


def reel(colours, style, cx, cy, turn, rp):
    flange, back, rim, screw = REELS[colours]
    wins = windows(style, cx, cy, turn)
    hole = 0.31 * R
    o = [f'<circle cx="{f(cx)}" cy="{f(cy)}" r="{f(R - 2)}" fill="{back}"/>',
         f'<circle cx="{f(cx)}" cy="{f(cy)}" r="{f(rp)}" fill="{PACK}"/>']
    for k in (8, 18, 30):
        if rp - k > 0.4 * R:
            o.append(f'<circle cx="{f(cx)}" cy="{f(cy)}" r="{f(rp - k)}" fill="none" stroke="{PACK_LINE}" stroke-width="1"/>')
    o.append(f'<path d="{circle_path(cx, cy, R)} {" ".join(wins)} {circle_path(cx, cy, hole)}" fill-rule="evenodd" fill="{flange}"/>')
    o.append(f'<circle cx="{f(cx)}" cy="{f(cy)}" r="{f(R - 1.25)}" fill="none" stroke="{rim}" stroke-width="2.5"/>')
    # NAB centre: the reel's drive keys, then the deck's white adapter
    o.append(f'<circle cx="{f(cx)}" cy="{f(cy)}" r="{f(hole)}" fill="{INK}"/>')
    for a in (0, 120, 240):
        p = [rot(q, turn + a) for q in ((-4, -hole), (4, -hole), (4, -0.255 * R), (-4, -0.255 * R))]
        o.append('<path d="M' + ' L'.join(f'{f(x + cx)},{f(y + cy)}' for x, y in p) + f' Z" fill="{flange}"/>')
    o.append(f'<circle cx="{f(cx)}" cy="{f(cy)}" r="{f(0.245 * R)}" fill="{WHITE}"/>')
    for a in (60, 180, 300):
        x, y = polar(cx, cy, 0.155 * R, a + turn)
        o.append(f'<circle cx="{f(x)}" cy="{f(y)}" r="2.6" fill="{INK}" fill-opacity="0.75"/>')
    o.append(f'<circle cx="{f(cx)}" cy="{f(cy)}" r="7" fill="{INK}"/><circle cx="{f(cx)}" cy="{f(cy)}" r="2.6" fill="{WHITE}"/>')
    for a in ((60, 180, 300) if style == 'three' else (45, 135, 225, 315)):
        x, y = rot((0, -0.36 * R), turn + a)
        o.append(f'<circle cx="{f(x + cx)}" cy="{f(y + cy)}" r="2.4" fill="{screw}"/>')
    return '\n'.join(o)


def toward(c, p, r):
    """The point at distance r from c, on the way to p."""
    d = math.hypot(p[0] - c[0], p[1] - c[1])
    return c[0] + (p[0] - c[0]) * r / d, c[1] + (p[1] - c[1]) * r / d


def tangent_from(point, c, r, side):
    """Tangent point on circle (c, r) for a line from an outside point; side=+1 picks the right-hand one."""
    dx, dy = point[0] - c[0], point[1] - c[1]
    d = math.hypot(dx, dy); phi = math.atan2(dy, dx); a = math.acos(r / d)
    pts = [(c[0] + r * math.cos(phi + s * a), c[1] + r * math.sin(phi + s * a)) for s in (1, -1)]
    return max(pts, key=lambda p: side * p[0])


# ---------------------------------------------------------------- VU meters (UVP's round meter)
def vu_angle(vu):
    return -50.0 + 100.0 * 10 ** ((vu - 3.0) / 20.0)


def vu_meter(cx, cy, vu, label):
    px, py = cx, cy + R * 0.6
    band = R * 1.04 - 6.0
    o = [f'<circle cx="{f(cx)}" cy="{f(cy)}" r="{f(R)}" fill="{WHITE}"/>',
         f'<path d="{arc(px, py, band, vu_angle(-20), vu_angle(0))}" fill="none" stroke="{YELLOW}" stroke-width="13" stroke-linecap="round"/>',
         f'<path d="{arc(px, py, band, vu_angle(0), vu_angle(3))}" fill="none" stroke="{SALMON}" stroke-width="13" stroke-linecap="round"/>']
    major = {-20: '20', -10: '10', -5: '5', 0: '0', 3: '+3'}
    for v in (-20, -10, -7, -5, -3, -2, -1, 0, 1, 2, 3):
        a = vu_angle(v); big = v in major
        p1, p2 = polar(px, py, band - 13, a), polar(px, py, band - (22 if big else 18), a)
        o.append(f'<line x1="{f(p1[0])}" y1="{f(p1[1])}" x2="{f(p2[0])}" y2="{f(p2[1])}" stroke="{INK}" '
                 f'stroke-opacity="{0.4 if big else 0.22}" stroke-width="{1.8 if big else 1.2}" stroke-linecap="round"/>')
        if big:
            tx, ty = polar(px, py, band - 36, a)
            o.append(text(tx, ty + 4.5, major[v], 13, 600, 0, INK, opacity=0.6))
    o.append(text(cx, py - 24, label, 12, 700, 2.4, INK, opacity=0.55))
    a = vu_angle(vu)
    n1, n2 = polar(px, py, 13, a), polar(px, py, band + 4, a)
    o.append(f'<line x1="{f(n1[0])}" y1="{f(n1[1])}" x2="{f(n2[0])}" y2="{f(n2[1])}" stroke="{INK}" stroke-opacity="0.85" stroke-width="2.4" stroke-linecap="round"/>')
    o.append(f'<circle cx="{f(px)}" cy="{f(py)}" r="9" fill="{SALMON}"/><circle cx="{f(px)}" cy="{f(py)}" r="3.4" fill="{WHITE}"/>')
    return '\n'.join(o)


# ---------------------------------------------------------------- knobs (UVP ParamKnob)
def knob(cx, cy, body, ring, value, label):
    rr = body + ring * 0.5 + max(4.0, body * 0.06)
    av = 135.0 * value
    o = [f'<path d="{arc(cx, cy, rr, -135, 135)}" fill="none" stroke="{SALMON_DK}" stroke-width="{f(ring)}" stroke-linecap="round"/>']
    if abs(av) > 0.5:
        o.append(f'<path d="{arc(cx, cy, rr, min(0.0, av), max(0.0, av))}" fill="none" stroke="{WHITE}" stroke-width="{f(ring)}" stroke-linecap="round"/>')
    else:
        x, y = polar(cx, cy, rr, 0)
        o.append(f'<circle cx="{f(x)}" cy="{f(y)}" r="{f(ring / 2)}" fill="{WHITE}"/>')
    o.append(f'<circle cx="{f(cx)}" cy="{f(cy)}" r="{f(body)}" fill="{WHITE}"/>')
    p1, p2 = polar(cx, cy, body * 0.42, av), polar(cx, cy, body * 0.8, av)
    o.append(f'<line x1="{f(p1[0])}" y1="{f(p1[1])}" x2="{f(p2[0])}" y2="{f(p2[1])}" stroke="{INK}" stroke-width="{f(max(2.6, body * 0.06))}" stroke-linecap="round"/>')
    o.append(text(cx, 634, label, 13, 700, 2.6))
    return '\n'.join(o)


# ---------------------------------------------------------------- static faceplate and the left-reel overlays
def base_svg(needles):
    o = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" role="img" '
         f'aria-label="Sound Aquarium Technical Tape Special 1500 faceplate">',
         f'<rect width="{W}" height="{H}" fill="{SALMON}"/>',
         f'<rect x="0" y="0" width="{W}" height="42" fill="{YELLOW}"/>',
         text(662, 26, 'OS', 13, 600, 1.6),
         f'<rect x="752" y="9" width="230" height="24" rx="12" fill="{WHITE}"/>',
         f'<line x1="900" y1="14" x2="900" y2="28" stroke="{INK}" stroke-opacity="0.3" stroke-width="1.5" stroke-linecap="round"/>']
    # tape: off the supply pack's right side and the take-up pack's left side, down into the tape select box
    e1, e2 = (SEL_X + BOX_W / 2 - 25, BOX_Y + 2), (SEL_X + BOX_W / 2 + 25, BOX_Y + 2)
    t1 = tangent_from(e1, REEL_L, RP_L, +1); t2 = tangent_from(e2, REEL_R, RP_R, -1)
    o.append(f'<path d="M{f(t1[0])},{f(t1[1])} L{f(e1[0])},{f(e1[1])} M{f(e2[0])},{f(e2[1])} L{f(t2[0])},{f(t2[1])}" '
             f'fill="none" stroke="{INK}" stroke-width="3" stroke-linecap="round"/>')
    o.append(reel('white', 'three', REEL_R[0], REEL_R[1], 47, RP_R))
    # the badge box in the middle of the faceplate
    o.append(f'<rect x="{f(LOGO_X + 1.5)}" y="{f(BOX_Y + 1.5)}" width="{f(BOX_W - 3)}" height="{f(BOX_H - 3)}" rx="11" fill="{SALMON}" stroke="{WHITE}" stroke-width="3"/>')
    # wordmark under the badge, set justified: every line exactly as wide as the box above it
    o.append(wordmark())
    o.append(vu_meter(*VU_IN, needles[0], 'INPUT VU'))
    o.append(vu_meter(*VU_OUT, needles[1], 'OUTPUT VU'))
    o.append(knob(REEL_L[0], ROW_Y, 46, 6, 1 / 3, 'INPUT'))     # under the supply reel
    o.append(knob(VU_OUT[0], ROW_Y, 40, 5, 0, 'OUTPUT'))      # under the output meter
    o.append('</svg>')
    return '\n'.join(o)


def left_reel_svg(tape):
    colours, style = ('atr', 'three') if tape == 'atr' else ('maxell', 'slots')
    label = 'ATR tape on a sky-blue reel' if tape == 'atr' else 'Maxell tape on a white reel'
    return '\n'.join([f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" role="img" aria-label="{label}">',
                      reel(colours, style, REEL_L[0], REEL_L[1], 12, RP_L), '</svg>'])


# ---------------------------------------------------------------- HTML controls (UVP behaviour)
def select_box(x, header, options, size, spacing, aria):
    """A white-bordered box like the badge's: a header over its choices, one row each. The chosen row is
    solid white with salmon letters. options: (label, state key, handler, aria label)."""
    row = (f'flex: 1 1 0; width: 100%; margin: 0; padding: 0 0 0 {f(spacing)}px; border: 0; font-family: {FONT}, sans-serif; '
           f'font-weight: 700; font-size: {f(size)}px; letter-spacing: {f(spacing)}px; cursor: pointer')
    o = [f'<div role="group" aria-label="{aria}" style="position: absolute; left: {f(x)}px; top: {f(BOX_Y)}px; '
         f'width: {f(BOX_W)}px; height: {f(BOX_H)}px; box-sizing: border-box; border: 3px solid {WHITE}; border-radius: 12.5px; '
         f'background: {SALMON}; overflow: hidden; display: flex; flex-direction: column">'
         f'<div style="flex: 0 0 {f(SEL_HEAD - 6)}px; display: flex; align-items: center; justify-content: center; '
         f'padding-left: 2.2px; border-bottom: 3px solid {WHITE}; color: {WHITE}; font-family: {FONT}, sans-serif; '
         f'font-weight: 700; font-size: 12.5px; letter-spacing: 2.2px">{header}</div>']
    for i, opt in enumerate(options):
        label, key, handler, a = opt[:4]
        tip = opt[4] if len(opt) > 4 else None
        rule = f'; border-bottom: 3px solid {WHITE}' if i < len(options) - 1 else ''
        if tip:
            # not available yet: faint, ignores the mouse, says why
            o.append(f'<button type="button" disabled aria-disabled="true" title="{tip}" aria-label="{a} ({tip})" '
                     f'style="{row}{rule}; background: {SALMON}; color: rgba(255, 255, 255, 0.4); cursor: default">{label}</button>')
            continue
        o.append(f'<button type="button" onClick="{{{{{handler}}}}}" aria-pressed="{{{{{key}On}}}}" aria-label="{a}" '
                 f'style="{row}{rule}; background: {{{{{key}Bg}}}}; color: {{{{{key}Fg}}}}">{label}</button>')
    o.append('</div>')
    return ''.join(o)


def round_button(key, cx, label, handler, aria):
    """AG and the speed buttons: a 56 px roller with a white ring. On: solid white, salmon letters."""
    return (f'<button type="button" onClick="{{{{{handler}}}}}" aria-pressed="{{{{{key}On}}}}" aria-label="{aria}" '
            f'style="position: absolute; left: {f(cx - SMALL_R)}px; top: {f(SMALL_Y - SMALL_R)}px; width: {f(2 * SMALL_R)}px; '
            f'height: {f(2 * SMALL_R)}px; margin: 0; padding: 0 0 0 1.2px; box-sizing: border-box; border: 3px solid {WHITE}; '
            f'border-radius: 50%; background: {{{{{key}Bg}}}}; color: {{{{{key}Fg}}}}; font-family: {FONT}, sans-serif; '
            f'font-weight: 700; font-size: 18px; letter-spacing: 1.2px; cursor: pointer">{label}</button>')


def pill(x, w, label, size, spacing, weight, extra=''):
    return (f'<button type="button" style="position: absolute; left: {x}px; top: 9px; width: {w}px; height: 24px; margin: 0; padding: 0; '
            f'border: 0; border-radius: 12px; font-family: {FONT}, sans-serif; font-weight: {weight}; font-size: {size}px; '
            f'letter-spacing: {spacing}px; cursor: pointer; {extra}">{label}</button>')


def artboard(title, tape, ips, needles, maxell_available=True):
    zone = f'font-family: {FONT}, sans-serif; font-weight: 600; font-size: 12px; letter-spacing: 1.6px; color: {INK}'
    body = '\n'.join([
        f'<div style="position: relative; width: {W}px; height: {H}px; overflow: hidden; background: {SALMON}; font-family: {FONT}, sans-serif">',
        f'<div style="position: absolute; left: 0; top: 0; width: {W}px; height: {H}px">', base_svg(needles), '</div>',
        f'<div style="position: absolute; left: 0; top: 0; width: {W}px; height: {H}px; display: {{{{atrShow}}}}">', left_reel_svg('atr'), '</div>',
        f'<div style="position: absolute; left: 0; top: 0; width: {W}px; height: {H}px; display: {{{{maxShow}}}}">', left_reel_svg('maxell'), '</div>',
        f'<img src="{PEGASUS_URL}" alt="Sound Aquarium Pegasus badge" style="position: absolute; left: {f(LOGO_X + BOX_W / 2 - 59 * 1218 / 1249)}px; '
        f'top: {f(BOX_Y + BOX_H / 2 - 59)}px; width: {f(118 * 1218 / 1249)}px; height: 118px">',
        pill(18, 84, 'BYPASS', 12, 1.4, 600, 'background: {{bypBg}}; color: {{bypFg}}').replace('<button type="button"', '<button type="button" onClick="{{toggleByp}}" aria-pressed="{{bypOn}}"'),
        pill(680, 56, '8X', 13, 0.6, 700, f'background: {WHITE}; color: {INK}').replace('<button type="button"', '<button type="button" aria-label="Oversampling: 8x"'),
        f'<button type="button" aria-label="Preset menu" style="position: absolute; left: 752px; top: 9px; width: 148px; height: 24px; margin: 0; '
        f'padding: 0; border: 0; border-radius: 12px 0 0 12px; background: transparent; cursor: pointer; text-align: left">'
        f'<svg width="14" height="12" viewBox="0 0 14 12" aria-hidden="true" style="position: absolute; left: 15px; top: 6px">'
        f'<path d="M1,1 H13 M1,6 H13 M1,11 H13" stroke="{INK}" stroke-width="2" stroke-linecap="round" fill="none"/></svg>'
        f'<span style="position: absolute; left: 40px; top: 5px; {zone}">MENU</span></button>',
        f'<button type="button" aria-label="Save preset" style="position: absolute; left: 900px; top: 9px; width: 82px; height: 24px; margin: 0; '
        f'padding: 0; border: 0; border-radius: 0 12px 12px 0; background: transparent; cursor: pointer; {zone}">SAVE</button>',
        select_box(SEL_X, 'TAPE SELECT', [('ATR', 'atr', 'pickAtr', 'ATR tape'),
                                          ('MAXELL', 'max', 'pickMax', 'Maxell tape') + (() if maxell_available else ('Maxell: not captured yet',))],
                   20, 2.5, 'Tape select'),
        select_box(SPD_X, 'SPEED SELECT', [('3¾ IPS', 'i375', 'pick375', '3¾ ips'), ('7½ IPS', 'i75', 'pick75', '7½ ips'),
                                           ('15 IPS', 'i15', 'pick15', '15 ips')], 18, 2.0, 'Speed select'),
        round_button('ag', AG_C[0], 'AG', 'toggleAg', 'Auto gain'),
        '</div>'])
    props = json.dumps({'tape': {'editor': 'enum', 'options': ['atr', 'maxell'], 'default': tape},
                        'ips': {'editor': 'enum', 'options': ['15', '7.5', '3.75'], 'default': ips},
                        'autoGain': {'editor': 'boolean', 'default': True},
                        '$preview': {'width': W, 'height': H}}, ensure_ascii=False)
    logic = f'''class Component extends DCLogic {{
renderVals() {{
const st = this.state || {{}};
const tape = st.tape || this.props.tape || '{tape}';
const ips = st.ips || this.props.ips || '{ips}';
const ag = st.ag !== undefined ? st.ag : (this.props.autoGain !== undefined ? this.props.autoGain : true);
const byp = st.byp === true;
const v = {{
atrShow: tape === 'atr' ? 'block' : 'none',
maxShow: tape === 'maxell' ? 'block' : 'none',
pickAtr: () => this.setState({{ tape: 'atr' }}),
pickMax: () => this.setState({{ tape: 'maxell' }}),
pick15: () => this.setState({{ ips: '15' }}),
pick75: () => this.setState({{ ips: '7.5' }}),
pick375: () => this.setState({{ ips: '3.75' }}),
toggleAg: () => this.setState({{ ag: !ag }}),
toggleByp: () => this.setState({{ byp: !byp }})
}};
const sel = {{ atr: tape === 'atr', max: tape === 'maxell', i15: ips === '15', i75: ips === '7.5', i375: ips === '3.75', ag: ag }};
Object.keys(sel).forEach((k) => {{ const on = sel[k]; v[k + 'Bg'] = on ? '{WHITE}' : '{SALMON}'; v[k + 'Fg'] = on ? '{SALMON}' : '{WHITE}'; v[k + 'On'] = on ? 'true' : 'false'; }});
const tog = {{ byp: byp }};
Object.keys(tog).forEach((k) => {{ const on = tog[k]; v[k + 'Bg'] = on ? '{INK}' : '{WHITE}'; v[k + 'Fg'] = on ? '{WHITE}' : '{INK}'; v[k + 'On'] = on ? 'true' : 'false'; }});
return v;
}}
}}'''
    return f'''<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>{title}</title>
<script src="./support.js"></script>
</head>
<body>
<x-dc>
<helmet>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=DM+Sans:ital,opsz,wght@0,9..40,100..1000;1,9..40,100..1000&amp;family=Yellowtail&amp;display=swap">
<style>
body{{margin:0}}
a{{color:#3A302C}}a:hover{{color:#000000}}
</style>
</helmet>
{body}
</x-dc>
<script type="text/x-dc" data-dc-script data-props='{props}'>
{logic}
</script>
</body>
</html>
'''


def update_canvas(root):
    p = f'{root}/project/canvas.json'
    c = json.load(open(p))
    c['boards']['TechnicalTape.dc.html'] = {'x': 0, 'y': 3600, 'w': W, 'h': H, 'title': 'Technical Tape · ATR (sky-blue reel) · 15 ips — MAXELL greyed out until it is captured', 'is_interactive': True}
    c['boards']['TechnicalTapeMaxell.dc.html'] = {'x': W + 80, 'y': 3600, 'w': W, 'h': H, 'title': 'Later: Maxell (white reel) · 7½ ips, once the Maxell captures are in', 'is_interactive': True}
    for b in ('TechnicalTape.dc.html', 'TechnicalTapeMaxell.dc.html'):
        if b not in c['order']:
            c['order'].append(b)
    c['notes']['titleV4'] = {'x': 0, 'y': 3340, 'text': 'Technical Tape — Special 1500: your sketch in UVP style', 'kind': 'title1', 'maxW': 2 * W + 80}
    open(p, 'w').write(json.dumps(c, indent=1, ensure_ascii=False))


if __name__ == '__main__':
    root, what = sys.argv[1], sys.argv[2]
    if what == 'canvas':
        update_canvas(root)
    elif what == 'atr':
        open(f'{root}/project/TechnicalTape.dc.html', 'w').write(artboard('Technical Tape Special 1500 — ATR, 15 ips', 'atr', '15', (-2.0, -4.0),
                                                                          maxell_available=False))
    elif what == 'maxell':
        open(f'{root}/project/TechnicalTapeMaxell.dc.html', 'w').write(artboard('Technical Tape Special 1500 — Maxell, 7½ ips', 'maxell', '7.5', (-5.0, -7.0)))
    print('ok', what)
