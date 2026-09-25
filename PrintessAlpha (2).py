#!/usr/bin/env python3
"""Hybrid bioprinter — single-file build (rebuild R2)."""
import sys

# ===================== compiler (matrix_to_gcode.py) =====================
#!/usr/bin/env python3
"""
matrix_to_gcode.py — hybrid bioprinter compiler (clean: ONE mapping, B-dispense).

Reads a Studio matrix (grid OR freeform) and emits a per-layer program via
emit_layers(M, cfg) -> ({'header':[...], 'layers':[{'i','z','lines':[...]}]}, msg).
The host app drives park / present / spray-pause around these per-layer bodies.

ONE coordinate mapping, fixed (no rotate / flip / home-corner knobs):
  HOME = the parked corner = max X/Y = canvas BOTTOM-LEFT. The design grows INTO
  the bed from home, so machineX = MAX - x_from_home and machineY = MAX - y_from_home
  with x/y_from_home >= 0 — nothing is ever commanded past the corner.
  Droplet = scaffold - DROPLET_DX in X (the bed rolls -X so the droplet head reaches
  where the scaffold printed); the result may be negative and that is fine — the blue
  droplet zone extends front-left of the corner.

Dispense is the B plunger only (G1 B). No solenoid / FAN / M106 anywhere.
"""

import json, argparse, sys
from dataclasses import dataclass


# ============================================================================ config
@dataclass
class Cfg:
    # scaffold extrusion (C squeeze)
    FLOW: float = 0.02              # was 0.05; rescaled so the old 0.40x is the new 1.0x
    PRIME: float = 0.6              # was 2.0; plunger pre-charge at stroke start (too big = start glob)
    RETRACT: float = 0.5            # was 1.0; plunger pull-back at stroke end (stops ooze on travel)
    CORNER_RETRACT: float = 0.0     # C pulled back + restored at each sharp turn to kill corner blobs
    PRINT_F: int = 240
    TRAVEL_F: int = 1500
    RETRACT_F: int = 800
    # machine frame — the one mapping
    COORD_MAX_X: float = 71.0       # home corner X (max)
    COORD_MAX_Y: float = 71.0       # home corner Y (max)
    SCAFFOLD_REACH: float = 60.0    # yellow square the scaffold head reaches from home
    DROPLET_REACH_X: float = 70.0   # blue width  the droplet head reaches
    DROPLET_REACH_Y: float = 60.0   # blue height the droplet head reaches
    DROPLET_DX: float = 31.0        # bed -X roll that puts the droplet head where scaffold printed
    DROPLET_DY: float = 0.0         # same in Y (the heads don't share a Y on the solenoid head)
    # scaffold lift (A)
    A_HOP: float = 2.0
    A_OFFSET: float = 0.0           # + digs the print plane down into the gel
    PRINT_LIFT: float = 0.0         # + raises the nozzle above the print plane every layer (anti-drag)
    AF: int = 300                   # scaffold-lift feed: matches the jog speed that returns cleanly
    LAYER_LIFT: float = None        # A raise per layer (default = bead height)
    DROPLET_RISE: float = None      # needle Z raise per layer (default = bead); set to real gel growth
    # droplet (Z dip + B push)
    DIP: float = 3.0                # relative needle dip at each stop
    DROP_HOLD: float = 1.0          # seconds the needle holds at depth AFTER pushing, before retract
    DROP_GAP: float = 0.5           # seconds between successive drops at the SAME spot (so each detaches)
    DROP_VOL: float = 0.30          # B mm per drop; rescaled so old 0.30x is the new 1.0x
    DROP_F: int = 150
    ZF: int = 600
    Z_SPOOF: float = 26.0           # G92 Z declared at home (needle up); app re-pins to its home Z
    # options
    raster: bool = False
    BRIDGE: bool = False            # strut isolated scaffold islands together
NB = [(1, 0), (-1, 0), (0, 1), (0, -1)]

def components(scaf, cols, rows):
    """4-connected components of filled cells -> list of sets of (c,r)."""
    seen = [[False] * cols for _ in range(rows)]
    comps = []
    for r in range(rows):
        for c in range(cols):
            if scaf[r][c] and not seen[r][c]:
                stack = [(c, r)]; seen[r][c] = True; comp = set()
                while stack:
                    x, y = stack.pop(); comp.add((x, y))
                    for dc, dr in NB:
                        nx, ny = x + dc, y + dr
                        if 0 <= nx < cols and 0 <= ny < rows and scaf[ny][nx] and not seen[ny][nx]:
                            seen[ny][nx] = True; stack.append((nx, ny))
                comps.append(comp)
    return comps

def _line_cells(c1, r1, c2, r2):
    """Bresenham cells from (c1,r1) to (c2,r2) inclusive."""
    cells = []; dx = abs(c2 - c1); dy = abs(r2 - r1)
    sx = 1 if c1 < c2 else -1; sy = 1 if r1 < r2 else -1
    err = dx - dy; c, r = c1, r1
    while True:
        cells.append((c, r))
        if c == c2 and r == r2: break
        e2 = 2 * err
        if e2 > -dy: err -= dy; c += sx
        if e2 < dx: err += dx; r += sy
    return cells

def bridge_islands(scaf, cols, rows):
    """Return a copy of scaf with thin struts added so every filled cell joins one
    4-connected body. Each island is linked to the main body along the shortest
    straight line between their nearest cells."""
    g = [row[:] for row in scaf]
    for _ in range(64):                       # safety cap
        comps = components(g, cols, rows)
        if len(comps) <= 1:
            break
        main = max(comps, key=len)
        best = None
        for comp in comps:
            if comp is main:
                continue
            for (c1, r1) in comp:
                for (c2, r2) in main:
                    d = abs(c1 - c2) + abs(r1 - r2)
                    if best is None or d < best[0]:
                        best = (d, (c1, r1), (c2, r2))
        if not best:
            break
        for (c, r) in _line_cells(best[1][0], best[1][1], best[2][0], best[2][1]):
            g[r][c] = 1
    return g

def trace_path(scaf, cols, rows):
    """Follow the connected line of cells, prefer straight, lift at dead-ends.
    Port of the studio's tracePath. Returns list of strokes; each stroke is [(c,r),...]."""
    on = lambda c, r: 0 <= c < cols and 0 <= r < rows and scaf[r][c]
    def deg(c, r): return sum(on(c + dc, r + dr) for dc, dr in NB)
    vis = [[False] * cols for _ in range(rows)]
    cells = [(c, r) for r in range(rows) for c in range(cols) if scaf[r][c]]
    if not cells:
        return []
    def un(c, r):
        return [(c + dc, r + dr, dc, dr) for dc, dr in NB
                if on(c + dc, r + dr) and not vis[r + dr][c + dc]]
    strokes, last = [], None
    while True:
        rem = [(c, r) for (c, r) in cells if not vis[r][c]]
        if not rem:
            break
        # start at a tip (deg<=1) else nearest to previous end
        best, start = 1e18, rem[0]
        for (c, r) in rem:
            tip = 0 if deg(c, r) <= 1 else 1
            dist = abs(c - last[0]) + abs(r - last[1]) if last else (r * cols + c) * 1e-6
            sc = tip * 1e6 + dist
            if sc < best:
                best, start = sc, (c, r)
        c, r = start
        stroke = [(c, r)]; vis[r][c] = True; d = None
        while True:
            opts = un(c, r)
            if not opts:
                break
            nxt = next((o for o in opts if d and o[2] == d[0] and o[3] == d[1]), None)
            if not nxt:                                  # finish dead-end branches first
                opts.sort(key=lambda o: len(un(o[0], o[1])))
                nxt = opts[0]
            c, r = nxt[0], nxt[1]; vis[r][c] = True
            stroke.append((c, r)); d = (nxt[2], nxt[3])
        strokes.append(stroke); last = stroke[-1]
    return optimize_order(strokes, (0.0, 0.0))     # order strokes to minimise travel (head starts at corner)

def _runs(scaf, cols, rows, axis):
    runs = []
    if axis == "h":
        for r in range(rows):
            c = 0
            while c < cols:
                if scaf[r][c]:
                    a = c
                    while c < cols and scaf[r][c]:
                        c += 1
                    runs.append((r, a, c - 1))
                else:
                    c += 1
    else:
        for c in range(cols):
            r = 0
            while r < rows:
                if scaf[r][c]:
                    a = r
                    while r < rows and scaf[r][c]:
                        r += 1
                    runs.append((c, a, r - 1))
                else:
                    r += 1
    return runs

def raster_path(scaf, cols, rows):
    """Scanline serpentine fill (for solid layers). Returns strokes as polylines."""
    rh, rv = _runs(scaf, cols, rows, "h"), _runs(scaf, cols, rows, "v")
    axis = "v" if len(rv) < len(rh) else "h"
    runs = rh if axis == "h" else rv
    if not runs:
        return []
    lines = {}
    for fixed, lo, hi in runs:
        lines.setdefault(fixed, []).append((lo, hi))
    strokes, cur, flip = [], [], False
    prev = None
    for f in sorted(lines):
        segs = sorted(lines[f])
        if flip:
            segs = list(reversed(segs))
        for lo, hi in segs:
            start, end = (hi, lo) if flip else (lo, hi)
            p = (start, f) if axis == "h" else (f, start)
            q = (end, f) if axis == "h" else (f, end)
            cont = (prev and abs(prev[0] - f) == 1 and
                    max(prev[1], lo) <= min(prev[2], hi))
            if cont:
                cur += [p, q]
            else:
                if cur:
                    strokes.append(cur)
                cur = [p, q]
            prev = (f, lo, hi)
        flip = not flip
    if cur:
        strokes.append(cur)
    return strokes

def optimize_order(strokes, start=(0.0, 0.0)):
    """Greedy nearest-neighbour ordering of polylines WITH reversal, to minimise the
    travel (pen-up) moves between strokes. Each stroke may be printed in either
    direction. start = where the head begins (the parked corner)."""
    rem = [list(s) for s in strokes if s]
    if len(rem) <= 1:
        return rem
    cur = start; out = []
    while rem:
        bi, brev, bd = 0, False, None
        for i, st in enumerate(rem):
            a, b = st[0], st[-1]
            da = (a[0] - cur[0]) ** 2 + (a[1] - cur[1]) ** 2
            db = (b[0] - cur[0]) ** 2 + (b[1] - cur[1]) ** 2
            if bd is None or da < bd: bd, bi, brev = da, i, False
            if db < bd: bd, bi, brev = db, i, True
        st = rem.pop(bi)
        if brev: st = st[::-1]
        out.append(st); cur = st[-1]
    return out

def optimize_points(pts, start=(0.0, 0.0)):
    """Greedy nearest-neighbour then 2-opt refinement over visit points (drops).
    Keeps a fixed virtual start (the parked corner). pts = [(x,y,...)]."""
    rem = list(pts)
    if len(rem) < 2:
        return rem
    cur = start; order = []
    while rem:                                      # greedy NN seed
        i = min(range(len(rem)), key=lambda k: (rem[k][0] - cur[0]) ** 2 + (rem[k][1] - cur[1]) ** 2)
        cur = rem[i]; order.append(rem.pop(i))
    D = lambda p, q: ((p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2) ** 0.5
    n = len(order); improved = True; guard = 0
    while improved and guard < 40:                  # 2-opt
        improved = False; guard += 1
        for i in range(n - 1):
            A = start if i == 0 else order[i - 1]
            B = order[i]
            for j in range(i + 1, n):
                C = order[j]; Dn = order[j + 1] if j + 1 < n else None
                before = D(A, B) + (D(C, Dn) if Dn else 0.0)
                after = D(A, C) + (D(B, Dn) if Dn else 0.0)
                if after + 1e-9 < before:
                    order[i:j + 1] = order[i:j + 1][::-1]
                    improved = True; B = order[i]
    return order

def droplet_path(drops, cols, rows):
    """Visit order over cells with drops, NN + 2-opt to minimise travel. Head starts
    over the parked corner (cell ~0,0)."""
    pts = [(c, r, drops[r][c]) for r in range(rows) for c in range(cols) if drops[r][c]]
    if len(pts) < 2:
        return pts
    return optimize_points(pts, start=(0.0, 0.0))


# ============================================================================ emit
def bridge_strokes(strokes, scaf, cols, rows, maxhop=6):
    """Join consecutive strokes when the end of one reaches the start of the next through a
    short chain of FILLED cells — so the connecting move stays on scaffold and prints as one
    continuous line instead of lifting mid-pattern. Strokes with no filled path within maxhop
    stay separate (they genuinely can't connect without crossing empty space)."""
    from collections import deque
    on = lambda c, r: 0 <= c < cols and 0 <= r < rows and scaf[r][c]
    def filled_path(a, b):
        if a == b: return []
        q = deque([a]); prev = {a: None}
        while q:
            c, r = q.popleft()
            for dc, dr in NB:
                n = (c + dc, r + dr)
                if on(*n) and n not in prev:
                    prev[n] = (c, r)
                    if n == b:
                        path, cur = [], n
                        while cur is not None: path.append(cur); cur = prev[cur]
                        path.reverse()
                        return path[1:] if (len(path) - 1) <= maxhop else None
                    q.append(n)
        return None
    if len(strokes) <= 1:
        return [list(s) for s in strokes]
    out = [list(strokes[0])]
    for st in strokes[1:]:
        bp = filled_path(tuple(out[-1][-1]), tuple(st[0]))
        if bp is not None:                       # bp ends at st[0]; drop it, then append the stroke
            out[-1].extend(bp[:-1]); out[-1].extend(st)
        else:
            out.append(list(st))
    return out

def emit_layers(M, cfg):
    """Per-layer scaffold+droplet bodies for the spray-pause workflow.
    Returns ({'header':[...], 'layers':[{'i','z','lines':[...]}]}, msg).
    Layer formats:  grid -> 'scaffold' bitmap + 'drops' grid ;
                    freeform -> 'paths' (polylines, bed mm) + 'points' [(x,y,n)]."""
    cols, rows = M["grid"]["cols"], M["grid"]["rows"]
    margin, pitch = M["margin"], M["pitch"]
    bead = M.get("bead_height", 1.0)
    mode = M.get("mode", "grid")
    layers = M["layers"]
    bedX = M.get("bed", {}).get("x", cfg.COORD_MAX_X)
    bedY = M.get("bed", {}).get("y", cfg.COORD_MAX_Y)
    XMAX, YMAX = cfg.COORD_MAX_X, cfg.COORD_MAX_Y
    pather = raster_path if cfg.raster else trace_path
    dip = cfg.DIP

    # ---- the ONE mapping: home = bottom-left = (XMAX,YMAX); design grows inward ----
    def MAP_grid(c, r):
        x_from_home = margin + (c + 0.5) * pitch              # rightward from home
        y_from_home = margin + (rows - 1 - r + 0.5) * pitch   # upward from home (row flipped)
        return (XMAX - x_from_home, YMAX - y_from_home)

    def MAP_free(x, y):
        # freeform design coords are canvas mm (y measured down from the top); home is the
        # bottom-left, so x_from_home = x and y_from_home = bedY - y.
        return (XMAX - x, YMAX - (bedY - y))

    # reachable boxes (machine coords) — warn (never clamp) if a design exceeds them
    sc_lo = XMAX - cfg.SCAFFOLD_REACH
    sc_loY = YMAX - cfg.SCAFFOLD_REACH
    dr_hiX = XMAX - cfg.DROPLET_DX
    dr_loX = dr_hiX - cfg.DROPLET_REACH_X
    dr_loY = YMAX - cfg.DROPLET_REACH_Y
    EPS = 0.05

    def layer_geometry(ly):
        # Emit in the SAME order the Studio preview draws, so the G-code matches the
        # on-screen path exactly. trace_path / droplet_path already order their output
        # from the canvas origin (and raster_path keeps its scan order) — the preview
        # uses those directly, so we must NOT re-order after mapping. For freeform we
        # order the drawn paths/points the same way the preview does (design space,
        # from the origin) BEFORE mapping to machine coords.
        if mode == "freeform":
            paths = optimize_order(ly.get("paths", []), (0.0, 0.0))
            pts = optimize_points(ly.get("points", []), (0.0, 0.0))
            strokes = [[MAP_free(x, y) for (x, y) in path] for path in paths if len(path) >= 1]
            stops = [(*MAP_free(x, y), n) for (x, y, n) in pts]
        else:
            scaf = bridge_islands(ly["scaffold"], cols, rows) if cfg.BRIDGE else ly["scaffold"]
            raw = pather(scaf, cols, rows)
            if not cfg.raster:                              # merge strokes that connect through scaffold
                raw = bridge_strokes(raw, scaf, cols, rows)
            strokes = [[MAP_grid(c, r) for (c, r) in st] for st in raw]
            stops = [(*MAP_grid(c, r), n) for (c, r, n) in droplet_path(ly["drops"], cols, rows)]
        return strokes, stops

    header = ["; ===== hybrid bioprinter — layered (spray-pause) =====",
              "; hybrid bioprinter",
              f"; mode = {mode}; HOME = bottom-left = machine max X/Y; design grows inward (machine = MAX - design)",
              "G21", "G90", "M82", "M211 S0", "M302 P1",
              f"G92 X{XMAX:.0f} Y{YMAX:.0f} A0 B0 C0 Z{cfg.Z_SPOOF:.0f}"]
    out, e, eb, ndrops, warns = [], 0.0, 0.0, 0, []
    lift_per_layer = cfg.LAYER_LIFT if cfg.LAYER_LIFT is not None else bead

    for li, ly in enumerate(layers):
        a_print = li * lift_per_layer - cfg.A_OFFSET + cfg.PRINT_LIFT
        clear = max(a_print + cfg.A_HOP, cfg.A_HOP)
        strokes, stops = layer_geometry(ly)
        g = [f"; --- layer {li} scaffold ---"]
        if strokes:
            for (X, Y) in (p for st in strokes for p in st):
                if not (sc_lo - EPS <= X <= XMAX + EPS and sc_loY - EPS <= Y <= YMAX + EPS):
                    warns.append(("scaffold", X, Y))
            g.append(f"G1 A{clear:.3f} F{cfg.AF}   ; raise clear, then travel HOME -> first point")
            for si, st in enumerate(strokes):
                x0, y0 = st[0]
                tag = "   ; travel HOME -> first scaffold point" if si == 0 else ""
                g.append(f"G0 X{x0:.2f} Y{y0:.2f} F{cfg.TRAVEL_F}{tag}")
                g.append(f"G1 A{a_print:.3f} F{cfg.AF}   ; lower to the print plane")
                if cfg.PRIME > 0:
                    e += cfg.PRIME; g.append(f"G1 C{e:.4f} F{cfg.RETRACT_F}   ; prime")
                px, py = x0, y0
                for k in range(1, len(st)):
                    x, y = st[k]
                    e += cfg.FLOW * ((x - px) ** 2 + (y - py) ** 2) ** 0.5
                    g.append(f"G1 X{x:.2f} Y{y:.2f} C{e:.4f} F{cfg.PRINT_F}")
                    if cfg.CORNER_RETRACT > 0 and k < len(st) - 1:
                        nx, ny = st[k + 1]
                        ax, ay = x - px, y - py                 # incoming direction
                        bx, by = nx - x, ny - y                 # outgoing direction
                        la = (ax * ax + ay * ay) ** 0.5; lb = (bx * bx + by * by) ** 0.5
                        if la > 1e-6 and lb > 1e-6 and (ax * bx + ay * by) / (la * lb) < 0.5:
                            # sharp turn (>60 deg): suck the plunger back then restore it, so gel
                            # doesn't pile up while the bed decelerates/re-accelerates through the corner
                            e -= cfg.CORNER_RETRACT; g.append(f"G1 C{e:.4f} F{cfg.RETRACT_F}   ; corner retract")
                            e += cfg.CORNER_RETRACT; g.append(f"G1 C{e:.4f} F{cfg.RETRACT_F}   ; corner prime")
                    px, py = x, y
                if cfg.RETRACT > 0:
                    e -= cfg.RETRACT; g.append(f"G1 C{e:.4f} F{cfg.RETRACT_F}   ; retract")
                g.append(f"G1 A{clear:.3f} F{cfg.AF}   ; raise between strokes")
        g.append(f"G1 A{max(a_print + bead + cfg.A_HOP, cfg.A_HOP):.3f} F{cfg.AF}   ; lift clear for travel / droplets")
        if stops:
            g.append(f"; --- layer {li} droplets ({len(stops)} stops) ---")
            g.append(f"G0 X{XMAX:.0f} Y{YMAX:.0f} F{cfg.TRAVEL_F}   ; bed home (scaffold over corner)")
            g.append(f"G0 X{dr_hiX:.2f} Y{YMAX - cfg.DROPLET_DY:.2f} F{cfg.TRAVEL_F}"
                     f"   ; bed -{cfg.DROPLET_DX:.0f}X/-{cfg.DROPLET_DY:.0f}Y -> droplet head over HOME")
            for si, (X, Y, n) in enumerate(stops):
                ndrops += int(n)
                xd, yd = X - cfg.DROPLET_DX, Y - cfg.DROPLET_DY   # droplet = scaffold - head offset (X and Y)
                if not (dr_loX - EPS <= xd <= dr_hiX + EPS and dr_loY - EPS <= yd <= YMAX + EPS):
                    warns.append(("droplet", xd, yd))
                tag = "   ; travel HOME -> first droplet point" if si == 0 else ""
                g.append(f"G0 X{xd:.2f} Y{yd:.2f} F{cfg.TRAVEL_F}{tag}")
                g.append("M400")
                g += ["G91", f"G1 Z-{dip:.2f} F{cfg.ZF}   ; lower needle", "G90"]
                # fire each droplet SEPARATELY so they detach one at a time instead of merging into
                # one blob. Between drops at the same spot, dwell DROP_GAP so the drop clears the tip.
                ndrop = int(n)
                for k in range(ndrop):
                    eb += cfg.DROP_VOL
                    g.append(f"G1 B{eb:.4f} F{cfg.DROP_F}   ; push plunger -> drop {k+1}/{ndrop}")
                    g.append("M400")
                    if cfg.DROP_HOLD > 0:            # needle-down dwell so the drop forms/detaches
                        g.append(f"; __HOLD__ {cfg.DROP_HOLD}")
                    if k < ndrop - 1 and cfg.DROP_GAP > 0:   # pause before the next drop at this spot
                        g.append(f"; __HOLD__ {cfg.DROP_GAP}")
                g += ["G91", f"G1 Z{dip:.2f} F{cfg.ZF}   ; retract needle", "G90"]
        if li < len(layers) - 1:            # between layers: step the droplet needle up too, but by
            # its OWN amount (the real gel growth per layer), NOT the scaffold bead — the scaffold
            # rises by bead for clearance, while the needle must track the actual gel surface height.
            drop_rise = cfg.DROPLET_RISE if cfg.DROPLET_RISE is not None else bead
            g += ["G91", f"G1 Z{drop_rise:.3f} F{cfg.ZF}   ; needle up one gel-layer for next layer", "G90"]
        out.append({"i": li, "z": ly.get("z", li * bead), "lines": g})

    msg = f"{mode}: {len(layers)} layer(s) · {ndrops} drops · B-dispense, no solenoid"
    allx = [float(t[1:]) for ly in out for L in ly["lines"] if L[:2] in ("G0", "G1")
            for t in L.split() if t[:1] == "X" and _isnum(t[1:])]
    ally = [float(t[1:]) for ly in out for L in ly["lines"] if L[:2] in ("G0", "G1")
            for t in L.split() if t[:1] == "Y" and _isnum(t[1:])]
    if allx and ally:
        msg += f"  | travel X {min(allx):.1f}..{max(allx):.1f}  Y {min(ally):.1f}..{max(ally):.1f} (home {XMAX:.0f},{YMAX:.0f})"
    if warns:
        wx = [w[1] for w in warns]; wy = [w[2] for w in warns]
        msg += (f"  ⚠ {len(warns)} point(s) outside the reachable zone "
                f"(X {min(wx):.1f}..{max(wx):.1f}, Y {min(wy):.1f}..{max(wy):.1f}); shrink the design.")
    return ({"header": header, "layers": out}, msg)


def _isnum(s):
    try: float(s); return True
    except ValueError: return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("matrix")
    ap.add_argument("--raster", action="store_true")
    a = ap.parse_args()
    M = json.load(open(a.matrix))
    payload, msg = emit_layers(M, Cfg(raster=a.raster))
    prefix = a.matrix.rsplit(".", 1)[0]
    for ly in payload["layers"]:
        open(f"{prefix}_layer{ly['i']}.gcode", "w").write(
            "\n".join(payload["header"] + ly["lines"]) + "\n")
    print(f"{prefix}: {msg}")

m2g = sys.modules[__name__]

# ===================== Studio canvas (studio_tab.py) =====================
#!/usr/bin/env python3
"""
studio_tab.py — native Tkinter "Studio" tab for printess_patched.py.

Replaces the separate Scaffold + Matrix tabs with one canvas where you:
  - paint scaffold cells and place droplet counts on a bed-accurate grid,
  - load an STL and slice it into layers by bead height,
  - flip through layers, preview the scaffold strokes + droplet visit order,
  - generate scaffold + droplet G-code in-app (with the 27mm reach guard).

Drop-in usage inside the main app:
    from studio_tab import StudioTab
    studio = StudioTab(notebook, on_generate=self.load_studio_gcode)
    notebook.add(studio, text="Studio")
where on_generate(scaffold_gcode:str|None, droplet_gcode:str|None, msg:str) hands
the strings back to the host app's print pipeline. If on_generate is None the tab
just writes <prefix>_scaffold.gcode / _droplet.gcode next to this file.

Toolpath + emission logic is reused from matrix_to_gcode.py so the on-screen
preview equals the printed path. Keep both files in the same folder.
"""

import os, struct, re, math, json, copy
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, simpledialog

# ----------------------------------------------------------------- STL slicing
def parse_stl(path):
    with open(path, "rb") as f:
        buf = f.read()
    tris = []
    if len(buf) > 84:
        (n,) = struct.unpack("<I", buf[80:84])
        if len(buf) == 84 + n * 50:                       # binary
            off = 84
            for _ in range(n):
                off += 12                                  # skip normal
                v = []
                for _ in range(3):
                    v.append(struct.unpack("<3f", buf[off:off + 12])); off += 12
                off += 2
                tris.append([tuple(p) for p in v])
            return tris
    txt = buf.decode("utf-8", "ignore")                    # ascii fallback
    vs = [(float(a), float(b), float(c)) for a, b, c in
          re.findall(r"vertex\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)", txt)]
    for i in range(0, len(vs) - 2, 3):
        tris.append([vs[i], vs[i + 1], vs[i + 2]])
    return tris

def _remap(p, up):
    x, y, z = p
    return (x, y, z) if up == "z" else (x, z, y) if up == "y" else (y, z, x)

def _chain_segments(segs, tol=1e-3):
    """Chain unordered 2D segments into polylines by endpoint proximity."""
    from collections import defaultdict
    key = lambda p: (round(p[0] / tol), round(p[1] / tol))
    adj = defaultdict(list)
    for i, (a, b) in enumerate(segs):
        adj[key(a)].append((i, a, b)); adj[key(b)].append((i, b, a))
    used = [False] * len(segs); polys = []
    for i in range(len(segs)):
        if used[i]: continue
        a, b = segs[i]; used[i] = True; poly = [a, b]; cur = b
        while True:
            nxt = next(((j, p, q) for (j, p, q) in adj[key(cur)] if not used[j]), None)
            if not nxt: break
            j, p, q = nxt; used[j] = True; poly.append(q); cur = q
            if key(cur) == key(poly[0]): break
        polys.append(poly)
    return polys

def _fill_to_snake(segs):
    """Connect horizontal scanline segments (from _poly_infill, sorted by y) into continuous
    zigzag snake paths. Breaks into a new path whenever connecting to the next row would leave
    the shape (multi-span row, or no x-overlap), so the deposited path stays on the solid."""
    paths, cur, flip = [], [], False
    for (a0, b0) in segs:
        a, b = (a0, b0) if not flip else (b0, a0)
        if cur:
            px, py = cur[-1]
            lo, hi = min(a[0], b[0]), max(a[0], b[0])
            if abs(a[1] - py) > 1e-9 and lo - 1e-6 <= px <= hi + 1e-6:
                cur.append(a); cur.append(b)          # step across to the next row, inside the shape
            else:
                paths.append(cur); cur = [a, b]
        else:
            cur = [a, b]
        flip = not flip
    if len(cur) >= 2: paths.append(cur)
    return paths

def _poly_infill(loops, spacing, horizontal):
    """Scanline infill segments inside closed contour loops (even-odd rule).
    Returns a list of 2-point paths. Direction alternates per layer for a lattice."""
    if not loops or spacing <= 0: return []
    edges = []
    for poly in loops:
        if len(poly) < 2: continue
        pts = poly if poly[0] == poly[-1] else poly + [poly[0]]
        for i in range(len(pts) - 1):
            edges.append((pts[i], pts[i + 1]))
    if not edges: return []
    xs = [p[0] for e in edges for p in e]; ys = [p[1] for e in edges for p in e]
    segs = []
    if horizontal:
        lo, hi = min(ys), max(ys); v = lo + spacing / 2
        while v < hi:
            xint = []
            for (a, b) in edges:
                y0, y1 = a[1], b[1]
                if (y0 <= v < y1) or (y1 <= v < y0):          # half-open avoids double-counting vertices
                    xint.append(a[0] + (v - y0) / (y1 - y0) * (b[0] - a[0]))
            xint.sort()
            for i in range(0, len(xint) - 1, 2):
                segs.append([(xint[i], v), (xint[i + 1], v)])
            v += spacing
    else:
        lo, hi = min(xs), max(xs); v = lo + spacing / 2
        while v < hi:
            yint = []
            for (a, b) in edges:
                x0, x1 = a[0], b[0]
                if (x0 <= v < x1) or (x1 <= v < x0):
                    yint.append(a[1] + (v - x0) / (x1 - x0) * (b[1] - a[1]))
            yint.sort()
            for i in range(0, len(yint) - 1, 2):
                segs.append([(v, yint[i]), (v, yint[i + 1])])
            v += spacing
    return segs

def slice_stl_contours(tris, bedX, bedY, margin, bead, scale=1.0, fit=False, up="z",
                       infill=True, infill_spacing=3.0):
    """Slice a mesh into per-layer paths (mm on the bed) for freeform mode: perimeter
    contours plus optional scanline infill (alternating direction per layer so the solid
    actually fills instead of showing only its outline)."""
    T = [[_remap(p, up) for p in t] for t in tris]
    pts = [p for t in T for p in t]
    if not pts: return [], "no triangles"
    minx = min(p[0] for p in pts); maxx = max(p[0] for p in pts)
    miny = min(p[1] for p in pts); maxy = max(p[1] for p in pts)
    minz = min(p[2] for p in pts); maxz = max(p[2] for p in pts)
    w, h, depth = maxx - minx, maxy - miny, maxz - minz
    availX, availY = bedX - 2 * margin, bedY - 2 * margin
    s = scale
    if fit and w > 0 and h > 0:
        s = min(availX / w, availY / h)
    nL = max(1, int(round(depth / bead)) if bead > 0 else 1)
    to_mm = lambda x, y: (margin + (x - minx) * s, margin + (y - miny) * s)
    sp_mm = max(0.1, infill_spacing)                         # spacing requested in bed mm
    layers = []; ninf = 0
    for li in range(nL):
        z = minz + (li + 0.5) * (depth / nL if nL else depth)
        segs = []
        for tri in T:
            cross = []
            for k in range(3):
                a, b = tri[k], tri[(k + 1) % 3]
                if (a[2] - z) * (b[2] - z) < 0:               # edge strictly crosses the plane
                    t = (z - a[2]) / (b[2] - a[2])
                    cross.append((a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])))
            if len(cross) == 2:
                segs.append((cross[0], cross[1]))
        loops = [[to_mm(x, y) for (x, y) in poly] for poly in _chain_segments(segs)]
        paths = list(loops)
        if infill:
            fill = _poly_infill(loops, sp_mm, horizontal=(li % 2 == 0))
            ninf += len(fill); paths += fill
        layers.append({"z": round((li + 0.5) * bead, 3), "paths": paths, "points": []})
    npaths = sum(len(L["paths"]) for L in layers)
    tag = f" · {ninf} infill" if infill else ""
    return layers, f"{len(tris)} tris · {nL} layers · {npaths} paths{tag} · {s*w:.0f}\u00d7{s*h:.0f}mm"

def slice_stl(tris, bedX, bedY, pitch, margin, bead, scale=1.0, fit=False, up="z"):
    tris = [[_remap(p, up) for p in t] for t in tris]
    xs = [p[0] for t in tris for p in t]; ys = [p[1] for t in tris for p in t]
    zs = [p[2] for t in tris for p in t]
    mnx, mxx, mny, mxy, mnz, mxz = min(xs), max(xs), min(ys), max(ys), min(zs), max(zs)
    usableW, usableH = bedX - 2 * margin, bedY - 2 * margin
    s = scale * (min(usableW / ((mxx - mnx) or 1), usableH / ((mxy - mny) or 1)) if fit else 1.0)
    fw, fh = (mxx - mnx) * s, (mxy - mny) * s
    ox, oy = margin + (usableW - fw) / 2, margin + (usableH - fh) / 2
    cols, rows = max(1, int(usableW // pitch)), max(1, int(usableH // pitch))
    nL = max(1, math.ceil((mxz - mnz) / bead))
    layers = []
    for k in range(nL):
        z = mnz + (k + 0.5) * bead
        segs = []
        for t in tris:
            pts = []
            for e in range(3):
                a, b = t[e], t[(e + 1) % 3]
                da, db = a[2] - z, b[2] - z
                if (da <= 0 < db) or (db <= 0 < da):
                    u = da / (da - db)
                    pts.append((a[0] + (b[0] - a[0]) * u, a[1] + (b[1] - a[1]) * u))
            if len(pts) == 2:
                A = (ox + (pts[0][0] - mnx) * s, oy + (pts[0][1] - mny) * s)
                B = (ox + (pts[1][0] - mnx) * s, oy + (pts[1][1] - mny) * s)
                segs.append((A, B))
        scaf = [[0] * cols for _ in range(rows)]
        for r in range(rows):
            yL = margin + (r + 0.5) * pitch
            xr = []
            for (A, B) in segs:
                y0, y1 = A[1], B[1]
                if (y0 <= yL < y1) or (y1 <= yL < y0):
                    u = (yL - y0) / (y1 - y0); xr.append(A[0] + (B[0] - A[0]) * u)
            xr.sort()
            for c in range(cols):
                cx = margin + (c + 0.5) * pitch
                if sum(1 for x in xr if x < cx) & 1:
                    scaf[r][c] = 1
        layers.append({"scaffold": scaf, "drops": [[0] * cols for _ in range(rows)],
                       "z": round((k + 0.5) * bead, 3)})
    info = f"{len(tris)} tris · {fw:.1f}×{fh:.1f}mm · {nL} layers"
    return cols, rows, layers, info

# ----------------------------------------------------------------- Studio tab
ALG, DROP, PATH, TRAVEL, REACH, BG, GRID = (
    "#3FB8AF", "#F2A03D", "#7AA2F7", "#3a4d6b", "#E5484D", "#10151b", "#26323c")

class StudioTab(ttk.Frame):
    def __init__(self, master, on_generate=None, app_params=None):
        super().__init__(master)
        self.on_generate = on_generate
        # optional (getter, setter) for app-side gel params (suction, extrude F) so
        # profiles capture them too. getter() -> dict ; setter(dict)
        self.app_get, self.app_set = (app_params or (None, None))
        # params
        self.bedX, self.bedY, self.pitch, self.margin, self.bead = 52, 60, 3.0, 0, 1.0
        self.offX, self.dead = 31, 0
        self.cols = self.rows = 0
        self.layers = []
        self._undo_stack = []          # snapshots of (active, layers) for Ctrl+Z / Undo
        self.active = 0
        self.mode = tk.StringVar(value="scaffold")
        self.fmode = tk.StringVar(value="grid")      # "grid" or "freeform"
        self._cur_path = None                         # in-progress freeform polyline (mm)
        self.ff_tool = tk.StringVar(value="freehand")  # freeform tool: freehand | line
        self._line_pts = []                            # in-progress straight-line polyline vertices (mm)
        self._ff = (0.0, 0.0, 1.0)                     # freeform transform: ox, oy, px-per-mm
        self.brush = tk.IntVar(value=1)
        self.up_axis = tk.StringVar(value="auto")   # STL slicing axis
        self._stl_tris = None                        # last-loaded mesh, for re-slicing
        self.v_path_s = tk.BooleanVar(value=True)
        self.v_path_d = tk.BooleanVar(value=False)
        self.v_raster = tk.BooleanVar(value=False)
        self.v_onion = tk.BooleanVar(value=False)
        self.o_bridge = tk.BooleanVar(value=False)   # auto-connect isolated islands
        self.v_infill = tk.BooleanVar(value=True)    # fill freeform STL solids (not just outline)
        self.v_noretract = tk.BooleanVar(value=True)  # disable ALL C-plunger retract/prime (Pluronic bubbles)
        # output params (combined mode is always native; these feed the compiler)
        self.o_flow = tk.DoubleVar(value=0.024)  # 0.02 x 1.2 baked in: old 1.20x scaffold flow is the new 1.0x
        self.o_corner = tk.DoubleVar(value=0.05) # C retract at sharp corners to kill blobs (0 = off)
        self.o_prime = tk.DoubleVar(value=0.6)   # C pre-charge at stroke start (lower = less start glob)
        self.o_retract = tk.DoubleVar(value=0.5) # C pull-back at stroke end (higher = less ooze/string)
        self.o_hx = tk.DoubleVar(value=51.5)     # head offset X (solenoid head; -2.5 from 51 for the +X miss)
        self.o_hy = tk.DoubleVar(value=7.0)      # head offset Y (solenoid head; calibrated)
        self.o_zprint = tk.DoubleVar(value=1.0)  # needle height at the gel (spoof frame)
        self.o_zsafe = tk.DoubleVar(value=5.0)   # needle clearance between drops
        self.o_zspoof = tk.DoubleVar(value=26.0) # G92 Z spoof = your Set-Home droplet Z
        self.o_hop = tk.DoubleVar(value=6.0)     # A lift between scaffold strokes (clear printed gel)
        self.o_depth = tk.DoubleVar(value=0.0)   # lowers scaffold print plane (+ = dig into gel)
        self.o_printlift = tk.DoubleVar(value=0.3) # raises nozzle above the print plane each layer (anti-drag)
        self.o_home = tk.StringVar(value="BL")   # canvas corner the head sits over at home
        self.o_rotate = tk.IntVar(value=0)       # orient design CW on bed: 0/90/180/270
        self.o_flipx = tk.BooleanVar(value=False)  # flip about X axis (top<->bottom)
        self.o_flipy = tk.BooleanVar(value=False)  # flip about Y axis (left<->right)
        self.o_dip = tk.DoubleVar(value=0.0)     # droplet needle dip at each spot (dispenses from hover)
        self.o_droprise = tk.DoubleVar(value=0.5) # needle Z rise per layer = REAL gel growth (< bead h)
        self.o_drophold = tk.DoubleVar(value=1.0) # seconds needle holds at depth after pushing a droplet
        self.o_dropgap = tk.DoubleVar(value=0.5)  # seconds between drops at the same spot (each detaches)
        self.o_dropmm = tk.DoubleVar(value=0.30) # B mm/drop; rescaled so old 0.30x droplet flow is the new 1.0x
        self.o_pulse = tk.IntVar(value=80)       # FAN0 open time per drop (ms)
        self.o_gap = tk.IntVar(value=120)        # closed dwell between drops (ms)
        self.cellpx, self.PAD = 24, 18

        self._build()
        self._rebuild_grid()

    # ---------------- UI ----------------
    def _build(self):
        self.columnconfigure(0, weight=1); self.rowconfigure(0, weight=1)
        self.canvas = tk.Canvas(self, bg=BG, highlightthickness=0)
        self.canvas.grid(row=0, column=0, sticky="nsew", padx=(6, 0), pady=6)
        self.canvas.bind("<Configure>", lambda e: self._draw())
        self.canvas.bind("<Button-1>", self._on_press)
        self.canvas.bind("<Double-Button-1>", self._line_finish)
        self.canvas.bind("<Motion>", self._on_hover)
        self.canvas.bind("<B1-Motion>", self._on_drag)
        self.canvas.bind("<ButtonRelease-1>", self._on_release)
        self.canvas.bind("<Button-3>", self._on_rclick)
        self.canvas.bind("<Control-z>", self._undo)
        self.canvas.bind("<Control-Z>", self._undo)
        self.canvas.bind("<B3-Motion>", self._on_rdrag)

        # --- scrollable side panel so every section (incl. Output) is reachable ---
        outer = ttk.Frame(self); outer.grid(row=0, column=1, sticky="ns", padx=6, pady=6)
        sc = tk.Canvas(outer, width=232, bg=BG, highlightthickness=0)
        vsb = ttk.Scrollbar(outer, orient="vertical", command=sc.yview)
        sc.configure(yscrollcommand=vsb.set)
        vsb.pack(side="right", fill="y"); sc.pack(side="left", fill="both", expand=True)
        side = ttk.Frame(sc)
        win = sc.create_window((0, 0), window=side, anchor="nw")
        side.bind("<Configure>", lambda e: sc.configure(scrollregion=sc.bbox("all")))
        sc.bind("<Configure>", lambda e: sc.itemconfig(win, width=e.width))
        def _wheel(e): sc.yview_scroll(int(-e.delta / 120), "units")
        outer.bind("<Enter>", lambda e: sc.bind_all("<MouseWheel>", _wheel))
        outer.bind("<Leave>", lambda e: sc.unbind_all("<MouseWheel>"))

        def hdr(t): ttk.Label(side, text=t, font=("Segoe UI", 8, "bold"),
                              foreground="#7b8893").pack(anchor="w", pady=(10, 2))

        hdr("MODE")
        modef = ttk.Frame(side); modef.pack(fill="x")
        ttk.Radiobutton(modef, text="Grid", value="grid", variable=self.fmode,
                        command=self._switch_mode).pack(side="left")
        ttk.Radiobutton(modef, text="Freeform", value="freeform", variable=self.fmode,
                        command=self._switch_mode).pack(side="left")

        hdr("PAINT")
        mf = ttk.Frame(side); mf.pack(fill="x")
        ttk.Radiobutton(mf, text="Scaffold", value="scaffold", variable=self.mode).pack(side="left")
        ttk.Radiobutton(mf, text="Droplets", value="droplet", variable=self.mode).pack(side="left")
        bf = ttk.Frame(side); bf.pack(fill="x", pady=2)
        ttk.Label(bf, text="drops/cell").pack(side="left")
        ttk.Spinbox(bf, from_=0, to=9, width=4, textvariable=self.brush).pack(side="left", padx=4)
        ttk.Button(side, text="Fill drawn shape (solid)", command=self._fill_shape).pack(fill="x", pady=(2, 0))
        ttk.Button(side, text="\u21b6 Undo  (Ctrl+Z)", command=self._undo).pack(fill="x", pady=(2, 0))
        tf = ttk.Frame(side); tf.pack(fill="x", pady=(2, 0))
        ttk.Label(tf, text="draw:").pack(side="left")
        ttk.Radiobutton(tf, text="Freehand", value="freehand", variable=self.ff_tool).pack(side="left")
        ttk.Radiobutton(tf, text="Line", value="line", variable=self.ff_tool).pack(side="left")
        ttk.Label(side, text="Line: click each corner, double-click to finish.",
                  font=("Segoe UI", 8), foreground="#888").pack(anchor="w")

        # legend (colour key like the standalone app)
        leg = ttk.Frame(side); leg.pack(fill="x", pady=(4, 0))
        for color, txt in [(ALG, "scaffold"), (DROP, "drops"), (PATH, "stroke path"), (REACH, "dead-zone")]:
            r = ttk.Frame(leg); r.pack(anchor="w")
            tk.Label(r, bg=color, width=2, height=1).pack(side="left", padx=(0, 5), pady=1)
            ttk.Label(r, text=txt, font=("Segoe UI", 8)).pack(side="left")

        hdr("STL")
        ttk.Button(side, text="Load STL \u2192 slice", command=self._load_stl).pack(fill="x")
        uf = ttk.Frame(side); uf.pack(fill="x", pady=(2, 0))
        ttk.Label(uf, text="up axis").pack(side="left")
        ttk.OptionMenu(uf, self.up_axis, "auto", "auto", "z", "y", "x",
                       command=lambda *_: self._reslice()).pack(side="left", padx=4)
        ttk.Checkbutton(side, text="Fill solids (100% infill)", variable=self.v_infill,
                        command=self._reslice).pack(anchor="w")
        self.stl_lbl = ttk.Label(side, text="no mesh", foreground="#888", wraplength=210); self.stl_lbl.pack(anchor="w")

        hdr("BED / GRID")
        self.e = {}
        for key, lab, val in [("bedX", "bed X", 52), ("bedY", "bed Y", 60),
                              ("pitch", "pitch", 3), ("margin", "margin", 0),
                              ("bead", "bead h", 1.0), ("dead", "dead cols", 0)]:
            row = ttk.Frame(side); row.pack(fill="x")
            ttk.Label(row, text=lab, width=10).pack(side="left")
            var = tk.StringVar(value=str(val)); self.e[key] = var
            ttk.Entry(row, textvariable=var, width=8).pack(side="left")
        ttk.Button(side, text="Apply grid", command=self._apply).pack(fill="x", pady=2)

        hdr("LAYERS")
        llf = ttk.Frame(side); llf.pack(fill="x")
        lsb = ttk.Scrollbar(llf, orient="vertical")
        self.layer_list = tk.Listbox(llf, height=6, yscrollcommand=lsb.set, bg="#161C22", fg=ALG,
                                     selectbackground=ALG, selectforeground="#0E1216",
                                     highlightthickness=0, font=("Consolas", 9), exportselection=False)
        lsb.config(command=self.layer_list.yview)
        lsb.pack(side="right", fill="y")
        self.layer_list.pack(side="left", fill="both", expand=True)
        self.layer_list.bind("<<ListboxSelect>>", self._on_layer_pick)
        lf = ttk.Frame(side); lf.pack(fill="x", pady=2)
        ttk.Button(lf, text="\u25b2", width=3, command=lambda: self._set_active(self.active + 1)).pack(side="left")
        ttk.Button(lf, text="\u25bc", width=3, command=lambda: self._set_active(self.active - 1)).pack(side="left")
        ttk.Button(lf, text="\uff0b", width=3, command=self._add_layer).pack(side="left", padx=(8, 0))
        ttk.Button(lf, text="\u29c9", width=3, command=self._dup_layer).pack(side="left")
        ttk.Button(lf, text="\uff0d", width=3, command=self._del_layer).pack(side="left")
        lf2 = ttk.Frame(side); lf2.pack(fill="x")
        ttk.Button(lf2, text="Clear layer", command=self._clear_layer).pack(side="left", fill="x", expand=True)
        ttk.Button(lf2, text="Reset all", command=self._reset_all).pack(side="left", fill="x", expand=True)

        hdr("VIEW")
        for txt, var in [("Scaffold path", self.v_path_s), ("Droplet path", self.v_path_d),
                         ("Raster fill", self.v_raster), ("Onion below", self.v_onion),
                         ("Bridge islands", self.o_bridge),
                         ("No retractions", self.v_noretract)]:
            ttk.Checkbutton(side, text=txt, variable=var, command=self._draw).pack(anchor="w")

        hdr("OUTPUT (per-layer)")
        for lab, var in [("flow", self.o_flow),
                         ("prime", self.o_prime),
                         ("retract", self.o_retract),
                         ("corner retract", self.o_corner),
                         ("droplet \u2212X", self.o_hx),
                         ("droplet \u2212Y", self.o_hy),
                         ("stroke lift", self.o_hop),
                         ("print depth", self.o_depth),
                         ("print lift", self.o_printlift),
                         ("droplet dip", self.o_dip), ("droplet rise", self.o_droprise),
                         ("drop hold s", self.o_drophold),
                         ("drop gap s", self.o_dropgap),
                         ("drop B mm", self.o_dropmm)]:
            row = ttk.Frame(side); row.pack(fill="x")
            ttk.Label(row, text=lab, width=10).pack(side="left")
            ttk.Entry(row, textvariable=var, width=8).pack(side="left")
        ttk.Button(side, text="Generate G-code", command=self._generate, style="Accent.TButton").pack(fill="x", pady=(6, 2))
        self.status = ttk.Label(side, text="", foreground=ALG, wraplength=210, justify="left")
        self.status.pack(anchor="w")

        hdr("GEL PROFILE")
        self.profile_var = tk.StringVar()
        self.profile_combo = ttk.Combobox(side, textvariable=self.profile_var, state="readonly")
        self.profile_combo.pack(fill="x")
        pbf = ttk.Frame(side); pbf.pack(fill="x", pady=2)
        ttk.Button(pbf, text="Save", command=self._profile_save).pack(side="left", fill="x", expand=True)
        ttk.Button(pbf, text="Load", command=self._profile_load).pack(side="left", fill="x", expand=True)
        self._refresh_profiles()

        hdr("IMPORT")
        imf = ttk.Frame(side); imf.pack(fill="x")
        ttk.Button(imf, text="STL\u2026", command=self._load_stl).pack(side="left", fill="x", expand=True)
        ttk.Button(imf, text="Job\u2026", command=self._job_import).pack(side="left", fill="x", expand=True)
        ttk.Button(imf, text="Profile\u2026", command=self._profile_import).pack(side="left", fill="x", expand=True)

        hdr("EXPORT")
        exf = ttk.Frame(side); exf.pack(fill="x")
        ttk.Button(exf, text="Job\u2026", command=self._job_export).pack(side="left", fill="x", expand=True)
        ttk.Button(exf, text="Profile\u2026", command=self._profile_export).pack(side="left", fill="x", expand=True)
        ttk.Button(exf, text="G-code\u2026", command=self._gcode_export).pack(side="left", fill="x", expand=True)

    def _on_layer_pick(self, e):
        sel = self.layer_list.curselection()
        if sel:
            idx = (len(self.layers) - 1) - sel[0]   # rail is top=highest
            if idx != self.active:
                self._set_active(idx)

    def _refresh_layer_list(self):
        self.layer_list.delete(0, tk.END)
        n = len(self.layers)
        for i in range(n - 1, -1, -1):              # highest layer on top
            ly = self.layers[i]
            ns = any(any(row) for row in ly["scaffold"])
            nd = any(any(row) for row in ly["drops"])
            tag = ("\u25cf" if ns else "\u00b7") + ("\u25cf" if nd else " ")
            self.layer_list.insert(tk.END, f"L{i}  z{ly['z']:.1f}  {tag}")
        if 0 <= self.active < n:
            row = (n - 1) - self.active
            self.layer_list.selection_clear(0, tk.END)
            self.layer_list.selection_set(row)
            self.layer_list.see(row)

    # ---------------- grid state ----------------
    def _blank(self, z=0.0):
        return {"scaffold": [[0] * self.cols for _ in range(self.rows)],
                "drops": [[0] * self.cols for _ in range(self.rows)],
                "paths": [], "points": [], "z": z}

    def _apply(self):
        try:
            self.bedX = float(self.e["bedX"].get()); self.bedY = float(self.e["bedY"].get())
            self.pitch = float(self.e["pitch"].get()); self.margin = float(self.e["margin"].get())
            self.bead = float(self.e["bead"].get())
            self.dead = max(0, int(float(self.e["dead"].get())))
        except ValueError:
            messagebox.showerror("Studio", "grid fields must be numbers"); return
        self._rebuild_grid()

    def _rebuild_grid(self):
        nc = max(1, int((self.bedX - 2 * self.margin) // self.pitch))
        nr = max(1, int((self.bedY - 2 * self.margin) // self.pitch))
        old = self.layers
        self.cols, self.rows = nc, nr
        if not old:
            self.layers = [self._blank(0.0)]
        else:  # resample preserving overlap
            new = []
            for ly in old:
                nl = self._blank(ly["z"])
                for r in range(min(len(ly["scaffold"]), nr)):
                    for c in range(min(len(ly["scaffold"][0]), nc)):
                        nl["scaffold"][r][c] = ly["scaffold"][r][c]
                        nl["drops"][r][c] = ly["drops"][r][c]
                nl["paths"] = [p[:] for p in ly.get("paths", [])]   # mm-based, grid-independent
                nl["points"] = list(ly.get("points", []))
                new.append(nl)
            self.layers = new
        self.active = min(self.active, len(self.layers) - 1)
        self._draw(); self._refresh_labels()

    # ---------------- layers ----------------
    def _set_active(self, i):
        self.active = max(0, min(len(self.layers) - 1, i)); self._draw(); self._refresh_labels()
    def _add_layer(self):
        self.layers.insert(self.active + 1, self._blank()); self.active += 1; self._renz(); self._draw(); self._refresh_labels()
    def _dup_layer(self):
        cur = self.layers[self.active]
        nl = self._blank()
        nl["scaffold"] = [row[:] for row in cur["scaffold"]]; nl["drops"] = [row[:] for row in cur["drops"]]
        nl["paths"] = [p[:] for p in cur.get("paths", [])]; nl["points"] = list(cur.get("points", []))
        self.layers.insert(self.active + 1, nl); self.active += 1; self._renz(); self._draw(); self._refresh_labels()
    def _del_layer(self):
        if len(self.layers) <= 1:
            self.layers = [self._blank()]
        else:
            self.layers.pop(self.active); self.active = min(self.active, len(self.layers) - 1)
        self._renz(); self._draw(); self._refresh_labels()
    def _renz(self):
        for i, ly in enumerate(self.layers): ly["z"] = round((i + 0.5) * self.bead, 3)

    def _refresh_labels(self):
        L = self.layers[self.active]
        self._refresh_layer_list()
        if self.fmode.get() == "freeform":
            paths = L.get("paths", []); pts = L.get("points", [])
            npts = sum(len(p) for p in paths)
            ndrop = sum(max(1, n) for (_x, _y, n) in pts)
            self.status.config(text=f"L{self.active}/{len(self.layers)-1}: {len(paths)} paths "
                                    f"({npts} pts) \u00b7 {len(pts)} drop spots ({ndrop} drops)")
            return
        ns = sum(sum(row) for row in L["scaffold"]); nd = sum(sum(row) for row in L["drops"])
        sp = (m2g.raster_path if self.v_raster.get() else m2g.trace_path)(L["scaffold"], self.cols, self.rows)
        strays = sum(1 for comp in m2g.components(L["scaffold"], self.cols, self.rows) if len(comp) == 1)
        warn = f" \u00b7 {strays} stray\u26a0" if strays and not self.o_bridge.get() else ""
        self.status.config(text=f"L{self.active}/{len(self.layers)-1}: {ns} cells \u00b7 {nd} drops \u00b7 "
                                f"{len(sp)} strokes/{max(0,len(sp)-1)} lifts{warn}")

    # ---------------- painting ----------------
    def _switch_mode(self):
        self._cur_path = None
        if self._stl_tris:                 # re-derive the right representation for the new mode
            self._slice_current()
        else:
            self._draw(); self._refresh_labels()

    def _on_press(self, e):
        self.canvas.focus_set()                        # so Ctrl+Z reaches the canvas
        if self.fmode.get() == "grid":
            rc = self._cell_at(e.x, e.y)
            self._grid_erase = False           # erase this stroke if the first cell is already filled
            if rc:
                c, r = rc; L = self.layers[self.active]
                filled = (L["scaffold"][r][c] if self.mode.get() == "scaffold" else L["drops"][r][c])
                self._grid_erase = bool(filled)
            self._push_undo()                          # snapshot before the paint stroke
            self._paint(e, not self._grid_erase); return
        x, y = self._px_to_mm(e.x, e.y)
        L = self.layers[self.active]
        if self.mode.get() == "scaffold" and self.ff_tool.get() == "line":
            self._line_pts.append((x, y))      # straight-line tool: each click drops a corner
            self._draw(); return               # (undo snapshot happens when the line is committed)
        if self.mode.get() == "scaffold":
            self._push_undo()                          # snapshot before a freehand stroke
            self._cur_path = [(x, y)]
        else:
            self._push_undo()                          # snapshot before adding a droplet
            L.setdefault("points", []).append((x, y, max(1, self.brush.get())))
            self._draw(); self._refresh_labels()

    def _push_undo(self):
        """Snapshot the whole drawing state before a change, so Undo / Ctrl+Z can restore it."""
        try:
            self._undo_stack.append((self.active, copy.deepcopy(self.layers)))
            if len(self._undo_stack) > 40:
                self._undo_stack.pop(0)
        except Exception:
            pass

    def _undo(self, e=None):
        """Revert the last drawing action (paint stroke, freehand, line, droplet, fill, clear)."""
        if not self._undo_stack:
            self.status.config(text="Nothing to undo.")
            return "break"
        self.active, self.layers = self._undo_stack.pop()
        self.active = max(0, min(self.active, len(self.layers) - 1))
        self._cur_path = None
        self._line_pts = []
        self._draw()
        self._refresh_labels()
        self.status.config(text="Undo.")
        return "break"

    def _on_drag(self, e):
        if self.fmode.get() == "grid":
            self._paint(e, not getattr(self, "_grid_erase", False)); return
        if self._cur_path is None: return
        x, y = self._px_to_mm(e.x, e.y)
        lx, ly = self._cur_path[-1]
        if (x - lx) ** 2 + (y - ly) ** 2 >= 0.49:     # >=0.7mm step, keeps it smooth but light
            self._cur_path.append((x, y)); self._draw()

    def _on_release(self, e):
        if self.fmode.get() == "grid": return
        if self.ff_tool.get() == "line": return        # line tool commits on double-click, not release
        if self._cur_path and len(self._cur_path) >= 2:
            self.layers[self.active].setdefault("paths", []).append(self._cur_path)
        self._cur_path = None
        self._draw(); self._refresh_labels()

    def _on_hover(self, e):
        """Rubber-band preview: while placing a straight line, draw a dashed guide from the last
        corner to the cursor so you can see the segment before clicking. Cheap (one canvas item)."""
        if self.fmode.get() != "freeform" or self.ff_tool.get() != "line" or not self._line_pts:
            return
        try:
            ox, oy, s = self._ff
        except Exception:
            return
        cv = self.canvas
        cv.delete("line_preview")
        lx, ly = self._line_pts[-1]
        cv.create_line(ox + lx * s, oy + ly * s, e.x, e.y, fill=PATH, dash=(4, 3), tags="line_preview")

    def _line_finish(self, e=None):
        """Finish the straight-line polyline (double-click). Removes the duplicate vertices the
        double-click adds, and auto-closes if the last corner is near the first so Fill sees a
        closed shape."""
        pts = list(self._line_pts); self._line_pts = []
        while len(pts) >= 2 and (pts[-1][0] - pts[-2][0]) ** 2 + (pts[-1][1] - pts[-2][1]) ** 2 < 0.25:
            pts.pop()                                   # drop the double-click's duplicate point(s)
        if len(pts) >= 3 and (pts[-1][0] - pts[0][0]) ** 2 + (pts[-1][1] - pts[0][1]) ** 2 < 4.0:
            pts.append(pts[0])                          # snap closed
        if len(pts) >= 2:
            self._push_undo()
            self.layers[self.active].setdefault("paths", []).append(pts)
        self._draw(); self._refresh_labels()

    def _fill_shape(self):
        """Freeform: fill the drawn outline(s) in the active layer with solid zigzag infill at
        the current pitch — so you get a solid without hand-scribbling. Keeps the outline."""
        if self.fmode.get() != "freeform":
            self.status.config(text="Fill is for Freeform mode (in Grid, just paint cells)."); return
        L = self.layers[self.active]
        loops = [p for p in L.get("paths", []) if len(p) >= 3]
        if not loops:
            self.status.config(text="Draw a closed outline first, then press Fill."); return
        try:
            segs = _poly_infill(loops, max(0.5, float(self.pitch)), horizontal=True)
        except Exception as ex:
            self.status.config(text=f"Fill failed: {ex}"); return
        snake = _fill_to_snake(segs)
        self._push_undo()
        L["paths"] = list(loops) + snake          # outline first, then the fill snake(s)
        self._draw(); self._refresh_labels()
        self.status.config(text=f"Filled with {len(snake)} infill path(s) at {float(self.pitch):.1f} mm spacing.")

    def _on_rclick(self, e):
        L = self.layers[self.active]
        if self.fmode.get() == "grid":
            self._paint(e, False); return
        x, y = self._px_to_mm(e.x, e.y)
        # delete nearest droplet point first (within 3mm), else nearest path
        pts = L.get("points", [])
        best_i, best_d = -1, 9.0
        for i, (px, py, n) in enumerate(pts):
            d = (px - x) ** 2 + (py - y) ** 2
            if d < best_d: best_d, best_i = d, i
        if best_i >= 0:
            pts.pop(best_i); self._draw(); self._refresh_labels(); return
        paths = L.get("paths", [])
        bi, bd = -1, 16.0
        for i, path in enumerate(paths):
            for (px, py) in path:
                d = (px - x) ** 2 + (py - y) ** 2
                if d < bd: bd, bi = d, i
        if bi >= 0:
            paths.pop(bi); self._draw(); self._refresh_labels()

    def _on_rdrag(self, e):
        # right-drag erases continuously in grid mode only; freeform deletes on click,
        # never on drag, so you can't accidentally wipe many paths at once.
        if self.fmode.get() == "grid":
            self._paint(e, False)

    def _px_to_mm(self, px, py):
        ox, oy, s = self._ff
        return ((px - ox) / s, (py - oy) / s) if s else (0.0, 0.0)

    def _cell_at(self, ex, ey):
        c = int((ex - self.PAD) // self.cellpx); r = int((ey - self.PAD) // self.cellpx)
        if 0 <= c < self.cols and 0 <= r < self.rows:
            return c, r
        return None
    def _paint(self, e, primary):
        rc = self._cell_at(e.x, e.y)
        if not rc: return
        c, r = rc; L = self.layers[self.active]
        if self.mode.get() == "scaffold":
            if primary: L["scaffold"][r][c] = 1
            else: L["scaffold"][r][c] = 0; L["drops"][r][c] = 0
        else:
            if primary:
                if L["scaffold"][r][c]: L["drops"][r][c] = self.brush.get()
            else: L["drops"][r][c] = 0
        self._draw(); self._refresh_labels()

    # ---------------- drawing ----------------
    def _draw(self):
        if self.fmode.get() == "freeform":
            self._draw_freeform(); return
        cv = self.canvas; cv.delete("all")
        cw, ch = cv.winfo_width(), cv.winfo_height()
        if cw < 10 or self.cols == 0: return
        self.cellpx = max(6, int(min((cw - 2 * self.PAD) / self.cols, (ch - 2 * self.PAD) / self.rows)))
        P, cp = self.PAD, self.cellpx
        x2c = lambda c: P + c * cp + cp / 2; y2r = lambda r: P + r * cp + cp / 2
        reach = self.dead > 0
        if reach:
            cv.create_rectangle(P + (self.cols - self.dead) * cp, P, P + self.cols * cp, P + self.rows * cp,
                                fill="#2a1416", outline="")
        L = self.layers[self.active]
        if self.v_onion.get() and self.active > 0:
            below = self.layers[self.active - 1]["scaffold"]
            for r in range(self.rows):
                for c in range(self.cols):
                    if below[r][c]:
                        cv.create_rectangle(P + c * cp + 1, P + r * cp + 1, P + (c + 1) * cp - 1, P + (r + 1) * cp - 1,
                                            fill="#203a38", outline="")
        for r in range(self.rows):
            for c in range(self.cols):
                x, y = P + c * cp, P + r * cp
                if L["scaffold"][r][c]:
                    dead = reach and c >= self.cols - self.dead
                    cv.create_rectangle(x + 1, y + 1, x + cp - 1, y + cp - 1,
                                        fill="#2f4f4d" if dead else ALG, outline="")
                d = L["drops"][r][c]
                if d:
                    cv.create_text(x + cp / 2, y + cp / 2, text=str(d), fill="#10151b",
                                   font=("TkDefaultFont", max(7, int(cp * 0.4)), "bold"))
        for c in range(self.cols + 1):
            cv.create_line(P + c * cp, P, P + c * cp, P + self.rows * cp, fill=GRID)
        for r in range(self.rows + 1):
            cv.create_line(P, P + r * cp, P + self.cols * cp, P + r * cp, fill=GRID)

        # HOME = bottom-left cell (= machine max X/Y = the pink dot on your bed)
        hx, hy = P + 0 * cp + cp / 2, P + (self.rows - 1) * cp + cp / 2
        cv.create_oval(hx - 6, hy - 6, hx + 6, hy + 6, fill="#FF4FD8", outline="#10151b")
        cv.create_text(hx, hy - cp * 0.7, text="HOME", fill="#FF4FD8",
                       font=("TkDefaultFont", max(7, int(cp * 0.32)), "bold"))
        cv.create_text(P + self.cols * cp / 2, P + self.rows * cp + 12,
                       text="FRONT  (home edge)", fill="#888", font=("TkDefaultFont", 9, "bold"))
        cv.create_text(P + self.cols * cp / 2, max(8, P - 8), text="BACK", fill="#888",
                       font=("TkDefaultFont", 9, "bold"))
        self._home_px = (hx, hy)


        if self.v_path_s.get():
            raw = L["scaffold"]
            scaf = m2g.bridge_islands(raw, self.cols, self.rows) if self.o_bridge.get() else raw
            sp = (m2g.raster_path if self.v_raster.get() else m2g.trace_path)(scaf, self.cols, self.rows)
            for i in range(len(sp) - 1):
                a, b = sp[i][-1], sp[i + 1][0]
                cv.create_line(x2c(a[0]), y2r(a[1]), x2c(b[0]), y2r(b[1]), fill=TRAVEL, dash=(2, 3))
            dot = max(3, int(cp * 0.18))
            for st in sp:
                if len(st) == 1:                         # isolated cell -> solid dot (a 1-pt line is invisible)
                    s = st[0]
                    cv.create_oval(x2c(s[0]) - dot, y2r(s[1]) - dot, x2c(s[0]) + dot, y2r(s[1]) + dot,
                                   fill=PATH, outline="")
                else:
                    pts = [coord for cell in st for coord in (x2c(cell[0]), y2r(cell[1]))]
                    cv.create_line(*pts, fill=PATH, width=max(2, int(cp * 0.12)))
            if sp:
                s = sp[0][0]; cv.create_oval(x2c(s[0]) - 4, y2r(s[1]) - 4, x2c(s[0]) + 4, y2r(s[1]) + 4, fill=PATH, outline="")
                hx, hy = self._home_px
                cv.create_line(hx, hy, x2c(s[0]), y2r(s[1]), fill="#FF4FD8", dash=(5, 3), width=2)
            if not self.o_bridge.get():                  # flag only stray single cells (loose dabs)
                rr = int(cp * 0.42)
                for comp in m2g.components(raw, self.cols, self.rows):
                    if len(comp) == 1:                   # a connected piece of 2+ cells is a valid separate part
                        c, r = next(iter(comp))
                        cv.create_oval(x2c(c) - rr, y2r(r) - rr, x2c(c) + rr, y2r(r) + rr,
                                       outline=REACH, width=2)
        if self.v_path_d.get():
            dp = m2g.droplet_path(L["drops"], self.cols, self.rows)
            if len(dp) > 1:
                pts = [coord for (c, r, n) in dp for coord in (x2c(c), y2r(r))]
                cv.create_line(*pts, fill=DROP, dash=(4, 3))
            for (c, r, n) in dp:
                cv.create_oval(x2c(c) - 3, y2r(r) - 3, x2c(c) + 3, y2r(r) + 3, fill=DROP, outline="")
            if dp:
                fc, fr, fn = dp[0]
                hx, hy = self._home_px
                cv.create_line(hx, hy, x2c(fc), y2r(fr), fill="#FF4FD8", dash=(5, 3), width=2)

    def _draw_freeform(self):
        cv = self.canvas; cv.delete("all")
        cw, ch = cv.winfo_width(), cv.winfo_height()
        if cw < 10: return
        P = self.PAD
        s = min((cw - 2 * P) / max(1e-6, self.bedX), (ch - 2 * P) / max(1e-6, self.bedY))
        ox, oy = P, P
        self._ff = (ox, oy, s)
        m2p = lambda x, y: (ox + x * s, oy + y * s)
        # bed outline
        cv.create_rectangle(ox, oy, ox + self.bedX * s, oy + self.bedY * s, outline=GRID)
        # HOME = bottom-left of the bed (= machine max X/Y = the pink dot)
        hpx, hpy = m2p(0, self.bedY)
        cv.create_oval(hpx - 6, hpy - 6, hpx + 6, hpy + 6, fill="#FF4FD8", outline="#10151b")
        cv.create_text(hpx + 4, hpy - 12, text="HOME", fill="#FF4FD8", anchor="w",
                       font=("TkDefaultFont", 9, "bold"))
        cv.create_text(ox + self.bedX * s / 2, oy + self.bedY * s + 12,
                       text="FRONT  (home edge)", fill="#888", font=("TkDefaultFont", 9, "bold"))
        cv.create_text(ox + self.bedX * s / 2, max(8, oy - 8), text="BACK", fill="#888",
                       font=("TkDefaultFont", 9, "bold"))
        self._home_px = (hpx, hpy)
        bw = max(2, int(self.bead * s))                    # brush width = bead
        L = self.layers[self.active]
        # onion: faint paths of the layer below
        if self.v_onion.get() and self.active > 0:
            for path in self.layers[self.active - 1].get("paths", []):
                if len(path) >= 2:
                    pts = [c for (x, y) in path for c in m2p(x, y)]
                    cv.create_line(*pts, fill="#24506a", width=bw, capstyle="round", joinstyle="round")
        # scaffold paths
        paths = L.get("paths", [])
        for path in paths:
            if len(path) >= 2:
                pts = [c for (x, y) in path for c in m2p(x, y)]
                cv.create_line(*pts, fill=ALG, width=bw, capstyle="round", joinstyle="round")
            elif path:
                px, py = m2p(*path[0]); cv.create_oval(px - bw, py - bw, px + bw, py + bw, fill=ALG, outline="")
        # in-progress straight-line polyline (Line tool): show placed corners + segments
        if self._line_pts:
            if len(self._line_pts) >= 2:
                lp = [c for (x, y) in self._line_pts for c in m2p(x, y)]
                cv.create_line(*lp, fill=PATH, width=2)
            for (x, y) in self._line_pts:
                px, py = m2p(x, y); cv.create_oval(px - 3, py - 3, px + 3, py + 3, fill=PATH, outline="")
        # print path: travel (lift) moves between separate scaffold strokes (travel-optimized order)
        if self.v_path_s.get() and len(paths) > 1:
            ordered = m2g.optimize_order(paths, (0.0, 0.0))
            for i in range(len(ordered) - 1):
                if ordered[i] and ordered[i + 1]:
                    ax, ay = m2p(*ordered[i][-1]); bx, by = m2p(*ordered[i + 1][0])
                    cv.create_line(ax, ay, bx, by, fill=TRAVEL, dash=(2, 3))
        if self.v_path_s.get() and paths:
            first = m2g.optimize_order(paths, (0.0, 0.0))[0]
            if first:
                fx, fy = m2p(*first[0]); hpx, hpy = self._home_px
                cv.create_line(hpx, hpy, fx, fy, fill="#FF4FD8", dash=(5, 3), width=2)
        # in-progress stroke
        if self._cur_path and len(self._cur_path) >= 2:
            pts = [c for (x, y) in self._cur_path for c in m2p(x, y)]
            cv.create_line(*pts, fill=PATH, width=bw, capstyle="round", joinstyle="round")
        # droplet points — small fixed dots (droplets are far smaller than the bead)
        points = L.get("points", [])
        if self.v_path_d.get() and len(points) > 1:        # droplet visit order (optimized)
            ordpts = m2g.optimize_points(points, (0.0, 0.0))
            pts = [c for (x, y, n) in ordpts for c in m2p(x, y)]
            cv.create_line(*pts, fill=DROP, dash=(4, 3))
        if self.v_path_d.get() and points:
            fp = m2g.optimize_points(points, (0.0, 0.0))[0]
            fx, fy = m2p(fp[0], fp[1]); hpx, hpy = self._home_px
            cv.create_line(hpx, hpy, fx, fy, fill="#FF4FD8", dash=(5, 3), width=2)
        for (x, y, n) in points:
            px, py = m2p(x, y)
            cv.create_oval(px - 2.5, py - 2.5, px + 2.5, py + 2.5, fill=DROP, outline="")
            if n > 1:
                cv.create_text(px + 6, py - 6, text=str(n), fill=DROP, font=("Segoe UI", 7))

    # ---------------- STL ----------------
    def _load_stl(self):
        path = filedialog.askopenfilename(filetypes=[("STL", "*.stl"), ("All", "*.*")])
        if not path: return
        try:
            tris = parse_stl(path)
            if not tris: raise ValueError("no triangles parsed")
            self._stl_tris = tris
            self._slice_current()
        except Exception as ex:
            messagebox.showerror("Studio", f"STL failed: {ex}")

    def _reslice(self):
        if self._stl_tris:
            self._slice_current()

    def _slice_current(self):
        tris = self._stl_tris
        if not tris: return
        self._apply()  # sync grid params first
        up = self.up_axis.get()
        if up == "auto":                      # height = shortest dimension (typical for scaffolds)
            xs = [p[0] for t in tris for p in t]
            ys = [p[1] for t in tris for p in t]
            zs = [p[2] for t in tris for p in t]
            spans = [max(xs) - min(xs), max(ys) - min(ys), max(zs) - min(zs)]
            up = ["x", "y", "z"][spans.index(min(spans))]
        if self.fmode.get() == "freeform":
            layers, info = slice_stl_contours(tris, self.bedX, self.bedY, self.margin,
                                              self.bead, fit=True, up=up,
                                              infill=self.v_infill.get(), infill_spacing=self.pitch)
            # give each contour layer the grid keys too, so grid ops never crash
            for ly in layers:
                ly.setdefault("scaffold", [[0] * self.cols for _ in range(self.rows)])
                ly.setdefault("drops", [[0] * self.cols for _ in range(self.rows)])
            self.layers, self.active = layers or [self._blank()], 0
            self.stl_lbl.config(text=f"{info} · up={up}")
        else:
            cols, rows, layers, info = slice_stl(tris, self.bedX, self.bedY, self.pitch,
                                                 self.margin, self.bead, fit=True, up=up)
            self.cols, self.rows, self.layers, self.active = cols, rows, layers, 0
            total = sum(sum(sum(row) for row in ly["scaffold"]) for ly in layers)
            self.stl_lbl.config(text=f"{info} · up={up} · {total} cells")
            if total == 0:
                messagebox.showwarning("Studio", "Sliced but no cells filled — try a smaller pitch "
                                       "or a different up axis.")
        self._draw(); self._refresh_labels()
        self._printability_warn("the imported mesh")

    def _min_clearance(self, paths, bead):
        """Smallest gap between non-adjacent contour points, via a bead-sized spatial
        hash. Only resolves gaps < bead (returns the smallest such, else None)."""
        from collections import defaultdict
        cs = max(bead, 0.1)
        pts = []
        for li_, poly in enumerate(paths):
            if len(poly) < 3:                 # skip 2-pt infill segments; check perimeters/curves
                continue
            n = len(poly)
            for vi, p in enumerate(poly):
                pts.append((p[0], p[1], li_, vi, n))
        grid = defaultdict(list)
        for idx, (x, y, _l, _v, _n) in enumerate(pts):
            grid[(int(x // cs), int(y // cs))].append(idx)
        mind = None
        for idx, (x, y, l, v, n) in enumerate(pts):
            gx, gy = int(x // cs), int(y // cs)
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    for j in grid.get((gx + dx, gy + dy), ()):
                        if j <= idx:
                            continue
                        xj, yj, lj, vj, nj = pts[j]
                        if l == lj and (abs(v - vj) <= 1 or abs(v - vj) >= n - 1):
                            continue          # adjacent vertices on the same loop
                        d = ((x - xj) ** 2 + (y - yj) ** 2) ** 0.5
                        if mind is None or d < mind:
                            mind = d
        return mind

    def _printability_warn(self, what="this design"):
        """Warn if features are finer than the gel/bead can resolve. Returns True if warned."""
        bead = self.bead; issues = []
        if self.fmode.get() == "grid" and self.pitch < bead:
            issues.append(f"\u2022 strand spacing (pitch {self.pitch:.2f} mm) is below the bead width "
                          f"{bead:.2f} mm \u2014 neighbouring strands will merge")
        if self.fmode.get() == "freeform":
            worst, worst_z = None, None
            for ly in self.layers:
                md = self._min_clearance(ly.get("paths", []), bead)
                if md is not None and (worst is None or md < worst):
                    worst, worst_z = md, ly.get("z")
            if worst is not None and worst < bead:
                issues.append(f"\u2022 narrowest feature/gap \u2248 {worst:.2f} mm (layer z={worst_z}) is below "
                              f"the bead width {bead:.2f} mm \u2014 it won't resolve and may blob")
        if issues:
            messagebox.showwarning("Printability",
                f"Some features of {what} are finer than the current gel parameters can print:\n\n"
                + "\n".join(issues)
                + f"\n\nThe bead width ({bead:.2f} mm) is the resolution limit. To fix: scale the model "
                "up, increase spacing, or use a finer nozzle / smaller bead.")
            return True
        return False

    def _clear_layer(self):
        self._push_undo()
        L = self.layers[self.active]
        L["scaffold"] = [[0] * self.cols for _ in range(self.rows)]
        L["drops"] = [[0] * self.cols for _ in range(self.rows)]
        L["paths"] = []; L["points"] = []
        self._draw(); self._refresh_labels()

    def _reset_all(self):
        if not messagebox.askokcancel("Reset", "Clear every layer and start with one empty layer?"):
            return
        self._stl_tris = None
        self.layers = [self._blank()]
        self.active = 0
        self.stl_lbl.config(text="no mesh")
        self._draw(); self._refresh_labels()

    # ---------------- gel profiles ----------------
    def _profiles_dir(self):
        d = os.path.join(os.path.dirname(os.path.abspath(__file__)), "profiles")
        os.makedirs(d, exist_ok=True)
        return d

    def _refresh_profiles(self):
        try:
            names = sorted(f[:-5] for f in os.listdir(self._profiles_dir()) if f.endswith(".json"))
        except Exception:
            names = []
        self.profile_combo["values"] = names
        if names and not self.profile_var.get():
            self.profile_var.set(names[0])

    def _profile_dict(self):
        d = {"flow": self.o_flow.get(), "hx": self.o_hx.get(), "hy": self.o_hy.get(),
             "zspoof": self.o_zspoof.get(), "hop": self.o_hop.get(), "depth": self.o_depth.get(),
             "dip": self.o_dip.get(), "dropmm": self.o_dropmm.get(),
             "pulse": self.o_pulse.get(), "gap": self.o_gap.get(),
             "bead": self.e["bead"].get(), "pitch": self.e["pitch"].get(),
             "margin": self.e["margin"].get(), "dead": self.e["dead"].get(),
             "rotate": self.o_rotate.get(), "flipx": self.o_flipx.get(), "flipy": self.o_flipy.get(),
             "home": self.o_home.get()}
        if self.app_get:
            try: d["app"] = self.app_get()
            except Exception: pass
        return d

    def _profile_apply(self, d):
        for k, var in [("flow", self.o_flow), ("hx", self.o_hx), ("hy", self.o_hy),
                       ("zspoof", self.o_zspoof), ("hop", self.o_hop), ("depth", self.o_depth),
                       ("dip", self.o_dip), ("dropmm", self.o_dropmm)]:
            if k in d:
                try: var.set(float(d[k]))
                except (ValueError, TypeError): pass
        for k, var in [("pulse", self.o_pulse), ("gap", self.o_gap)]:
            if k in d:
                try: var.set(int(float(d[k])))
                except (ValueError, TypeError): pass
        for k in ("bead", "pitch", "margin", "dead"):
            if k in d: self.e[k].set(str(d[k]))
        if "rotate" in d:
            try: self.o_rotate.set(int(d["rotate"]))
            except (ValueError, TypeError): pass
        if "flipx" in d:
            self.o_flipx.set(bool(d["flipx"]))
        if "flipy" in d:
            self.o_flipy.set(bool(d["flipy"]))
        if d.get("home") in ("TL", "TR", "BL", "BR"):
            self.o_home.set(d["home"])
        if self.app_set and "app" in d:
            try: self.app_set(d["app"])
            except Exception: pass
        self._apply(); self._draw(); self._refresh_labels()

    def _profile_save(self):
        name = simpledialog.askstring("Save profile", "Profile name (e.g. Alginate 3% / CaCl2 3%):",
                                      parent=self)
        if not name: return
        safe = "".join(c for c in name if c.isalnum() or c in " -_").strip() or "profile"
        path = os.path.join(self._profiles_dir(), safe + ".json")
        try:
            json.dump({"name": name, **self._profile_dict()}, open(path, "w"), indent=2)
        except Exception as ex:
            messagebox.showerror("Profile", f"Save failed: {ex}"); return
        self._refresh_profiles(); self.profile_var.set(safe)
        self.status.config(text=f"saved profile: {name}")

    def _profile_load(self):
        name = self.profile_var.get()
        if not name:
            messagebox.showinfo("Profile", "No profile selected."); return
        path = os.path.join(self._profiles_dir(), name + ".json")
        try:
            self._profile_apply(json.load(open(path)))
        except Exception as ex:
            messagebox.showerror("Profile", f"Load failed: {ex}"); return
        self.status.config(text=f"loaded profile: {name}")

    def _profile_export(self):
        path = filedialog.asksaveasfilename(defaultextension=".json",
                                            filetypes=[("Gel profile", "*.json")],
                                            initialfile=(self.profile_var.get() or "gel_profile") + ".json")
        if not path: return
        try:
            json.dump({"name": self.profile_var.get() or "gel_profile", **self._profile_dict()},
                      open(path, "w"), indent=2)
        except Exception as ex:
            messagebox.showerror("Profile", f"Export failed: {ex}"); return
        self.status.config(text=f"exported \u2192 {os.path.basename(path)}")

    def _profile_import(self):
        path = filedialog.askopenfilename(filetypes=[("Gel profile", "*.json"), ("All", "*.*")])
        if not path: return
        try:
            d = json.load(open(path))
            self._profile_apply(d)
            base = d.get("name") or os.path.splitext(os.path.basename(path))[0]
            safe = "".join(c for c in base if c.isalnum() or c in " -_").strip() or "imported"
            json.dump(d, open(os.path.join(self._profiles_dir(), safe + ".json"), "w"), indent=2)
            self._refresh_profiles(); self.profile_var.set(safe)
        except Exception as ex:
            messagebox.showerror("Profile", f"Import failed: {ex}"); return
        self.status.config(text=f"imported: {os.path.basename(path)}")

    # ---------------- design geometry (scaffold + droplets) ----------------
    def _design_dict(self):
        return {"mode": self.fmode.get(), "bed": {"x": self.bedX, "y": self.bedY},
                "pitch": self.pitch, "margin": self.margin, "bead": self.bead, "dead": self.dead,
                "cols": self.cols, "rows": self.rows,
                "layers": [{"z": ly.get("z", 0.0),
                            "scaffold": ly.get("scaffold", []), "drops": ly.get("drops", []),
                            "paths": ly.get("paths", []), "points": ly.get("points", [])}
                           for ly in self.layers]}

    def _design_apply(self, d):
        if "layers" not in d:
            raise ValueError("no geometry in file")
        self.fmode.set(d.get("mode", "grid"))
        bed = d.get("bed", {})
        for key, val in [("bedX", bed.get("x", self.bedX)), ("bedY", bed.get("y", self.bedY)),
                         ("pitch", d.get("pitch", self.pitch)), ("margin", d.get("margin", self.margin)),
                         ("bead", d.get("bead", self.bead)), ("dead", d.get("dead", self.dead))]:
            self.e[key].set(str(val))
        self.bedX = float(bed.get("x", self.bedX)); self.bedY = float(bed.get("y", self.bedY))
        self.pitch = float(d.get("pitch", self.pitch)); self.margin = float(d.get("margin", self.margin))
        self.bead = float(d.get("bead", self.bead)); self.dead = int(d.get("dead", self.dead))
        self.cols = int(d.get("cols", self.cols)); self.rows = int(d.get("rows", self.rows))
        layers = []
        for ly in d["layers"]:
            layers.append({"z": ly.get("z", 0.0),
                           "scaffold": ly.get("scaffold") or [[0] * self.cols for _ in range(self.rows)],
                           "drops": ly.get("drops") or [[0] * self.cols for _ in range(self.rows)],
                           "paths": ly.get("paths", []), "points": ly.get("points", [])})
        self.layers = layers or [self._blank()]
        self.active = 0; self._stl_tris = None

    # ---------------- jobs (design + gel profile in one file) ----------------
    def _job_export(self):
        path = filedialog.asksaveasfilename(defaultextension=".json",
                                            filetypes=[("HBL job", "*.json")], initialfile="job.json")
        if not path: return
        job = {"format": "hbl-job/v1", "design": self._design_dict(),
               "profile": {"name": self.profile_var.get() or "job", **self._profile_dict()}}
        try:
            json.dump(job, open(path, "w"))
        except Exception as ex:
            messagebox.showerror("Job", f"Export failed: {ex}"); return
        self.status.config(text=f"exported job \u2192 {os.path.basename(path)} ({len(self.layers)} layers)")

    def _job_import(self):
        path = filedialog.askopenfilename(filetypes=[("HBL job", "*.json"), ("All", "*.*")])
        if not path: return
        try:
            job = json.load(open(path))
            if "design" in job:                       # full job: geometry + profile
                self._design_apply(job["design"])
                if "profile" in job: self._profile_apply(job["profile"])
            elif "layers" in job:                     # bare design file
                self._design_apply(job)
            else:
                raise ValueError("not a job/design file")
        except Exception as ex:
            messagebox.showerror("Job", f"Import failed: {ex}"); return
        self._draw(); self._refresh_labels()
        self.status.config(text=f"imported job: {os.path.basename(path)} ({len(self.layers)} layers, {self.fmode.get()})")
        self._printability_warn("the imported job")

    # ---------------- G-code export ----------------
    def _gcode_export(self):
        path = filedialog.asksaveasfilename(defaultextension=".gcode",
                                            filetypes=[("G-code", "*.gcode"), ("All", "*.*")],
                                            initialfile="print.gcode")
        if not path: return
        try:
            payload, msg = m2g.emit_layers(self._matrix(), self._cfg())
            lines = list(payload["header"])
            for ly in payload["layers"]:
                lines.append(f"; ===== LAYER {ly['i']} (z={ly['z']}) \u2014 present + spray + continue here =====")
                lines += ly["lines"]
            open(path, "w").write("\n".join(lines) + "\n")
        except Exception as ex:
            messagebox.showerror("G-code", f"Export failed: {ex}"); return
        self.status.config(text=f"exported G-code \u2192 {os.path.basename(path)} ({msg})")

    # ---------------- generate ----------------
    def _cfg(self):
        nr = self.v_noretract.get()               # No retractions: kill prime + retract + corner retract
        return m2g.Cfg(raster=self.v_raster.get(), FLOW=self.o_flow.get(),
                       CORNER_RETRACT=0.0 if nr else self.o_corner.get(),
                       PRIME=0.0 if nr else self.o_prime.get(),
                       RETRACT=0.0 if nr else self.o_retract.get(),
                       DROPLET_DX=self.o_hx.get(), DROPLET_DY=self.o_hy.get(),
                       A_HOP=self.o_hop.get(), A_OFFSET=self.o_depth.get(),
                       PRINT_LIFT=self.o_printlift.get(),
                       DIP=self.o_dip.get(), DROP_VOL=self.o_dropmm.get(),
                       DROPLET_RISE=self.o_droprise.get(), DROP_HOLD=self.o_drophold.get(),
                       DROP_GAP=self.o_dropgap.get(),
                       BRIDGE=self.o_bridge.get())

    def _matrix(self):
        base = {"format": "scaffold-grid/v2", "mode": self.fmode.get(),
                "bed": {"x": self.bedX, "y": self.bedY},
                "pitch": self.pitch, "margin": self.margin, "bead_height": self.bead,
                "grid": {"cols": self.cols, "rows": self.rows},
                "droplet_head": {"offset_x": self.o_hx.get(), "dead_cols": self.dead}}
        if self.fmode.get() == "freeform":
            base["layers"] = [{"z": ly["z"], "paths": ly.get("paths", []),
                               "points": ly.get("points", [])} for ly in self.layers]
        else:
            base["layers"] = [{"z": ly["z"], "scaffold": ly["scaffold"], "drops": ly["drops"]}
                              for ly in self.layers]
        return base

    def _generate(self):
        cfg = self._cfg()
        M = self._matrix()
        if self.fmode.get() == "grid" and not self.o_bridge.get():
            bad = [ly["z"] for ly in self.layers
                   if any(len(c) == 1 for c in m2g.components(ly["scaffold"], self.cols, self.rows))]
            if bad:
                if not messagebox.askokcancel("Stray cells",
                        f"{len(bad)} layer(s) have stray single cells with no neighbours \u2014 each prints "
                        "as a loose dab.\n\n(Separate multi-cell pieces are fine and not flagged.)\n\n"
                        "Enable 'Bridge islands' to strut everything together, or remove the stray "
                        "cells.\n\nGenerate anyway?"):
                    return
        payload, msg = m2g.emit_layers(M, cfg)
        if self.on_generate:
            self.on_generate(payload, msg)            # ({'header','layers'}, message)
        else:
            base = os.path.dirname(os.path.abspath(__file__))
            for ly in payload["layers"]:
                open(os.path.join(base, f"studio_layer{ly['i']}.gcode"), "w").write(
                    "\n".join(payload["header"] + ly["lines"]) + "\n")
            msg += f"\nwrote {len(payload['layers'])} layer files"
        self.status.config(text=msg)


# ---------------- standalone test harness ----------------

# ===================== host app (printess_patched.py) ====================
import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import serial
import serial.tools.list_ports
import time
import csv
import sv_ttk 
import webbrowser
import threading
import math
import os
import re
from PIL import Image, ImageTk
import sys, faulthandler, datetime
# voron imports
import urllib.request
import urllib.error
import urllib.parse
import mimetypes
import uuid

# --- crash capture -------------------------------------------------------------
# Writes a full traceback to bioprinter_crash.log on ANY failure, including hard
# native crashes (tkinter-off-thread, serial, etc.) that don't raise a catchable
# Python exception. If the app dies, send me that file.
try:
    _crashlog_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "bioprinter_crash.log")
except Exception:
    _crashlog_path = "bioprinter_crash.log"
try:
    _crashf = open(_crashlog_path, "a", buffering=1)
    _crashf.write(f"\n===== session {datetime.datetime.now():%Y-%m-%d %H:%M:%S} =====\n")
    faulthandler.enable(file=_crashf, all_threads=True)
except Exception:
    _crashf = None
    faulthandler.enable(all_threads=True)

def _log_crash(kind, et, ev, tb):
    import traceback
    for stream in (sys.stderr, _crashf):
        if stream is None: continue
        try:
            stream.write(f"\n--- {kind} {datetime.datetime.now():%H:%M:%S} ---\n")
            traceback.print_exception(et, ev, tb, file=stream); stream.flush()
        except Exception: pass

sys.excepthook = lambda t, e, tb: _log_crash("uncaught (main thread)", t, e, tb)
try:
    threading.excepthook = lambda a: _log_crash("uncaught (worker thread)", a.exc_type, a.exc_value, a.exc_traceback)
except Exception: pass
# -------------------------------------------------------------------------------

def _status(msg):
    """Thread-safe status update. tkinter must only be touched from the main thread,
    so background (print) threads marshal their status text here via root.after."""
    try: root.after(0, lambda: scaffold_status_var.set(msg))
    except Exception: pass

# Config snapshot: tkinter vars can only be safely read on the main thread, so we copy
# every value the print worker needs into plain Python here, once, at print start.
_cfg = {"slow_mode": False, "slow_cap": 300, "droplet_flow": 1.0, "scaffold_flow": 1.0, "dry_run": False,
        "extrude_f": 0.0, "speed": "300", "auto_park": False, "suction": False,
        "load_mm": 30.0, "pull_f": 150, "manual_offset": False, "present": False, "xlink_dwell": 0.0,
        "continuous": False, "solenoid": False, "solenoid_us": 2000,
        "pneumatic": False, "valve_time_ms": 2.0}

def _snapshot_print_cfg():
    def _f(fn, d):
        try: return fn()
        except Exception: return d
    _cfg["slow_mode"]     = _f(lambda: bool(slow_mode_var.get()), False)
    _cfg["slow_cap"]      = _f(lambda: int(float(slow_cap_var.get())), 300)
    _cfg["droplet_flow"]  = _f(lambda: float(droplet_flow_var.get()), 1.0)
    _cfg["scaffold_flow"] = _f(lambda: float(scaffold_flow_var.get()), 1.0)
    _cfg["dry_run"]       = _f(lambda: bool(dry_run_var.get()), False)
    _cfg["extrude_f"]     = _f(lambda: float(extrude_speed_entry.get()), 0.0)
    _cfg["speed"]         = _f(lambda: idle_speed_entry.get(), "300")
    _cfg["auto_park"]     = _f(lambda: bool(auto_park_var.get()), False)
    _cfg["suction"]       = _f(lambda: bool(suction_load_var.get()), False)
    _cfg["load_mm"]       = _f(lambda: float(load_volume_entry.get()), 30.0)
    _cfg["pull_f"]        = _f(lambda: int(float(pullup_speed_entry.get())), 150)
    _cfg["manual_offset"] = _f(lambda: bool(manual_offset_var.get()), False)
    _cfg["present"]       = _f(lambda: bool(present_var.get()), False)
    _cfg["raise_heads"]   = _f(lambda: bool(raise_heads_var.get()), False)
    _cfg["xlink_dwell"]   = _f(lambda: float(xlink_dwell_var.get()), 0.0)
    _cfg["continuous"]    = _f(lambda: bool(continuous_var.get()), False)
    _cfg["solenoid"]      = _f(lambda: bool(solenoid_var.get()), False)
    _cfg["solenoid_us"]   = _f(lambda: int(float(solenoid_us_var.get())), 2000)
    _cfg["solenoid_snap_ms"] = _f(lambda: int(float(solenoid_snap_var.get())), 0)
    _cfg["drop_f"]        = _f(lambda: int(float(drop_speed_var.get())), 3000)
    _cfg["drop_accel"]    = _f(lambda: int(float(drop_accel_var.get())), 30000)
    _cfg["pneumatic"]     = _f(lambda: bool(pneumatic_var.get()), False)
    _cfg["valve_time_ms"] = _f(lambda: float(valve_time_ms_var.get()), 2.0)

# Bump this string whenever the code changes so you can confirm at a glance which build is
# actually running. Shown in the window title, the startup status line, and the terminal.
BUILD_VERSION = "build 2026-06-29 R7 — B-dispense + offset-locked + present-last"
print(f"\n=== Hybrid bioprinter {BUILD_VERSION} ===\n")

# Studio tab (replaces the old Matrix tab). Needs studio_tab.py + matrix_to_gcode.py
# sitting in the same folder as this file.
_STUDIO_OK = True      # StudioTab is defined above in this single-file build
_STUDIO_ERR = None

# --- HARDWARE OFFSETS & GLOBALS ---
SOLENOID_OFFSET_X = 40.0
SOLENOID_OFFSET_Y = 50.0

# --- HOME / PARK CONFIGURATION ---
# At HOME: bed at the corner (max X/Y), scaffold (right/A) lift at its lowest,
# droplet (left/Z) lift 26 mm up from its bottom. In that pose the scaffold
# nozzle tip and the droplet needle tip sit at the SAME height. Moving the bed
# HEAD_OFFSET_X in -X then puts the needle over the scaffold nozzle's point.
HOME_X = 71.0            # max X (corner)
HOME_Y = 71.0            # max Y (corner)
# The Y motor was physically flipped, which reverses Y for BOTH jogging and prints. Rather than
# reflash firmware (INVERT_Y_DIR), reflect Y at the single send point below. Set _Y_FLIPPED=False
# if you ever undo the motor flip or fix it in firmware (don't do both -> double flip).
_Y_FLIPPED = False
_y_abs_mode = True       # tracks G90(absolute)/G91(relative) so the reflection uses the right form
HOME_SCAFFOLD_A = 0.0    # scaffold (right) lift at its lowest
HOME_DROPLET_Z  = 26.0   # droplet (left) lift, 26 mm up from its bottom
HEAD_OFFSET_X   = 31.0   # droplet needle sits 31 mm in -X of the scaffold nozzle

# --- AUTO-CONVERTER SETTINGS (normal G-code -> this printer) ---
CONV_XMAX, CONV_YMAX = 71.0, 71.0   # parked corner = max X/Y on the bed
CONV_MARGIN     = 5.0               # keep the print this far off the maxed corner
CONV_TARGET_MAX = 45.0              # largest drawing dimension after scaling (<= reach)
CONV_FLOW       = 0.05              # C squeeze mm per mm of travel  (tune)
CONV_PRIME      = 2.0               # squeeze advance at start so stroke 1 isn't dry
CONV_HOP        = 1.5               # A lift height for pen-up between strokes
CONV_USE_HOPS   = True              # False = flat & fast, but joins strokes into a scribble
CONV_PRINT_F, CONV_TRAVEL_F, CONV_ZF, CONV_PRIME_F = 300, 1500, 600, 100

motor_serial = None
valve_serial = None
is_printing = False
_abort_evt = threading.Event()   # set by E-stop to break any serial 'ok' wait
_rx_lines = []                   # board reply lines from the most recent send (for the console)
_serial_busy = False             # (legacy flag, unused) — serialization is via _serial_lock
_serial_lock = threading.Lock()  # only one command touches the port at a time; others are dropped
_serial_dead = False             # set when a write fails (board killed/unplugged); cleared on reconnect
_valve_ok = False                # True once the valve Arduino answers on connect
_connecting = False              # True while connect is (re)opening the port — poller stands down
drop_matrix = [] 
loaded_gcode = [] # legacy single store (kept for compatibility)
loaded_scaffold_gcode = []   # scaffold stream (from Studio or disk)
loaded_droplet_gcode  = []   # droplet stream (from Studio or disk)

# --- COM PORT DETECTION ---
def get_active_ports():
    ports = serial.tools.list_ports.comports()
    port_list = [port.device for port in ports]
    return port_list if port_list else ["No Ports Found"]

def refresh_com_ports():
    new_ports = get_active_ports()
    motor_com_combo['values'] = new_ports
    valve_com_combo['values'] = new_ports
    if new_ports and new_ports[0] != "No Ports Found":
        motor_com_combo.current(0)
        if len(new_ports) > 1:
            valve_com_combo.current(1)
        else:
            valve_com_combo.current(0)
    else:
        motor_com_combo.set("No Ports Found")
        valve_com_combo.set("No Ports Found")

def connect_all_hardware():
    """Read the port names on the main thread, then do ALL serial work on a background
    thread so the UI can never freeze during connect, no matter how the ports behave."""
    m_port = motor_com_combo.get()
    v_port = valve_com_combo.get()
    global _connecting
    _connecting = True
    try: connect_btn.config(state="disabled", text="CONNECTING\u2026")
    except Exception: pass
    sys_status_var.set("Hardware: connecting\u2026")
    threading.Thread(target=_do_connect, args=(m_port, v_port), daemon=True).start()

def _do_connect(m_port, v_port):
    """Runs on a worker thread. No tkinter here — results are marshaled back via root.after.
    Everything is wrapped so the button ALWAYS returns from 'CONNECTING…', and the port has a
    write timeout so a board that isn't reading can't block the connect forever."""
    global motor_serial, valve_serial, _serial_dead, _valve_ok
    _valve_ok = False
    try:
        try:
            if motor_serial and motor_serial.is_open: motor_serial.close()
            motor_serial = serial.Serial(m_port, 115200, timeout=1, write_timeout=1)
            _serial_dead = False           # fresh connection -> clear the dead flag
            time.sleep(2)  # Give STM32 time to reboot
            try:
                motor_serial.reset_input_buffer()
                motor_serial.reset_output_buffer()
                motor_serial.write(b"\nM410\n")   # flush any motion left queued from a prior/crashed session
                motor_serial.flush()
                time.sleep(0.2); motor_serial.reset_input_buffer()
            except Exception as e:
                print(f"connect flush warning: {e}")
        except Exception as e:
            motor_serial = None
            print(f"Warning: Could not connect Motor to {m_port}: {e}")
        if not v_port or v_port in ("", "No Ports Found"):
            valve_serial = None
            print("Valve: no port selected (solenoid will not fire)")
        else:
            try:
                if valve_serial and valve_serial.is_open: valve_serial.close()
                valve_serial = serial.Serial(v_port, 115200, timeout=1, write_timeout=2)
                time.sleep(2.5)                       # the Arduino RESETS when its port opens; let it boot
                for ping in ("C", "F1"):              # ping: works with either firmware
                    try:
                        valve_serial.reset_input_buffer()
                        valve_serial.write((ping + "\n").encode()); valve_serial.flush()
                    except Exception:
                        break
                    t0 = time.time()
                    while time.time() - t0 < 2.0:
                        if valve_serial.in_waiting > 0:
                            resp = valve_serial.readline().decode("utf-8", "ignore").strip()
                            if resp: print(f"VALVE RX: {resp}")
                            if "ok" in resp.lower() or "ready" in resp.lower():
                                _valve_ok = True; break
                        time.sleep(0.01)
                    if _valve_ok: break
                if not _valve_ok:
                    print(f"Valve on {v_port} opened but did NOT answer -- wrong port, or the "
                          "Nano isn't running the valve firmware.")
            except Exception as e:
                valve_serial = None
                print(f"Warning: Could not connect Valve to {v_port}: {e}")
    finally:
        root.after(0, _connect_done)          # ALWAYS re-enable the button, even on error/hang-recovery

def _connect_done():
    """Back on the main thread: update the UI with the connect result."""
    global _connecting
    _connecting = False
    try: connect_btn.config(state="normal", text="CONNECT ALL HARDWARE")
    except Exception: pass
    if motor_serial or valve_serial:
        v = "valve OK" if _valve_ok else ("valve NO REPLY" if valve_serial else "no valve")
        sys_status_var.set(f"Hardware: ONLINE \u00b7 {v}")
        sys_status_label.config(foreground="#A8E6CF" if (not valve_serial or _valve_ok) else "#E6C08F")
        msg = "Motor: " + ("connected" if motor_serial else "NOT connected")
        if _valve_ok:
            msg += "\nValve: connected and answering."
        elif valve_serial:
            msg += ("\nValve: port opened but the Arduino did NOT answer.\n"
                    "Wrong port, or the Nano isn't running the valve firmware.")
        else:
            msg += "\nValve: no port selected \u2014 the solenoid will not fire.\n" \
                   "Set the VALVE COM dropdown to the Nano (CH340) and reconnect."
        messagebox.showinfo("Connected", msg)
    else:
        sys_status_var.set("Hardware: OFFLINE")
        sys_status_label.config(foreground="#E6A8A8")
        messagebox.showwarning("Not connected",
                               "Could not open the motor port. Check it's powered, plugged in, "
                               "and that the COM port is right and not held by another program.")

# --- HARDWARE COMMUNICATION ---
def send_valve_command(cmd_str):
    """Send a command to the Arduino valve controller and wait for its 'ok' (with a safety
    timeout). Matches the firmware author's protocol: 'F<us>' fires a drop for <us> microseconds."""
    global valve_serial
    full_cmd = cmd_str + '\n'
    print(f"VALVE TX: {cmd_str}")
    if valve_serial and valve_serial.is_open:
        try:
            valve_serial.reset_input_buffer()
            valve_serial.write(full_cmd.encode('utf-8'))
            t0 = time.time()
            while time.time() - t0 < 3.0:            # wait for the drop to actually fire
                if valve_serial.in_waiting > 0:
                    resp = valve_serial.readline().decode('utf-8', errors='ignore').strip()
                    if "ok" in resp.lower():
                        break
                time.sleep(0.005)
        except Exception as e:
            print(f"valve write failed: {e}")

def _flip_y(cmd_str):
    """Reflect Y for a physically flipped Y motor, so jog AND print move the correct way.
    G91 (relative): negate Y.  G90 (absolute): mirror around HOME_Y (the parked corner, a fixed
    point so home stays put).  G92 is left untouched so Set-Home still declares the corner."""
    global _y_abs_mode
    if not _Y_FLIPPED:
        return cmd_str
    u = cmd_str.strip().upper()
    if u == "G90": _y_abs_mode = True;  return cmd_str
    if u == "G91": _y_abs_mode = False; return cmd_str
    if u[:3] == "G92" or u[:2] not in ("G0", "G1"):
        return cmd_str
    parts = cmd_str.split()
    if not any(p[:1] in ("Y", "y") for p in parts):
        return cmd_str
    out = []
    for p in parts:
        if p[:1] in ("Y", "y") and _mo_isfloat(p[1:]):
            y = float(p[1:])
            y = (-y) if not _y_abs_mode else (2 * HOME_Y - y)
            out.append(f"Y{y:.4f}")
        else:
            out.append(p)
    return " ".join(out)

def send_motor_command(cmd_str, timeout=20.0):
    """Send a line and wait for the board's 'ok'. Serialized by _serial_lock: if another
    command is already mid-flight (a re-entrant jog via root.update, a worker thread, or the
    poller), this one is dropped rather than interleaving bytes on the port (which corrupts
    the stream and kills the port). Times out fast so a dead board can't freeze the UI."""
    global motor_serial, _serial_dead
    if _serial_dead:                       # connection already flagged dead -> stop instantly, no writes,
        return False                       # so a killed/unplugged board can't cause a write-timeout flood
    cmd_str = _flip_y(cmd_str)             # correct for the physically flipped Y motor (jog + print)
    if not _serial_lock.acquire(blocking=False):
        print(f"MOTOR BUSY, dropped: {cmd_str}")
        return False
    try:
        print(f"MOTOR TX: {cmd_str}")
        _rx_lines.clear()
        if not motor_serial:
            return True
        try:
            motor_serial.reset_input_buffer()
            motor_serial.write((cmd_str + '\n').encode('utf-8'))
        except Exception as e:
            print(f"serial write failed: {e}")
            _serial_dead = True            # board not accepting writes -> stop hammering it
            root.after(0, lambda: sys_status_var.set("Hardware: LOST \u2014 reconnect"))
            return False
        to = 180.0 if cmd_str.upper().startswith(("M400", "G4", "M84")) else timeout
        t0 = time.time()
        _on_main = threading.current_thread() is threading.main_thread()
        while True:
            if _abort_evt.is_set():
                return False
            # root.update() is safe on the MAIN thread only. With the lock held, a jog fired
            # during this update can't interleave — it hits the busy-drop above.
            if _on_main:
                try: root.update()
                except Exception: pass
            try:
                if motor_serial.in_waiting > 0:
                    response = motor_serial.readline().decode('utf-8', 'ignore').strip()
                    if response:
                        _rx_lines.append(response)          # keep the reply text (M114/M122/M503)
                        print(f"MOTOR RX: {response}")
                    if "ok" in response.lower():
                        return True
            except Exception as e:
                print(f"serial read failed: {e}")
                return False
            if time.time() - t0 > to:
                print(f"TIMEOUT waiting for ok on: {cmd_str}")
                return False
            time.sleep(0.002)
    finally:
        _serial_lock.release()

def jog(axis, direction):
    """
    Jogs the translation stages (X, Y) and vertical/extruder stages (Z, A, B, C)
    """
    if is_printing: return 
    _abort_evt.clear()
    distance = step_size_var.get() * direction
    speed = idle_speed_entry.get()
    send_motor_command("G91")
    send_motor_command(f"G0 {axis}{distance} F{speed}")
    send_motor_command("G90")

SAFE_CLEAR = 8.0    # lift before Go-Home travel (smaller = less distance for the A lift to slip over)
# "present for spraying" pose: bed (X/Y only) comes forward from home between layers
PRESENT_DX = 0.0    # forward offset from home in X
PRESENT_DY = 70.0   # forward offset from home in Y when presenting (clamped to bed travel below)
PAUSE_LIFT = 10.0   # needle Z raise while paused, so the droplet syringe can be swapped cleanly
# Mechanical travel of the bed from the home corner (effective bed). The present move is
# clamped to this minus a safety gap so the bed can't keep driving into the end stop.
BED_TRAVEL_X = 52.0
BED_TRAVEL_Y = 60.0
TRAVEL_SAFETY = 2.0   # stop this far short of the mechanical limit

def set_home():
    """Declare the CURRENT pose as home. No motion - jog to the corner with the
    scaffold lift at its lowest and the droplet needle at the matching height,
    then press this. Only X/Y and the two lifts (A/Z) are referenced; the extrude
    axes (B/C) are deliberately left untouched."""
    if is_printing: return
    _abort_evt.clear()
    send_motor_command("G90")
    send_motor_command("M82")
    send_motor_command("M211 S0")
    send_motor_command(f"G92 X{HOME_X} Y{HOME_Y} A{HOME_SCAFFOLD_A} Z{HOME_DROPLET_Z}")
    sys_status_var.set("Home set at current position")

def go_home():
    """Return to the home set by Set Home: lift both heads clear, travel to the
    corner, then lower to the home heights. Only moves X/Y/A/Z - never B/C.
    Requires Set Home to have been done first (no endstops to find it otherwise)."""
    if is_printing: return
    _abort_evt.clear()
    speed = idle_speed_entry.get()
    send_motor_command("G90")
    send_motor_command("M82")
    send_motor_command("M211 S0")
    # lift both heads clear of the block before travelling
    send_motor_command(f"G0 A{HOME_SCAFFOLD_A + SAFE_CLEAR} Z{HOME_DROPLET_Z + SAFE_CLEAR} F{speed}")
    send_motor_command(f"G0 X{HOME_X} Y{HOME_Y} F{speed}")            # bed to the corner
    send_motor_command(f"G0 A{HOME_SCAFFOLD_A} Z{HOME_DROPLET_Z} F{speed}")  # down to home heights
    send_motor_command("M400")

# --- DYNAMIC GRID & EXCEL ---
def generate_dynamic_grid():
    for widget in matrix_frame.winfo_children():
        widget.destroy()
    drop_matrix.clear()
    try:
        r_count = int(rows_entry.get())
        c_count = int(cols_entry.get())
    except ValueError:
        return
    for r in range(r_count):
        row_entries = []
        for c in range(c_count):
            box = ttk.Entry(matrix_frame, width=4, justify="center", font=("Segoe UI", 12))
            box.grid(row=r, column=c, padx=2, pady=2)
            box.insert(0, "0") 
            row_entries.append(box)
        drop_matrix.append(row_entries)

def import_from_excel():
    filepath = filedialog.askopenfilename(filetypes=[("CSV Files", "*.csv")])
    if not filepath: return 
    try:
        with open(filepath, newline='') as f:
            reader = csv.reader(f)
            data = list(reader)
        if not data: return
        rows_entry.delete(0, tk.END); rows_entry.insert(0, str(len(data)))
        cols_entry.delete(0, tk.END); cols_entry.insert(0, str(len(data[0])))
        generate_dynamic_grid()
        for r in range(len(data)):
            for c in range(len(data[r])):
                drop_matrix[r][c].delete(0, tk.END)
                drop_matrix[r][c].insert(0, str(data[r][c]).strip())
    except Exception as e:
        messagebox.showerror("Import Error", f"Could not read file.\n{e}")

def emergency_halt():
    """Actually stop the machine NOW: abort all queued/in-progress motion and
    shut the valve, written straight to the port so it jumps the command queue.
    (Requires EMERGENCY_PARSER enabled in Marlin so M410 acts immediately; if it
    isn't, swap M410 for M112 - a hard stop that needs a reconnect afterward.)"""
    global is_printing, _serial_dead
    is_printing = False
    _abort_evt.set()                          # break any in-flight serial 'ok' wait
    hard = False
    try: hard = hard_stop_var.get()           # defined in the UI; guard in case of early call
    except Exception: hard = False
    # Wait briefly for the print thread's current command to release the port (the abort flag
    # above makes it return fast), THEN write the stop — so we never collide bytes with an
    # in-flight command, which is what corrupts and kills the port.
    got = _serial_lock.acquire(timeout=2.0)
    try:
        if motor_serial and motor_serial.is_open:
            if hard:
                motor_serial.write(b"\nM112\n")   # emergency kill: stops buffered moves NOW (needs reconnect)
                _serial_dead = True               # board is now halted -> stop sending until reconnect
                root.after(0, lambda: sys_status_var.set("Hardware: HALTED (M112) \u2014 reconnect to continue"))
            else:
                motor_serial.write(b"\nM410\n")   # quickstop: only instant if EMERGENCY_PARSER is enabled
            motor_serial.write(b"M106 P0 S0\n")   # valve OFF (FAN0)
            motor_serial.flush()
    except Exception as e:
        print(f"E-STOP write error (motor): {e}")
    finally:
        if got: _serial_lock.release()
    try:
        if valve_serial and valve_serial.is_open:
            valve_serial.write(b"M106 P0 S0\n")   # valve off on a separate controller too
            valve_serial.flush()
    except Exception as e:
        print(f"E-STOP write error (valve): {e}")

def stop_print_job():
    emergency_halt()
    status_var.set("Status: STOPPED")


# --- AUTO-CONVERTER (normal slicer G-code -> this printer's format) ---
def is_printer_format(lines):
    """True if the file is already in this printer's format (don't convert)."""
    for ln in lines:
        if "hybrid bioprinter" in ln.lower():          # our converter signature
            return True
        code = ln.split(";")[0].upper().split()
        if not code:
            continue
        if "G92" in code and any(p == "X71" for p in code):   # inverted corner self-zero
            return True
        if code[0] in ("G0", "G1") and any(p.startswith("C") for p in code[1:]):
            return True
    return False


def convert_gcode_lines(raw_lines):
    """Plain slicer G-code (X/Y/Z/E from a normal origin) -> this printer:
    flipped/scaled into the corner, E->C flow, Z->A lift, primed, stripped."""
    cx = cy = 0.0
    segs = []
    for line in raw_lines:
        line = line.split(";")[0].strip()
        if not line:
            continue
        p = line.split()
        if p[0] in ("G0", "G1"):
            x, y, e, hasxy = cx, cy, 0.0, False
            for t in p[1:]:
                k = t[:1]
                try:
                    if k == "X": x = float(t[1:]); hasxy = True
                    elif k == "Y": y = float(t[1:]); hasxy = True
                    elif k == "E": e = float(t[1:])
                except ValueError:
                    pass
            if hasxy:
                if e > 0:
                    segs.append((cx, cy, x, y))
                cx, cy = x, y
    if not segs:
        return list(raw_lines)

    polys, cur = [], None
    for x0, y0, x1, y1 in segs:
        if cur and abs(cur[-1][0]-x0) < 1e-6 and abs(cur[-1][1]-y0) < 1e-6:
            cur.append((x1, y1))
        else:
            if cur: polys.append(cur)
            cur = [(x0, y0), (x1, y1)]
    if cur: polys.append(cur)

    pts = [pt for pl in polys for pt in pl]
    minx = min(p[0] for p in pts); maxx = max(p[0] for p in pts)
    miny = min(p[1] for p in pts); maxy = max(p[1] for p in pts)
    w, h = maxx-minx, maxy-miny
    s = CONV_TARGET_MAX / max(w, h, 1e-9)

    def MX(px): return CONV_XMAX - CONV_MARGIN - (px-minx)*s
    def MY(py): return CONV_YMAX - CONV_MARGIN - (py-miny)*s

    out = ["; converted for hybrid bioprinter (auto-converted on load)",
           "G21", "G90",
           f"G92 X{CONV_XMAX:.0f} Y{CONV_YMAX:.0f}", "G92 A0", "G92 C0"]
    e = 0.0
    for i, pl in enumerate(polys):
        sx, sy = MX(pl[0][0]), MY(pl[0][1])
        if CONV_USE_HOPS: out.append(f"G0 A{CONV_HOP:.2f} F{CONV_ZF}")
        out.append(f"G0 X{sx:.2f} Y{sy:.2f} F{CONV_TRAVEL_F}")
        if CONV_USE_HOPS: out.append(f"G0 A0 F{CONV_ZF}")
        if i == 0:
            e += CONV_PRIME; out.append(f"G1 C{e:.4f} F{CONV_PRIME_F}")
        px, py = sx, sy
        for raw in pl[1:]:
            qx, qy = MX(raw[0]), MY(raw[1])
            e += CONV_FLOW * math.hypot(qx-px, qy-py)
            out.append(f"G1 X{qx:.2f} Y{qy:.2f} C{e:.4f} F{CONV_PRINT_F}")
            px, py = qx, qy
    if CONV_USE_HOPS: out.append(f"G0 A{CONV_HOP:.2f} F{CONV_ZF}")
    out += [f"G0 X{CONV_XMAX:.0f} Y{CONV_YMAX:.0f} F{CONV_TRAVEL_F}", "M400", "M84", "; done"]
    return out

# --- PRINT-CONTROL LOGIC (per-layer program with a spray pause between layers) ---
print_header = []        # one-time setup lines (G21/G90/M211/G92)
print_layers = []        # [{'i','z','lines':[exec lines]}]
loaded_view_gcode = []   # full preview text
_continue_evt = threading.Event()
_run_evt = threading.Event(); _run_evt.set()   # cleared = paused

def get_app_gel_params():
    """App-side gel params for the Studio profile system (called at runtime)."""
    return {"suction_load": bool(suction_load_var.get()),
            "load_vol": load_volume_entry.get(),
            "pullup_f": pullup_speed_entry.get(),
            "extrude_f": extrude_speed_entry.get(),
            "purge": purge_entry.get()}

def set_app_gel_params(d):
    try:
        if "suction_load" in d: suction_load_var.set(bool(d["suction_load"]))
        for key, ent in (("load_vol", load_volume_entry), ("pullup_f", pullup_speed_entry),
                         ("extrude_f", extrude_speed_entry), ("purge", purge_entry)):
            if key in d:
                ent.delete(0, tk.END); ent.insert(0, str(d[key]))
    except Exception:
        pass

def _strip_exec(text_lines):
    out = []
    for ln in text_lines:
        s = ln.strip()
        if not s or s.startswith(";"):
            continue
        s = s.split(";")[0].strip()
        if s:
            out.append(s)
    return out

def _show_preview():
    gcode_listbox.delete(0, tk.END)
    for ln in loaded_view_gcode:
        gcode_listbox.insert(tk.END, ln)
    scaffold_status_var.set(f"{len(print_layers)} layers ready")

def load_studio_gcode(payload, msg):
    """Receive the per-layer payload {'header','layers'} from the Studio tab."""
    global print_header, print_layers, loaded_view_gcode
    print_header = _strip_exec(payload.get("header", []))
    print_layers = []
    view = list(payload.get("header", []))
    for ly in payload.get("layers", []):
        print_layers.append({"i": ly["i"], "z": ly["z"], "lines": _strip_exec(ly["lines"])})
        view.append(f"; ===== LAYER {ly['i']}  z={ly['z']:.2f} =====")
        view += ly["lines"]
    loaded_view_gcode = [l for l in view if l.strip()]
    _show_preview()
    notebook.select(tab_scaffold)
    messagebox.showinfo("Studio \u2192 Print", msg)

def load_gcode_file():
    """Load a flat .gcode file from disk as a single layer (manual/testing path)."""
    global print_header, print_layers, loaded_view_gcode
    filepath = filedialog.askopenfilename(filetypes=[("G-Code Files", "*.gcode *.txt")])
    if not filepath: return
    try:
        with open(filepath, 'r') as f:
            raw = [ln.rstrip('\n') for ln in f]
    except Exception as e:
        messagebox.showerror("Error", f"Failed to load file:\n{e}"); return
    src = raw if is_printer_format(raw) else convert_gcode_lines(raw)
    print_header = []
    print_layers = [{"i": 0, "z": 0.0, "lines": _strip_exec(src)}]
    loaded_view_gcode = [l for l in src if l.strip()]
    _show_preview()

def _scale_C(line, factor, st):
    """Delta-scale the C (extrude) value so a live slider change is smooth. G92 C resets."""
    toks = line.split()
    if toks and toks[0].upper() == "G92":
        for p in toks:
            if p[:1] == "C":
                try: v = float(p[1:]); st['oc'] = v; st['sc'] = v
                except ValueError: pass
        return line
    out = []
    for p in toks:
        if p[:1] == "C":
            try:
                c = float(p[1:])
            except ValueError:
                out.append(p); continue
            d = c - st['oc']; delta = factor * d
            if delta > 12.0: delta = 12.0        # anti-dump: no single command shoves the plunger
            elif delta < -12.0: delta = -12.0     #   more than 12mm (a real stroke is < ~3mm)
            sc = st['sc'] + delta
            st['oc'] = c; st['sc'] = sc
            out.append(f"C{sc:.4f}")
        else:
            out.append(p)
    return " ".join(out)

def _cap_feedrate(line):
    """Bad-wire workaround: when Slow Mode is on, clamp any Fxxxx in a move down to the
    cap, so fast travels run at the slow, reliable jog speed instead of skipping steps
    on a marginal connection. Small/slow jog moves already pass through unchanged."""
    try:
        if not _cfg["slow_mode"]:
            return line
        cap = int(float(_cfg["slow_cap"]))
    except Exception:
        return line
    out = []
    for t in line.split():
        if t[:1] == "F":
            try: out.append(f"F{min(int(float(t[1:])), cap)}")
            except ValueError: out.append(t)
        else:
            out.append(t)
    return " ".join(out)

def _send_scaled(line, st, vstate):
    """Send one line with live flow scaling. Returns send_motor_command's result
    (False = timeout/abort/serial error) so the caller can stop a dead print.
    In dry-run, the valve never opens and plunger (C) moves are stripped."""
    line = _cap_feedrate(line)                       # slow-mode travel cap (bad-wire workaround)
    if line.startswith("; __HOLD__"):                # needle-at-depth dwell (host-side; no board G4 needed)
        try: time.sleep(float(line.split()[-1]))
        except Exception: pass
        return True
    # droplet push speed: rewrite the compiled droplet feedrate (a B-only move) with the live
    # "drop speed" so a fast push flicks the drop off the needle instead of letting it bead.
    _u0 = line.upper()
    if _u0[:2] == "G1" and "B" in _u0 and not any(a in _u0 for a in ("X", "Y", "C", "Z", "A")):
        _dspd = int(_cfg.get("drop_f", 0) or 0)
        if _dspd > 0:
            line = " ".join(p for p in line.split() if p[:1].upper() != "F") + f" F{_dspd}"
    # Solenoid dispense: a droplet's B-plunger push (G1 B..., no X/Y/C) becomes a valve pulse to the
    # Arduino instead of a plunger move. The needle dip/hold/retract still happen normally.
    if _cfg.get("pneumatic"):
        u0 = line.upper()
        if u0[:2] == "G1" and "B" in u0 and "X" not in u0 and "Y" not in u0 and "C" not in u0:
            # PNEUMATIC: constant air pressure, plunger OUT. Valve OPEN TIME is the only thing that
            # sets droplet size. One droplet = one valve pulse; no stepper/plunger move at all.
            us = int(round(_live_valve_ms * 1000))   # LIVE: reflects the field even mid-print
            send_valve_command(f"F{max(0, us)}")
            return True
    if _cfg.get("solenoid"):
        u0 = line.upper()
        if u0[:2] == "G1" and "B" in u0 and "X" not in u0 and "Y" not in u0 and "C" not in u0:
            # MODEL B: the plunger meters the dose; the valve just gates the flow. Open the valve,
            # let the plunger do its fast push, then CLOSE it -- the close breaks the fluid column
            # at the tip (an air gap) so the next drop ejects discretely instead of a bead growing.
            send_valve_command("O")               # open the gate
            send_motor_command(line)              # fast plunger push (the actual dose) THROUGH the valve
            send_motor_command("M400")            # wait for the push to finish
            try: _sh = float(_cfg.get("solenoid_snap_ms", 0)) / 1000.0
            except Exception: _sh = 0.0
            if _sh > 0: time.sleep(_sh)            # brief settle so the drop clears the tip
            send_valve_command("C")               # close -> snap the column, break surface tension
            return True
    # droplet dose: the compiler emits the B plunger push at DROP_VOL mm/drop; scale it live
    # so the dispense can go as small as you want (no floor). B is cumulative-absolute, so
    # scaling every B by the same factor keeps the per-drop amount consistent.
    if "B" in line.upper() and line[:2].upper() in ("G0", "G1"):
        try: _df = float(_cfg["droplet_flow"])
        except Exception: _df = 1.0
        if _df != 1.0:
            line = " ".join((f"B{float(p[1:]) * _df:.4f}" if p[:1] in ("B", "b") and _mo_isfloat(p[1:]) else p)
                            for p in line.split())
    u = line.upper()
    dry = _cfg["dry_run"]
    if u.startswith("M106") and "S255" in u:
        if not dry: return send_motor_command(line)
        vstate[0] = True; return True
    elif u.startswith("M106") or u.startswith("M107"):
        vstate[0] = False; return send_motor_command("M106 P0 S0")
    elif u.startswith("G4") and vstate[0]:
        f = _cfg["droplet_flow"]; newP = None
        for p in line.split():
            if p[:1] == "P":
                try: newP = max(1, int(round(float(p[1:]) * f)))
                except ValueError: newP = None
        return send_motor_command(f"G4 P{newP}" if newP is not None else line)
    else:
        f = _cfg["scaffold_flow"]
        has_c = any(p[:1] == "C" for p in line.split())
        has_b = any(p[:1] == "B" for p in line.split())
        if (has_c or has_b) and dry:                 # dry-run: suppress ALL dispensing (C scaffold + B droplet)
            for p in line.split():                   # keep tracking C so a dry-run toggle can't desync/dump
                if p[:1] == "C":
                    try: st['oc'] = float(p[1:])
                    except ValueError: pass
            toks = [t for t in line.split() if t[:1] not in ("C", "B")]
            if len(toks) <= 1: return True
            return send_motor_command(" ".join(toks))
        if has_c:
            line = _scale_C(line, f, st)
            try: ef = float(_cfg["extrude_f"])
            except (ValueError, NameError): ef = 0
            if ef > 0 and (line.startswith("G0") or line.startswith("G1")):
                line = " ".join(t for t in line.split() if t[:1] != "F") + f" F{int(ef)}"
        return send_motor_command(line)

def _mo_isfloat(s):
    try: float(s); return True
    except ValueError: return False

def _manual_offset_transform(lines):
    """Manual-offset stopgap. Removes the automatic -X roll, drops in a pause marker
    where you hand-jog the droplet head into place, and shifts every droplet XY move
    back into scaffold coordinates (+offset). After your manual jog and a re-zero, the
    bed then only has to cover the reachable scaffold range instead of rolling the full
    offset it currently can't."""
    out, offset, in_drop = [], None, False
    for ln in lines:
        low = ln.lower()
        if "droplets (" in low:
            in_drop = True; out.append(ln); continue
        if not in_drop:
            out.append(ln); continue
        if "bed home (scaffold over corner)" in low:
            continue                                        # skip auto move to the corner
        if "droplet head over home" in low:                 # parse offset, swap for a pause
            try:
                xtok = next(t for t in ln.split() if t[:1] == "X")
                offset = round(HOME_X - float(xtok[1:]), 2)
            except Exception:
                offset = HEAD_OFFSET_X
            out.append("; __MANUAL_PAUSE__")
            continue
        if offset and ln[:3] in ("G0 ", "G1 ") and "X" in ln:
            out.append(" ".join(
                (f"X{float(t[1:]) + offset:.2f}" if t[:1] == "X" and _mo_isfloat(t[1:]) else t)
                for t in ln.split()))
            continue
        out.append(ln)
    return out

def execute_layered():
    global is_printing
    speed = _cfg["speed"]
    if not print_layers:
        _status("Status: nothing to print"); is_printing = False; return
    if _cfg["auto_park"]:                          # drive to corner via Set-Home reference
        _status("Status: parking to corner...")
        send_motor_command("G90"); send_motor_command("M82"); send_motor_command("M211 S0")
        send_motor_command(f"G0 A{HOME_SCAFFOLD_A + SAFE_CLEAR} Z{HOME_DROPLET_Z + SAFE_CLEAR} F{speed}")
        send_motor_command(f"G0 X{HOME_X} Y{HOME_Y} F{speed}")
        send_motor_command(f"G0 A{HOME_SCAFFOLD_A} Z{HOME_DROPLET_Z} F{speed}")
        send_motor_command("M400")
    st = {'oc': 0.0, 'sc': 0.0}; vstate = [False]
    send_motor_command("M302 P1")                     # allow cold extrusion (C/B are extruder axes)
    send_motor_command("M82")                         # ABSOLUTE extruder/lift moves (no drift on A0)
    _dpf = int(_cfg.get("drop_f", 3000) or 3000)      # let the fast droplet push actually reach speed:
    send_motor_command(f"M203 B{max(50, _dpf // 60 + 5)}")  # else Marlin clamps F to B's max feedrate
    send_motor_command(f"M201 B{int(_cfg.get('drop_accel', 30000))}")  # accel: a tiny push is accel-limited
    if _cfg["suction"] and not _cfg["dry_run"]:   # pull the plunger up to draw gel in
        _status("Status: loading gel (suction)...")
        try: load_mm = float(_cfg["load_mm"])
        except ValueError: load_mm = 30.0
        try: pull_f = int(float(_cfg["pull_f"]))
        except ValueError: pull_f = 150
        send_motor_command("G91")
        send_motor_command(f"G0 C-{load_mm} F{pull_f}")   # -C = plunger up = suction load
        send_motor_command("M400")
        send_motor_command("G90")
        send_motor_command("G92 C0")                  # reset plunger zero for the print
    for ln in print_header:
        if not is_printing: break
        _send_scaled(ln, st, vstate)
    # the header re-zeros Z to its spoof value; pin it back to the SAME reference Set Home uses
    # so the droplet Z (and therefore Go Home) returns to the exact spot you set.
    send_motor_command(f"G92 Z{HOME_DROPLET_Z}")
    N = len(print_layers)
    for li, layer in enumerate(print_layers):
        if not is_printing: break
        tag = "DRY-RUN " if _cfg["dry_run"] else ""
        _status(f"Status: {tag}Layer {li+1}/{N} printing...")
        lines = layer["lines"]
        if _cfg["manual_offset"]:
            lines = _manual_offset_transform(lines)
        total = max(1, len(lines))
        for i, ln in enumerate(lines):
            if not is_printing: break
            if ln.strip() == "; __MANUAL_PAUSE__":       # manual-offset stopgap: hand-jog here
                _status("Status: SCAFFOLD DONE \u2014 jog the droplet needle over the "
                                        "HOME corner, then press CONTINUE")
                _continue_evt.clear()
                root.after(0, lambda: continue_btn.config(state="normal"))
                while is_printing and not _continue_evt.is_set():
                    time.sleep(0.1)
                root.after(0, lambda: continue_btn.config(state="disabled"))
                if not is_printing: break
                send_motor_command(f"G92 X{HOME_X} Y{HOME_Y}")   # re-zero at the hand-jogged spot
                continue
            _paused_lifted = False
            while is_printing and not _run_evt.is_set():    # pause gate
                if not _paused_lifted:                      # entering pause: raise needle clear of the gel
                    send_motor_command("G91")               # so you can swap the syringe without dripping
                    send_motor_command(f"G0 Z{PAUSE_LIFT:.1f} F300")
                    send_motor_command("G90")
                    _paused_lifted = True
                _status(f"Status: PAUSED \u2014 swap syringe, then RESUME  (Layer {li+1}/{N})")
                time.sleep(0.1)
            if _paused_lifted and is_printing:              # resuming: drop the needle back to the exact spot
                send_motor_command("G91")
                send_motor_command(f"G0 Z-{PAUSE_LIFT:.1f} F300")
                send_motor_command("G90")
            if not is_printing: break
            if not _send_scaled(ln, st, vstate) and is_printing:
                _status("Status: serial timeout/error — print aborted")
                emergency_halt(); break
            if i % 4 == 0:
                root.after(0, lambda v=(i / total) * 100: scaffold_progress.config(value=v))
        if not is_printing: break
        root.after(0, lambda: scaffold_progress.config(value=100))
        last = (li == N - 1)
        # lift both heads clear (always, so nothing drags during the transition)
        if is_printing:
            send_motor_command(f"G0 A{HOME_SCAFFOLD_A + SAFE_CLEAR} Z{HOME_DROPLET_Z + SAFE_CLEAR} F{speed}")
        if _cfg.get("continuous") and not last:
            _status(f"Status: Layer {li+1}/{N} done \u2014 continuing (one-go)")
            continue                                  # straight into the next layer: no home/present/pause/dwell
        # return to the home corner
        if is_printing:
            send_motor_command(f"G0 X{HOME_X} Y{HOME_Y} F{speed}")
        # present forward (-Y) only for a spray pause (intermediate layer) or if the user
        # explicitly asked to present the finished print. Default = no forward move at the end,
        # so the bed just sits at home. This was the source of the big -Y lurch.
        if is_printing and (not last or _cfg["present"]):
            fwd = min(PRESENT_DY, BED_TRAVEL_Y - TRAVEL_SAFETY)
            py = max(0.0, HOME_Y - fwd)
            for _pcmd in (f"G0 Y{py:.1f} F{speed}", "M400"):
                if not is_printing: break
                send_motor_command(_pcmd)
        else:
            send_motor_command("M400")
        if last:
            if is_printing and _cfg.get("raise_heads"):   # lift both toolheads 2cm clear for access
                send_motor_command("G91")
                send_motor_command(f"G0 A20 Z20 F{speed}")  # +A raises scaffold, +Z raises droplet
                send_motor_command("M400")
                send_motor_command("G90")
            _status("Status: Print complete" +
                                    (" \u2014 presented" if _cfg["present"] else " (bed at home)"))
            break                                    # done
        # intermediate layer only: pause so you can spray, then return home for the next layer
        _status(f"Status: Layer {li+1}/{N} done — SPRAY, then press Continue")
        _continue_evt.clear()
        root.after(0, lambda: continue_btn.config(state="normal"))
        while is_printing and not _continue_evt.is_set():
            time.sleep(0.1)
        root.after(0, lambda: continue_btn.config(state="disabled"))
        if not is_printing: break
        dwell = _cfg.get("xlink_dwell", 0)                 # let the sprayed layer set before printing on it
        if dwell > 0:
            t_end = time.time() + dwell
            while is_printing and time.time() < t_end:
                _status(f"Status: crosslink dwell {max(0, t_end - time.time()):.0f}s\u2026")
                time.sleep(0.2)
            if not is_printing: break
        send_motor_command(f"G0 X{HOME_X} Y{HOME_Y} F{speed}")   # home before the next layer prints
        if send_motor_command("M114"):                    # Z-verify: confirm the heads stepped up a bead
            zrep = next((r for r in _rx_lines if ("A:" in r or "Z:" in r)), "")
            m = dict(re.findall(r"([AZ]):(-?\d+\.?\d*)", zrep)) if zrep else {}
            if m:
                _status(f"Status: Layer {li+2}/{N} \u2014 heads at A{m.get('A','?')} Z{m.get('Z','?')}")
    _status("Status: Print complete" if is_printing else "Status: Print aborted")
    is_printing = False

def continue_layer():
    _continue_evt.set()

def toggle_pause():
    if not is_printing: return
    if _run_evt.is_set():
        _run_evt.clear(); pause_btn.config(text="\u23f5 RESUME")
    else:
        _run_evt.set(); pause_btn.config(text="\u23f8 PAUSE")

def purge_gel():
    if is_printing:
        messagebox.showinfo("Purge", "Can't purge while printing."); return
    _abort_evt.clear()
    try: amt = float(purge_entry.get())
    except (ValueError, NameError): amt = 2.0
    try: ef = float(extrude_speed_entry.get())
    except (ValueError, NameError): ef = 0
    f = int(ef) if ef > 0 else 150
    send_motor_command("M302 P1"); send_motor_command("G91")
    send_motor_command(f"G1 C{amt} F{f}"); send_motor_command("M400"); send_motor_command("G90")
    scaffold_status_var.set(f"Status: purged {amt} mm")

def test_valve():
    """Dispense ONE test drop the way a print actually does it: by advancing the B plunger
    (the 'Extrude Left (B+)' axis). The old FAN0/M106 solenoid command does nothing on this
    machine, so this now pushes B by a small fixed amount instead."""
    if is_printing:
        messagebox.showinfo("Test dispense", "Can't test while printing."); return
    if not motor_serial:
        messagebox.showwarning("Not connected", "Connect the motor controller on the Control tab first."); return
    amt = 1.0                                     # mm of B to advance for one test drop
    if not send_motor_command("M302 P1"):         # allow cold extrusion; also probes for a response
        messagebox.showwarning("No response",
                               "The board didn't answer. Is it powered on and connected?")
        return
    send_motor_command("G91")                     # relative, so this one B push doesn't depend on B's zero
    send_motor_command(f"G1 B{amt:.3f} F{cfg_drop_f()}")
    send_motor_command("M400")
    send_motor_command("G90")
    scaffold_status_var.set(f"Test dispense: pushed B +{amt:.1f} mm — did gel come out the needle?")

def cfg_drop_f():
    """Feed for the B test push; mirrors the print's droplet feed. The 'extrude F' field
    uses 0 = auto, and anything too slow (< 10) would make a 1 mm push take ages and stall
    the port, so those fall back to a sane 150 mm/min."""
    try:
        v = int(float(extrude_speed_entry.get()))
        return v if v >= 10 else 150
    except (ValueError, NameError):
        return 150

def _line_time(line, pos, rel):
    """Rough seconds for one gcode line; updates pos and the [relative] flag."""
    if line.startswith("; __HOLD__"):                # needle-at-depth dwell (host-side)
        try: return float(line.split()[-1])
        except Exception: return 0.0
    toks = line.split()
    if not toks: return 0.0
    cmd = toks[0].upper()
    if cmd == "G90": rel[0] = False; return 0.0
    if cmd == "G91": rel[0] = True; return 0.0
    if cmd == "G92":
        for p in toks[1:]:
            if p[0] in "XYZABC":
                try: pos[p[0]] = float(p[1:])
                except ValueError: pass
        return 0.0
    if cmd in ("G0", "G1"):
        f = None; d2 = 0.0
        for p in toks[1:]:
            a = p[0]
            if a == "F":
                try: f = float(p[1:])
                except ValueError: pass
            elif a in "XYZABC":
                try: v = float(p[1:])
                except ValueError: continue
                if rel[0]:
                    d2 += v * v; pos[a] = pos.get(a, 0.0) + v
                else:
                    d2 += (v - pos.get(a, 0.0)) ** 2; pos[a] = v
        feed_mm_min = f or 1500
        try:                                          # honor Slow Mode's travel cap
            if slow_mode_var.get():
                feed_mm_min = min(feed_mm_min, float(slow_cap_var.get()))
        except Exception:
            pass
        feed = feed_mm_min / 60.0
        return (d2 ** 0.5) / feed if feed > 0 else 0.0
    if cmd == "G4":
        for p in toks[1:]:
            if p[0] == "P":
                try: return float(p[1:]) / 1000.0
                except ValueError: return 0.0
            if p[0] == "S":
                try: return float(p[1:])
                except ValueError: return 0.0
    return 0.02

def estimate_print():
    if not print_layers:
        messagebox.showinfo("Estimate", "Generate or load a print first."); return
    pos = {a: 0.0 for a in "XYZABC"}; rel = [False]; move_t = 0.0; rows = []
    for ln in print_header: move_t += _line_time(ln, pos, rel)
    for layer in print_layers:
        drops = sum(1 for l in layer["lines"] if "S255" in l.upper())
        t = sum(_line_time(l, pos, rel) for l in layer["lines"])
        rows.append((layer["i"], drops, t)); move_t += t
    # crosslink dwell happens in each gap between layers (N-1 of them)
    try: xlink = float(xlink_dwell_var.get())
    except Exception: xlink = 0.0
    gaps = max(0, len(print_layers) - 1)
    xlink_total = xlink * gaps
    total = move_t + xlink_total
    try: slow_on = slow_mode_var.get(); cap = int(float(slow_cap_var.get()))
    except Exception: slow_on, cap = False, 0
    msg = f"Total \u2248 {total/60:.1f} min  ·  {len(print_layers)} layers\n\n"
    msg += "\n".join(f"L{i}:  {d} drops  ·  ~{t/60:.1f} min" for i, d, t in rows)
    if xlink_total > 0:
        msg += f"\n\n+ crosslink dwell: {xlink_total/60:.1f} min  ({xlink:.0f}s x {gaps} gaps)"
    slow_note = f"  ·  Slow Mode capping travels at F{cap}" if slow_on else ""
    msg += f"\n\n(approximate{slow_note}; excludes manual spray pauses and accel)"
    messagebox.showinfo("Print estimate", msg)

def query_position():
    """Poll M114 when idle and update the read-out. Acquires the serial lock so it can
    never collide with a jog/print command; if the lock is busy it just skips this tick.
    Skipped entirely during prints and while (re)connecting."""
    try:
        if motor_serial and not is_printing and not _connecting and _serial_lock.acquire(blocking=False):
            try:
                _abort_evt.clear()                  # idle: clear any leftover stop flag
                motor_serial.reset_input_buffer()
                motor_serial.write(b"M114\n")
                buf = ""; t0 = time.time()
                while time.time() - t0 < 0.3:
                    if motor_serial.in_waiting > 0:
                        resp = motor_serial.readline().decode("utf-8", "ignore").strip()
                        if any(a + ":" in resp for a in "XYZABC"):
                            buf = resp
                        if "ok" in resp.lower():
                            break
                    else:
                        time.sleep(0.01)
                if buf:
                    vals = dict(re.findall(r"([XYZABC]):(-?\d+\.?\d*)", buf))
                    pos_var.set("  ".join(f"{a}{vals[a]}" for a in "XYZABC" if a in vals) or "\u2014")
            finally:
                _serial_lock.release()
    except Exception:
        pass
    root.after(1200, query_position)

def start_layered_print():
    global is_printing
    if is_printing: return
    if not print_layers:
        messagebox.showwarning("No G-code", "Generate in the Studio tab, or load a file first.")
        return
    if not motor_serial:
        if not messagebox.askokcancel("Not connected",
            "No motor controller is connected — this will run as a no-op (no motion).\n\n"
            "Connect hardware on the Control tab first, or continue for a UI test?"):
            return
    parked = ("Auto-park is ON: the bed drives to the corner first "
              "(needs Set Home done once this session).\n\n"
              if auto_park_var.get() else
              "Auto-park is OFF: printing starts from the current position.\n\n")
    dry = "DRY-RUN: valve and plunger are suppressed (motion only).\n\n" if dry_run_var.get() else ""
    if not messagebox.askokcancel("Ready to print?",
        "Per-layer print with a spray pause between layers. NO homing \u2014 it trusts "
        "the current position.\n\n" + dry + parked +
        "Confirm:\n"
        "\u2022 Bed at the corner, scaffold nozzle at the gel\n"
        "\u2022 Droplet needle at its Set-Home standoff\n"
        "\u2022 Syringes primed\n\n"
        "After each layer the bed presents forward; spray, then press Continue. Proceed?"):
        return
    _run_evt.set(); pause_btn.config(text="\u23f8 PAUSE")
    _abort_evt.clear()
    _snapshot_print_cfg()          # read all UI vars here on the MAIN thread; worker uses the snapshot
    is_printing = True; scaffold_progress['value'] = 0
    def _run_print():
        global is_printing
        try:
            execute_layered()
        except Exception as ex:                      # never let a crash leave the board moving
            _log_crash("execute_layered", type(ex), ex, ex.__traceback__)
            _status(f"Status: print error \u2014 {ex}")
            emergency_halt()
        finally:
            is_printing = False
    threading.Thread(target=_run_print, daemon=True).start()

def stop_scaffold_print():
    global is_printing
    is_printing = False
    _continue_evt.set()      # release any spray-pause wait
    _run_evt.set()           # release any pause wait
    emergency_halt()
    scaffold_status_var.set("Status: STOPPED")

def open_lab_website(event): webbrowser.open_new("https://labs.utdallas.edu/hbl/")
def jump_to_settings(): notebook.select(tab_control)

# Voron UI layout (added by ricky, check here for inaccuracies)
# ============================================================
# VORON / SLA PRINTER — MOONRAKER CONTROL
# Active printer.cfg target:
#   - Z motion only
#   - build plate lift
#   - manual vat rotation
#   - projector/SLA integration handled by printer-side config
# ============================================================

_voron_connected = False
_voron_polling = False
_voron_busy = False

# These Tk variables are created later when the Voron tab is built.
voron_host_var = None
voron_state_var = None
voron_connection_var = None
voron_z_var = None
voron_homed_var = None
voron_job_var = None
voron_progress_var = None
voron_step_var = None
voron_file_var = None


def _voron_base_url():
    """Return normalized Moonraker base URL from the UI field."""
    if voron_host_var is None:
        return ""

    host = voron_host_var.get().strip()

    if not host:
        return ""

    if not host.startswith(("http://", "https://")):
        host = "http://" + host

    # printer's uploaded moonraker.conf uses the standard 7125 port
    parsed = urllib.parse.urlsplit(host)

    if parsed.port is None:
        host = host.rstrip("/") + ":7125"

    return host.rstrip("/")


def _voron_log(msg):
    """Append one message to the Voron console safely."""
    try:
        stamp = time.strftime("%H:%M:%S")
        voron_log.configure(state="normal")
        voron_log.insert(tk.END, f"{stamp}  {msg}\n")
        voron_log.see(tk.END)
        voron_log.configure(state="disabled")
    except Exception:
        print("VORON:", msg)


def _voron_request(path, method="GET", data=None,
                   content_type="application/json", timeout=5.0):
    """
    Basic Moonraker HTTP request.

    Returns decoded JSON dictionary.
    Raises on connection / HTTP errors.
    """

    base = _voron_base_url()
    if not base:
        raise RuntimeError("No Voron host configured")

    url = base + path

    body = data
    if isinstance(data, str):
        body = data.encode("utf-8")

    req = urllib.request.Request(
        url,
        data=body,
        method=method
    )

    if body is not None and content_type:
        req.add_header("Content-Type", content_type)

    with urllib.request.urlopen(req, timeout=timeout) as resp:
        raw = resp.read()

    if not raw:
        return {}

    return json.loads(raw.decode("utf-8", "ignore"))


def _voron_post_json(path, payload=None, timeout=5.0):
    return _voron_request(
        path,
        method="POST",
        data=json.dumps(payload or {}),
        content_type="application/json",
        timeout=timeout
    )


def _voron_send_gcode(script):
    """
    Send ordinary G-code through Moonraker.
    Used only for controls supported by the active printer.cfg.
    """
    payload = {"script": script}
    return _voron_post_json(
        "/printer/gcode/script",
        payload,
        timeout=10.0
    )


def voron_connect():
    """
    Test Moonraker and Klipper state without freezing Tkinter.
    """
    global _voron_connected

    if _voron_busy:
        return

    host = _voron_base_url()

    if not host:
        messagebox.showwarning(
            "Voron Connection",
            "Enter the Voron hostname or IP address first."
        )
        return

    voron_connection_var.set("Connecting…")
    voron_connect_btn.config(state="disabled")

    def worker():
        global _voron_connected

        try:
            server = _voron_request("/server/info", timeout=4.0)
            printer = _voron_request("/printer/info", timeout=4.0)

            p_result = printer.get("result", {})
            state = p_result.get("state", "unknown")
            state_message = p_result.get("state_message", "")

            _voron_connected = True

            def done():
                voron_connection_var.set("Moonraker: Connected")
                voron_state_var.set(f"State: {state.upper()}")
                voron_connect_btn.config(
                    state="normal",
                    text="RECONNECT"
                )

                _voron_log("Connected to Moonraker")
                _voron_log(f"Printer state: {state}")

                if state_message:
                    _voron_log(state_message.replace("\n", " "))

                _start_voron_polling()

            root.after(0, done)

        except Exception as ex:
            _voron_connected = False

            def failed():
                voron_connection_var.set("Moonraker: OFFLINE")
                voron_state_var.set("State: unavailable")
                voron_connect_btn.config(
                    state="normal",
                    text="CONNECT"
                )
                _voron_log(f"Connection failed: {ex}")

            root.after(0, failed)

    threading.Thread(target=worker, daemon=True).start()


def _query_voron_status():
    """
    Read only objects relevant to the active Z-only SLA configuration.

    toolhead.position  -> live XYZ array (we display Z only)
    toolhead.homed_axes
    print_stats
    virtual_sdcard.progress
    """

    global _voron_busy

    if not _voron_connected or _voron_busy:
        return

    _voron_busy = True

    try:
        path = (
            "/printer/objects/query"
            "?toolhead=position,homed_axes"
            "&print_stats=state,filename,message"
            "&virtual_sdcard=progress"
        )

        result = _voron_request(path, timeout=3.0)
        status = result.get("result", {}).get("status", {})

        toolhead = status.get("toolhead", {})
        print_stats = status.get("print_stats", {})
        virtual_sd = status.get("virtual_sdcard", {})

        position = toolhead.get("position", [])
        if isinstance(position, list) and len(position) >= 3:
            try:
                z = float(position[2])
                voron_z_var.set(f"{z:.3f} mm")
            except Exception:
                pass

        homed = str(toolhead.get("homed_axes", ""))
        voron_homed_var.set(
            "Z" if "z" in homed.lower() else "Not homed"
        )

        job_state = print_stats.get("state", "standby")
        filename = print_stats.get("filename", "")

        if filename:
            voron_job_var.set(f"{filename}  ({job_state})")
        else:
            voron_job_var.set(job_state)

        try:
            progress = float(virtual_sd.get("progress", 0.0))
            pct = max(0.0, min(100.0, progress * 100.0))
            voron_progress_var.set(pct)
            voron_progress_label.config(text=f"{pct:.1f}%")
        except Exception:
            pass

        voron_state_var.set(f"State: {job_state.upper()}")

    except Exception as ex:
        _voron_log(f"Status update failed: {ex}")

    finally:
        _voron_busy = False


def _voron_poll_tick():
    if _voron_connected:
        threading.Thread(
            target=_query_voron_status,
            daemon=True
        ).start()

    root.after(1500, _voron_poll_tick)


def _start_voron_polling():
    global _voron_polling

    if _voron_polling:
        return

    _voron_polling = True
    root.after(100, _voron_poll_tick)


def voron_jog_z(direction):
    """
    Jog the actual configured build-plate Z axis.

    Uses relative positioning only for the requested jog, then restores
    absolute positioning.
    """

    if not _voron_connected:
        messagebox.showwarning(
            "Voron",
            "Connect to the Voron first."
        )
        return

    try:
        step = float(voron_step_var.get())
    except Exception:
        step = 1.0

    dz = step * direction

    script = (
        "G91\n"
        f"G1 Z{dz:.3f} F300\n"
        "G90"
    )

    def worker():
        try:
            _voron_send_gcode(script)
            root.after(
                0,
                lambda: _voron_log(
                    f"Jog Z {dz:+.3f} mm"
                )
            )
        except Exception as ex:
            root.after(
                0,
                lambda: _voron_log(
                    f"Z jog failed: {ex}"
                )
            )

    threading.Thread(target=worker, daemon=True).start()


def voron_home_z():
    if not _voron_connected:
        messagebox.showwarning(
            "Voron",
            "Connect to the Voron first."
        )
        return

    if not messagebox.askokcancel(
        "Home Z",
        "Home the Voron build-plate Z axis?"
    ):
        return

    def worker():
        try:
            _voron_send_gcode("G28 Z")
            root.after(
                0,
                lambda: _voron_log("Z homing command sent")
            )
        except Exception as ex:
            root.after(
                0,
                lambda: _voron_log(
                    f"Z homing failed: {ex}"
                )
            )

    threading.Thread(target=worker, daemon=True).start()


def voron_choose_gcode():
    path = filedialog.askopenfilename(
        title="Choose Voron G-code",
        filetypes=[
            ("G-code Files", "*.gcode"),
            ("All Files", "*.*")
        ]
    )

    if not path:
        return

    voron_file_var.set(path)


def _multipart_file_body(field_name, filepath):
    """
    Build a multipart/form-data upload body using only stdlib.
    Moonraker's /server/files/upload endpoint accepts multipart uploads.
    """

    boundary = "----PrintessVoron" + uuid.uuid4().hex

    filename = os.path.basename(filepath)

    mime = (
        mimetypes.guess_type(filename)[0]
        or "application/octet-stream"
    )

    with open(filepath, "rb") as f:
        payload = f.read()

    chunks = []

    chunks.append(
        f"--{boundary}\r\n"
        f'Content-Disposition: form-data; name="{field_name}"; '
        f'filename="{filename}"\r\n'
        f"Content-Type: {mime}\r\n\r\n"
    .encode("utf-8"))

    chunks.append(payload)
    chunks.append(f"\r\n--{boundary}--\r\n".encode("utf-8"))

    return (
        b"".join(chunks),
        f"multipart/form-data; boundary={boundary}"
    )


def voron_upload_gcode():
    """
    Upload selected G-code to Moonraker without automatically starting it.
    """

    if not _voron_connected:
        messagebox.showwarning(
            "Voron",
            "Connect to the Voron first."
        )
        return

    filepath = voron_file_var.get().strip()

    if not filepath:
        messagebox.showwarning(
            "Voron",
            "Choose a G-code file first."
        )
        return

    if not os.path.isfile(filepath):
        messagebox.showwarning(
            "Voron",
            "The selected G-code file no longer exists."
        )
        return

    voron_upload_btn.config(
        state="disabled",
        text="UPLOADING…"
    )

    def worker():
        try:
            body, ctype = _multipart_file_body(
                "file",
                filepath
            )

            result = _voron_request(
                "/server/files/upload",
                method="POST",
                data=body,
                content_type=ctype,
                timeout=30.0
            )

            item = result.get("result", {}).get("item", {})
            uploaded = (
                item.get("path")
                or os.path.basename(filepath)
            )

            def done():
                # Store printer-side filename after upload.
                voron_file_var.set(uploaded)

                voron_upload_btn.config(
                    state="normal",
                    text="UPLOAD"
                )

                _voron_log(
                    f"Uploaded G-code: {uploaded}"
                )

            root.after(0, done)

        except Exception as ex:

            def failed():
                voron_upload_btn.config(
                    state="normal",
                    text="UPLOAD"
                )

                _voron_log(
                    f"Upload failed: {ex}"
                )

            root.after(0, failed)

    threading.Thread(target=worker, daemon=True).start()


def voron_start_print():
    if not _voron_connected:
        messagebox.showwarning(
            "Voron",
            "Connect to the Voron first."
        )
        return

    filename = voron_file_var.get().strip()

    # If the field still contains a local absolute path, upload it first.
    if os.path.isfile(filename):
        messagebox.showinfo(
            "Voron",
            "Upload the selected file first, then press START PRINT."
        )
        return

    if not filename:
        messagebox.showwarning(
            "Voron",
            "Choose and upload a G-code file first."
        )
        return

    if not messagebox.askokcancel(
        "Start Voron Print",
        f"Start this job?\n\n{filename}"
    ):
        return

    quoted = urllib.parse.quote(filename, safe="/")

    def worker():
        try:
            _voron_post_json(
                "/printer/print/start"
                f"?filename={quoted}"
            )

            root.after(
                0,
                lambda: _voron_log(
                    f"Started print: {filename}"
                )
            )

        except Exception as ex:
            root.after(
                0,
                lambda: _voron_log(
                    f"Start failed: {ex}"
                )
            )

    threading.Thread(target=worker, daemon=True).start()


def voron_pause_print():
    if not _voron_connected:
        return

    def worker():
        try:
            _voron_post_json("/printer/print/pause")
            root.after(
                0,
                lambda: _voron_log("Pause requested")
            )
        except Exception as ex:
            root.after(
                0,
                lambda: _voron_log(
                    f"Pause failed: {ex}"
                )
            )

    threading.Thread(target=worker, daemon=True).start()


def voron_resume_print():
    if not _voron_connected:
        return

    def worker():
        try:
            _voron_post_json("/printer/print/resume")
            root.after(
                0,
                lambda: _voron_log("Resume requested")
            )
        except Exception as ex:
            root.after(
                0,
                lambda: _voron_log(
                    f"Resume failed: {ex}"
                )
            )

    threading.Thread(target=worker, daemon=True).start()


def voron_cancel_print():
    if not _voron_connected:
        return

    if not messagebox.askokcancel(
        "Cancel Voron Print",
        "Cancel the current Voron print?"
    ):
        return

    def worker():
        try:
            _voron_post_json("/printer/print/cancel")
            root.after(
                0,
                lambda: _voron_log(
                    "Print cancel requested"
                )
            )
        except Exception as ex:
            root.after(
                0,
                lambda: _voron_log(
                    f"Cancel failed: {ex}"
                )
            )

    threading.Thread(target=worker, daemon=True).start()


def voron_emergency_stop():
    """
    Moonraker emergency stop / Klipper shutdown.

    This intentionally requires confirmation because recovery generally
    requires firmware restart/reconnection.
    """

    if not _voron_connected:
        return

    if not messagebox.askyesno(
        "VORON EMERGENCY STOP",
        "Emergency-stop the Voron?\n\n"
        "This will shut down the Klipper printer controller and "
        "will require recovery/restart before printing again."
    ):
        return

    def worker():
        global _voron_connected

        try:
            _voron_post_json(
                "/printer/emergency_stop",
                timeout=4.0
            )

        except Exception:
            # Connection may disappear immediately after an E-stop.
            pass

        _voron_connected = False

        def done():
            voron_connection_var.set(
                "Moonraker: controller stopped"
            )
            voron_state_var.set(
                "State: EMERGENCY STOP"
            )
            _voron_log(
                "EMERGENCY STOP sent"
            )

        root.after(0, done)

    threading.Thread(target=worker, daemon=True).start()

#end of ricky voron commit

# --- THE VISUALS (Layout) ---
root = tk.Tk()
root.title(f"Bioprinting Technologies Demo  —  {BUILD_VERSION}")
root.geometry("1100x800") 

notebook = ttk.Notebook(root)
notebook.pack(fill="both", expand=True, padx=15, pady=15)

# THE NEW 4-TAB STRUCTURE
tab_home = ttk.Frame(notebook)
tab_matrix = ttk.Frame(notebook)
tab_control = ttk.Frame(notebook)
tab_scaffold = ttk.Frame(notebook)
tab_voron = ttk.Frame(notebook)

notebook.add(tab_home, text=" Home ")
# Studio replaces the old Matrix Generate tab
if _STUDIO_OK:
    tab_studio = StudioTab(notebook, on_generate=load_studio_gcode,
                           app_params=(get_app_gel_params, set_app_gel_params))
else:
    tab_studio = ttk.Frame(notebook)
    ttk.Label(tab_studio,
              text="Studio failed to load.\nKeep studio_tab.py and matrix_to_gcode.py\n"
                   "in the same folder as this file.\n\n" + str(_STUDIO_ERR),
              foreground="#FF6B6B", justify="center").pack(pady=60)
notebook.add(tab_studio, text=" Studio ")
notebook.add(tab_scaffold, text=" Print Control ")
notebook.add(tab_control, text=" Printess Control and Set Up ")
notebook.add(tab_voron, text=" Voron Control ")
# tab_matrix is intentionally NOT added to the notebook; Studio supersedes it.

# ==========================================
# --- 1. HOME TAB ---
# ==========================================
notebook.select(tab_home)
try:
    bg_image = Image.open("bg.png").resize((1100, 800))
    bg_photo = ImageTk.PhotoImage(bg_image)
    bg_label = tk.Label(tab_home, image=bg_photo)
    bg_label.image = bg_photo 
    bg_label.place(x=0, y=0, relwidth=1, relheight=1)
except: pass

home_wrapper = ttk.Frame(tab_home)
home_wrapper.place(relx=0.5, rely=0.5, anchor="center")
ttk.Label(home_wrapper, text="Bioprinting Technologies", font=("Segoe UI", 42, "bold"), foreground="#4DA8DA").pack(pady=(20, 5))
ttk.Label(home_wrapper, text="Hybrid Bioprinting Lab Microfluidic Control System", font=("Segoe UI", 16, "italic")).pack(pady=(0, 10))
link_label = ttk.Label(home_wrapper, text="🔗 labs.utdallas.edu/hbl", font=("Segoe UI", 11, "underline"), cursor="hand2", foreground="#4DA8DA")
link_label.pack(pady=(0, 20)); link_label.bind("<Button-1>", open_lab_website) 

dashboard_frame = ttk.Frame(home_wrapper); dashboard_frame.pack(pady=20)
card1 = ttk.LabelFrame(dashboard_frame, text="System Status"); card1.grid(row=0, column=0, padx=10, ipadx=20, ipady=10)
sys_status_var = tk.StringVar(value="Hardware: OFFLINE")
sys_status_label = ttk.Label(card1, textvariable=sys_status_var, font=("Segoe UI", 12, "bold"), foreground="#FF6B6B"); sys_status_label.pack(pady=10)

ttk.Button(home_wrapper, text="INITIALIZE SYSTEM", command=jump_to_settings, style="Accent.TButton", width=30).pack(pady=30, ipady=10)

# ==========================================
# --- 2. MATRIX GENERATE TAB ---
# ==========================================
# Left Side: All settings
mat_left = ttk.Frame(tab_matrix)
mat_left.pack(side=tk.LEFT, fill="y", padx=20, pady=15)

# Grid Settings
control_frame = ttk.LabelFrame(mat_left, text="Grid Setup", padding=10)
control_frame.pack(fill="x", pady=10)
ttk.Label(control_frame, text="Rows:").grid(row=0, column=0, padx=5, pady=5)
rows_entry = ttk.Entry(control_frame, width=5); rows_entry.insert(0, "3"); rows_entry.grid(row=0, column=1, padx=5, pady=5)
ttk.Label(control_frame, text="Cols:").grid(row=0, column=2, padx=5, pady=5)
cols_entry = ttk.Entry(control_frame, width=5); cols_entry.insert(0, "3"); cols_entry.grid(row=0, column=3, padx=5, pady=5)
ttk.Button(control_frame, text="Generate Blank Grid", command=generate_dynamic_grid).grid(row=1, column=0, columnspan=2, padx=5, pady=5, sticky="we")
ttk.Button(control_frame, text="Import .CSV", command=import_from_excel, style="Accent.TButton").grid(row=1, column=2, columnspan=2, padx=5, pady=5, sticky="we")

# Coordinates Settings
dim_frame = ttk.LabelFrame(mat_left, text="Spatial Coordinates", padding=10)
dim_frame.pack(fill="x", pady=10)
ttk.Label(dim_frame, text="X Distance:").grid(row=0, column=0, sticky="e", pady=5)
x_dist_entry = ttk.Entry(dim_frame, width=10); x_dist_entry.insert(0, "2"); x_dist_entry.grid(row=0, column=1, padx=5)
ttk.Label(dim_frame, text="Y Distance:").grid(row=0, column=2, sticky="e", pady=5)
y_dist_entry = ttk.Entry(dim_frame, width=10); y_dist_entry.insert(0, "2"); y_dist_entry.grid(row=0, column=3, padx=5)
ttk.Label(dim_frame, text="X Start:").grid(row=1, column=0, sticky="e", pady=5)
x_start_entry = ttk.Entry(dim_frame, width=10); x_start_entry.insert(0, "0"); x_start_entry.grid(row=1, column=1, padx=5)
ttk.Label(dim_frame, text="Y Start:").grid(row=1, column=2, sticky="e", pady=5)
y_start_entry = ttk.Entry(dim_frame, width=10); y_start_entry.insert(0, "0"); y_start_entry.grid(row=1, column=3, padx=5)

# Z-Heights
z_frame = ttk.LabelFrame(mat_left, text="Nozzle Z-Heights (Spoofing Z=40)", padding=10)
z_frame.pack(fill="x", pady=10)
ttk.Label(z_frame, text="Safe Travel (mm):").grid(row=0, column=0, sticky="e", pady=5)
z_safe_entry = ttk.Entry(z_frame, width=10); z_safe_entry.insert(0, "5.0"); z_safe_entry.grid(row=0, column=1, padx=5)
ttk.Label(z_frame, text="Print Drop (mm):").grid(row=0, column=2, sticky="e", pady=5)
z_print_entry = ttk.Entry(z_frame, width=10); z_print_entry.insert(0, "1.0"); z_print_entry.grid(row=0, column=3, padx=5)

# Valve Settings
valve_frame = ttk.LabelFrame(mat_left, text="Valve Parameters", padding=10)
valve_frame.pack(fill="x", pady=10)
ttk.Label(valve_frame, text="Channel:").grid(row=0, column=0, sticky="e", pady=5)
channel_combo = ttk.Combobox(valve_frame, values=["Channel 1", "Channel 2", "Channel 3", "Channel 4"], width=10); channel_combo.current(0); channel_combo.grid(row=0, column=1, padx=5)
ttk.Label(valve_frame, text="Pulse (us):").grid(row=0, column=2, sticky="e", pady=5)
pulse_entry = ttk.Entry(valve_frame, width=10); pulse_entry.insert(0, "2000"); pulse_entry.grid(row=0, column=3, padx=5)

# Right Side: The Map
matrix_frame = ttk.LabelFrame(tab_matrix, text="Drop Pattern Map (Drops per Well)", padding=20)
matrix_frame.pack(side=tk.RIGHT, fill="both", expand=True, padx=20, pady=15)
generate_dynamic_grid()

# ==========================================
# --- 3. PRINTESS CONTROL AND SET UP TAB ---
# ==========================================
# Top Panel: Connection
conn_frame = ttk.LabelFrame(tab_control, text="1. Hardware Connection", padding=15)
conn_frame.pack(fill="x", padx=20, pady=10)

ttk.Label(conn_frame, text="MOTOR COM:").grid(row=0, column=0, padx=5)
initial_ports = get_active_ports()
motor_com_combo = ttk.Combobox(conn_frame, values=initial_ports, width=12, state="readonly")
motor_com_combo.grid(row=0, column=1, padx=5)

ttk.Label(conn_frame, text="VALVE COM:").grid(row=0, column=2, padx=15)
valve_com_combo = ttk.Combobox(conn_frame, values=initial_ports, width=12, state="readonly")
valve_com_combo.grid(row=0, column=3, padx=5)

if initial_ports and initial_ports[0] != "No Ports Found":
    motor_com_combo.current(0)
    if len(initial_ports) > 1: valve_com_combo.current(1)
    else: valve_com_combo.current(0)

ttk.Button(conn_frame, text="↻ Refresh", command=refresh_com_ports).grid(row=0, column=4, padx=15)
connect_btn = ttk.Button(conn_frame, text="CONNECT ALL HARDWARE", style="Accent.TButton", command=connect_all_hardware)
connect_btn.grid(row=0, column=5, padx=15)

# Middle Section Wrapper
ctrl_wrapper = ttk.Frame(tab_control)
ctrl_wrapper.pack(fill="both", expand=True, padx=10, pady=10)

# Left Side of Middle: Jogging
jog_frame = ttk.LabelFrame(ctrl_wrapper, text="2. Manual Stage Jogging (Full Toolhead Control)", padding=15)
jog_frame.pack(side=tk.LEFT, fill="both", expand=True, padx=10)

speed_frame = ttk.Frame(jog_frame)
speed_frame.pack(fill="x", pady=(0,10))
ttk.Label(speed_frame, text="Jog Speed:").grid(row=0, column=0, padx=5)
idle_speed_entry = ttk.Entry(speed_frame, width=8); idle_speed_entry.insert(0, "300"); idle_speed_entry.grid(row=0, column=1, padx=5)
ttk.Label(speed_frame, text="Print Speed:").grid(row=0, column=2, padx=15)
work_speed_entry = ttk.Entry(speed_frame, width=8); work_speed_entry.insert(0, "300"); work_speed_entry.grid(row=0, column=3, padx=5)

# live position read-out (polled M114 when idle)
pos_var = tk.StringVar(value="—")
pos_frame = ttk.Frame(jog_frame); pos_frame.pack(fill="x", pady=(0, 6))
ttk.Label(pos_frame, text="Position:", font=("Segoe UI", 9, "bold")).pack(side=tk.LEFT, padx=5)
ttk.Label(pos_frame, textvariable=pos_var, font=("Consolas", 10), foreground="#A8E6CF").pack(side=tk.LEFT)

step_size_var = tk.DoubleVar(value=1.0)
step_frame = ttk.Frame(jog_frame)
step_frame.pack(fill="x", pady=10)
ttk.Label(step_frame, text="Step Size (mm):", font=("Segoe UI", 10, "bold")).pack(side=tk.LEFT, padx=5)
ttk.Radiobutton(step_frame, text="0.1", variable=step_size_var, value=0.1).pack(side=tk.LEFT, padx=5)
ttk.Radiobutton(step_frame, text="1.0", variable=step_size_var, value=1.0).pack(side=tk.LEFT, padx=5)
ttk.Radiobutton(step_frame, text="10.0", variable=step_size_var, value=10.0).pack(side=tk.LEFT, padx=5)

dpad_frame = ttk.Frame(jog_frame)
dpad_frame.pack(pady=10)
# X / Y Translation Pad
ttk.Button(dpad_frame, text="Y +", width=6, command=lambda: jog('Y', 1)).grid(row=0, column=1, pady=2)
ttk.Button(dpad_frame, text="X -", width=6, command=lambda: jog('X', -1)).grid(row=1, column=0, padx=2)
ttk.Button(dpad_frame, text="Go Home", width=6, style="Accent.TButton", command=go_home).grid(row=1, column=1)
ttk.Button(dpad_frame, text="X +", width=6, command=lambda: jog('X', 1)).grid(row=1, column=2, padx=2)
ttk.Button(dpad_frame, text="Y -", width=6, command=lambda: jog('Y', -1)).grid(row=2, column=1, pady=2)

# LEFT TOOL (Droplet Matrix: Z and B Axes)
ttk.Button(dpad_frame, text="Lift Left (Z+)", width=14, command=lambda: jog('Z', 1)).grid(row=0, column=3, padx=(30,5), pady=2)
ttk.Button(dpad_frame, text="Drop Left (Z-)", width=14, command=lambda: jog('Z', -1)).grid(row=1, column=3, padx=(30,5), pady=2)
ttk.Button(dpad_frame, text="Extrude Left (B+)", width=14, command=lambda: jog('B', 1)).grid(row=2, column=3, padx=(30,5), pady=2)
ttk.Button(dpad_frame, text="Retract Left (B-)", width=14, command=lambda: jog('B', -1)).grid(row=3, column=3, padx=(30,5), pady=2)

# RIGHT TOOL (Scaffold Extruder: A and C Axes)
ttk.Button(dpad_frame, text="Lift Right (A+)", width=14, command=lambda: jog('A', 1)).grid(row=0, column=4, padx=5, pady=2)
ttk.Button(dpad_frame, text="Drop Right (A-)", width=14, command=lambda: jog('A', -1)).grid(row=1, column=4, padx=5, pady=2)
ttk.Button(dpad_frame, text="Extrude Right (C+)", width=14, command=lambda: jog('C', 1)).grid(row=2, column=4, padx=5, pady=2)
ttk.Button(dpad_frame, text="Retract Right (C-)", width=14, command=lambda: jog('C', -1)).grid(row=3, column=4, padx=5, pady=2)

# Set Home: jog to the corner pose first, then press this to declare it home
home_set_frame = ttk.Frame(jog_frame)
home_set_frame.pack(pady=(4, 0))
ttk.Button(home_set_frame, text="Set Home (declare current pose)", command=set_home).pack(fill="x", ipady=4)
ttk.Label(home_set_frame, text="Jog to corner + heights, press Set Home, then Go Home returns to it.",
          font=("Segoe UI", 8, "italic")).pack(pady=(2, 0))


# Right Side of Middle: Execution Panel
exec_frame = ttk.LabelFrame(ctrl_wrapper, text="3. Execution Panel", padding=20)
exec_frame.pack(side=tk.RIGHT, fill="both", expand=True, padx=10)

ttk.Label(exec_frame, text="Design + slice in Studio, then run\nthe streams from the Print Control tab.", font=("Segoe UI", 10, "bold"), foreground="#A8E6CF").pack(pady=10)

ttk.Button(exec_frame, text="\u25b6 GO TO PRINT CONTROL", style="Accent.TButton", command=lambda: notebook.select(tab_scaffold)).pack(fill="x", pady=20, ipady=15)
ttk.Button(exec_frame, text="\u23f9 EMERGENCY STOP", command=stop_print_job).pack(fill="x", pady=5, ipady=5)

status_var = tk.StringVar(value="Status: Ready to Print")
ttk.Label(exec_frame, textvariable=status_var, font=("Segoe UI", 12, "bold")).pack(pady=20)
progress_bar = ttk.Progressbar(exec_frame, orient="horizontal", mode="determinate", value=0)
progress_bar.pack(fill="x", ipady=5)

# --- Manual G-code console (diagnostics): type a line, press Send or Enter ---
manual_frame = ttk.LabelFrame(tab_control, text="Manual G-code (diagnostics)", padding=12)
manual_frame.pack(fill="x", padx=20, pady=(0, 12))
_mrow = ttk.Frame(manual_frame); _mrow.pack(fill="x")
manual_entry = ttk.Entry(_mrow)
manual_entry.pack(side=tk.LEFT, fill="x", expand=True, padx=(0, 8))
manual_out = tk.StringVar(value="type one line, e.g.  G0 X33.5 F800   \u2014  Enter sends. Each G0 X is absolute.")

# scrollable reply log — M114/M122/M503 come back as multi-line text, shown here
_logrow = ttk.Frame(manual_frame); _logrow.pack(fill="x", pady=(6, 0))
manual_log = tk.Text(_logrow, height=8, wrap="none", font=("Consolas", 8),
                     bg="#12161c", fg="#A8E6CF", insertbackground="#A8E6CF")
_logsb = ttk.Scrollbar(_logrow, command=manual_log.yview)
manual_log.configure(yscrollcommand=_logsb.set, state="disabled")
manual_log.pack(side=tk.LEFT, fill="both", expand=True)
_logsb.pack(side=tk.RIGHT, fill="y")

def _log(text):
    manual_log.configure(state="normal")
    manual_log.insert(tk.END, text + "\n")
    manual_log.see(tk.END)
    manual_log.configure(state="disabled")

def send_manual(_evt=None):
    if is_printing:
        manual_out.set("busy: a print is running \u2014 stop it first"); return
    line = manual_entry.get().strip()
    if not line:
        return
    manual_entry.delete(0, tk.END)
    ok = send_motor_command(line)
    _log(f">> {line}")
    for r in _rx_lines:                       # echo whatever the board sent back
        _log(f"   {r}")
    if not _rx_lines:
        _log(f"   ({'ok' if ok else 'no reply / timeout'})")
    manual_out.set(("ok    " if ok else "ERR/timeout   ") + line)

manual_entry.bind("<Return>", send_manual)
ttk.Button(_mrow, text="Send", command=send_manual).pack(side=tk.LEFT)
for _lbl, _cmd in (("M114", "M114"), ("M503", "M503"), ("M122", "M122")):
    ttk.Button(_mrow, text=_lbl, width=6,
               command=lambda c=_cmd: (manual_entry.delete(0, tk.END),
                                       manual_entry.insert(0, c), send_manual())).pack(side=tk.LEFT, padx=(6, 0))
ttk.Label(manual_frame, textvariable=manual_out, foreground="#A8E6CF",
          font=("Consolas", 9)).pack(anchor="w", pady=(6, 0))
ttk.Label(manual_frame,
          text="M503 = read steps/mm (M92 X__) & limits    |    M122 = TMC driver status (X fault flags)\n"
               "Scale/direction test:  G92 X71  \u2192  G0 X40 F800  \u2192  M114  \u2192  G0 X71 F800  \u2192  M114",
          font=("Consolas", 8), foreground="#888", justify="left").pack(anchor="w", pady=(4, 0))

# ==========================================
# --- 4. PRINT CONTROL TAB (view combined G-code + run the print) ---
# ==========================================
scaffold_flow_var = tk.DoubleVar(value=1.0)
droplet_flow_var  = tk.DoubleVar(value=1.0)
auto_park_var     = tk.BooleanVar(value=False)
manual_offset_var = tk.BooleanVar(value=False)
hard_stop_var     = tk.BooleanVar(value=False)    # M410 by default; M112 (needs reconnect) only if ticked
present_var       = tk.BooleanVar(value=False)    # move bed forward (-Y) to present when done
raise_heads_var   = tk.BooleanVar(value=True)    # raise both toolheads 2cm after a finished print
xlink_dwell_var   = tk.StringVar(value="10")       # seconds to wait after spray/Continue before next layer
continuous_var    = tk.BooleanVar(value=False)    # print all layers in one go (no spray pause/present)
solenoid_var      = tk.BooleanVar(value=True)     # gate the plunger push with the valve (model B)
pneumatic_var     = tk.BooleanVar(value=True)    # constant air pressure; valve time alone meters drops
valve_time_ms_var = tk.StringVar(value="1.0")     # pneumatic: valve open time per droplet (ms)
_live_valve_ms = 1.0     # live mirror of valve_time_ms_var, safe to read from the print thread
def _update_live_valve_ms(*_a):
    """Keep _live_valve_ms in sync with the field so valve open time can be tuned MID-PRINT."""
    global _live_valve_ms
    try: _live_valve_ms = max(0.0, float(valve_time_ms_var.get()))
    except (ValueError, tk.TclError): pass
valve_time_ms_var.trace_add("write", _update_live_valve_ms)
solenoid_us_var   = tk.StringVar(value="2000")     # (model A only) fire pulse us when valve meters
solenoid_snap_var = tk.StringVar(value="0")        # ms to wait after the push before closing the valve
slow_mode_var     = tk.BooleanVar(value=True)     # cap travel feedrate (bad-wire workaround)
slow_cap_var      = tk.StringVar(value="300")

_scaf_outer = ttk.Frame(tab_scaffold)
_scaf_outer.pack(side=tk.LEFT, fill="y", padx=20, pady=15)
_scaf_canvas = tk.Canvas(_scaf_outer, width=300, highlightthickness=0)
_scaf_vsb = ttk.Scrollbar(_scaf_outer, orient="vertical", command=_scaf_canvas.yview)
_scaf_canvas.configure(yscrollcommand=_scaf_vsb.set)
_scaf_vsb.pack(side="right", fill="y"); _scaf_canvas.pack(side="left", fill="both", expand=True)
scaf_left = ttk.Frame(_scaf_canvas)
_scaf_win = _scaf_canvas.create_window((0, 0), window=scaf_left, anchor="nw")
scaf_left.bind("<Configure>", lambda e: _scaf_canvas.configure(scrollregion=_scaf_canvas.bbox("all")))
_scaf_canvas.bind("<Configure>", lambda e: _scaf_canvas.itemconfig(_scaf_win, width=e.width))
def _scaf_wheel(e): _scaf_canvas.yview_scroll(int(-e.delta / 120), "units")
_scaf_outer.bind("<Enter>", lambda e: _scaf_canvas.bind_all("<MouseWheel>", _scaf_wheel))
_scaf_outer.bind("<Leave>", lambda e: _scaf_canvas.unbind_all("<MouseWheel>"))

# 1. source
src_frame = ttk.LabelFrame(scaf_left, text="1. G-Code Source", padding=15)
src_frame.pack(fill="x", pady=8)
ttk.Label(src_frame, text="Generate the combined print in the\nStudio tab, or load a file here.",
          foreground="#A8E6CF").pack(pady=(0, 8))
ttk.Button(src_frame, text="Load G-Code File", command=load_gcode_file, style="Accent.TButton").pack(fill="x", ipady=4)

# 2. live flow sliders
flow_frame = ttk.LabelFrame(scaf_left, text="2. Flow Rate (live)", padding=15)
flow_frame.pack(fill="x", pady=8)
ttk.Label(flow_frame, text="Scaffold head  (C squeeze)").pack(anchor="w")
sflab = ttk.Label(flow_frame, text="1.00\u00d7", font=("Segoe UI", 9, "bold")); sflab.pack(anchor="e")
ttk.Scale(flow_frame, from_=0.2, to=3.0, variable=scaffold_flow_var, orient="horizontal",
          command=lambda v: sflab.config(text=f"{float(v):.2f}\u00d7")).pack(fill="x")
ttk.Label(flow_frame, text="Droplet head  (B dose \u00d7 per drop)").pack(anchor="w", pady=(8, 0))
dflab = ttk.Label(flow_frame, text="1.00\u00d7", font=("Segoe UI", 9, "bold")); dflab.pack(anchor="e")
ttk.Scale(flow_frame, from_=0.02, to=3.0, variable=droplet_flow_var, orient="horizontal",
          command=lambda v: dflab.config(text=f"{float(v):.2f}\u00d7")).pack(fill="x")
_dfrow = ttk.Frame(flow_frame); _dfrow.pack(fill="x", pady=(2, 0))
ttk.Label(_dfrow, text="exact \u00d7", width=7).pack(side="left")
_dfent = ttk.Entry(_dfrow, width=7); _dfent.pack(side="left")
_dfent.insert(0, "0.10")
def _set_dose(_e=None):
    try:
        v = float(_dfent.get()); droplet_flow_var.set(v); dflab.config(text=f"{v:.2f}\u00d7")
    except ValueError: pass
_dfent.bind("<Return>", _set_dose)
ttk.Button(_dfrow, text="set", width=4, command=_set_dose).pack(side="left", padx=(4, 0))

# --- PNEUMATIC MODE: constant air pressure, valve open time is the only droplet control ---
pneu_frame = ttk.LabelFrame(scaf_left, text="Pneumatic Dispense", padding=15)
pneu_frame.pack(fill="x", pady=8)
ttk.Checkbutton(pneu_frame,
                text="Pneumatic mode (plunger OUT \u2014 valve time sets the droplet)",
                variable=pneumatic_var).pack(anchor="w")
_pt = ttk.Frame(pneu_frame); _pt.pack(fill="x", pady=(6, 0))
ttk.Label(_pt, text="valve open time (ms)", width=18).pack(side="left")
ttk.Entry(_pt, width=8, textvariable=valve_time_ms_var).pack(side="left")
ttk.Label(pneu_frame,
          text="Droplet size = this time \u00d7 regulator pressure. Bigger time = bigger drop.\n"
               "Set pressure at the regulator (~3-5 psi); dial the drop with this number.\n"
               "This can be changed DURING a print \u2014 it takes effect on the next droplet.",
          foreground="#8FA5B8", font=("Segoe UI", 8), justify="left").pack(anchor="w", pady=(4, 0))

def _apply_pneumatic(*_a):
    """When pneumatic mode is on, disable everything tied to stepper/plunger droplet dispensing."""
    on = pneumatic_var.get()
    st_stepper = "disabled" if on else "normal"
    for w in (_dd_ent, _ds_ent, _da_ent):
        try: w.config(state=st_stepper)
        except Exception: pass
    try: _solchk.config(state=st_stepper)      # the plunger-gate checkbox is meaningless in pneumatic
    except Exception: pass
    try: _release_btn.config(text="\u25cf Release 1 droplet (valve pulse)" if on
                             else "\u25cf Release 1 droplet (dial in size)")
    except Exception: pass
pneumatic_var.trace_add("write", _apply_pneumatic)

# 3. plunger / gel load (ported from the original app)
suction_load_var = tk.BooleanVar(value=False)
plunger_frame = ttk.LabelFrame(scaf_left, text="3. Plunger / Gel Load", padding=15)
plunger_frame.pack(fill="x", pady=8)
ttk.Checkbutton(plunger_frame, text="Suction-load gel before print", variable=suction_load_var).pack(anchor="w")
pr = ttk.Frame(plunger_frame); pr.pack(fill="x", pady=2)
ttk.Label(pr, text="load vol (mm)", width=12).pack(side="left")
load_volume_entry = ttk.Entry(pr, width=8); load_volume_entry.insert(0, "30"); load_volume_entry.pack(side="left")
pr2 = ttk.Frame(plunger_frame); pr2.pack(fill="x", pady=2)
ttk.Label(pr2, text="pull-up F", width=12).pack(side="left")
pullup_speed_entry = ttk.Entry(pr2, width=8); pullup_speed_entry.insert(0, "150"); pullup_speed_entry.pack(side="left")
pr3 = ttk.Frame(plunger_frame); pr3.pack(fill="x", pady=2)
ttk.Label(pr3, text="extrude F (0=auto)", width=12).pack(side="left")
extrude_speed_entry = ttk.Entry(pr3, width=8); extrude_speed_entry.insert(0, "0"); extrude_speed_entry.pack(side="left")
pr4 = ttk.Frame(plunger_frame); pr4.pack(fill="x", pady=2)
ttk.Label(pr4, text="purge (mm)", width=12).pack(side="left")
purge_entry = ttk.Entry(pr4, width=8); purge_entry.insert(0, "2"); purge_entry.pack(side="left")
ttk.Button(plunger_frame, text="Purge / Prime now", command=purge_gel).pack(fill="x", pady=(4, 0))

# --- single droplet dial-in: dose + speed, and a button that fires exactly one ---
_dd = ttk.Frame(plunger_frame); _dd.pack(fill="x", pady=(8, 2))
ttk.Label(_dd, text="drop dose (mm)", width=13).pack(side="left")
drop_dose_var = tk.StringVar(value="0.10")
_dd_ent=ttk.Entry(_dd, width=8, textvariable=drop_dose_var); _dd_ent.pack(side="left")
_ds = ttk.Frame(plunger_frame); _ds.pack(fill="x", pady=2)
ttk.Label(_ds, text="drop speed (F)", width=13).pack(side="left")
drop_speed_var = tk.StringVar(value="3000")
_ds_ent=ttk.Entry(_ds, width=8, textvariable=drop_speed_var); _ds_ent.pack(side="left")
_da = ttk.Frame(plunger_frame); _da.pack(fill="x", pady=2)
ttk.Label(_da, text="drop accel", width=13).pack(side="left")
drop_accel_var = tk.StringVar(value="30000")
_da_ent=ttk.Entry(_da, width=8, textvariable=drop_accel_var); _da_ent.pack(side="left")
_dcap = ttk.Label(plunger_frame, text="", foreground="#8FA5B8", font=("Segoe UI", 8), justify="left")
_dcap.pack(anchor="w")

def _refresh_drop_cap(*_a):
    """A tiny push is ACCEL-limited: over 0.05-0.08mm the plunger never reaches 'drop speed',
    so acceleration is what sets how hard the drop is flicked. Show the real achievable speed."""
    try:
        d = float(drop_dose_var.get()); a = float(drop_accel_var.get()); f = float(drop_speed_var.get())
        if d <= 0 or a <= 0: _dcap.config(text=""); return
        v_cap = (a * d) ** 0.5              # peak speed of a triangular accel/decel move, mm/s
        if f / 60.0 <= v_cap:
            _dcap.config(text=f"{d:.3f}mm push reaches F{int(f)} (speed-limited). "
                              f"Raising accel won't help here.")
        else:
            _dcap.config(text=f"{d:.3f}mm push only reaches F{int(v_cap*60)} of the F{int(f)} asked "
                              f"\u2014 ACCEL-LIMITED.\nRaise 'drop accel' to flick smaller drops off.")
    except Exception:
        _dcap.config(text="")
for _v in (drop_dose_var, drop_speed_var, drop_accel_var):
    _v.trace_add("write", _refresh_drop_cap)
_refresh_drop_cap()

ttk.Label(plunger_frame,
          text="Fast push flicks the drop off before it beads. Smaller dose = smaller drop.",
          foreground="#8FA5B8", font=("Segoe UI", 8), justify="left").pack(anchor="w")

def release_one_droplet():
    """Fire exactly one droplet at the current settings, so you can dial in a small drop.
    Pneumatic mode: one valve pulse (no plunger). Otherwise: a fast high-accel plunger push,
    optionally gated by the valve."""
    if is_printing:
        messagebox.showinfo("Release droplet", "Can't dispense while printing."); return

    if pneumatic_var.get():                              # PNEUMATIC: one valve pulse, no plunger/motor
        if not (valve_serial and valve_serial.is_open):
            messagebox.showwarning("Valve not connected",
                                   "Pneumatic mode dispenses through the valve only. Set the VALVE "
                                   "COM dropdown to the Nano and connect first.")
            return
        try:    ms = float(valve_time_ms_var.get())
        except ValueError: messagebox.showwarning("Release droplet", "Valve time isn't a number."); return
        us = max(0, int(round(ms * 1000)))
        threading.Thread(target=lambda: send_valve_command(f"F{us}"), daemon=True).start()
        scaffold_status_var.set(f"Released 1 droplet: valve open {ms:.2f} ms (pneumatic)")
        return

    if not motor_serial:
        messagebox.showwarning("Not connected", "Connect the motor controller first."); return
    try:    dose = float(drop_dose_var.get())
    except ValueError: messagebox.showwarning("Release droplet", "Drop dose isn't a number."); return
    if dose <= 0:
        messagebox.showinfo("Release droplet", "Drop dose is 0 \u2014 raise it above 0."); return
    try:    spd = int(float(drop_speed_var.get()))
    except ValueError: spd = 3000
    gate = solenoid_var.get()
    if gate and not (valve_serial and valve_serial.is_open):
        messagebox.showwarning("Valve gate on, but no valve",
                               "'Solenoid valve gates plunger' is checked, but the valve isn't "
                               "connected. Set the VALVE COM dropdown to the Nano and reconnect, "
                               "or uncheck the gate to dispense with the plunger only.")
        return
    gate = gate and (valve_serial and valve_serial.is_open)
    send_motor_command("M302 P1")                        # allow cold extrusion
    send_motor_command(f"M203 B{max(50, spd // 60 + 5)}")# raise B speed limit, else F is clamped
    try:    dacc = int(float(drop_accel_var.get()))
    except Exception: dacc = 30000
    send_motor_command(f"M201 B{dacc}")                  # accel is the real lever on a tiny push
    if gate: send_valve_command("O")                     # open the gate
    send_motor_command("G91")                            # relative push
    send_motor_command(f"G1 B{dose:.4f} F{spd}")         # the metered dose, fast
    send_motor_command("M400")
    send_motor_command("G90")
    if gate:
        try: _sh = int(float(solenoid_snap_var.get())) / 1000.0
        except Exception: _sh = 0.0
        if _sh > 0: time.sleep(_sh)
        send_valve_command("C")                          # close -> snap the column
    scaffold_status_var.set(f"Released 1 droplet: B +{dose:.4f} mm at F{spd}"
                            + ("  \u00b7 valve-gated" if gate else ""))

_release_btn = ttk.Button(plunger_frame, text="\u25cf Release 1 droplet (dial in size)",
           command=release_one_droplet, style="Accent.TButton"); _release_btn.pack(fill="x", pady=(4, 0))

pr5 = ttk.Frame(plunger_frame); pr5.pack(fill="x", pady=(6, 2))
ttk.Label(pr5, text="valve test ms", width=12).pack(side="left")
valve_test_ms = tk.StringVar(value="150")
ttk.Entry(pr5, width=8, textvariable=valve_test_ms).pack(side="left")
ttk.Button(plunger_frame, text="Test droplet valve (fire once)", command=test_valve).pack(fill="x", pady=(2, 0))

# 4. run (per-layer print with spray pause)
run_frame = ttk.LabelFrame(scaf_left, text="4. Run", padding=15)
run_frame.pack(fill="x", pady=8)
dry_run_var = tk.BooleanVar(value=False)
ttk.Checkbutton(run_frame, text="Dry run (no valve / no plunger)", variable=dry_run_var).pack(anchor="w")
ttk.Checkbutton(run_frame, text="Auto-park to corner first (needs Set Home)", variable=auto_park_var).pack(anchor="w", pady=(0, 6))
ttk.Checkbutton(run_frame, text="Manual droplet offset (pause after scaffold, hand-jog the 31mm)",
                variable=manual_offset_var).pack(anchor="w", pady=(0, 6))
_slowrow = ttk.Frame(run_frame); _slowrow.pack(fill="x", pady=(0, 2))
ttk.Checkbutton(_slowrow, text="Slow mode: cap all moves to F", variable=slow_mode_var).pack(side="left")
ttk.Entry(_slowrow, width=6, textvariable=slow_cap_var).pack(side="left", padx=(2, 6))
ttk.Label(_slowrow, text="(bad-wire workaround \u2014 runs travels at jog speed)",
          foreground="#888").pack(side="left")
ttk.Checkbutton(run_frame, text="Hard stop on abort \u2014 M112 (stops buffered moves, needs reconnect after)",
                variable=hard_stop_var).pack(anchor="w", pady=(0, 6))
ttk.Checkbutton(run_frame, text="Present bed forward when done (moves bed -Y; off = stay at home)",
                variable=present_var).pack(anchor="w", pady=(0, 6))
ttk.Checkbutton(run_frame, text="Raise heads for bed removal",
                variable=raise_heads_var).pack(anchor="w", pady=(0, 6))
_xdrow = ttk.Frame(run_frame); _xdrow.pack(fill="x", pady=(0, 6))
ttk.Label(_xdrow, text="Crosslink dwell after spray (s)").pack(side="left")
ttk.Entry(_xdrow, width=6, textvariable=xlink_dwell_var).pack(side="left", padx=(6, 0))
ttk.Checkbutton(run_frame, text="Print in one go (no spray pause / present between layers)",
                variable=continuous_var).pack(anchor="w", pady=(0, 6))
_solrow = ttk.Frame(run_frame); _solrow.pack(fill="x", pady=(0, 6))
_solchk = ttk.Checkbutton(_solrow, text="Solenoid valve gates plunger (close snaps the drop off)",
                variable=solenoid_var); _solchk.pack(side="left")
_solrow2 = ttk.Frame(run_frame); _solrow2.pack(fill="x", pady=(0, 6))
ttk.Label(_solrow2, text="snap settle (ms)", width=16).pack(side="left")
ttk.Entry(_solrow2, width=6, textvariable=solenoid_snap_var).pack(side="left")
ttk.Label(_solrow2, text="wait after push before the valve closes (0 = close immediately)",
          foreground="#8FA5B8", font=("Segoe UI", 8)).pack(side="left", padx=(6, 0))
_svtest = ttk.Frame(run_frame); _svtest.pack(fill="x", pady=(0, 6))
ttk.Label(_svtest, text="Valve test:").pack(side="left")
ttk.Button(_svtest, text="Open", width=6, command=lambda: send_valve_command("O")).pack(side="left", padx=2)
ttk.Button(_svtest, text="Pulse", width=6,
           command=lambda: send_valve_command(f"F{solenoid_us_var.get()}")).pack(side="left", padx=2)
ttk.Button(_svtest, text="Close", width=6, command=lambda: send_valve_command("C")).pack(side="left", padx=2)
ttk.Button(run_frame, text="\u25b6 PRINT (layer by layer)", command=start_layered_print, style="Accent.TButton").pack(fill="x", ipady=10, pady=3)
pause_btn = ttk.Button(run_frame, text="\u23f8 PAUSE", command=toggle_pause)
pause_btn.pack(fill="x", ipady=6, pady=3)
continue_btn = ttk.Button(run_frame, text="\u23ed CONTINUE (after spraying)", command=continue_layer, state="disabled")
continue_btn.pack(fill="x", ipady=8, pady=3)
ttk.Button(run_frame, text="\u2139 Estimate time / drops", command=estimate_print).pack(fill="x", ipady=4, pady=3)
try:
    ttk.Style().configure("Estop.TButton", foreground="#FF6B6B", font=("Segoe UI", 10, "bold"))
except Exception:
    pass
ttk.Button(run_frame, text="\u23f9 EMERGENCY STOP", command=stop_scaffold_print, style="Estop.TButton").pack(fill="x", ipady=8, pady=(10, 3))

scaffold_status_var = tk.StringVar(value="Status: waiting for G-code...")
ttk.Label(run_frame, textvariable=scaffold_status_var, font=("Segoe UI", 10, "bold")).pack(pady=8)
scaffold_progress = ttk.Progressbar(run_frame, orient="horizontal", mode="determinate", value=0)
scaffold_progress.pack(fill="x")

# preview
scaf_right = ttk.LabelFrame(tab_scaffold, text="G-Code Preview", padding=10)
scaf_right.pack(side=tk.RIGHT, fill="both", expand=True, padx=20, pady=15)
gcode_scroll = ttk.Scrollbar(scaf_right)
gcode_scroll.pack(side=tk.RIGHT, fill="y")
gcode_listbox = tk.Listbox(scaf_right, yscrollcommand=gcode_scroll.set, bg="#2E2E2E", fg="#FFFFFF", font=("Consolas", 10))
gcode_listbox.pack(side=tk.LEFT, fill="both", expand=True)
gcode_scroll.config(command=gcode_listbox.yview)

# voron control tab (start of ricky voron commit)

# ==========================================
# --- VORON / SLA CONTROL TAB ---
# ==========================================

# ----- variables -----

voron_host_var = tk.StringVar(value="voron.local")
voron_connection_var = tk.StringVar(
    value="Moonraker: Not connected"
)
voron_state_var = tk.StringVar(
    value="State: unknown"
)

voron_z_var = tk.StringVar(value="—")
voron_homed_var = tk.StringVar(value="Not homed")
voron_job_var = tk.StringVar(value="—")

voron_progress_var = tk.DoubleVar(value=0.0)
voron_step_var = tk.DoubleVar(value=1.0)
voron_file_var = tk.StringVar(value="")


# ----- main wrapper -----

voron_wrapper = ttk.Frame(tab_voron)
voron_wrapper.pack(
    fill="both",
    expand=True,
    padx=20,
    pady=15
)


# =========================================================
# CONNECTION
# =========================================================

voron_conn = ttk.LabelFrame(
    voron_wrapper,
    text="1. Voron / SLA Connection",
    padding=15
)

voron_conn.pack(fill="x", pady=(0, 10))

ttk.Label(
    voron_conn,
    text="Host / IP:"
).grid(
    row=0,
    column=0,
    padx=5,
    pady=5,
    sticky="e"
)

voron_host_entry = ttk.Entry(
    voron_conn,
    textvariable=voron_host_var,
    width=32
)

voron_host_entry.grid(
    row=0,
    column=1,
    padx=5,
    pady=5,
    sticky="we"
)

voron_connect_btn = ttk.Button(
    voron_conn,
    text="CONNECT",
    style="Accent.TButton",
    command=voron_connect
)

voron_connect_btn.grid(
    row=0,
    column=2,
    padx=10,
    pady=5
)

ttk.Label(
    voron_conn,
    textvariable=voron_connection_var,
    font=("Segoe UI", 10, "bold"),
    foreground="#A8E6CF"
).grid(
    row=0,
    column=3,
    padx=15,
    pady=5
)

ttk.Label(
    voron_conn,
    textvariable=voron_state_var,
    font=("Segoe UI", 10, "bold")
).grid(
    row=1,
    column=1,
    padx=5,
    sticky="w"
)

ttk.Label(
    voron_conn,
    text=(
        "Active target: Z-axis SLA build-plate motion. "
        "Vat rotation is manual."
    ),
    foreground="#8FA5B8",
    font=("Segoe UI", 8)
).grid(
    row=1,
    column=2,
    columnspan=2,
    padx=10,
    sticky="w"
)

voron_conn.columnconfigure(1, weight=1)


# =========================================================
# LOWER TWO-COLUMN AREA
# =========================================================

voron_body = ttk.Frame(voron_wrapper)
voron_body.pack(fill="both", expand=True)

voron_left = ttk.Frame(voron_body)
voron_left.pack(
    side=tk.LEFT,
    fill="both",
    expand=True,
    padx=(0, 8)
)

voron_right = ttk.Frame(voron_body)
voron_right.pack(
    side=tk.RIGHT,
    fill="both",
    expand=True,
    padx=(8, 0)
)


# =========================================================
# LEFT — BUILD PLATE / Z CONTROL
# =========================================================

voron_motion = ttk.LabelFrame(
    voron_left,
    text="2. Build Plate / Z Control",
    padding=15
)

voron_motion.pack(fill="x", pady=(0, 10))


# position

zpos = ttk.Frame(voron_motion)
zpos.pack(fill="x", pady=(0, 8))

ttk.Label(
    zpos,
    text="Current Z:",
    font=("Segoe UI", 10, "bold")
).pack(side=tk.LEFT)

ttk.Label(
    zpos,
    textvariable=voron_z_var,
    font=("Consolas", 12, "bold"),
    foreground="#A8E6CF"
).pack(side=tk.LEFT, padx=(8, 20))

ttk.Label(
    zpos,
    text="Homed:"
).pack(side=tk.LEFT)

ttk.Label(
    zpos,
    textvariable=voron_homed_var,
    font=("Segoe UI", 10, "bold")
).pack(side=tk.LEFT, padx=6)


# step size

step_row = ttk.Frame(voron_motion)
step_row.pack(fill="x", pady=8)

ttk.Label(
    step_row,
    text="Step size:",
    font=("Segoe UI", 10, "bold")
).pack(side=tk.LEFT, padx=(0, 8))

for value, label in (
    (0.1, "0.1 mm"),
    (1.0, "1 mm"),
    (10.0, "10 mm"),
):
    ttk.Radiobutton(
        step_row,
        text=label,
        variable=voron_step_var,
        value=value
    ).pack(side=tk.LEFT, padx=6)


# Z movement

z_btns = ttk.Frame(voron_motion)
z_btns.pack(pady=12)

ttk.Button(
    z_btns,
    text="▲  Z +",
    width=16,
    command=lambda: voron_jog_z(+1)
).grid(
    row=0,
    column=0,
    padx=5,
    pady=4
)

ttk.Button(
    z_btns,
    text="▼  Z −",
    width=16,
    command=lambda: voron_jog_z(-1)
).grid(
    row=1,
    column=0,
    padx=5,
    pady=4
)

ttk.Button(
    z_btns,
    text="HOME Z",
    width=16,
    style="Accent.TButton",
    command=voron_home_z
).grid(
    row=2,
    column=0,
    padx=5,
    pady=(10, 4)
)

ttk.Label(
    voron_motion,
    text=(
        "Only the configured Z build-plate axis is exposed here. "
        "No X/Y motion is assumed."
    ),
    foreground="#888",
    font=("Segoe UI", 8),
    wraplength=400,
    justify="left"
).pack(anchor="w", pady=(8, 0))


# =========================================================
# LEFT — G-CODE / JOB
# =========================================================

voron_job_frame = ttk.LabelFrame(
    voron_left,
    text="3. G-Code / Job",
    padding=15
)

voron_job_frame.pack(fill="x", pady=(0, 10))

ttk.Label(
    voron_job_frame,
    text="Selected / uploaded file:"
).pack(anchor="w")

voron_file_entry = ttk.Entry(
    voron_job_frame,
    textvariable=voron_file_var
)

voron_file_entry.pack(
    fill="x",
    pady=(4, 8)
)

job_buttons = ttk.Frame(voron_job_frame)
job_buttons.pack(fill="x")

ttk.Button(
    job_buttons,
    text="CHOOSE FILE",
    command=voron_choose_gcode
).pack(
    side=tk.LEFT,
    fill="x",
    expand=True,
    padx=(0, 4)
)

voron_upload_btn = ttk.Button(
    job_buttons,
    text="UPLOAD",
    style="Accent.TButton",
    command=voron_upload_gcode
)

voron_upload_btn.pack(
    side=tk.LEFT,
    fill="x",
    expand=True,
    padx=(4, 0)
)

ttk.Label(
    voron_job_frame,
    text=(
        "Upload stages the file on the Voron. "
        "Starting the print is a separate action."
    ),
    foreground="#888",
    font=("Segoe UI", 8),
    wraplength=400
).pack(anchor="w", pady=(8, 0))


# =========================================================
# RIGHT — PRINTER STATUS
# =========================================================

voron_status = ttk.LabelFrame(
    voron_right,
    text="4. Printer Status",
    padding=15
)

voron_status.pack(fill="x", pady=(0, 10))

status_grid = ttk.Frame(voron_status)
status_grid.pack(fill="x")

ttk.Label(
    status_grid,
    text="Firmware:"
).grid(
    row=0,
    column=0,
    sticky="w",
    pady=3
)

ttk.Label(
    status_grid,
    textvariable=voron_state_var,
    font=("Segoe UI", 10, "bold")
).grid(
    row=0,
    column=1,
    sticky="w",
    padx=10
)

ttk.Label(
    status_grid,
    text="Homed axes:"
).grid(
    row=1,
    column=0,
    sticky="w",
    pady=3
)

ttk.Label(
    status_grid,
    textvariable=voron_homed_var
).grid(
    row=1,
    column=1,
    sticky="w",
    padx=10
)

ttk.Label(
    status_grid,
    text="Current job:"
).grid(
    row=2,
    column=0,
    sticky="w",
    pady=3
)

ttk.Label(
    status_grid,
    textvariable=voron_job_var
).grid(
    row=2,
    column=1,
    sticky="w",
    padx=10
)

ttk.Label(
    status_grid,
    text="Vat:"
).grid(
    row=3,
    column=0,
    sticky="w",
    pady=3
)

ttk.Label(
    status_grid,
    text="Manual rotation"
).grid(
    row=3,
    column=1,
    sticky="w",
    padx=10
)

ttk.Label(
    voron_status,
    text="Print progress:"
).pack(anchor="w", pady=(12, 3))

voron_progress = ttk.Progressbar(
    voron_status,
    variable=voron_progress_var,
    maximum=100.0
)

voron_progress.pack(fill="x")

voron_progress_label = ttk.Label(
    voron_status,
    text="0.0%",
    font=("Segoe UI", 9, "bold")
)

voron_progress_label.pack(anchor="e", pady=(2, 0))


# =========================================================
# RIGHT — PRINT CONTROL
# =========================================================

voron_print_controls = ttk.LabelFrame(
    voron_right,
    text="5. Print Control",
    padding=15
)

voron_print_controls.pack(fill="x", pady=(0, 10))

ttk.Button(
    voron_print_controls,
    text="▶ START PRINT",
    style="Accent.TButton",
    command=voron_start_print
).pack(
    fill="x",
    ipady=8,
    pady=3
)

pause_resume = ttk.Frame(voron_print_controls)
pause_resume.pack(fill="x", pady=3)

ttk.Button(
    pause_resume,
    text="⏸ PAUSE",
    command=voron_pause_print
).pack(
    side=tk.LEFT,
    fill="x",
    expand=True,
    padx=(0, 3)
)

ttk.Button(
    pause_resume,
    text="▶ RESUME",
    command=voron_resume_print
).pack(
    side=tk.LEFT,
    fill="x",
    expand=True,
    padx=(3, 0)
)

ttk.Button(
    voron_print_controls,
    text="CANCEL PRINT",
    command=voron_cancel_print
).pack(
    fill="x",
    ipady=5,
    pady=3
)

ttk.Button(
    voron_print_controls,
    text="⏹ EMERGENCY STOP",
    style="Estop.TButton",
    command=voron_emergency_stop
).pack(
    fill="x",
    ipady=8,
    pady=(14, 3)
)


# =========================================================
# RIGHT — CONSOLE
# =========================================================

voron_console = ttk.LabelFrame(
    voron_right,
    text="6. Voron / Moonraker Messages",
    padding=10
)

voron_console.pack(
    fill="both",
    expand=True
)

voron_log = tk.Text(
    voron_console,
    height=12,
    wrap="word",
    font=("Consolas", 9),
    bg="#12161c",
    fg="#A8E6CF",
    insertbackground="#A8E6CF",
    state="disabled"
)

voron_log_scroll = ttk.Scrollbar(
    voron_console,
    command=voron_log.yview
)

voron_log.configure(
    yscrollcommand=voron_log_scroll.set
)

voron_log.pack(
    side=tk.LEFT,
    fill="both",
    expand=True
)

voron_log_scroll.pack(
    side=tk.RIGHT,
    fill="y"
)

_voron_log(
    "Voron SLA control ready. Enter printer host/IP and connect."
)

#end of ricky voron ui commit ^

sv_ttk.set_theme("dark") 
try: scaffold_status_var.set(f"Status: ready — {BUILD_VERSION}")
except Exception: pass
try: _apply_pneumatic()            # pneumatic on by default -> grey out stepper-droplet controls now
except Exception: pass
root.after(1500, query_position)   # begin live position polling
root.mainloop()