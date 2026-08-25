#!/usr/bin/env python3
"""Draw residue-number labels onto a finished pestifer-snapshot image.

VMD's own text does not survive ray tracing -- `graphics text` is a GL stroke font, and
TachyonInternal drops it silently, producing an image with no labels and no warning.  So the
Tcl side emits a screen position per highlighted residue (`-labels 1`) and the labels are
composited here, which also buys real typography and a halo that reads over both the pale
ribbon and the dark sticks underneath.

Reads "label <resid> <resname> <x> <y> <colorid>" lines on stdin; edits the image in place.

An optional second argument nudges individual labels, for the one or two that automatic
placement puts somewhere awkward: "15:1.2,1.2;52:0,-0.5", in units of the label's own type size,
x positive right and y positive down.  Nudges are applied after the labels have been spread
apart, so a nudge is not undone by the spreading.
"""
import sys, re, os
from PIL import Image, ImageDraw, ImageFont

# VMD ColorID -> RGB, for the ids pestifer-snapshot uses as highlight colors
VMD_COLORS = {1: (255, 0, 0), 3: (255, 128, 0), 4: (255, 255, 0), 5: (210, 155, 100),
              7: (0, 210, 0), 9: (255, 130, 165), 10: (80, 200, 200)}
FONTS = ["/usr/share/fonts/truetype/DejaVuSans-Bold.ttf",
         "/usr/share/fonts/texlive-dejavu/DejaVuSans-Bold.ttf",
         "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf"]


def _font(px):
    for f in FONTS:
        if os.path.exists(f):
            return ImageFont.truetype(f, px)
    return ImageFont.load_default()


def main():
    if len(sys.argv) < 2:
        sys.exit("usage: pestifer_label_overlay.py IMAGE < vmd-output")
    path = sys.argv[1]
    nudges = {}
    if len(sys.argv) > 2 and sys.argv[2].strip():
        for part in sys.argv[2].split(";"):
            part = part.strip()
            if not part:
                continue
            try:
                who, delta = part.split(":")
                dx, dy = delta.split(",")
                nudges[who.strip()] = (float(dx), float(dy))
            except ValueError:
                print(f"pestifer_snapshot: ignoring bad -labelnudge term '{part}'",
                      file=sys.stderr)
    labels = []
    for line in sys.stdin:
        m = re.search(r'label (\S+) (\S+) ([-\d.]+) ([-\d.]+) (\d+)', line)
        if m:
            labels.append((m.group(1), m.group(2), float(m.group(3)), float(m.group(4)),
                           int(m.group(5))))
    if not labels:
        return
    im = Image.open(path).convert("RGB")
    W, H = im.size
    d = ImageDraw.Draw(im)
    size = max(18, int(H * 0.030))
    font = _font(size)
    # push each label away from the frame centre so it clears the sticks it names
    off = H * 0.038
    placed = []
    for resid, _resname, x, y, cid in labels:
        vx, vy = x - W / 2.0, y - H / 2.0
        n = (vx * vx + vy * vy) ** 0.5 or 1.0
        placed.append([resid, cid, x + off * vx / n, y + off * vy / n, x, y])

    # Residues that sit close together in the structure get pushed to nearly the same spot, and
    # two numbers on top of each other are worse than no numbers at all.  Nudge overlapping
    # labels apart, then let each keep a leader back to the atom it names.
    # Measure the real glyph box per label: a fixed width treats "5" and "132" alike, which
    # leaves two adjacent numbers just far enough apart not to overlap and close enough to read
    # as one number -- "13" beside "55" becomes "1355".
    gap = size * 0.55
    tall = size * 1.15
    for lab in placed:
        bb = font.getbbox(lab[0])
        lab.append(bb[2] - bb[0])          # text width
    for _ in range(200):
        moved = False
        for i in range(len(placed)):
            for j in range(i + 1, len(placed)):
                dx = placed[j][2] - placed[i][2]
                dy = placed[j][3] - placed[i][3]
                need = (placed[i][6] + placed[j][6]) / 2.0 + gap
                if abs(dx) < need and abs(dy) < tall:
                    dist = (dx * dx + dy * dy) ** 0.5 or 0.01
                    push = (need - abs(dx)) * 0.25 + 1.0
                    ux, uy = dx / dist, dy / dist
                    placed[i][2] -= ux * push; placed[i][3] -= uy * push
                    placed[j][2] += ux * push; placed[j][3] += uy * push
                    moved = True
        if not moved:
            break

    for lab in placed:
        if lab[0] in nudges:
            dx, dy = nudges[lab[0]]
            lab[2] += dx * size
            lab[3] += dy * size

    for resid, cid, tx, ty, ax, ay, _w in placed:
        col = VMD_COLORS.get(cid, (0, 0, 0))
        # leader line when the label had to travel far enough that the pairing is not obvious
        if ((tx - ax) ** 2 + (ty - ay) ** 2) ** 0.5 > off * 1.6:
            d.line([(ax, ay), (tx, ty)], fill=(255, 255, 255), width=max(4, size // 7))
            d.line([(ax, ay), (tx, ty)], fill=col, width=max(2, size // 14))
        halo = max(2, size // 9)
        # white halo first so the number reads over dark ribbon as well as pale background
        for dx in range(-halo, halo + 1):
            for dy in range(-halo, halo + 1):
                if dx * dx + dy * dy <= halo * halo:
                    d.text((tx + dx, ty + dy), resid, font=font, fill=(255, 255, 255),
                           anchor="mm")
        d.text((tx, ty), resid, font=font, fill=col, anchor="mm")
    im.save(path)
    print(f"pestifer_snapshot: overlaid {len(labels)} label(s) on {path}")


if __name__ == "__main__":
    main()
