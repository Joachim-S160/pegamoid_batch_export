# Tutorial 3 — Make a video from orbital montages

This tutorial shows how to loop `pegamoid_batch_export.py` over multiple HDF5 files
(e.g. one per geometry along a dissociation curve) and assemble the resulting montage
PNGs into an MP4 video.

## Additional requirements

Beyond the base `requirements.txt`, the video script needs:

```bash
pip install opencv-python
```

`Pillow` is already included in `requirements.txt` and is used here for frame labelling.

## Overview

The workflow is:

1. Loop over all `.rasscf.h5` files (one per geometry / run directory)
2. Run `pegamoid_batch_export.py` on each → produces a montage PNG
3. Load each montage PNG, add a text label, convert to a video frame
4. Write all frames to an MP4 with OpenCV

## Step 1 — Export montages for each geometry

Run `pegamoid_batch_export.py` in a bash loop. Adjust the glob pattern and output
directory naming to match your directory structure.

```bash
# Example: Mertens CAS(8,7) runs, active space + 10 neighbours
for h5 in path/to/CAS_8_7/run_*/PbO_cas87_mertens/PbO_cas87_mertens.rasscf.h5; do
    run=$(echo "$h5" | grep -oP 'run_\K[0-9]+')
    python pegamoid_batch_export.py "$h5" \
        --active \
        --neighbors 10 \
        --camera 45,30 \
        --output-dir "./montages/run${run}"
done
```

For a specific orbital range instead of the active space, replace `--active --neighbors 10`
with `--orbitals 19,20,21,...,31`.

## Step 2 — Assemble into an MP4

Save the following as `make_video.py` and run it after Step 1.

```python
import cv2
import numpy as np
from pathlib import Path
from PIL import Image, ImageDraw

# ── Configuration ──────────────────────────────────────────────────────────────
MONTAGE_DIRS = sorted(
    Path("./montages").glob("run*"),
    key=lambda p: int(p.name.replace("run", ""))
)
MONTAGE_GLOB = "*_montage.png"   # filename pattern produced by pegamoid_batch_export.py
OUTPUT_MP4   = "orbital_scan.mp4"
FPS          = 2                 # frames per second (slow enough to inspect each frame)
LABEL_FMT    = "run_{run:02d}"   # text overlaid on each frame
# ───────────────────────────────────────────────────────────────────────────────

frames = []
for run_dir in MONTAGE_DIRS:
    pngs = list(run_dir.glob(MONTAGE_GLOB))
    if not pngs:
        print(f"WARNING: no montage found in {run_dir}")
        continue
    png = pngs[0]
    run_num = int(run_dir.name.replace("run", ""))

    img = Image.open(png).convert("RGB")
    draw = ImageDraw.Draw(img)
    label = LABEL_FMT.format(run=run_num)

    # Black bar + yellow label in the top-left corner
    draw.rectangle([0, 0, 300, 28], fill=(0, 0, 0))
    draw.text((8, 6), label, fill=(255, 255, 0))

    frames.append(np.array(img))

if not frames:
    raise RuntimeError("No frames found — check MONTAGE_DIRS and MONTAGE_GLOB.")

h, w = frames[0].shape[:2]
writer = cv2.VideoWriter(
    OUTPUT_MP4,
    cv2.VideoWriter_fourcc(*"mp4v"),
    FPS,
    (w, h),
)
for frame in frames:
    writer.write(cv2.cvtColor(frame, cv2.COLOR_RGB2BGR))
writer.release()

print(f"Written: {OUTPUT_MP4}  ({len(frames)} frames, {w}x{h}, {FPS} fps)")
```

Run it:

```bash
python make_video.py
```

## Customising the label

Replace `LABEL_FMT` to show the bond length, system number, or any other per-frame
annotation. You can retrieve it from the directory name or a separate list:

```python
DISTANCES = [1.55, 1.65, 1.75, 1.80, 1.85, 1.90, 1.95, 2.00, 2.10, ...]

for i, run_dir in enumerate(MONTAGE_DIRS):
    ...
    label = f"R = {DISTANCES[i]:.2f} Å"
```

## Adjusting playback speed

Change `FPS` in `make_video.py`:
- `FPS = 1` — one frame per second (useful for detailed inspection)
- `FPS = 2` — default; comfortable pace for ~27-frame dissociation curves
- `FPS = 5` — fast scan for long series (50+ geometries)

## Full working example

```bash
# 1. Export montages for all 27 Mertens CAS(8,7) geometries
for h5 in CAS_8_7/run_*/PbO_cas87_mertens/PbO_cas87_mertens.rasscf.h5; do
    run=$(echo "$h5" | grep -oP 'run_\K[0-9]+')
    python pegamoid_batch_export.py "$h5" \
        --active --neighbors 10 --camera 45,30 \
        --output-dir "montages/run${run}"
done

# 2. Assemble video
python make_video.py
# → orbital_scan.mp4  (27 frames, 2 fps)
```

Open the resulting MP4 in any video player. Step through frame by frame to compare
orbital character across geometries and identify problematic geometries where d/f or
Rydberg orbitals have contaminated the active space.
