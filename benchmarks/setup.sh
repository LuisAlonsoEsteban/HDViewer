#!/usr/bin/env bash
# Create the benchmark environments (Linux x86_64). Re-running is safe: finished steps are skipped.
#
#   ENVS=/scratch/me/envs bash benchmarks/setup.sh
#
# Needs: python3.11 (PY), conda/mamba with conda-forge (MAMBA), and node/npm (NODE_ENV conda env).
# Paths written here must match benchmarks/config.yaml (python:, env:, xvfb:).
set -euo pipefail

REPO="$(cd "$(dirname "$0")/.." && pwd)"
ENVS="${ENVS:-/scratch/lalonsoeste/PhD/envs}"
PY="${PY:-python3.11}"
MAMBA="${MAMBA:-mamba}"
NODE_ENV="${NODE_ENV:-hdv-node}"         # conda env providing node and npm
export PLAYWRIGHT_BROWSERS_PATH="${PLAYWRIGHT_BROWSERS_PATH:-$ENVS/playwright-browsers}"

# 1. HDViewer itself: the project's .venv (see the main README)
[ -x "$REPO/.venv/bin/python" ] || { echo "Create $REPO/.venv first (README: Install)"; exit 1; }

# 2. SpatialData + Vitessce + Playwright driver
if [ ! -x "$ENVS/bench-sd/bin/python" ]; then
  "$PY" -m venv "$ENVS/bench-sd"
  "$ENVS/bench-sd/bin/pip" install -q --upgrade pip
  "$ENVS/bench-sd/bin/pip" install -r "$REPO/benchmarks/envs/bench-sd.txt"
fi
"$ENVS/bench-sd/bin/python" -m playwright install chromium

# 3. napari (PySide6/Qt 6)
if [ ! -x "$ENVS/bench-napari/bin/python" ]; then
  "$PY" -m venv "$ENVS/bench-napari"
  "$ENVS/bench-napari/bin/pip" install -q --upgrade pip
  "$ENVS/bench-napari/bin/pip" install -r "$REPO/benchmarks/envs/bench-napari.txt"
fi

# 4. An X server with OpenGL for napari: conda-forge Xvfb, the system's Mesa DRI drivers,
#    and the xcb libraries Qt's xcb plugin needs (linked alone, so conda's libGL does not
#    shadow the system's).
X="$ENVS/bench-xvfb"
if [ ! -x "$X/bin/Xvfb" ]; then
  "$MAMBA" create -y -q -p "$X" -c conda-forge --override-channels \
    xorg-xvfb-server libxkbcommon xcb-util-cursor xcb-util-wm xcb-util-image xcb-util-keysyms xcb-util-renderutil
fi
mkdir -p "$X/qtlibs"
for f in "$X"/lib/libxkbcommon*.so* "$X"/lib/libxcb-{cursor,icccm,util,image,keysyms,render-util}.so*; do
  ln -sf "$f" "$X/qtlibs/"
done

# 5. The page hosting Vitessce (React + vitessce from npm, built once)
cd "$REPO/benchmarks/adapters/vitessce_app"
if [ ! -f dist/index.html ]; then
  conda run -n "$NODE_ENV" npm install --no-audit --no-fund
  conda run -n "$NODE_ENV" npm run build
fi

# 6. Smoke check: napari draws through Xvfb
cd "$REPO"
LIBGL_DRIVERS_PATH=/usr/lib64/dri "$X/bin/Xvfb" :99 -screen 0 1280x1024x24 +extension GLX -nolisten tcp \
  >/dev/null 2>&1 & XPID=$!
sleep 2
DISPLAY=:99 QT_API=pyside6 QT_QPA_PLATFORM=xcb LD_LIBRARY_PATH="$X/qtlibs" "$ENVS/bench-napari/bin/python" - <<'PY'
import numpy as np, napari
v = napari.Viewer(show=False)
v.add_image(np.random.randint(0, 255, (256, 256, 3), dtype=np.uint8), rgb=True)
shot = v.screenshot(canvas_only=True)
assert shot.mean() > 0, "napari rendered an empty canvas"
print("napari OpenGL rendering: ok")
PY
kill $XPID
echo "Environments ready. Generate data and run: .venv/bin/python benchmarks/run.py --list"
