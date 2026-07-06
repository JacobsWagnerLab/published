# -*- coding: utf-8 -*-
"""
Export Imaris 3D-rendered biofilm images as grouped montage videos (PNG / MP4).

Two layout modes, selected per position via CHANNEL_MODE:

  "sytox"  (xy01–xy14): GFP + SYTOX Orange
      Panels: GFP side view | SYTOX side view | GFP/SYTOX overlay top view
      Folder suffixes expected: _gfp, _sytox, _comp_top

  "ebba"   (xy15–xy16): GFP + EbbaBiolight680
      Panels: GFP side view | EbbaBiolight680 side view | GFP top view | EbbaBiolight680 top view
      Folder suffixes expected: _gfp, _ebba_glow, _gfp_top, _ebba_glow_top

Each position must have a subfolder per rendered channel/view inside ROOT_DIR.
Outputs go to ROOT_DIR/png_grouped and ROOT_DIR/mp4_grouped.

Requires: imageio-ffmpeg, tifffile, Pillow, matplotlib

@author: alessio fragasso
"""

import re
import shutil
import subprocess
import tempfile
from pathlib import Path

import matplotlib.cm as mcm
import matplotlib.pyplot as plt
import numpy as np
import tifffile as tiff
import imageio_ffmpeg as iio_ffmpeg
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.colors import LinearSegmentedColormap, ListedColormap, Normalize
from matplotlib.lines import Line2D
from PIL import Image


'''Settings'''
ROOT_DIR = #Path(r"path_to_imaris_export_folder")

MAKE_PNG = False
MAKE_MP4 = True


'''Condition and channel mode per position'''
# xy01–xy14: GFP + SYTOX Orange (3-panel layout)
# xy15–xy16: GFP + EbbaBiolight680 (4-panel layout)
CONDITION_MAP = {
    "xy01": "Ctrl",  "xy02": "Ctrl",
    "xy03": "LL-37", "xy04": "LL-37",
    "xy05": "CecA",  "xy06": "CecA",
    "xy07": "Indo",  "xy08": "Indo",
    "xy09": "Tur1A", "xy10": "Tur1A",
    "xy11": "Bac7",  "xy12": "Bac7",
    "xy13": "PR-39", "xy14": "PR-39",
    "xy15": "Ctrl",  "xy16": "Ctrl",
}

CHANNEL_MODE = {
    "xy01": "sytox", "xy02": "sytox",
    "xy03": "sytox", "xy04": "sytox",
    "xy05": "sytox", "xy06": "sytox",
    "xy07": "sytox", "xy08": "sytox",
    "xy09": "sytox", "xy10": "sytox",
    "xy11": "sytox", "xy12": "sytox",
    "xy13": "sytox", "xy14": "sytox",
    "xy15": "ebba",  "xy16": "ebba",
}


'''Video parameters'''
FPS            = 8
REPEAT_FACTOR  = 5
MP4_CRF        = 18
MP4_PRESET     = "slow"

# per-mode output width (EBBA 4-panel needs more pixels)
MP4_SCALE_WIDTH = {
    "sytox": 1920,
    "ebba":  3000,
}

MINUTES_PER_FRAME  = 30
AMP_START_FRAME_NUM = 3   # frame number (1-indexed) from which AMP label appears


'''Crop parameters'''
TOP_CROP_SIDE_PX = 250   # px to remove from top of side-view images
CROP_LR_SIDE_PX  = 180   # px to remove from each side of side-view images
CROP_LR_TOP_PX   = 555   # px to remove from each side of top-view images


'''Typography'''
FONT_FAMILY          = "Arial"
TITLE_FONTSIZE       = 80
TIMESTAMP_FONTSIZE   = 80
SCALE_TEXT_FONTSIZE  = 80
CBAR_LABEL_FONTSIZE  = 80
CBAR_TICK_FONTSIZE   = 80

plt.rcParams["font.family"] = FONT_FAMILY


'''Layout'''
DPI                = 100
TOP_TITLE_BAND_PX  = 140
AMP_LINE_GAP_PX    = 110
TOP_MARGIN_PX      = 120
BOTTOM_MARGIN_PX   = 120
ROW_GAP_PX         = 45
COL_GAP_PX         = 34

CBAR_PAD_LEFT_PX   = 60
CBAR_W_PX          = 24
CBAR_GAP_PX        = 500
CBAR_RIGHT_MARGIN_PX = 600
CBAR_HEIGHT_FRAC   = 0.56
SIDEBAR_PX = CBAR_PAD_LEFT_PX + 2 * CBAR_W_PX + CBAR_GAP_PX + CBAR_RIGHT_MARGIN_PX


'''Scale bar'''
SCALE_LABEL                    = "20 µm"
SCALE_LINEWIDTH                = 8
SCALE_TEXT_OFFSET_PX           = 10
SCALE_MARGIN_X_PX              = 22
SCALE_CLEAR_PAD_PX             = 5
SCALE_BAR_Y_FROM_IMAGE_BOTTOM_PX = -65


'''Colormaps and intensity ranges'''
CMAP_GFP   = LinearSegmentedColormap.from_list("imaris_green",   ["black", "#00FF00"])
CMAP_SYTOX = LinearSegmentedColormap.from_list("imaris_magenta", ["black", "#FF00FF"])

GFP_RANGE   = (500,  5000)
SYTOX_RANGE = (2000, 20000)
GLOW_RANGE  = (250,  2500)


def make_glow_cmap():
    """
    Imaris/Fiji-style Glow LUT.
    Loads glow.lut from the script directory if present; otherwise uses a
    built-in approximation.
    """
    script_dir = Path(__file__).resolve().parent if "__file__" in globals() else Path.cwd()
    lut_path = script_dir / "glow.lut"

    if lut_path.exists():
        rgb = np.loadtxt(lut_path, skiprows=1, usecols=(1, 2, 3)) / 255.0
        if rgb.shape[0] >= 256:
            rgb = rgb[1:-1]
        return ListedColormap(rgb, name="imaris_fiji_glow")

    glow_rgb_255 = [
        (0,   0,   0),   (21,  1,   1),   (46,  2,   4),   (80,  6,   6),
        (110, 13,  6),   (141, 21,  10),  (168, 33,  9),   (192, 48,  8),
        (206, 70,  8),   (219, 96,  7),   (227, 123, 12),  (234, 148, 23),
        (238, 172, 41),  (246, 192, 67),  (249, 213, 99),  (249, 228, 132),
        (253, 241, 168), (253, 249, 203), (253, 252, 238), (253, 253, 252),
    ]

    return LinearSegmentedColormap.from_list(
        "imaris_fiji_glow_approx",
        [(r / 255, g / 255, b / 255) for r, g, b in glow_rgb_255],
        N=256,
    )


CMAP_GLOW = make_glow_cmap()


'''Image helpers'''
def tif_to_rgb(frame: np.ndarray) -> np.ndarray:
    frame = np.squeeze(frame)
    if frame.ndim == 3 and frame.shape[-1] == 3:
        pass
    elif frame.ndim == 2:
        frame = np.stack([frame] * 3, axis=-1)
    else:
        raise ValueError(f"Unexpected frame shape: {frame.shape}")
    if frame.dtype != np.uint8:
        vmin, vmax = frame.min(), frame.max()
        if vmax > vmin:
            frame = ((frame - vmin) / (vmax - vmin) * 255).astype(np.uint8)
        else:
            frame = np.zeros_like(frame, dtype=np.uint8)
    return frame


def rgb_to_gray(img: np.ndarray) -> np.ndarray:
    return (0.2126 * img[..., 0] + 0.7152 * img[..., 1] + 0.0722 * img[..., 2]).astype(np.float32)


def crop_top(img: np.ndarray, top_px: int) -> np.ndarray:
    if top_px <= 0:
        return img
    return img[min(top_px, img.shape[0] - 1):, :, :]


def crop_lr(img: np.ndarray, crop_px: int) -> np.ndarray:
    h, w = img.shape[:2]
    x0 = max(0, crop_px)
    x1 = max(x0 + 1, w - crop_px)
    return img[:, x0:x1, :]


def resize_to_height(img: np.ndarray, target_h: int):
    h, w = img.shape[:2]
    scale = target_h / h
    target_w = max(1, int(round(w * scale)))
    resampling = getattr(Image, "Resampling", Image).BILINEAR
    out = Image.fromarray(img).resize((target_w, target_h), resample=resampling)
    return np.asarray(out), scale


'''Scale-bar detection'''
def _find_long_runs(binary_row: np.ndarray):
    idx = np.flatnonzero(binary_row)
    if idx.size == 0:
        return []
    breaks = np.where(np.diff(idx) > 1)[0]
    starts = np.r_[0, breaks + 1]
    ends   = np.r_[breaks, len(idx) - 1]
    return [(idx[s], idx[e]) for s, e in zip(starts, ends)]


def detect_scale_bar(frame_rgb: np.ndarray):
    """Detect the white Imaris scale bar in the bottom-left corner."""
    h, w = frame_rgb.shape[:2]
    gray = rgb_to_gray(frame_rgb)
    x_max = int(0.42 * w)
    y_min = int(0.70 * h)
    roi   = gray[y_min:h, 0:x_max]
    thr   = max(235, np.percentile(roi, 99.85))
    bright = roi >= thr

    candidates = []
    for ry in range(bright.shape[0]):
        for x0, x1 in _find_long_runs(bright[ry]):
            run_len = x1 - x0 + 1
            if run_len < 20:
                continue
            thickness = sum(
                1 for off in range(1, 6)
                if (ry + off < bright.shape[0]) and bright[ry + off, x0:x1 + 1].mean() > 0.85
            ) + 1
            y_score = ry / max(1, bright.shape[0] - 1)
            score = run_len + 18 * y_score - 2.5 * x0 + 2 * thickness
            candidates.append({"score": score, "x0": int(x0), "x1": int(x1),
                                "y": int(ry), "thickness": int(thickness)})

    if not candidates:
        return None

    best = max(candidates, key=lambda d: d["score"])
    full_x0 = best["x0"]
    full_x1 = best["x1"]
    full_y0 = y_min + best["y"]
    full_y1 = full_y0 + best["thickness"] - 1

    return {
        "x0": full_x0, "x1": full_x1,
        "y0": full_y0, "y1": full_y1,
        "px_len": full_x1 - full_x0 + 1,
        "clear_x0": max(0, full_x0 - SCALE_CLEAR_PAD_PX),
        "clear_x1": min(w, full_x1 + 120),
        "clear_y0": max(0, full_y0 - 52),
        "clear_y1": min(h, full_y1 + 16),
    }


def erase_detected_scale_bar(frame_rgb: np.ndarray, sb_info):
    if sb_info is None:
        return frame_rgb
    out = frame_rgb.copy()
    out[sb_info["clear_y0"]:sb_info["clear_y1"],
        sb_info["clear_x0"]:sb_info["clear_x1"], :] = 0
    return out


'''Figure helpers'''
def add_scale_bar_fig(fig, x_left_px, y_image_bottom_px, px_len, total_w, total_h):
    x0 = x_left_px + SCALE_MARGIN_X_PX
    x1 = x0 + px_len
    y_bar  = y_image_bottom_px + SCALE_BAR_Y_FROM_IMAGE_BOTTOM_PX
    y_text = y_bar + SCALE_TEXT_OFFSET_PX + 5

    fig.add_artist(Line2D(
        [x0 / total_w, x1 / total_w], [y_bar / total_h, y_bar / total_h],
        transform=fig.transFigure, color="white",
        lw=SCALE_LINEWIDTH, solid_capstyle="butt", zorder=3000,
    ))
    fig.text(x0 / total_w, y_text / total_h, SCALE_LABEL,
             color="white", fontsize=SCALE_TEXT_FONTSIZE, fontfamily=FONT_FAMILY,
             fontweight="regular", va="bottom", ha="left", zorder=3001)


def add_fig_title(fig, x_center_px, y_center_px, text, total_w, total_h):
    fig.text(x_center_px / total_w, y_center_px / total_h, text,
             color="white", fontsize=TITLE_FONTSIZE, fontfamily=FONT_FAMILY,
             fontweight="regular", ha="center", va="center")


def style_colorbar(cbar, label: str):
    cbar.set_label(label, fontsize=CBAR_LABEL_FONTSIZE, color="white", labelpad=16)
    cbar.ax.yaxis.label.set_fontfamily(FONT_FAMILY)
    cbar.ax.yaxis.label.set_fontweight("regular")
    cbar.ax.tick_params(colors="white", labelsize=CBAR_TICK_FONTSIZE, width=2.0, length=7)
    for tick in cbar.ax.get_yticklabels():
        tick.set_fontfamily(FONT_FAMILY)
        tick.set_fontsize(CBAR_TICK_FONTSIZE)
    cbar.outline.set_edgecolor("white")
    cbar.outline.set_linewidth(1.8)
    cbar.ax.yaxis.tick_right()
    cbar.ax.yaxis.set_label_position("right")


def _add_timestamp_artists(fig, total_w, total_h, image_area_h):
    """Create and return the four text artists for timestamp + AMP label (with shadows)."""
    x_px = 30
    ts_y  = BOTTOM_MARGIN_PX + image_area_h + TOP_TITLE_BAND_PX * 1.2
    amp_y = ts_y - AMP_LINE_GAP_PX

    kw = dict(fontfamily=FONT_FAMILY, fontweight="regular", va="center", ha="left", clip_on=False)

    ts_shadow  = fig.text((x_px+3)/total_w, (ts_y-3)/total_h,  "", color="black",  fontsize=TIMESTAMP_FONTSIZE, zorder=1000, **kw)
    ts_text    = fig.text( x_px/total_w,     ts_y/total_h,      "", color="white",  fontsize=TIMESTAMP_FONTSIZE, zorder=1001, **kw)
    amp_shadow = fig.text((x_px+3)/total_w, (amp_y-3)/total_h, "", color="black",  fontsize=TIMESTAMP_FONTSIZE, zorder=1000, **kw)
    amp_text   = fig.text( x_px/total_w,     amp_y/total_h,     "", color="yellow", fontsize=TIMESTAMP_FONTSIZE, fontweight="bold", zorder=1001, **kw)

    return ts_shadow, ts_text, amp_shadow, amp_text


def _update_labels(ts_shadow, ts_text, amp_shadow, amp_text, frame_i, amp_name):
    total_mins = frame_i * MINUTES_PER_FRAME
    time_label = f"{total_mins // 60:02d}h {total_mins % 60:02d}m"
    ts_shadow.set_text(time_label)
    ts_text.set_text(time_label)

    show_amp = (frame_i + 1) >= AMP_START_FRAME_NUM and amp_name and "Ctrl" not in amp_name
    label = f"+ {amp_name}" if show_amp else ""
    amp_shadow.set_text(label)
    amp_text.set_text(label)


'''Frame preparation'''
def prepare_side_frame(path: Path, sb_side):
    img = tif_to_rgb(tiff.imread(path))
    img = crop_top(img, TOP_CROP_SIDE_PX)
    img = erase_detected_scale_bar(img, sb_side)
    img = crop_lr(img, CROP_LR_SIDE_PX)
    return img


def prepare_top_frame(path: Path, sb_top, target_h: int):
    img = tif_to_rgb(tiff.imread(path))
    img = erase_detected_scale_bar(img, sb_top)
    img = crop_lr(img, CROP_LR_TOP_PX)
    img, scale = resize_to_height(img, target_h)
    return img, scale


'''MP4 encoding'''
def encode_png_sequence_to_mp4(png_dir, out_mp4, fps, repeat_factor, scale_width):
    def extract_number(p):
        m = re.search(r"frame_(\d+(?:\.\d+)?)", p.stem, re.IGNORECASE) or \
            re.search(r"(\d+(?:\.\d+)?)", p.stem)
        return float(m.group(1)) if m else float("inf")

    files = sorted([p for p in png_dir.iterdir() if p.suffix.lower() == ".png"],
                   key=extract_number)

    if not files:
        raise RuntimeError(f"No PNG files found in:\n{png_dir}")

    out_mp4.parent.mkdir(parents=True, exist_ok=True)
    print(f"  N original PNG files: {len(files)}")

    ffmpeg_exe = iio_ffmpeg.get_ffmpeg_exe()

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        frame_i = 0
        for p in files:
            for _ in range(repeat_factor):
                shutil.copy2(p, tmpdir / f"img_{frame_i:06d}.png")
                frame_i += 1

        print(f"  N video frames: {frame_i}")

        cmd = [
            ffmpeg_exe, "-y",
            "-framerate", str(fps),
            "-i", str(tmpdir / "img_%06d.png"),
            "-vf", f"scale={scale_width}:-2,format=yuv420p",
            "-c:v", "libx264", "-profile:v", "main", "-level", "4.0",
            "-pix_fmt", "yuv420p", "-crf", str(MP4_CRF),
            "-preset", MP4_PRESET, "-movflags", "+faststart",
            str(out_mp4),
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print("STDERR:", result.stderr)
        raise RuntimeError("ffmpeg failed.")

    print(f"  Saved MP4: {out_mp4} ({out_mp4.stat().st_size / 1e6:.1f} MB)")


def _save_frames_and_encode(canvas, out_dir, tmp_png_dir, out_mp4,
                            frames_gfp, make_png, make_mp4,
                            render_fn, n_frames, mode):
    try:
        for i in range(n_frames):
            render_fn(i)
            canvas.draw()
            rendered_rgb = np.asarray(canvas.buffer_rgba())[:, :, :3]

            fname = frames_gfp[i].stem + ".png"

            if make_png:
                Image.fromarray(rendered_rgb).save(out_dir / fname)

            if make_mp4:
                Image.fromarray(rendered_rgb).save(tmp_png_dir / fname)

            if (i + 1) % 10 == 0 or i == n_frames - 1:
                print(f"    rendered {i + 1}/{n_frames}")

    finally:
        plt.close("all")

    if make_mp4 and out_mp4 is not None:
        encode_png_sequence_to_mp4(
            tmp_png_dir, out_mp4,
            fps=FPS, repeat_factor=REPEAT_FACTOR,
            scale_width=MP4_SCALE_WIDTH[mode],
        )


'''SYTOX export: 3-panel layout (xy01–xy14)'''
def export_sytox_position(
    xy_prefix, folder_gfp, folder_sytox, folder_top,
    out_dir, out_mp4, amp_name, make_png, make_mp4,
):
    frames_gfp   = sorted(folder_gfp.glob("*.tif"))
    frames_sytox = sorted(folder_sytox.glob("*.tif"))
    frames_top   = sorted(folder_top.glob("*.tif"))

    if not (frames_gfp and frames_sytox and frames_top):
        print(f"  [Skip] Missing frames for {xy_prefix}")
        return

    n_frames = min(len(frames_gfp), len(frames_sytox), len(frames_top))

    if make_png:
        out_dir.mkdir(parents=True, exist_ok=True)

    tmp_png_dir = None
    if make_mp4:
        tmp_png_dir = out_dir if make_png else out_dir.parent / f"{out_dir.name}_tmp_for_mp4"
        tmp_png_dir.mkdir(parents=True, exist_ok=True)

    print(f"  {xy_prefix} [GFP+SYTOX]: {n_frames} frames")

    # Detect scale bars
    sb_side = detect_scale_bar(crop_top(tif_to_rgb(tiff.imread(frames_gfp[0])), TOP_CROP_SIDE_PX)) \
        or detect_scale_bar(crop_top(tif_to_rgb(tiff.imread(frames_sytox[0])), TOP_CROP_SIDE_PX))
    sb_top  = detect_scale_bar(tif_to_rgb(tiff.imread(frames_top[0])))

    side_bar_px     = sb_side["px_len"] if sb_side else 80
    top_bar_px_raw  = sb_top["px_len"]  if sb_top  else side_bar_px

    first_gfp   = prepare_side_frame(frames_gfp[0],   sb_side)
    first_sytox = prepare_side_frame(frames_sytox[0], sb_side)
    side_h, side_w = first_gfp.shape[:2]

    left_h = 2 * side_h + ROW_GAP_PX

    first_top, top_scale = prepare_top_frame(frames_top[0], sb_top, left_h)
    top_h, top_w = first_top.shape[:2]
    top_bar_px = max(1, int(round(top_bar_px_raw * top_scale)))

    total_h = BOTTOM_MARGIN_PX + left_h + TOP_TITLE_BAND_PX + TOP_MARGIN_PX
    total_w = side_w + COL_GAP_PX + top_w + SIDEBAR_PX

    def fx(px): return px / total_w
    def fy(px): return px / total_h

    fig = plt.figure(figsize=(total_w / DPI, total_h / DPI), dpi=DPI)
    fig.patch.set_facecolor("black")

    ax_sytox = fig.add_axes([0,            fy(BOTTOM_MARGIN_PX),                  fx(side_w), fy(side_h)])
    ax_gfp   = fig.add_axes([0,            fy(BOTTOM_MARGIN_PX + side_h + ROW_GAP_PX), fx(side_w), fy(side_h)])
    ax_top   = fig.add_axes([fx(side_w + COL_GAP_PX), fy(BOTTOM_MARGIN_PX), fx(top_w), fy(left_h)])

    for ax in (ax_gfp, ax_sytox, ax_top):
        ax.set_facecolor("black")
        ax.axis("off")

    im_gfp   = ax_gfp.imshow(first_gfp,   interpolation="nearest", aspect="auto")
    im_sytox = ax_sytox.imshow(first_sytox, interpolation="nearest", aspect="auto")
    im_top   = ax_top.imshow(first_top,   interpolation="nearest", aspect="auto")

    title_y_top  = BOTTOM_MARGIN_PX * 0.7 + left_h + TOP_TITLE_BAND_PX * 0.1
    title_y_sytox = BOTTOM_MARGIN_PX * 0.7 + side_h + ROW_GAP_PX * 0.1

    add_fig_title(fig, side_w / 2,                   title_y_top,   "GFP side view",                   total_w, total_h)
    add_fig_title(fig, side_w / 2,                   title_y_sytox, "SYTOX Orange side view",           total_w, total_h)
    add_fig_title(fig, side_w + COL_GAP_PX + top_w/2, title_y_top,   "GFP / SYTOX Orange overlay top view", total_w, total_h)

    add_scale_bar_fig(fig, 0,                  BOTTOM_MARGIN_PX, side_bar_px, total_w, total_h)
    add_scale_bar_fig(fig, side_w + COL_GAP_PX, BOTTOM_MARGIN_PX, top_bar_px,  total_w, total_h)

    ts_shadow, ts_text, amp_shadow, amp_text = _add_timestamp_artists(fig, total_w, total_h, left_h)

    # Colorbars
    cbar_h_px  = int(round(left_h * CBAR_HEIGHT_FRAC))
    cbar_bot   = BOTTOM_MARGIN_PX + int(round((left_h - cbar_h_px) / 2))
    cax1_left  = side_w + COL_GAP_PX + top_w + CBAR_PAD_LEFT_PX
    cax2_left  = cax1_left + CBAR_W_PX + CBAR_GAP_PX

    cax_gfp   = fig.add_axes([fx(cax1_left), fy(cbar_bot), fx(CBAR_W_PX), fy(cbar_h_px)])
    cax_sytox = fig.add_axes([fx(cax2_left), fy(cbar_bot), fx(CBAR_W_PX), fy(cbar_h_px)])

    sm_gfp = mcm.ScalarMappable(cmap=CMAP_GFP,   norm=Normalize(*GFP_RANGE));   sm_gfp.set_array([])
    sm_syt = mcm.ScalarMappable(cmap=CMAP_SYTOX, norm=Normalize(*SYTOX_RANGE)); sm_syt.set_array([])
    style_colorbar(fig.colorbar(sm_gfp, cax=cax_gfp),   "GFP fluor. intensity (a.u.)")
    style_colorbar(fig.colorbar(sm_syt, cax=cax_sytox), "SYTOX Orange fluor. intensity (a.u.)")

    canvas = FigureCanvasAgg(fig)
    canvas.draw()

    def render_fn(i):
        im_gfp.set_data(prepare_side_frame(frames_gfp[i], sb_side))
        im_sytox.set_data(prepare_side_frame(frames_sytox[i], sb_side))
        im_top.set_data(prepare_top_frame(frames_top[i], sb_top, left_h)[0])
        _update_labels(ts_shadow, ts_text, amp_shadow, amp_text, i, amp_name)

    _save_frames_and_encode(canvas, out_dir, tmp_png_dir, out_mp4,
                            frames_gfp, make_png, make_mp4, render_fn, n_frames, mode="sytox")


'''EBBA export: 4-panel layout (xy15–xy16)'''
def export_ebba_position(
    xy_prefix, folder_gfp, folder_ebba, folder_gfp_top, folder_ebba_top,
    out_dir, out_mp4, amp_name, make_png, make_mp4,
):
    frames_gfp      = sorted(folder_gfp.glob("*.tif"))
    frames_ebba     = sorted(folder_ebba.glob("*.tif"))
    frames_gfp_top  = sorted(folder_gfp_top.glob("*.tif"))
    frames_ebba_top = sorted(folder_ebba_top.glob("*.tif"))

    if not (frames_gfp and frames_ebba and frames_gfp_top and frames_ebba_top):
        print(f"  [Skip] Missing frames for {xy_prefix}")
        return

    n_frames = min(len(frames_gfp), len(frames_ebba), len(frames_gfp_top), len(frames_ebba_top))

    if make_png:
        out_dir.mkdir(parents=True, exist_ok=True)

    tmp_png_dir = None
    if make_mp4:
        tmp_png_dir = out_dir if make_png else out_dir.parent / f"{out_dir.name}_tmp_for_mp4"
        tmp_png_dir.mkdir(parents=True, exist_ok=True)

    print(f"  {xy_prefix} [GFP+EBBA]: {n_frames} frames")

    # Detect scale bars
    sb_side     = detect_scale_bar(crop_top(tif_to_rgb(tiff.imread(frames_gfp[0])), TOP_CROP_SIDE_PX)) \
        or detect_scale_bar(crop_top(tif_to_rgb(tiff.imread(frames_ebba[0])), TOP_CROP_SIDE_PX))
    sb_gfp_top  = detect_scale_bar(tif_to_rgb(tiff.imread(frames_gfp_top[0])))
    sb_ebba_top = detect_scale_bar(tif_to_rgb(tiff.imread(frames_ebba_top[0])))

    side_bar_px     = sb_side["px_len"]     if sb_side     else 80
    gfp_top_bar_raw = sb_gfp_top["px_len"]  if sb_gfp_top  else side_bar_px
    ebba_top_bar_raw = sb_ebba_top["px_len"] if sb_ebba_top else side_bar_px

    first_gfp  = prepare_side_frame(frames_gfp[0],  sb_side)
    first_ebba = prepare_side_frame(frames_ebba[0], sb_side)
    side_h, side_w = first_gfp.shape[:2]
    left_h = 2 * side_h + ROW_GAP_PX

    first_gfp_top,  gfp_scale  = prepare_top_frame(frames_gfp_top[0],  sb_gfp_top,  left_h)
    first_ebba_top, ebba_scale = prepare_top_frame(frames_ebba_top[0], sb_ebba_top, left_h)

    _, gfp_top_w  = first_gfp_top.shape[:2]
    _, ebba_top_w = first_ebba_top.shape[:2]

    gfp_top_bar  = max(1, int(round(gfp_top_bar_raw  * gfp_scale)))
    ebba_top_bar = max(1, int(round(ebba_top_bar_raw * ebba_scale)))

    total_h = BOTTOM_MARGIN_PX + left_h + TOP_TITLE_BAND_PX + TOP_MARGIN_PX
    total_w = side_w + COL_GAP_PX + gfp_top_w + COL_GAP_PX + ebba_top_w + SIDEBAR_PX

    gfp_top_left  = side_w + COL_GAP_PX
    ebba_top_left = side_w + COL_GAP_PX + gfp_top_w + COL_GAP_PX

    def fx(px): return px / total_w
    def fy(px): return px / total_h

    fig = plt.figure(figsize=(total_w / DPI, total_h / DPI), dpi=DPI)
    fig.patch.set_facecolor("black")

    ax_ebba     = fig.add_axes([0,                  fy(BOTTOM_MARGIN_PX),                   fx(side_w),    fy(side_h)])
    ax_gfp      = fig.add_axes([0,                  fy(BOTTOM_MARGIN_PX + side_h + ROW_GAP_PX), fx(side_w), fy(side_h)])
    ax_gfp_top  = fig.add_axes([fx(gfp_top_left),   fy(BOTTOM_MARGIN_PX), fx(gfp_top_w),   fy(left_h)])
    ax_ebba_top = fig.add_axes([fx(ebba_top_left),  fy(BOTTOM_MARGIN_PX), fx(ebba_top_w),  fy(left_h)])

    for ax in (ax_gfp, ax_ebba, ax_gfp_top, ax_ebba_top):
        ax.set_facecolor("black")
        ax.axis("off")

    im_gfp      = ax_gfp.imshow(first_gfp,      interpolation="nearest", aspect="auto")
    im_ebba     = ax_ebba.imshow(first_ebba,     interpolation="nearest", aspect="auto")
    im_gfp_top  = ax_gfp_top.imshow(first_gfp_top,  interpolation="nearest", aspect="auto")
    im_ebba_top = ax_ebba_top.imshow(first_ebba_top, interpolation="nearest", aspect="auto")

    title_y_top   = BOTTOM_MARGIN_PX * 0.7 + left_h + TOP_TITLE_BAND_PX * 0.1
    title_y_ebba  = BOTTOM_MARGIN_PX * 0.7 + side_h + ROW_GAP_PX * 0.1

    add_fig_title(fig, side_w / 2,                         title_y_top,  "GFP side view",            total_w, total_h)
    add_fig_title(fig, side_w / 2,                         title_y_ebba, "EbbaBiolight680 side view", total_w, total_h)
    add_fig_title(fig, gfp_top_left  + gfp_top_w  / 2,    title_y_top,  "GFP top view",              total_w, total_h)
    add_fig_title(fig, ebba_top_left + ebba_top_w / 2,     title_y_top,  "EbbaBiolight680 top view",  total_w, total_h)

    add_scale_bar_fig(fig, 0,            BOTTOM_MARGIN_PX, side_bar_px,  total_w, total_h)
    add_scale_bar_fig(fig, gfp_top_left, BOTTOM_MARGIN_PX, gfp_top_bar,  total_w, total_h)

    ts_shadow, ts_text, amp_shadow, amp_text = _add_timestamp_artists(fig, total_w, total_h, left_h)

    # Colorbars
    cbar_h_px = int(round(left_h * CBAR_HEIGHT_FRAC))
    cbar_bot  = BOTTOM_MARGIN_PX + int(round((left_h - cbar_h_px) / 2))
    cbar_base = ebba_top_left + ebba_top_w + CBAR_PAD_LEFT_PX
    cax_gfp_left  = cbar_base
    cax_ebba_left = cbar_base + CBAR_W_PX + CBAR_GAP_PX

    cax_gfp  = fig.add_axes([fx(cax_gfp_left),  fy(cbar_bot), fx(CBAR_W_PX), fy(cbar_h_px)])
    cax_ebba = fig.add_axes([fx(cax_ebba_left), fy(cbar_bot), fx(CBAR_W_PX), fy(cbar_h_px)])

    sm_gfp  = mcm.ScalarMappable(cmap=CMAP_GFP,  norm=Normalize(*GFP_RANGE));  sm_gfp.set_array([])
    sm_ebba = mcm.ScalarMappable(cmap=CMAP_GLOW, norm=Normalize(*GLOW_RANGE)); sm_ebba.set_array([])
    style_colorbar(fig.colorbar(sm_gfp,  cax=cax_gfp),  "GFP fluor. intensity (a.u.)")
    style_colorbar(fig.colorbar(sm_ebba, cax=cax_ebba), "EbbaBiolight680 fluor. intensity (a.u.)")

    canvas = FigureCanvasAgg(fig)
    canvas.draw()

    def render_fn(i):
        im_gfp.set_data(prepare_side_frame(frames_gfp[i], sb_side))
        im_ebba.set_data(prepare_side_frame(frames_ebba[i], sb_side))
        im_gfp_top.set_data(prepare_top_frame(frames_gfp_top[i], sb_gfp_top, left_h)[0])
        im_ebba_top.set_data(prepare_top_frame(frames_ebba_top[i], sb_ebba_top, left_h)[0])
        _update_labels(ts_shadow, ts_text, amp_shadow, amp_text, i, amp_name)

    _save_frames_and_encode(canvas, out_dir, tmp_png_dir, out_mp4,
                            frames_gfp, make_png, make_mp4, render_fn, n_frames, mode="ebba")


'''Batch driver'''
def batch_export(root_dir: Path):
    png_out = root_dir / "png_grouped"
    mp4_out = root_dir / "mp4_grouped"
    png_out.mkdir(exist_ok=True)
    mp4_out.mkdir(exist_ok=True)

    all_dirs = [d for d in root_dir.iterdir() if d.is_dir()
                and d not in (png_out, mp4_out)]

    by_prefix = {}
    for d in all_dirs:
        name = d.name
        # longest suffix first to avoid partial matches
        for suffix, key in [
            ("_ebba_glow_top", "ebba_top"),
            ("_ebba_glow",     "ebba"),
            ("_gfp_top",       "gfp_top"),
            ("_comp_top",      "top"),
            ("_sytox",         "sytox"),
            ("_gfp",           "gfp"),
        ]:
            if name.endswith(suffix):
                prefix = name[:-len(suffix)]
                by_prefix.setdefault(prefix, {})[key] = d
                break

    prefixes = sorted(p for p in by_prefix if p.startswith("xy"))
    print(f"Found {len(prefixes)} positions.\n")

    for prefix in prefixes:
        group = by_prefix[prefix]
        mode  = CHANNEL_MODE.get(prefix, "sytox")
        amp   = CONDITION_MAP.get(prefix, "")

        if mode == "sytox":
            required = {"gfp", "sytox", "top"}
            if not required.issubset(group):
                print(f"[Skip] {prefix}: missing {sorted(required - set(group))}")
                continue

            export_sytox_position(
                xy_prefix=prefix,
                folder_gfp=group["gfp"],
                folder_sytox=group["sytox"],
                folder_top=group["top"],
                out_dir=png_out / f"{prefix}_grouped",
                out_mp4=mp4_out / f"{prefix}_grouped.mp4",
                amp_name=amp,
                make_png=MAKE_PNG,
                make_mp4=MAKE_MP4,
            )

        elif mode == "ebba":
            required = {"gfp", "ebba", "gfp_top", "ebba_top"}
            if not required.issubset(group):
                print(f"[Skip] {prefix}: missing {sorted(required - set(group))}")
                continue

            export_ebba_position(
                xy_prefix=prefix,
                folder_gfp=group["gfp"],
                folder_ebba=group["ebba"],
                folder_gfp_top=group["gfp_top"],
                folder_ebba_top=group["ebba_top"],
                out_dir=png_out / f"{prefix}_grouped_ebba",
                out_mp4=mp4_out / f"{prefix}_grouped_ebba.mp4",
                amp_name=amp,
                make_png=MAKE_PNG,
                make_mp4=MAKE_MP4,
            )

        else:
            print(f"[Skip] {prefix}: unknown mode '{mode}'")


'''Run'''
batch_export(ROOT_DIR)
