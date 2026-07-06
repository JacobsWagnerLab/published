# -*- coding: utf-8 -*-
"""
Assemble PNG frames produced by plot_single_trench_growth_rate.py into an MP4 video.

Steps:
1) Read PNG frames from input_dir, sorted by frame number
2) Optionally repeat each frame repeat_factor times to slow down playback
3) Encode with ffmpeg (libx264, crf=18) and write to output_path

Requires: imageio-ffmpeg (pip install imageio-ffmpeg)

@author: alessio fragasso
"""

import os
import re
import shutil
import subprocess
import tempfile
import imageio_ffmpeg as iio_ffmpeg


'''Settings'''
input_dir  = #r"path_to_png_frames_folder"     # folder containing frame_0.png, frame_1.png, ...
output_dir = #r"path_to_output_folder"
output_file = #"output_filename.mp4"

fps = 10           # playback speed (frames per second)
repeat_factor = 1  # set > 1 to duplicate each frame and slow down the video


'''Assemble and encode'''
os.makedirs(output_dir, exist_ok=True)
output_path = os.path.join(output_dir, output_file)

def extract_number(fname):
    m = re.search(r"frame_(\d+(?:\.\d+)?)\.png$", fname, re.IGNORECASE)
    return float(m.group(1)) if m else float("inf")

files = [f for f in os.listdir(input_dir) if f.lower().endswith(".png")]
files.sort(key=extract_number)

if not files:
    raise RuntimeError(f"No PNG files found in:\n{input_dir}")

print("First files:", files[:10])
print("Last files:", files[-10:])
print("N original files:", len(files))

ffmpeg_exe = iio_ffmpeg.get_ffmpeg_exe()

with tempfile.TemporaryDirectory() as tmpdir:
    frame_i = 0
    for f in files:
        src = os.path.join(input_dir, f)
        for _ in range(repeat_factor):
            dst = os.path.join(tmpdir, f"img_{frame_i:06d}.png")
            shutil.copy2(src, dst)
            frame_i += 1

    print("N video frames:", frame_i)

    input_pattern = os.path.join(tmpdir, "img_%06d.png")
    cmd = [
        ffmpeg_exe,
        "-y",
        "-framerate", str(fps),
        "-i", input_pattern,
        "-vf", "scale=1920:-2,format=yuv420p",   # rescale to 1920-wide for player compatibility
        "-c:v", "libx264",
        "-profile:v", "main",
        "-level", "4.0",
        "-pix_fmt", "yuv420p",
        "-crf", "18",
        "-preset", "slow",
        "-movflags", "+faststart",
        output_path
    ]

    print("Running ffmpeg...")
    print(" ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True)

print("RETURN CODE:", result.returncode)
print("STDOUT:"); print(result.stdout)
print("STDERR:"); print(result.stderr)

if result.returncode != 0:
    raise RuntimeError("ffmpeg failed.")

print("Output exists:", os.path.exists(output_path))
print("Output size MB:", os.path.getsize(output_path) / 1e6)
print(f"\nDone -> {output_path}")
