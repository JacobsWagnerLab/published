# -*- coding: utf-8 -*-
"""
Created on Sun Jun 15 14:11:43 2025

@author: fragasso
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import shift
from skimage.registration import phase_cross_correlation
from scipy.signal import savgol_filter

'''This function prompts overlay of first and last frame of the phase contrast image sequence, 
and allows user to define ROI to crop and plot a single wide trench over time'''
def get_user_roi(image1, image2, title = None):
    """Overlay first and last image with a grid and prompt user for ROI coordinates, allowing corrections."""
    while True:
        overlay = np.zeros((*image1.shape, 3), dtype=np.uint8)
        overlay[..., 1] = (image1 / image1.max() * 255).astype(np.uint8)  # Green channel
        overlay[..., 2] = (image2 / image2.max() * 255).astype(np.uint8)  # Magenta channel
        
        plt.figure(figsize=(12, 8))
        plt.imshow(overlay)
        plt.title(title if title else 'Enter x,y coordinates')
        
        # Add grid and increase tick density
        plt.grid(color='white', linestyle='-', linewidth=0.5)
        plt.xticks(np.arange(0, image1.shape[1], step=50), rotation=45)  # Increase tick density
        plt.yticks(np.arange(0, image1.shape[0], step=50), rotation=45)
        
        plt.show()
        
        x_min = int(input("Enter x_min: "))
        x_max = int(input("Enter x_max: "))
        y_min = int(input("Enter y_min: "))
        y_max = int(input("Enter y_max: "))
        
        # Replot with selected ROI
        fig = plt.figure(figsize=(12, 12))
        ax1 = fig.add_subplot(2,1,1)
        plt.imshow(image1)
        # plt.title(title if title else 'Chosen region of interest')
        plt.grid(color='white', linestyle='--', linewidth=0.5)
        plt.xticks(np.arange(0, image1.shape[1], step=50), rotation=45)
        plt.yticks(np.arange(0, image1.shape[0], step=50), rotation=45)
        plt.gca().add_patch(plt.Rectangle((x_min, y_min), x_max-x_min, y_max-y_min, edgecolor='red', facecolor='none', linewidth=2))
        
        ax1 = fig.add_subplot(2,1,2)
        plt.imshow(image2)
        # plt.title(title if title else 'Chosen region of interest')
        plt.grid(color='white', linestyle='--', linewidth=0.5)
        plt.xticks(np.arange(0, image1.shape[1], step=50), rotation=45)
        plt.yticks(np.arange(0, image1.shape[0], step=50), rotation=45)
        plt.gca().add_patch(plt.Rectangle((x_min, y_min), x_max-x_min, y_max-y_min, edgecolor='red', facecolor='none', linewidth=2))
        
        plt.show()
        
        confirm = input("Are the coordinates correct? (y/n): ").strip().lower()
        if confirm == 'y':
            break
        else:
            print("Re-enter ROI coordinates.")
    
    return y_min, y_max, x_min, x_max

"""
This function corrects Y-axis drift by aligning each frame to the first frame using phase cross-correlation.
Handles top-padding artifacts, removes sharp outliers based on first derivative (remove_spike_outliers), then applies Savitzky-Golay-smoothed Y-axis shifts to all channels (phase, fluor1, and fluor2)
"""
def drift_correction_wide_trenches(image_paths, fluor1_paths=None, fluor2_paths=None, show=False):

    def remove_spike_outliers(y, threshold=5.0):
        """Detects spikes based on derivative and replaces with interpolated values."""
        y = y.copy()
        dy = np.diff(y, prepend=y[0])
        spike_idx = np.where(np.abs(dy) > threshold)[0]
        for i in spike_idx:
            if 0 < i < len(y) - 1:
                y[i] = 0.5 * (y[i - 1] + y[i + 1])
        return y, spike_idx

    print("Performing Y-axis drift correction...")

    images = image_paths
    min_height = min(img.shape[0] for img in images)
    min_width = min(img.shape[1] for img in images)
    images = [img[:min_height, :min_width] for img in images]

    y_min, y_max, x_min, x_max = get_user_roi(images[0], images[-1], title="Select region with a fixed feature")

    ref_roi = images[0][y_min:y_max, x_min:x_max]
    y_shifts = [0.0]  # First frame is reference

    for i in range(1, len(images)):
        curr_roi_full = images[i][y_min:y_max, x_min:x_max]

        # Identify top padding
        top_row = curr_roi_full[0, :]
        nonpad_start = np.argmax([not np.all(row == top_row) for row in curr_roi_full])
        if np.all(curr_roi_full == top_row):
            nonpad_start = 0  # fallback if the whole ROI is flat

        curr_roi = curr_roi_full[nonpad_start:, :]
        ref_roi_cropped = ref_roi[:curr_roi.shape[0], :]

        if nonpad_start > 1 and show:
            fig, axs = plt.subplots(1, 2, figsize=(10, 6))
            axs[0].imshow(ref_roi_cropped, cmap='gray')
            axs[0].set_title("Reference cropped")
            axs[1].imshow(curr_roi, cmap='gray')
            axs[1].set_title("Current unpadded")
            plt.show()

        shift_vector, _, _ = phase_cross_correlation(ref_roi_cropped, curr_roi, upsample_factor=100)
        y_shift = shift_vector[0] - nonpad_start  # account for removed rows
        y_shifts.append(y_shift)
        print(f"Frame {i}: raw y-shift = {y_shift:.2f}; nonpad_start = {nonpad_start}")

    raw_shifts = np.array(y_shifts)
    cleaned_y, spike_idx = remove_spike_outliers(raw_shifts, threshold=5.0)

    # Interactive smoothing window selection
    window = 35
    while True:
        try:
            smoothed_y = savgol_filter(cleaned_y, window_length=window, polyorder=2)
        except ValueError:
            print(f"Invalid window length ({window}). Must be odd and <= number of frames.")
            window = int(input("Enter an odd Savitzky-Golay window size (<= number of frames): "))
            continue

        smoothed_shifts = np.column_stack((smoothed_y, np.zeros_like(smoothed_y)))

        plt.figure()
        plt.plot(smoothed_y, label='Smoothed Y-shift')
        plt.scatter(range(len(cleaned_y)), cleaned_y, s=6, marker='.', label='Cleaned Y-shift')
        if len(spike_idx) > 0:
            plt.scatter(spike_idx, raw_shifts[spike_idx], c='r', s=10, label='Spike outliers')
        plt.title("Drift magnitude per frame, window = " + str(window))
        plt.xlabel("Frame")
        plt.ylabel("Shift (pixels)")
        plt.legend()
        plt.show()

        user_input = input(f"Is window = {window} ok? (y/n): ").strip().lower()
        if user_input == 'y':
            break
        elif user_input == 'n':
            window = int(input("Enter a new odd window size (>=3 and <= number of frames): "))
        else:
            print("Please enter 'y' or 'n'.")

    print("Applying Y-drift correction...")

    def apply_shift(images, shifts):
        if images is None:
            return None
        return [shift(img, shift=sh, mode='reflect', order=3) for img, sh in zip(images, shifts)]

    corrected_images = apply_shift(images, smoothed_shifts)
    corrected_fluor1 = apply_shift(fluor1_paths, smoothed_shifts)
    corrected_fluor2 = apply_shift(fluor2_paths, smoothed_shifts)

    print("Y-axis drift correction complete.")
    return corrected_images, corrected_fluor1, corrected_fluor2, smoothed_shifts, raw_shifts
