# -*- coding: utf-8 -*-
"""
Code for generating a wide-trench training dataset for Omnipose retraining.

Imports curated wide-trench images and their masks. For each image-mask pair it
creates four matched training examples:
    1. original image and mask
    2. horizontally flipped image and mask
    3. photometrically / defocus-augmented image and original mask
    4. horizontally flipped augmented image and flipped mask

Unlike the microwell training-data script, this script does not rotate images,
because preserving the vertical wide-trench geometry is important.

The augmented individual trenches can optionally be concatenated horizontally
into larger tile images for training.

@author: Alessio Fragasso
July 2026
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import tifffile
import imageio.v3 as iio
from skimage import measure, morphology


''' Inputs:
 - folder_path = folder containing one subfolder per condition. Each condition subfolder must contain paired TIFF files, e.g. image_name.tif and image_name_masks.tif.
 - save_folder = destination folder where to store the curated images
'''

# Script containing the augment_phase() function used below.
aug_functions_path = Path(
    r"N:\Alessio_Fragasso\Retrain_omnipose\omnipose_retraining_wide_trenches\code\aug_functions.py"
)

cond_list = [
    "1_PR-39_4uM",
    "14_CecA_10uM",
    "15_LL-37_20uM",
    "11_TM1_5uM",
]

# Set True to inspect the original, flipped, augmented, and augmented-flipped
# images with their corresponding mask contours.
plot_augmentations = False

# Reproducible, but different, augmentation draws for consecutive images.
random_seed = 123

# Run the two stages independently when needed.
run_augmentation = True
run_tiling = True

# %% Inputs for creating tiled training images
# Folder that receives the final tile image-mask pairs.
tile_save_path = Path(
    r"N:\Alessio_Fragasso\Retrain_omnipose\omnipose_retraining_wide_trenches_fluor\training_data"
)

# Conditions to tile. This may be a subset of cond_list.
tile_cond_list = ["11_TM1_5uM"]

# Number of wide-trench images to place side by side in each tile. The last
# tile is padded with zeros if fewer than this number of images remain.
images_per_tile = 3

# Retain the original cleaning criterion before tiling.
min_mask_size = 200

# The original code had a special crop for this condition. Keep it configurable.
# Use None to disable the crop.
special_crop_condition = "1_PR-39_4uM"
special_crop_height = 800

# Set True to inspect each final tile before it is written.
plot_tiles = False


# %% Helper functions
def load_augmentation_functions(script_path):
    """Load augment_phase() from the external augmentation helper script."""
    if not script_path.is_file():
        raise FileNotFoundError(
            f"Could not find aug_functions.py:\n{script_path}\n"
            "Update aug_functions_path before running this script."
        )

    with open(script_path, "r", encoding="utf-8") as file:
        exec(compile(file.read(), str(script_path), "exec"), globals())

    if "augment_phase" not in globals():
        raise NameError(
            "aug_functions.py was loaded, but it did not define augment_phase()."
        )


def get_training_pairs(condition_path):
    """Return matched image-mask TIFF pairs for one experimental condition."""
    tif_paths = sorted(
        [*condition_path.glob("*.tif"), *condition_path.glob("*.tiff")]
    )

    mask_paths = [path for path in tif_paths if "masks" in path.stem]
    image_paths = [path for path in tif_paths if "masks" not in path.stem]
    image_by_stem = {path.stem: path for path in image_paths}

    pairs = []
    missing_images = []

    for mask_path in mask_paths:
        image_stem = mask_path.stem.replace("_masks", "")
        image_path = image_by_stem.get(image_stem)

        if image_path is None:
            missing_images.append(mask_path.name)
        else:
            pairs.append((image_path, mask_path))

    if missing_images:
        missing_text = "\n  ".join(missing_images)
        raise FileNotFoundError(
            "Could not find a matching image for the following mask file(s):\n"
            f"  {missing_text}"
        )

    if not pairs:
        raise FileNotFoundError(
            f"No matched image-mask TIFF pairs were found in:\n{condition_path}"
        )

    return pairs


def get_mask_contours(mask):
    """Extract per-label contours for optional visual inspection."""
    contours_all = []

    for label_id in np.unique(mask):
        if label_id == 0:
            continue

        single_mask = mask == label_id
        contours = measure.find_contours(single_mask, level=0.5)
        contours_all.extend((label_id, contour) for contour in contours)

    unique_labels = np.unique(mask)
    unique_labels = unique_labels[unique_labels > 0]

    if len(unique_labels) == 0:
        return contours_all, {}

    cmap = plt.colormaps.get_cmap("tab20").resampled(len(unique_labels))
    color_for = {
        label_id: cmap(index)
        for index, label_id in enumerate(unique_labels)
    }

    return contours_all, color_for


def plot_image_with_contours(ax, image, contours_all, color_for, title):
    """Plot an image with mask contours for visual quality control."""
    ax.imshow(image, cmap="gray")
    ax.set_title(title)
    ax.axis("off")

    for label_id, contour in contours_all:
        ax.plot(
            contour[:, 1],
            contour[:, 0],
            linewidth=1,
            color=color_for[label_id],
        )


def make_augmented_examples(image, mask, rng):
    """Create the four image-mask variants used for wide-trench training."""
    image_augmented = augment_phase(image, rng)

    # Restore the original image dynamic range and datatype after augmentation.
    scale = np.percentile(image, 99.8) + 1e-8
    image_augmented = np.clip(image_augmented * scale, 0, np.iinfo(np.uint16).max)
    image_augmented = image_augmented.astype(np.uint16)

    image_flipped = np.fliplr(image)
    mask_flipped = np.fliplr(mask)
    image_augmented_flipped = np.fliplr(image_augmented)

    return {
        "_0": (image, mask),
        "_1": (image_flipped, mask_flipped),
        "_0_aug": (image_augmented, mask),
        "_1_aug": (image_augmented_flipped, mask_flipped),
    }


def save_augmented_condition(condition, pairs, destination_path, rng):
    """Create and save all individual augmented examples for one condition."""
    condition_save_path = destination_path / condition
    condition_save_path.mkdir(parents=True, exist_ok=True)

    for index, (image_path, mask_path) in enumerate(pairs, start=1):
        image = iio.imread(image_path)
        mask = iio.imread(mask_path)

        if image.shape != mask.shape:
            raise ValueError(
                f"Image and mask shape mismatch for:\n"
                f"  image: {image_path.name} {image.shape}\n"
                f"  mask:  {mask_path.name} {mask.shape}"
            )

        examples = make_augmented_examples(image, mask, rng)
        base_name = f"{condition}_train_{index:05d}"

        if plot_augmentations:
            original_contours, original_colors = get_mask_contours(mask)
            flipped_mask = examples["_1"][1]
            flipped_contours, flipped_colors = get_mask_contours(flipped_mask)

            fig, axes = plt.subplots(1, 4, figsize=(16, 5), constrained_layout=True)
            fig.suptitle(image_path.stem)

            plot_image_with_contours(
                axes[0], image, original_contours, original_colors, "original"
            )
            plot_image_with_contours(
                axes[1], examples["_1"][0], flipped_contours, flipped_colors, "flipped"
            )
            plot_image_with_contours(
                axes[2], examples["_0_aug"][0], original_contours, original_colors,
                "augmented",
            )
            plot_image_with_contours(
                axes[3], examples["_1_aug"][0], flipped_contours, flipped_colors,
                "augmented + flipped",
            )
            plt.show()

        for suffix, (example_image, example_mask) in examples.items():
            image_name = f"{base_name}{suffix}.tif"
            mask_name = f"{base_name}{suffix}_masks.tif"

            tifffile.imwrite(condition_save_path / image_name, example_image)
            tifffile.imwrite(condition_save_path / mask_name, example_mask)

        if index % 50 == 0 or index == len(pairs):
            print(
                f"Saved {index}/{len(pairs)} augmented image-mask pairs "
                f"for condition {condition}."
            )


def apply_optional_crop(image, condition):
    """Apply the condition-specific crop retained from the original script."""
    if (
        condition == special_crop_condition
        and special_crop_height is not None
        and image.shape[0] > special_crop_height
    ):
        return image[:special_crop_height, :]

    return image


def get_augmented_pairs(condition_path):
    """Return matched augmented image-mask pairs, preserving filename pairing."""
    return get_training_pairs(condition_path)


def create_tiles_for_condition(condition, augmented_path, destination_path, rng):
    """Concatenate wide-trench image-mask pairs horizontally into training tiles."""
    condition_input_path = augmented_path / condition
    condition_output_path = destination_path / condition
    condition_output_path.mkdir(parents=True, exist_ok=True)

    pairs = get_augmented_pairs(condition_input_path)
    order = rng.permutation(len(pairs))
    pairs = [pairs[index] for index in order]

    tile_index = 0
    for start in range(0, len(pairs), images_per_tile):
        selected_pairs = pairs[start:start + images_per_tile]

        images = []
        masks = []
        for image_path, mask_path in selected_pairs:
            image = apply_optional_crop(iio.imread(image_path), condition)
            mask = apply_optional_crop(iio.imread(mask_path), condition)
            mask = morphology.remove_small_objects(mask, min_size=min_mask_size)

            if image.shape != mask.shape:
                raise ValueError(
                    f"Image and mask shape mismatch after processing:\n"
                    f"  image: {image_path.name} {image.shape}\n"
                    f"  mask:  {mask_path.name} {mask.shape}"
                )

            images.append(image)
            masks.append(mask)

        image_heights = {image.shape[0] for image in images}
        image_widths = {image.shape[1] for image in images}
        if len(image_heights) != 1 or len(image_widths) != 1:
            raise ValueError(
                "All images within one tiled condition must have the same dimensions. "
                f"Found shapes: {[image.shape for image in images]}"
            )

        image_height, image_width = images[0].shape
        tile_width = images_per_tile * image_width

        tile_image = np.zeros(
            (image_height, tile_width),
            dtype=images[0].dtype,
        )
        tile_mask = np.zeros(
            (image_height, tile_width),
            dtype=masks[0].dtype,
        )

        for column, (image, mask) in enumerate(zip(images, masks)):
            x_start = column * image_width
            x_end = x_start + image_width
            tile_image[:, x_start:x_end] = image
            tile_mask[:, x_start:x_end] = mask

        tile_index += 1
        tile_image_name = f"train_{tile_index:05d}.tif"
        tile_mask_name = f"train_{tile_index:05d}_masks.tif"

        tifffile.imwrite(condition_output_path / tile_image_name, tile_image)
        tifffile.imwrite(condition_output_path / tile_mask_name, tile_mask)

        if plot_tiles:
            fig, axes = plt.subplots(1, 2, figsize=(16, 6), constrained_layout=True)
            axes[0].imshow(tile_image, cmap="gray")
            axes[0].set_title("Tile image")
            axes[0].axis("off")
            axes[1].imshow(tile_mask, cmap="inferno")
            axes[1].set_title("Tile mask")
            axes[1].axis("off")
            plt.show()

        if tile_index % 50 == 0 or start + images_per_tile >= len(pairs):
            print(
                f"Saved {tile_index} tile(s) for condition {condition}."
            )


# %% Run the workflow
def main():
    """Run the selected wide-trench dataset-generation stages."""
    if run_augmentation:
        load_augmentation_functions(aug_functions_path)
        augmentation_rng = np.random.default_rng(random_seed)

        for condition in cond_list:
            condition_path = folder_path / condition
            pairs = get_training_pairs(condition_path)
            save_augmented_condition(
                condition=condition,
                pairs=pairs,
                destination_path=save_path,
                rng=augmentation_rng,
            )

    if run_tiling:
        tiling_rng = np.random.default_rng(random_seed)

        for condition in tile_cond_list:
            create_tiles_for_condition(
                condition=condition,
                augmented_path=save_path,
                destination_path=tile_save_path,
                rng=tiling_rng,
            )


if __name__ == "__main__":
    main()
