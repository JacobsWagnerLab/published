# 5 — Omnipose retraining & segmentation models

Tools to build the curated training datasets and the **retrained Omnipose
models** used for cell segmentation throughout this study. 

## Scripts

| Script | Purpose |
|---|---|
| `screening_good_bad.py` | Screen good vs. bad segmentation masks; attempts to recover poorly segmented cells with `binary_fill_holes`; |
| `generate_training_dataset.py` | Build the training dataset from curated images: create 3 augmented replicates per image by 90/180/270° rotation; save the full dataset into one folder. |
| `generate_training_dataset_wide-trench.py` | Build the wide-trench training dataset. Unlike `generate_training_dataset.py`, it does **not** rotate images, since that would break the vertical trench geometry. Instead, each curated image-mask pair yields 4 variants: original, horizontally flipped, photometrically/defocus-augmented, and augmented+flipped. Can optionally tile several wide-trench images side by side into larger training images. |

## Retrained models (`retrained_omnipose_models/`)

| Model | Used for | Figures |
|---|---|---|
| `omnipose_phase_microwells_model` | Phase-contrast segmentation of single cells in **PDMS-microwell** experiments (folder 2). Trained from scratch on bact_phase_omni-derived, hand-curated masks. | all microwell single-cell data |
| `omnipose_ribo_wide-trench_model` | **Ribosome-channel** segmentation of wide-trench cells (folder 3, growth rate). Fine-tuned from bact_fluor_omni (N = 592 curated 4-trench images). Preferred for single-cell segmentation. | 5C–F, 5I, S8A–B; Videos 5–6 |
| `omnipose_phase_wide-trench_model` | **Phase-contrast** segmentation of wide-trench cells; used to define cell-containing regions for the SYTOX Green / LL-37 kymographs (folder 3b). Trained from scratch (N = 868 curated images). | 5G–H, S8C–G |

## Workflow

Every retraining round shares the same curation loop: an existing model is
applied to new data, its masks are curated, the curated set is augmented, and a
model is (re)trained.

```
apply current model        screening_good_bad.py        generate_training_dataset.py           Omnipose training
to new data            ──▶  keep good / recover via  ──▶ (90/180/270° rotational           ──▶  (external; https://
(segment images)            hole-fill / discard bad      augmentation) if single-cell,           github.com/kevinjohncutler/omnipose)
                                                          or generate_training_dataset_wide-trench.py
                                                          (horizontal flip + photometric/defocus
                                                          augmentation, no rotation, optional tiling)
                                                          if wide-trenches → training folder
```

The **wide-trench** models were not trained in a single pass — they were built
over three bootstrapped rounds, each generating training data with the model
from the previous round (see Methods, *"Training of Omnipose model for
segmentation of wide-trench data"*):

```
Round 1 — preliminary phase model (intermediate, not used)
    bact_fluor_omni ──apply to ribosome channel──▶ curate + pair with phase
        └──▶ train phase model from scratch   (N = 353 curated images)

Round 2 — omnipose_ribo_wide-trench_model  (final ribosome model)
    preliminary phase model ──apply to multi-AMP data, curate vs ribosome signal──▶
        └──▶ fine-tune bact_fluor_omni        (N = 592 curated 4-trench images)

Round 3 — omnipose_phase_wide-trench_model  (final phase model)
    omnipose_ribo_wide-trench_model ──apply to subset, screen vs phase──▶ curate
        └──▶ train phase model from scratch   (N = 868 curated images)

                              ▼
        retrained_omnipose_models/   (used by OmniSegger in folders 2 & 3)
```

The **microwell** model (`omnipose_phase_microwells_model`, folder 2) was also retrained in multiple
rounds, see details in published paper Thappeta et al. 2025.

Training parameters (per Methods): trained from scratch or fine-tuned, 4000
epochs, learning rate 0.1, GPU (NVIDIA A40), 224×224-px tiles.

> The training datasets and final models are also archived at the BioImage Archive
> (accession [S-BIAD1823](https://www.ebi.ac.uk/biostudies/bioimages/studies/S-BIAD1823),
> dataset folder `training_data_omnipose`) and the Jacobs-Wagner lab Github
> repository, respectively.
