# Preprocessing extracted features

In this module, we apply preprocessing scripts to:

1. Harmonize the extracted features from SQLite to single-cell parquet profile with CytoTable.
2. Clean the single-cell profile to remove mis-segmentations and blurry cells with coSMicQC.
3. Annotate, normalize, and feature select the profile with Pycytominer.

## Perform preprocessing on data

To perform preprocessing on the profile, run the bash script [perform_preprocessing.sh](./perform_preprocessing.sh) using the command below:

```bash
source perform_preprocessing.sh
```
