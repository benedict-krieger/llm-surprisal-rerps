# About

Implementation of the rERP analysis reported in

[Krieger, B., Brouwer, H., Aurnhammer, C., & Crocker, M. W. (2025). On the limits of LLM surprisal as a functional explanation of the N400 and P600. Manuscript submitted for publication.]()

Evaluating ERP data from

[Aurnhammer, C., Delogu, F., Brouwer, H., & Crocker, M. (2023). The P600 as a continuous index of integration effort. Psychophysiology](https://doi.org/10.1111/psyp.14302)

[Aurnhammer, C., Delogu, F., Schulz, M., Brouwer, H., & Crocker, M. (2021). Retrieval (N400) and integration (P600) in expectation-based comprehension. PLOS ONE](https://doi.org/10.1371/journal.pone.0257430)

[Delogu, F., Brouwer, H., & Crocker, M. (2019). Event-related potentials index lexical retrieval (N400) and integration (P600) during language comprehension. Brain and Cognition](https://doi.org/https://doi.org/10.1016/j.bandc.2019.05.007)

The code for the rERP analysis (Julia & R) is based on the implementation by Christoph Aurnhammer, found [here](https://github.com/caurnhammer/psyp23rerps).

# Requirements

- 15 GB free disk space

- `data` (see release) must be downloaded and extracted in main project folder

- Python (tested on 3.11.9)
    - numpy
    - pandas
    - torch
    - transformers
    - seaborn
    - matplotlib

Using `conda`, you can create a dedicated environment:

```
conda env create -f llm-surprisal-rerps.yml
```

- Julia (tested on 1.9)
    - CategoricalArrays
    - Combinatorics
    - CSV
    - DataFrames
    - Distributions
    - LinearAlgebra
    - MixedModels
    - StatsBase

- R (tested on 4.3.1)
    - here
    - dplyr
    - glue
    - data.table
    - ggplot2
    - grid
    - gridExtra
    - stringr
    - viridisLite

- GNU Make (optional)


# Usage

From within `code` directory, to reproduce all results:

```
make analysis
```

This will take some time to run. Collection of surprisal values is excluded, since specifically the Llama-2 LLM (Leo13b) requires a lot of computational resources. Surprisal values are included in the stimulus data files (e.g. `adbc23.csv`), and can optionally be reproduced by running `python 01_collect_surprisal -m [model_id]` (which is not recommended to be run locally for Leo13b).

# File structure

## code
- numbered files indicate subsequent steps of analysis (with 01 being optional)

## data
- organized by study ids and surprisal ids:
    - adsbc21 (Study 1): Aurnhammer et al. (2021)
    - dbc19 (Study 2): Delogu et al. (2019)
    - dbc19_corrected (Study 2): Delogu et al. (2019), corrected for component overlap
    - adbc23 (Study3): Aurnhammer et al. (2023)
    - leo13b_surp : Leo-13b surprisal
    - gerpt2large_surp : GerPT-2 large surprisal
    - gerpt2_surp : GerPT-2 surprisal


## results
- organized by study ids and surprisal ids (see above)