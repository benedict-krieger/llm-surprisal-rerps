# About

Implementation of the rERP analysis reported in

[Krieger, B., Brouwer, H., Aurnhammer, C., & Crocker, M. W. (2024). On the limits of LLM surprisal as functional Explanation of ERPs. Proceedings of the Annual Meeting of the Cognitive Science Society, 46.](https://escholarship.org/uc/item/2m53k85t#main)

Evaluating ERP data from

[Aurnhammer, C., Delogu, F., Brouwer, H., & Crocker, M. (2023). The P600 as a continuous index of integration effort. Psychophysiology](https://doi.org/10.1111/psyp.14302)

[Aurnhammer, C., Delogu, F., Schulz, M., Brouwer, H., & Crocker, M. (2021). Retrieval (N400) and integration (P600) in expectation-based comprehension. PLOS ONE](https://doi.org/10.1371/journal.pone.0257430)

[Delogu, F., Brouwer, H., & Crocker, M. (2019). Event-related potentials index lexical retrieval (N400) and integration (P600) during language comprehension. Brain and Cognition](https://doi.org/https://doi.org/10.1016/j.bandc.2019.05.007)

The code for the rERP analysis (Julia & R) is based on the implementation by Christoph Aurnhammer, found [here](https://github.com/caurnhammer/psyp23rerps).

# Requirements

- 10 GB free disk space

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
    - StatsBase

- R (tested on 4.3.1)
    - here
    - glue
    - data.table
    - ggplot2
    - grid
    - gridExtra
    - stringr

- GNU Make (optional)


# Usage

From within `code` directory, to reproduce all results:

```
make analysis
```

This may take some time to run. Collection of surprisal values is excluded, since specifically the Llama-2 LLM (Leo13b) requires a lot of computational resources. Surprisal values are included in the stimulus data files (e.g. `adbc23.csv`), and can be reproduced by running `make collect_secretgpt2_surp` and `make collect_leo13b_surp` (the latter is not recommended to be run locally).

The individual parts of the analysis can also be run separately:

```
    make paths # create required directory structure for results
    make eval_surp # merge surprisal and original ERP data, create density plots and check BPE splits
    make rERP # run rERP analysis (Julia)
    make plot_rERP # create plots of rERP results (R)
    make clean # remove all results and rERP data files
```


If not using GNU Make, the code can also be run directly in the order indicated by the naming scheme, with `01_collect_surprisal.py` being optional.

# File structure

## code
- numbered files indicate subsequent steps of analysis (with 01 being optional)

## data
- organized by study ids and surprisal ids:
    - adbc23: Aurnhammer et al. (2023)
    - adsbc21: Aurnhammer et al. (2021)
    - dbc19: Delogu et al. (2019)
    - leo13b_surp : Leo-13b surprisal
    - secretgpt2_surp : GPT-2 surprisal

- each study's directory initially contains
    - `{study_id}.csv` : stimulus data, including human ratings and llm surprisal values
    - `{study_id}_erp.csv` : erp data

- `02_eval_surprisal.py` will merge ERP data and surprisal values, resulting in `{study_id}_surp_erp.csv`

- this file will be used in `03_do_rERP.jl` to create an LLM-specific rERP file: `{study_id}_{surp_id}_rERP.csv`

- this file will be used to then create two further rERP files `{study_id}_{surp_id}_rERP_data.csv` and `{study_id}_{surp_id}_rERP_models.csv`

- the latter two files will also have a version used for inferential statistics, indicated by `{across_subj}`

## results
- organized by study ids and surprisal ids (see above)

- `{study_id}_corr.csv` contains Kendall correlations between human judgements and surprisal values

- `plots` contains density plots of surprisal values per condition (one for each LLM)

- `plots` contains further subdirectories with plots of the rERP results (observed, coefficients, forward estimates, residuals)

- the `{across_subj}` subdirectories contain t-values