# MFDMA (Multifractal Detrended Moving Average)
A Python library for performing 1D and 2D multifractal analysis using the Detrended Moving Average (DMA) method. Supports multifractal spectrum estimation with customizable parameters for in-depth data analysis.

## Citation

If you use the `mfdma` library in your research, please cite the following paper:

Gu, G.-F., & Zhou, W.-X. (2010). Detrended moving average algorithm for multifractals.  
*Physical Review E, 82*(1), 011136. [https://doi.org/10.1103/PhysRevE.82.011136](https://doi.org/10.1103/PhysRevE.82.011136)

## Features

- **MFDMA_1D**: Perform multifractal analysis on 1D time series data.
- **MFDMA_2D**: Conduct multifractal analysis on 2D data sets.

## Installation

You can install the package directly from the GitHub repository using pip:

```bash
pip install git+https://github.com/rhimonsouza/mfdma.git

#private

    ```bash
    pip install git+https://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@github.com/rhimonsouza/mfdma.git
    ```

## Usage

### 1D MF-DMA Analysis

To perform multifractal analysis on a 1D time series, use the `MFDMA_1D` function.

#### Parameters

- `x`: 1D array-like input data.
- `n_min`: Minimum window size for the log-spaced window sizes.
- `n_max`: Maximum window size for the log-spaced window sizes.
- `N`: Number of scales to consider between `n_min` and `n_max`.
- `theta`: Offset parameter controlling detrending (0 ≤ theta ≤ 1).
- `q`: Array of moments for multifractal analysis (can include positive and negative values).

#### Example

```python
from mfdma import MFDMA_1D

# Example data
x = [your_1D_data]

# Parameters
n_min = 10
n_max = 100
N = 50
theta = 0.5
q = [-2, -1, 0, 1, 2]

# Perform 1D MF-DMA
n, Fq, tau, alpha, f = MFDMA_1D(x, n_min, n_max, N, theta, q)


## Usage

### 2D MF-DMA Analysis

To perform multifractal analysis on a 2D dataset, use the `MFDMA_2D` function.

#### Parameters

- `X`: 2D numpy array representing the spatial data to analyze.
- `n_min`: Minimum window size for the moving average.
- `n_max`: Maximum window size for the moving average.
- `N`: Number of window sizes to consider between `n_min` and `n_max`.
- `theta`: Overlap parameter, controlling the degree of overlap between windows (0 ≤ theta ≤ 1).
- `q`: Array of q-order moments used in the multifractal analysis (including positive and negative values).

#### Example

```python
from mfdma import MFDMA_2D
import numpy as np

# Example 2D data, such as an image
X = np.random.rand(100, 100)  # Replace with your 2D data

# Parameters for 2D MF-DMA analysis
n_min = 5
n_max = 50
N = 30
theta = 0.3
q = [-2, -1, 0, 1, 2]

# Perform 2D MF-DMA
n, Fq, tau, alpha, f = MFDMA_2D(X, n_min, n_max, N, theta, q)



