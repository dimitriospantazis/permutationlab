# PermutationLab Toolbox

PermutationLab is a MATLAB-based toolbox designed for statistical inference on time-series and multidimensional data using non-parametric permutation-based tests. This approach controls for multiple comparisons and provides robust statistical significance testing. The toolbox supports both max-statistic-based inference and cluster-based inference, making it a powerful tool for analyzing complex datasets.

## Key Features

- **One-sample and Two-sample Tests:** Perform statistical tests on one or two datasets.
- **Max-Statistic-Based Inference:** Control for familywise error rate (FWER) using the maximum statistic.
- **Cluster-Based Inference:** Identify and test significant clusters in the data using various cluster statistics.
- **Multidimensional Data Support:** Analyze time-series and higher-dimensional data such as temporal generalization maps.
- **Flexible Statistical Parameters:** Customize tests with options for significance level, tail of the hypothesis, number of permutations, and more.
- **Visualization Tools:** Generate plots for statistical maps and significant clusters.
- **User-Friendly Interface:** Simple function calls with customizable parameters.

## Getting Started

### Prerequisites

- MATLAB (version 2020 or later recommended)
- Basic knowledge of MATLAB scripting

### Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/yourusername/permutationlab.git
   ```

2. Add the toolbox to your MATLAB path:
    ```matlab
    addpath(genpath('path_to/permutationlab'))
    ```

### Example Usage

Max-Statistic-Based Permutation Tests

    ```matlab
    % Perform a one-sample test
    stat = pl_permtest(data, 'tail', 'both', 'statistic', 'tstat', 'numpermutation', 1000, 'alpha', 0.05);    
    ```

Cluster-Based Permutation Tests

    ```matlab
    % Perform a one-sample cluster test
    stat = pl_permtestcluster(data, 'verbose', 10, 'clusteralpha', 0.05);
    ```

Visualization

    ```matlab
    figure; hold on; set(gcf, 'color', 'white');
    plot(time, stat.statmap, 'linewidth', 2, 'color', [0 0 0]);
    xlabel('Time (s)');
    ylabel('Decoding Accuracy (%)');
    ```

### Functions

- `pl_permtest`: One-sample permutation test  
- `pl_permtest2`: Two-sample permutation test  
- `pl_permtestcluster`: One-sample cluster-based permutation test  
- `pl_permtestcluster2`: Two-sample cluster-based permutation test  
- `pl_permtestclustertess`: Cluster-based tests with tessellated data  
- `pl_permtestclustertess2`: Two-sample cluster tests with tessellated data  

### Output

The primary output is a structure stat containing:

- `.statmap`: Statistical map of the test
- `.statmappv`: p-value map
- `.FDR.criticalmap`: False discovery rate (FDR) significance map
- `.FWER.criticalmap`: Familywise error rate (FWER) significance map
- `.clusters`: Significant clusters and their properties (for cluster tests)

### Licence

This toolbox is provided "as is," without any guarantees or warranties, and is available for unrestricted use. If you use this toolbox in your work, please consider citing the author.

