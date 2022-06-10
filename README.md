# Chaterjee's xi correlation for MATLAB
[![View Chaterjee's xi correlation on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/112530-chaterjee-s-xi-correlation)

What is Chaterjee's xi correlation:
- Measure of the dependency between 2 variables.
- Able to capture **non-linear** dependencies. Pearson/Spearman only work for linear trends (see below).
- Tool for inference (if offers p-values or confidence intervals).

How to use it:
````
x = linspace(0, 20, 1000);
y = cos(x) + 0.3*randn(1, 1000);

xi = xicor(x, y);
````

Post any bugs/ideas/comments as [issues!](https://github.com/drombas/xicor-mat/issues)

The difference between Pearson's ($r$) and xi-correlation ($\xi$): 

![Image description](fig_example.jpg)

References:
1. Sourav Chatterjee, A New Coefficient of Correlation, Journal of the American Statistical Association, 116:536, 2009-2022, 2021. [DOI:10.1080/01621459.2020.1758115](https://doi.org/10.1080/01621459.2020.1758115)
2. XICOR R package: https://cran.r-project.org/web/packages/XICOR/index.html
