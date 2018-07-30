# UnobservableQueues.jl
Order-based estimation algorithm for unobservable queues

## Installation
`Pkg.clone("https://github.com/ajkeith/UnobservableQueues.jl")`

## Algorithm
The order-based estimation algorithm is implemented in src/estimators.jl. This file also includes our implementation of the variance estimator due to Park, Kim & Willemain (2011). The LCFS versions are in src/lcfs.jl.

## Testing
Tests for the correctness of each algorithm are available in test/runtests.jl.

## Micro-Benchmarks
Micro-benchmarks are available in the benchmarks directory.

## Analysis
The experimental design, regression analysis, and plots are located in the analysis directory. All analysis was conducted in JMP.

## Output Data

### FCFS
- *err11.csv*: first replicate of estimation error with column order of order-based, variance, uninformed
- *err12.csv*: second replicate of estimation error with column order of order-based, variance, uninformed
- *err_meas11.csv*: first replicate of estimation error under measurement error with column order of order-based, variance, uninformed
- *err_meas12.csv*: second replicate of estimation error under measurement error with column order of order-based, variance, uninformed
- *conv11.csv*: first replicate of convergence with column order of order-based, variance, uninformed
- *conv12.csv*: second replicate of convergence with column order of order-based, variance, uninformed
- *conv_meas11.csv*: first replicate of convergence with measurement error and with column order of order-based, variance, uninformed
- *conv_meas12.csv*: second replicate of convergence with measurement error and with column order of order-based, variance, uninformed

### LCFS
- *err13.csv*: first replicate of LCFS estimation error with column order of order-based, variance, uninformed
- *err_meas13.csv*: first replicate of LCFS estimation error under measurement error with column order of order-based, variance, uninformed
- *conv13.csv*: first replicate of LCFS convergence with column order of order-based, variance, uninformed
- *conv_meas13.csv*: first replicate of LCFS convergence with measurement error and with column order of order-based, variance, uninformed

### Archived Data
These data files use older versions of the code.
- *err.csv*: first replicate of estimation error with column order of order-based, variance, uninformed
- *err2.csv*: second replicate of estimation error with column order of order-based, variance, uninformed
- *err_meas.csv*: first replicate of estimation error under measurement error with column order of order-based, variance, uninformed
- *err_meas2.csv*: second replicate of estimation error under measurement error with column order of order-based, variance, uninformed
- *conv.csv*: convergence with 400 obs requirement, with column order of order-based, variance, uninformed
- *conv_meas.csv*: convergence with 400 obs requirement and measurement error, with column order of order-based, variance, uninformed
- *conv2.csv*: convergence with 200 obs requirement, with column order of order-based, variance, uninformed
- *conv_meas2.csv*: convergence with 200 obs requirement and measurement error, with column order of order-based, variance, uninformed

## Primary References
Park, J., Kim, Y. B. Y. B., & Willemain, T. R. T. R. (2011). Analysis of an unobservable queue using arrival and departure times. Computers and Industrial Engineering, 61(3), 842–847. https://doi.org/10.1016/j.cie.2011.05.017
