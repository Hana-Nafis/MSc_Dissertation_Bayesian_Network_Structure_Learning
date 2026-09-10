# MSc_Dissertation_Bayesian_Network_Structure_Learning
Repository for evaluating Bayesian network structure learning algorithms under dimensional scaling and causal do-calculus. Features two core scripts: data simulation and end-to-end pipeline execution for structural and interventional validation.

## Empirical Evaluation of Bayesian Network Structure Learning Algorithms under Dimensional Scaling and Causal do-calculus

Official repository for evaluating Bayesian network structure learning algorithms under dimensional scaling and causal $do$-calculus. This codebase implements multi-paradigm structural benchmarks across sub-critical, intermediate, and super-critical network topologies (ALARM, Hailfinder, and Hepar II) to assess how structural errors propagate into downstream interventional predictions.

### Repository Structure

* `data_generation.R`: Simulates and exports tiered datasets across sample sizes $n \in \{500, 1000, 2500, 5000, 10000\}$ using distinct, controlled random seeds for each benchmark network.
* `Pipeline_Execution_.R`: Executes end-to-end constraint-based (PC), score-based (Hill Climbing), and hybrid (MMHC) learning algorithms, computes structural hamming distances (SHD), evaluates localized Markov blanket recoveries, runs interventional simulations under Pearl's $do$-calculus, and generates comparative validation plots.

### Prerequisites and Dependencies

The pipeline requires R along with the following packages:

```r
install.packages(c("bnlearn", "ggplot2", "gRain"))
