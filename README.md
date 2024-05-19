


Two successful introductions can satisfy the two mutation separation constraint either as:
- a single clade descending on a two-mutation branch from a basal clade (here denoted $(\tau_p,\tau_p)_{1c}$); or
- two clades descending on one-mutation branches from the MRCA (here denoted $(\tau_p,\tau_p)_{2c}$).

These correspond, respectively, to the single introduction topologies $\tau_{1c}$ and $\tau_{2c}$. They have the same compatibility:
- $\tau = \tau_{1c}$ is compatible with $S_{MRCA} \in \{S_{A}, S_{B}\}$
- $\tau = \tau_{2c}$ is compatible with $S_{MRCA} \in \{S_{C/C}, S_{T/T}\}$
- $\tau = (\tau_p,\tau_p)\_{1c}$ is compatible with $S_{MRCA} \in \{S_{A}, S_{B}\}$
- $\tau = (\tau_p,\tau_p)\_{2c}$ is compatible with $S_{MRCA} \in \{S_{C/C}, S_{T/T}\}$

With these more specific two-introduction topologies, the compatibility vectors $P(S_{MRCA}|\tau)$ become:
- $P(S_{MRCA}|\tau = \tau_{1c}) = (0.5, 0.5, 0, 0)$
- $P(S_{MRCA}|\tau = \tau_{2c}) = (0, 0, 0.5, 0.5)$
- $P(S_{MRCA}|\tau = (\tau_p,\tau_p)\_{1c}) = (0.5, 0.5, 0, 0)$
- $P(S_{MRCA}|\tau = (\tau_p,\tau_p)\_{2c}) = (0, 0, 0.5, 0.5)$

The Bayes factor equations for the unconstrained and recCA rootings then become:

$$
BF_{uncon.} = \frac{ 1.76 \cdot P((\tau_p,\tau_p)\_{1c}|I_2) + 0.24\cdot  P(\tau_p,\tau_p)\_{2c}|I_2)}
{ 1.76 \cdot P(\tau_{1c}|I_1) + 0.24\cdot  P(\tau_{2c}|I_1)},
\quad
BF_{recCA} = \frac{ 1.72 \cdot P((\tau_p,\tau_p)\_{1c}|I_2) + 0.28\cdot  P(\tau_p,\tau_p)\_{2c}|I_2)}
{ 1.72 \cdot P(\tau_{1c}|I_1) + 0.28\cdot  P(\tau_{2c}|I_1)}
$$

To determine if two simulations have the $(\tau_p,\tau_p)\_{1c}$ or $(\tau_p,\tau_p)\_{2c}$ topologies, the numbers of mutations between the MRCA and each clade root can be sampled using the molecular clock and the times between the MRCA and each clade root. Pekar et al.'s  simulations provide the times between each introduction and subsequent clade root. Absent information on upstream population structures and dynamics, a reasonable model for the time between the MRCA and each introduction, and the time between the first introduction and the second introduction, samples from exponential distributions with expected values of $t_1$ and $t_2$, respectively. This allows possibilities to be explored for a range of parameters $t_1$ (representing the upstream effective population size) and $t_2$ (representing the introduction intensity).


The time between the first and second introductions $t_2$ also provides a lag that can be incorporated into the relative size test.


There is some ambiguity about the separation constraint. Readers may assume that two clades must be separated by two mutations. In the [code](https://github.com/niemasd/multi-introduction/blob/main/FAVITES-COVID-Lite/scripts/stableCoalescence_cladeAnalysis.py) the clades must separated by two or more mutations. For completeness, results for both constraints are presented here. The stricter constraint reduces the single introduction likelihoods: $P(\tau_{1c})$ from 0.2 to 0.1; $P(\tau_{2c})$ from 3.3 to 2.5.


Code and instructions for reproducing these results are available at [this branch](https://github.com/nizzaneela/multi-introduction/tree/corrected_with_relative_size_and_separation_conditions) of the authors' repository.


**Two or more mutation separation and recCA rooting:**

| $t_1$ \ $t_2$ | 0 days | 1 day | 2 days | 4 days | 7 days | 14 days | 28 days |
|---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| **0 days** |0.21 | 0.22 | 0.23 | 0.23 | 0.23 | 0.21 | 0.16 |
| **1 day** |0.23 | 0.24 | 0.25 | 0.25 | 0.25 | 0.21 | 0.17 |
| **2 days** |0.24 | 0.24 | 0.24 | 0.25 | 0.24 | 0.22 | 0.17 |
| **4 days** |0.25 | 0.26 | 0.25 | 0.26 | 0.26 | 0.22 | 0.16 |
| **7 days** |0.26 | 0.27 | 0.26 | 0.26 | 0.26 | 0.21 | 0.16 |
| **14 days** |0.27 | 0.27 | 0.26 | 0.26 | 0.24 | 0.21 | 0.16 |
| **28 days** |0.24 | 0.24 | 0.24 | 0.24 | 0.23 | 0.19 | 0.14 |

**Two or more mutation separation and unconstrained rooting:**
| $t_1$ \ $t_2$ | 0 days | 1 day | 2 days | 4 days | 7 days | 14 days | 28 days |
|---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| **0 days** |0.21 | 0.21 | 0.22 | 0.23 | 0.22 | 0.20 | 0.16 |
| **1 day** |0.22 | 0.23 | 0.24 | 0.24 | 0.24 | 0.20 | 0.16 |
| **2 days** |0.23 | 0.23 | 0.23 | 0.24 | 0.23 | 0.21 | 0.16 |
| **4 days** |0.24 | 0.25 | 0.24 | 0.25 | 0.25 | 0.21 | 0.15 |
| **7 days** |0.25 | 0.26 | 0.25 | 0.24 | 0.25 | 0.20 | 0.15 |
| **14 days** |0.25 | 0.25 | 0.24 | 0.24 | 0.22 | 0.19 | 0.14 |
| **28 days** |0.22 | 0.22 | 0.22 | 0.22 | 0.20 | 0.17 | 0.13 |

**Exactly two mutation separation and recCA rooting:**
| $t_1$ \ $t_2$ | 0 days | 1 day | 2 days | 4 days | 7 days | 14 days | 28 days |
|---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| **0 days** |0.16 | 0.16 | 0.17 | 0.17 | 0.16 | 0.13 | 0.10 |
| **1 day** |0.17 | 0.17 | 0.18 | 0.17 | 0.17 | 0.13 | 0.10 |
| **2 days** |0.18 | 0.17 | 0.17 | 0.17 | 0.16 | 0.13 | 0.09 |
| **4 days** |0.17 | 0.17 | 0.17 | 0.17 | 0.16 | 0.13 | 0.08 |
| **7 days** |0.16 | 0.17 | 0.16 | 0.15 | 0.15 | 0.11 | 0.07 |
| **14 days** |0.14 | 0.14 | 0.13 | 0.12 | 0.11 | 0.09 | 0.06 |
| **28 days** |0.09 | 0.09 | 0.09 | 0.09 | 0.08 | 0.06 | 0.04 |

**Exactly two mutation separation and unconstrained rooting:**
| $t_1$ \ $t_2$ | 0 days | 1 day | 2 days | 4 days | 7 days | 14 days | 28 days |
|---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| **0 days** |0.16 | 0.16 | 0.16 | 0.16 | 0.15 | 0.13 | 0.09 |
| **1 day** |0.17 | 0.17 | 0.18 | 0.17 | 0.16 | 0.13 | 0.10 |
| **2 days** |0.17 | 0.17 | 0.17 | 0.17 | 0.15 | 0.13 | 0.09 |
| **4 days** |0.17 | 0.17 | 0.16 | 0.17 | 0.16 | 0.13 | 0.08 |
| **7 days** |0.16 | 0.17 | 0.16 | 0.15 | 0.14 | 0.11 | 0.07 |
| **14 days** |0.14 | 0.14 | 0.13 | 0.12 | 0.11 | 0.09 | 0.06 |
| **28 days** |0.09 | 0.09 | 0.09 | 0.08 | 0.07 | 0.06 | 0.04 |


Irresepctive of the constraints and parameters, the Bayes factors provide at least moderate support (($<0.3$) for the single introduction hypothesis.

