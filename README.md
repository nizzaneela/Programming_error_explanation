


Two successful introductions have a chance of satisfying the two mutation separation constraint as:
- a single clade descending on a two-mutation branch from a basal clade (here denoted $(\tau_p,\tau_p)_{1c}$); and
- two clades descending on one-mutation branches from the MRCA (here denoted $(\tau_p,\tau_p)_{2c}$).

These correspond, respectively, to the single introduction topologies $\tau_{1c}$ and $\tau_{2c}$, and have the same compatibility:
- $\tau = \tau_{1c}$ is compatible with $S_{MRCA} \in \{S_{A}, S_{B}\}$
- $\tau = \tau_{2c}$ is compatible with $S_{MRCA} \in \{S_{C/C}, S_{T/T}\}$
- $\tau = (\tau_p,\tau_p)\_{1c}$ is compatible with $S_{MRCA} \in \{S_{A}, S_{B}\}$
- $\tau = (\tau_p,\tau_p)\_{2c}$ is compatible with $S_{MRCA} \in \{S_{C/C}, S_{T/T}\}$

The compatibility vectors $P(S_{MRCA}|\tau)$ for the relevant topologies then become:
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

The likelihoods $P((\tau_p,\tau_p)\_{1c}|I_2)$ and $P(\tau_p,\tau_p)\_{2c}|I_2)$ can be measured by drawing pairs of successful simulations and counting the frequency with which they produce basal polytomies, satisfy the relative size constraint and conform to the relevant topology. The topology is determined by the combined mutations between the MRCA and each introduction, and between each introduction and subsequent clade root. The number of mutations between each introduction and subsequent clade root can be sampled using the molecular clock and the time between each introduction and subsequent clade root, which is provided by the simulations. The numbers of mutations between the MRCA and each introduction can be sampled using the molecular clock and the times between the MRCA and each introduction, where:
- the time between the MRCA and the first introduction is sampled from an exponential distribution with expected value $t_1$; and
- the time between the first and second introductions is sampled from an exponential distribution with expected value $t_2$.

This model allows exploration across different values for the parameters $t_1$ (representing the upstream effective population size) and $t_2$ (representing the introduction intensity).

The time between the first and second introductions can also be incorporated into the relative size test.






can be generated using the corrected phylogenies from #11 and the molecular clock used to simulate mutations over the rest of the phylogenies. The [results](https://github.com/nizzaneela/multi-introduction/blob/corrected_with_relative_size_and_separation_conditions/notebooks/cladeAnalysis.ipynb) are roughly geometrically distributed: 50.5% have no mutations, 26.0% have one mutation, 12.4% have two mutations and 11.0% have three or more mutations.

Mutations between the MRCA and each introduction are not mentioned by the authors. In particular, there is no information in the paper that could inform a strong prior for these mutations. Absent such information, a model with minimal assumptions at least allow possibilities to be explored across a range of parameters that a physically meaningful.

The time between each introduction and subsequent clade root is usually the sum of a small number of Poisson processes. The simulated mutations are also a Poisson process and have a similar expected waiting time (~13 days). Therefore, the number of mutations between each introduction and subsequent clade root should be roughly geometrically distributed, with an expected value around 1. [Results](https://github.com/nizzaneela/multi-introduction/blob/corrected_with_relative_size_and_separation_conditions/notebooks/cladeAnalysis.ipynb) from simulating mutations between each introduction and subsequent clade root in the corrected phylogenies from #11 agree; 50.5% have no mutations, 26.0% have one mutation, 12.4% have two mutations and 11.0% have three or more mutations.

The authors do not mention the number of mutations between the MRCA and each introduction. However, it cannot be assumed that the number of mutations between the MRCA and an introduction always combines with the number of mutations between the introduction and clade root and the number of mutations between the MRCA and the other clade root in a way that somehow ensures two clades always satisfy the two mutation separation constraint. There may be historical information or expert opinion that can inform strong priors for these mutations, but it is not in the paper and beyond the scope of this comment. Absent such information, a model is necessary.  


The average time between each introduction and subsequent clade root

### recCA - relaxed
|  |$t_2$ = 0 days | $t_2$ = 1 day | $t_2$ = 2 days | $t_2$ = 4 days | $t_2$ = 7 days | $t_2$ = 14 days | $t_2$ = 28 days |
|---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| $t_1$ **= 0 days** |0.21 | 0.22 | 0.23 | 0.23 | 0.23 | 0.21 | 0.16 |
| $t_1$ **= 1 day** |0.23 | 0.24 | 0.25 | 0.25 | 0.25 | 0.21 | 0.17 |
| $t_1$ **= 2 days** |0.24 | 0.24 | 0.24 | 0.25 | 0.24 | 0.22 | 0.17 |
| $t_1$ **= 4 days** |0.25 | 0.26 | 0.25 | 0.26 | 0.26 | 0.22 | 0.16 |
| $t_1$ **= 7 days** |0.26 | 0.27 | 0.26 | 0.26 | 0.26 | 0.21 | 0.16 |
| $t_1$ **= 14 days** |0.27 | 0.27 | 0.26 | 0.26 | 0.24 | 0.21 | 0.16 |
| $t_1$ **= 28 days** |0.24 | 0.24 | 0.24 | 0.24 | 0.23 | 0.19 | 0.14 |

### unconstrained - relaxed
|  |$t_2$ = 0 days | $t_2$ = 1 day | $t_2$ = 2 days | $t_2$ = 4 days | $t_2$ = 7 days | $t_2$ = 14 days | $t_2$ = 28 days |
|---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| $t_1$ **= 0 days** |0.21 | 0.21 | 0.22 | 0.23 | 0.22 | 0.20 | 0.16 |
| $t_1$ **= 1 day** |0.22 | 0.23 | 0.24 | 0.24 | 0.24 | 0.20 | 0.16 |
| $t_1$ **= 2 days** |0.23 | 0.23 | 0.23 | 0.24 | 0.23 | 0.21 | 0.16 |
| $t_1$ **= 4 days** |0.24 | 0.25 | 0.24 | 0.25 | 0.25 | 0.21 | 0.15 |
| $t_1$ **= 7 days** |0.25 | 0.26 | 0.25 | 0.24 | 0.25 | 0.20 | 0.15 |
| $t_1$ **= 14 days** |0.25 | 0.25 | 0.24 | 0.24 | 0.22 | 0.19 | 0.14 |
| $t_1$ **= 28 days** |0.22 | 0.22 | 0.22 | 0.22 | 0.20 | 0.17 | 0.13 |

### recCA - strict
|  |$t_2$ = 0 days | $t_2$ = 1 day | $t_2$ = 2 days | $t_2$ = 4 days | $t_2$ = 7 days | $t_2$ = 14 days | $t_2$ = 28 days |
|---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| $t_1$ **= 0 days** |0.16 | 0.16 | 0.17 | 0.17 | 0.16 | 0.13 | 0.10 |
| $t_1$ **= 1 day** |0.17 | 0.17 | 0.18 | 0.17 | 0.17 | 0.13 | 0.10 |
| $t_1$ **= 2 days** |0.18 | 0.17 | 0.17 | 0.17 | 0.16 | 0.13 | 0.09 |
| $t_1$ **= 4 days** |0.17 | 0.17 | 0.17 | 0.17 | 0.16 | 0.13 | 0.08 |
| $t_1$ **= 7 days** |0.16 | 0.17 | 0.16 | 0.15 | 0.15 | 0.11 | 0.07 |
| $t_1$ **= 14 days** |0.14 | 0.14 | 0.13 | 0.12 | 0.11 | 0.09 | 0.06 |
| $t_1$ **= 28 days** |0.09 | 0.09 | 0.09 | 0.09 | 0.08 | 0.06 | 0.04 |

### unconstrained - strict
|  |$t_2$ = 0 days | $t_2$ = 1 day | $t_2$ = 2 days | $t_2$ = 4 days | $t_2$ = 7 days | $t_2$ = 14 days | $t_2$ = 28 days |
|---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| $t_1$ **= 0 days** |0.16 | 0.16 | 0.16 | 0.16 | 0.15 | 0.13 | 0.09 |
| $t_1$ **= 1 day** |0.17 | 0.17 | 0.18 | 0.17 | 0.16 | 0.13 | 0.10 |
| $t_1$ **= 2 days** |0.17 | 0.17 | 0.17 | 0.17 | 0.15 | 0.13 | 0.09 |
| $t_1$ **= 4 days** |0.17 | 0.17 | 0.16 | 0.17 | 0.16 | 0.13 | 0.08 |
| $t_1$ **= 7 days** |0.16 | 0.17 | 0.16 | 0.15 | 0.14 | 0.11 | 0.07 |
| $t_1$ **= 14 days** |0.14 | 0.14 | 0.13 | 0.12 | 0.11 | 0.09 | 0.06 |
| $t_1$ **= 28 days** |0.09 | 0.09 | 0.09 | 0.08 | 0.07 | 0.06 | 0.04 |





The published results look suspicious to me, but not egregious. I'll get p-values for them eventually and doubt they'll be extreme. I assume UCSD investigated the errors that were corrected, and didn't find anything. Maybe FOIA could turn up more, but 

The derivation of the Bayes factors is described over pages 11 to 14 of the [Supplementary Materials](https://www.science.org/doi/suppl/10.1126/science.abp8337/suppl_file/science.abp8337_sm.v2.pdf). The approach is set out on page 11:

![](approach.png)

That is, the Bayes factors are the ratio of the likelihood of the sequence data $\bf{Y}$ arising from two introductions $I_2$ vs the likelihood of that sequence data $\bf{Y}$ arising from one introduction $I_1$:   

$$
BF = \frac{P(\mathbf{Y}|I_2)}{P(\mathbf{Y} |I_1)}
$$

Page 14 describes how the Bayes factor is taken from the combination of the posterior and prior odds:

$$
BF = \frac{P(I_2|\mathbf{Y})}{P(I_1|\mathbf{Y})} \cdot \frac{P(I_1)}{P(I_2)}
$$

Page 12 rewrites the posterior probabilities as joint probabilities and breaks them down over the relevant MRCA haplotypes $S_{MRCA}$:

$$
BF = \frac{\sum_{S_{MRCA}} P(S_{MRCA} | Y)P(I_2, S_{MRCA})}{\sum_{S_{MRCA}} P(S_{MRCA} | Y)P(I_1, S_{MRCA})} \cdot \frac{P(I_1)}{P(I_2)}
$$

Page 13 breaks down the joint probabilities $P(I_2, S_{MRCA})$ and $P(I_1, S_{MRCA})$ over the relevant topologies $\tau$:

$$
BF = \frac{\sum_{S_{MRCA}} \left( P(S_{MRCA} | Y)\left( \sum_{\tau} P(S_{MRCA} | \tau)P(\tau | I_2)P(I_2) \right) \right)}
{\sum_{S_{MRCA}}\left(  P(S_{MRCA} | Y)\left( \sum_{\tau} P(S_{MRCA} | \tau)P(\tau | I_1)P(I_1) \right) \right)} 
\cdot \frac{P(I_1)}{P(I_2)}
$$

Page 13 then describes the compatibility between the different topologies and MRCA haploytpes:

![](compatibility.png)

As noted in #7, the topology $(\tau_p, \tau_p)$ encompasses the sub-topologies $(\tau_p, \tau_{1c})$, $(\tau_{1c}, \tau_p)$ and $(\tau_{1c}, \tau_{1c})$. Therefore, it is sufficient to state:
- $\tau = \tau_{2c}$ is compatible with $S_{MRCA} \in \{S_{C/C}, S_{T/T}\}$
- $\tau = \tau_{1c}$ is compatible with $S_{MRCA} \in \{S_{A}, S_{B}\}$
- $\tau = (\tau_{P}, \tau_{P})$ is compatible with $S_{MRCA} \in \{S_{A}, S_{B}, S_{C/C}, S_{T/T}\}$

Page 13 then describes how the compatbility statements define $P(S_{MRCA}|\tau)$:

![](renorm.png)

This description is wrong. The renormalization is applied to the vector of multiple MRCA haplotypes $(A, B, C/C, T/T)$ that are compatible with a topology. The renormalised vectors for the three relevant topologies are:
- $P(S_{MRCA}|\tau = \tau_{2c}) = (0, 0, 0.5, 0.5)$
- $P(S_{MRCA}|\tau = \tau_{1c}) = (0.5, 0.5, 0, 0)$
- $P(S_{MRCA}|\tau = (\tau_p,\tau_p)) = (0.25, 0.25, 0.25, 0.25)$

These allow the equation to be expanded out, with the priors $P(I_2)/P(I_1)$ cancelled by their inverse $P(I_1)/P(I_2)$:

$$
BF = 
\frac{0.25 \cdot P((\tau_P,\tau_P)|I_2) \cdot \left( P(S_A|Y) + P(S_B|Y) + P(S_{C/C}|Y) + P(S_{T/T}|Y) \right)} 
{ 0.5 \cdot P(\tau_{1c}|I_1) \cdot \left( P(S_A|Y) + P(S_B|Y) \right) + 0.5 \cdot P(\tau_{2c}|I_1)\cdot \left( P(S_{C/C}|Y) + P(S_{T/T}|Y) \right)}
$$

This is mathematically identical to the authors' equations, but simpler because it avoids the superfluous marginalization over the sub-topologies $(\tau_p, \tau_{1c})$, $(\tau_{1c}, \tau_p)$ and $(\tau_{1c}, \tau_{1c})$.

The authors inferred the posterior probabilities $P(S_A|Y)$, $P(S_B|Y)$, $P(S_{C/C}|Y)$ and $P(S_{T/T}|Y)$ using the phylodynamic software [BEAST](https://beast.community/), producing results shown in [Table 1](https://www.science.org/doi/10.1126/science.abp8337#T1). The results for the recCA and unconstrained rootings produce the equations:

$$
BF_{uncon.} = \frac{P((\tau_P,\tau_P)|I_2)} 
{ 1.76 \cdot P(\tau_{1c}|I_1) + 0.24 P(\tau_{2c}|I_1)},
\quad
BF_{recCA} = \frac{P((\tau_P,\tau_P)|I_2)} 
{ 1.72 \cdot P(\tau_{1c}|I_1) + 0.28 P(\tau_{2c}|I_1)}
$$

The authors measured the likelihoods $P(\tau_P|I_1)$, $P(\tau_{1c}|I_1)$ and $P(\tau_{2c}|I_1)$ from the frequencies with which successful simulated introductions produced the topologies $\tau_p$, $\tau_{1c}$ and $\tau_{2c}$. 

 $\tau_p$ is defined on page 10:

![](tau_p.png)

$P((\tau_P,\tau_P)|I_2)$ was then calculated as $P(\tau_P|I_1)^2$. Therefore, $P((\tau_p, \tau_p)|I_2)$ is the likelihood of two successful introductions producing two clades where there is a polytomy at the root of each clade.

$\tau_{1c}$ and  $\tau_{2c}$ are defined on page 11:

![](tau_12c.png)

Therefore, $P(\tau_{1c}|I_1)$ and $P(\tau_{2c}|I_1)$ are the likelihoods of a single successful introduction producing two clades where there is a polytomy at the root of each clade, and
1. each clade includes 30-70% of the total taxa (taxa being sampled from the first 50,000 infections); and
2. the clades are separated by two mutations (with one being basal for $P(\tau_{1c}|I_1)$ and both being derived for $P(\tau_{2c}|I_1)$ ).

The additional conditions 1 and 2 reduce the single introduction likelihoods. This makes the published Bayes factors (4.2 and 4.3) meaningless, since they may be caused by this reduction rather than a higher plausiblity of the two introduction hypothesis.

This can be corrected by applying conditions 1 and 2 to the two introduction likelihood $P((\tau_P,\tau_P)|I_2)$.

This comment examines the effects of condition 1. A subsequent comment will examine condition 2.

Two simulations can be be tested against condition 1 by going through the first 50,000 infections amongst both, and comparing the number of samples from each, e.g.:
```
def test_sizes(run_1,run_2):
    # read samples into seperate sets for each run
    sample_dict = {run_1: set(), run_2: set()}
    # read transmissions into one list for both runs
    transmission_list = []
    for run in [run_1, run_2]:
        with open(f'simulations/{run:04d}/subsample_times.corrected.txt', 'r') as file:
            for line in file:
                u, t = line.strip().split('\t')
                sample_dict[run].add(u)
        with open(f'simulations/{run:04d}/transmission_network.subsampled.corrected.txt', 'r') as file:
            for line in file:
                u,v,t = line.strip().split()
                if u != v:
                    transmission_list.append((v,float(t),run))
    # sort the transmissions in order of time
    transmission_list = sorted(transmission_list, key=lambda x: x[1])
    # start sample counts from -1 because the primary case is not to be included 
    n_valid_samples = {run_1: -1, run_2: -1}
    # go through the first 50,000 transmissions and count those that are sampled from each run
    for event_no, transmission in enumerate(transmission_list):
        if event_no == 50000:
            break
        if transmission[0] in sample_dict[transmission[2]]:
            n_valid_samples[transmission[2]] += 1
    # compare the numbers of sampled transmissions to the relative size condition 
    if 0.3 <= n_valid_samples[run_1]/(n_valid_samples[run_1] + n_valid_samples[run_2]) <= 0.7:
        return True
    else:
        return False
```

The likelihood of two successful introductions producing basal polytomies and satisfying condition 1 can be measured by repeatedly drawing random pairs of simulations and testing them against both conditions, e.g.:
```
# sample 1100 topologies, as for the others
for i in range(1100):
    # draw two simulations
    run_1, run_2 = np.random.choice(range(1, 1101), 2, replace=True)
    # get the clade sizes
    clade_sizes_1 = clade_analyses_CC_d[run_1]['clade_sizes'] 
    clade_sizes_2 = clade_analyses_CC_d[run_2]['clade_sizes']
    # test them for basal polytomies
    if len(clade_sizes_1) >= min_polytomy_descendants and len(clade_sizes_2) >= min_polytomy_descendants: 
        # test their sizes
        if test_sizes(run_1, run_2):
            count_atLeastMinDescendants += 1
```

The function `clade_analysis_updated` in the notebook [stableCoalescence_cladeAnalysis.py](https://github.com/nizzaneela/multi-introduction/blob/corrected/FAVITES-COVID-Lite/scripts/stableCoalescence_cladeAnalysis.py) can be adapted to do this instead of counting the number of simulations with basal polytomies. The `function calculate_bf` can then be adapted to use the resulting likelihood as $P((\tau_p,\tau_p)|I_2)$. $(\tau_p,\tau_p)$ encompasses all topologies compatible with the two introduction hypothesis, so the likelihoods of the other two introduction topologies ( $(\tau_p,\tau_{1c})$, $(\tau_{1c},\tau_p)$, $(\tau_{1c},\tau_{1c})$ ) are set to zero to avoid double counting.
```
def calculate_bf(asr_results, simulation_results):
    # Let t_p be a polytomy, t_1C be one clade, and t_2c be two clades. Note that t_p includes t_1c. t_p equals all topologies with a basal polytomy (Fig. 2a). 
    # trees are in the order
    # (t_p, t_1C, t_2C, (t_p,(t_p,t_1C,t_2C)), (t_1C,(t_p,t_1C,t_2C)), (t_2c,(t_p,t_1C,t_2C)))
    ...
    pr_trees_given_I2 = np.concatenate([np.array([0]*3), np.array([simulation_results[0]]), np.array([0]*8)])
```

The notebook can then be rerun to produce the following results.

![](results2.png)

Reducing inconsistency in the treatment of the different hyoptheses reduces the Bayes factors by a factor of six.

Note that the authors measured the two week doubling times of the simulated epidemics at the 1000th infection, and found a 95% highest density interval (HDI) of 1.35 to 5.44 days. This range of early growth rates suggests that two introductions are unlikely to grow to similar sizes within the first 50,000 infections, and the results confirm this.

Code and instructions for reproducing these results are available at [this branch](https://github.com/nizzaneela/multi-introduction/tree/corrected_with_relative_size_condition) of the authors' repository.

#### Limitations
1. `test_sizes` assumes both introductions are simultaneous. This maximises the likelihood of the two introductions having similar sizes;
2. Each introduction has a different delay before the first samples due to the different times of first hospitalization. Realistically, the earlier time should apply to both. Therefore, this analysis misses some early samples in the simulation with the later time of first hospitalization. In some cases, missing these samples will drop branch counts below the polytomy threshold, or clade sizes below the relative size threshold, reducing the likelihood of a basal polytomy and satisfying condition 1;
3. Each introduction has a stable coalescence based on when it reaches 50,000 infections. Assuming the stable coalescence should instead be based on when simulations reach 50,000 combined infections, this analysis removes more basal lineages than it should, increasing the likelihood of basal polytomies.

Overall, the limitations increase the likelihood of two simulations producing two basal polytomies and satisfying conditions 1, so the corrected Bayes factors may be considered a rough upper bound.

It is theoretically possible that applying condition 3 to the two introduction likelihoods could increase the Bayes factors. This would require an increase in the likelihood of two introductions producing topologies compatible with the more probable MRCA haplotypes (i.e. A or B). This is considered unlikely, but will be confirmed either way in the next comment.
