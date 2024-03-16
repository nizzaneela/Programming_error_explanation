The Bayes factors weigh the single introduction likelihoods according to the posterior probabilites of their compatible MRCA haplotypes. Defining weights $\alpha_{1c}$ and $\alpha_{2c}$ as $\frac{2 \cdot (P(S_A|Y) + P(S_B|Y))}{P(S_A|Y) + P(S_B|Y) + P(S_{C/C}|Y) + P(S_{T/T}|Y)}$ and $\frac{2 \cdot (P(S_{C/C}|Y) + P(S_{T/T}|Y)}{P(S_A|Y) + P(S_B|Y) + P(S_{C/C}|Y) + P(S_{T/T}|Y)}$ , respectively, the Bayes factor to be written as

$$
BF = 
\frac{P((\tau_P,\tau_P)|I_2)} 
{ \alpha_{1c} \cdot P(\tau_{1c}|I_1) + \alpha_{2c} \cdot P(\tau_{2c}|I_1)}
$$

$P((\tau_P,\tau_P)|I_2)$ is the likelihood of two successful introductions producing two clades where:
1. there is a polytomy at the root of each clade.

$P(\tau_{1c}|I_1)$ and $P(\tau_{2c}|I_1)$ are the likelihoods of a single successful introduction producing two clades where:
1. there is a polytomy at the root of each clade;
2. each clade includes 30-70% of the total taxa (taxa being sampled from the first 50,000 infections); and
3. the clades are separated by two mutations (with one being basal for $P(\tau_{1c}|I_1)$ and both being derived for $P(\tau_{2c}|I_1)$).

The published Bayes factors (4.2 and 4.3) might reflect a higher plausiblity of the two introduction hypothesis, or the additional conditions 2 and 3 applied to single introduction likelihoods. This makes the Bayes factors meaningless.

This can be corrected by applying conditions 2 and 3 to the two introduction likelihoods.

Two successful introductions can be be tested against condition 2 by going through the first 50,000 infections amongst both, and comparing the number taxa (samples) from each, e.g.:
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

The likelihood of two successful introductions satisfying conditions 1 and 2 can be measured by repeatedly drawing random pairs of simulations and testing them against both conditions, e.g.:
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

The function `clade_analysis_updated` in the notebook [stableCoalescence_cladeAnalysis.py](https://github.com/nizzaneela/multi-introduction/blob/corrected/FAVITES-COVID-Lite/scripts/stableCoalescence_cladeAnalysis.py) can be adapted to do this instead of counting the number of simulations with basal polytomies. The `function calculate_bf` can then be adapted to use the resulting likelihood as $P((\tau_p,\tau_p)|I_2)$, where $(\tau_p,\tau_p)$ encompasses all topologies compatible with the two introduction hypothesis, so that the other two introduction likelihoods are zero (as in #7).
```
def calculate_bf(asr_results, simulation_results):
    # Let t_p be a polytomy, t_1C be one clade, and t_2c be two clades. Note that t_p includes t_1c. t_p equals all topologies with a basal polytomy (Fig. 2a). 
    # trees are in the order
    # (t_p, t_1C, t_2C, (t_p,(t_p,t_1C,t_2C)), (t_1C,(t_p,t_1C,t_2C)), (t_2c,(t_p,t_1C,t_2C)))
    ...
    pr_trees_given_I2 = np.concatenate([np.array([0]*3), np.array([simulation_results[0]]), np.array([0]*8)])
```

The compatibility matrix can also be updated to negate compatiblity with the topologies that are made redundant, but this isn't really necessary because their likelihoods are zero.
```
    compatibility_matrix = np.array([np.array([0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0]), # S_A
                                     np.array([0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0]), # S_B
                                     np.array([0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0]), # S_CC
                                     np.array([0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0])] # S_TT
                                    )
```

The notebook can then be rerun to produce the following results.

![](results2.png)

The paper measures the doubling times of the simulated epidemics at the 1000th infection, and found a 95% highest density interval (HDI) of 1.35 to 5.44 days. This range of early growth rates suggests that two introductions are unlikely to grow to similar sizes within the first 50,000 infections, and the results confirm this.

Code and instructions for reproducing these results are available at [this branch](https://github.com/nizzaneela/multi-introduction/tree/corrected_with_relative_size_condition) of the author's repository.

#### Limitations
1. `test_sizes` assumes both introductions are simultaneous. This maximises the likelihood of the two introductions having similar sizes;
2. Each introduction has a different delay before the first samples due to the different times of first hospitalization. Realistically, the earlier time should apply to both. Therefore, this analysis misses some early samples in the simulation with the later time of first hospitalization. In some cases, missing these samples will drop branch counts below the polytomy threshold, or clade sizes below the relative size threshold, reducing the likelihood of satisfying conditions 1 and 2;
3. Each introduction has a stable coalescence based on when it reaches 50,000 infections. Realistically, these should be based on when they reach 50,000 combined infections. Therefore, this analysis removes more basal lineages than it should, increasing the likelihood of satisfying condition 1.

Overall, the limitations increase the likelihood of two simulations satisfying conditions 1 and 2, so the resulting Bayes factors may be considered a  rough upper bound.

