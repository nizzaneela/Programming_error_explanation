Two successful introductions are deemed consistent with the evidence if they produce two clades where:
1. there is a polytomy at the root of each clade.

A single successful introduction is deemed consistent with the evidence if it produces two clades where:
1. there is a polytomy at the root of each clade;
2. each clade includes 30-70% of the total samples; and
3. the clades are separated by two mutations.

This means the published Bayes factors (4.2 and 4.3) cannot support the two introduction hypothesis, because their values could reflect the more stringent conditions applied to single introduction hypothesis, rather than a higher plausiblity of the two introduction hypothesis.

This can be corrected by applying conditions 2 and 3 to two introduction hypothesis.

Two successful introductions can be be tested against condition 2 by going through the first 50,000 infections amongst both, counting the taxa (samples) in each, and comparing the results, e.g.:
```
def test_sizes(run_1,run_2):
    # read samples into seperate sets for each run
    sample_dict = defaultdict(set)
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
    # start sample counts from -1 to compensate for the primary sample
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

The function `clade_analysis_updated` in the notebook [stableCoalescence_cladeAnalysis.py](https://github.com/nizzaneela/multi-introduction/blob/corrected/FAVITES-COVID-Lite/scripts/stableCoalescence_cladeAnalysis.py) can be adapted to do this instead of counting the number of simulations with a basal polytomies. The `function calculate_bf` must then be adapted to receive the likelihood $P(\tau_p,\tau_p|I_2)$ where it previously received $P(\tau_p|I_1)$:

- the compatiblity matix must be reset to ignore the previous marginalization of the $(\tau_p,\tau_p)$ topology, i.e. to reflect $\tau_p$ including $\tau_{1C}$:
```
    # Let t_p be a polytomy, t_1C be one clade, and t_2c be two clades. Note that t_p is a topology with a polytomy including t_1c. t_p equals all topologies with a basal polytomy (Fig. 2a). 
    # trees are in the order
    # (t_p, t_1C, t_2C, (t_p,(t_p,t_1C,t_2C)), (t_1C,(t_p,t_1C,t_2C)), (t_2c,(t_p,t_1C,t_2C)))
    compatibility_matrix = np.array([np.array([0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0]), # S_A
                                     np.array([0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0]), # S_B
                                     np.array([0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0]), # S_CC
                                     np.array([0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0])] # S_TT
```

- the likelihood vector for the two introductions must receive the new likelihood $P(\tau_p,\tau_p|I_2)$ directly and place zeros elsewhere:
```
     pr_trees_given_I2 = np.concatenate([np.array([0]*3), np.array([simulation_results[0]]), np.array([0]*8)])
```

The notebook can then be rerun to produce Bayes factors that at least apply two of the three conditions equally.

#### Limitations
1. `test_sizes` assumes both introductions are simultaneous. This maximises the likelihood of the two introductions having similar sizes;
2. Each introduction has a different delay before the first samples due to the different times of first hospitalization. Realistically, the earlier time should apply to both. Therefore, this analysis misses some early samples in the simulation with the later time of first hospitalization. In some cases, missing these samples will drop branch counts below the polytomy threshold, or clade sizes below the relative size threshold, reducing the likelihood of satisfying conditions 1 and 2;
3. Each introduction has a stable coalescence based on when it reaches 50,000 infections. Realistically, these should be based on when they reach 50,000 combined infections. Therefore, this analysis removes more basal lineages than it should, increasing the likelihood of satisfying condition 1.

Overall, the limitations increase the likelihood of two simulations satisfying conditions 1 and 2, so that the analysis is conservative.
