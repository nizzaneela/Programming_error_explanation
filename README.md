There are discrepencies between the methods described in the text and those carried out by the code.

The effects are difficult to quantify because the stochastic simulations are not reproducible and the results are sensitive to sampling noise. However, resampling the final stochastic phase of the simulations, using corrected code, reduces the Bayes factors of the primary analysis by ~15%.

## Discrepencies

### The stable coalescence

The simulated phylogenies are constructed by coalescing lineages sampled from the first 50,000 simulated cases. Partial and delayed sampling is simulated by sampling only a portion of simulated cases, and by ignoring samples that precede the first simulated hospitalization (page 8 of the [Supplementary Materials](https://www.science.org/doi/suppl/10.1126/science.abp8337/suppl_file/science.abp8337_sm.v2.pdf)).

![Excerpt from page 8 of the Supplementary Materials](https://github.com/nizzaneela/Programming_error_explanation/blob/4b653347fb1b4642c98d82c50fcea29200c4add1/sample.png)

Samples are also ignored if they are on branches ancestral to a stable coalescence (page 10 of the [Supplementary Materials](https://www.science.org/doi/suppl/10.1126/science.abp8337/suppl_file/science.abp8337_sm.v2.pdf)).  

![Excerpt from page 10 of the Supplementary Materials](https://github.com/nizzaneela/Programming_error_explanation/blob/b988d5b5b507d88619c9b9fb9fcaceb5349ff771/sctext.png)

Contrary to the described method, the function `coalescent_timing` in the script [stableCoalescence_cladeAnalysis.py](https://github.com/sars-cov-2-origins/multi-introduction/blob/78ec9e3b90215267b45ed34be2720566b7398b77/FAVITES-COVID-Lite/scripts/stableCoalescence_cladeAnalysis.py) does not stop the tMRCA calculations when 50,000 individuals have been infected. Calculations continue until the last day of the simulation.

```
def coalescent_timing(time_inf_dict, current_inf_dict, total_inf_dict, tree, num_days=100):
    # determine the stable coalescence of the tree
    ...
    times = list(time_inf_dict.keys()) # the number of days in the sim
    ...
    for index, time in enumerate(times):
        if time > num_days: # sometimes gemf goes past the limit but we don't always know when
            break
        labels = time_inf_dict[time] # currently circulating infections
        if time == 0:
            ...        
        else: # get height of the subtree; the tMRCA
            subtree = tree.extract_tree_with(labels)
            height_subtree = distance_from_zero_dict[subtree.root.label]
            heights.append(height_subtree)
    ...
    coalescent_timing_results = [used_times, heights, total_inf, current_inf, current_samples]
    return coalescent_timing_results
```

For example, the first successful simulation run `0001` reaches 50,000 cases on day 39, when the tMRCA of active sampled cases is 0.000333 years (~3 hours). The calculations continue until the end of the simulation, after another 61 days and 1.3 million cases. By then, only 710 sampled cases are still active and their tMRCA is 0.016277 years (~6 days). This can be seen in the file `simulations_01/0001/coalData_parameterized.txt` in the published [data](https://zenodo.org/records/6899613/files/simulations_01.zip?download=1) on [Zenodo](https://zenodo.org/records/6899613). 

```
time	coalescence time	total infected	currently infected	current samples
...
38	0.000333	45444	33079	5731
39	0.000333	53606	38638	6175
...
99	0.016277	1363477	149166	743
100	0.016277	1371985	144107	710
```

The tMRCA from the last day of the simulation is used as the time of stable coalescence. This can be seen by comparing the coalescence times in each simulation's `coalData_parameterized.txt` to those recorded in the [summary data](https://github.com/sars-cov-2-origins/multi-introduction/raw/main/FAVITES-COVID-Lite/cumulative_results/FAVITES_results.zip) on GitHub.
```
Run	Coalescent time	First ascertained	First unascertained	First hospitalized	infections
0001	0.016277	0.003101	0.018548	0.038307	3
...
```

This means the code removes basal lineages that do not have active sampled cases at the end of simulation (day 100), even if they do have active sampled cases at the end of sampling period (infection 50,000th case). This does not agree with the method described in the text

As described in the text, the stable coalescence is based on the MRCA of lineages that have a sampled case that is active when the 50,000th case is infected. Lineages are more likely to have a sampled case that is active at this time if they undergo more early growth, so they are more likely to be associated with an early polytomy. By rooting the tree at the stable coalescence, basal lineages that are less likely to be associated with an early polytomy are removed. If the next derived structure is a polytomy, it is made basal.

As implemented in the code, the stable coalescence is the MRCA of lineages that have a sampled case that is active at the end of the simulation, even though the last sampled case might have been infected months earlier. Lineages are more likely to have a sampled case that is active at this later time if they undergo a lot more early growth, so they are much more likely to be associated with an early polytomy. This means more basal lineages are removed, and more early polytomies are made basal.

Lineages are not removed if they descend from the stable coalescence.

### The primary case state

The epidemic simulation script [FAVITES-COVID-Lite_noSeqgen.py](https://github.com/sars-cov-2-origins/multi-introduction/blob/78ec9e3b90215267b45ed34be2720566b7398b77/FAVITES-COVID-Lite/scripts/FAVITES-COVID-Lite_noSeqgen.py) accepts commands to initialise the primary case in arbitrary states. However, the script is hard-coded to initialise the primary case in an infectious state (`P1`), so that transmissions being immediately.
```
    # write GEMF status file
    out_file = open(out_fn, 'w')
    status_fn = "%s/%s" % (gemf_tmp_dir, GEMF_STATUS_FN)
    status_file = open(status_fn, 'w')
    seeds = set(sample(range(max_node_label+1), k=num_seeds))
    gemf_state = list()
    for node in range(max_node_label+1):
        if node in seeds:
            status_file.write("%d\n" % gemf_state_to_num['P1'])
            gemf_state.append(gemf_state_to_num['P1'])
            out_file.write("None\t%d\t0\n" % node)
        else:
            status_file.write("%d\n" % gemf_state_to_num['S'])
            gemf_state.append(gemf_state_to_num['S'])
    status_file.close()
    print_log("Wrote GEMF '%s' file: %s" % (GEMF_STATUS_FN, status_fn))
```
This is despite the published [command](https://github.com/sars-cov-2-origins/multi-introduction/blob/main/FAVITES-COVID-Lite/commands/command.0.28TF_0.15r.txt) indicating that the primary case should start out in the exposed but non-infectious state (`--tn_freq_e 0.00000020`).
```
~/scripts/FAVITES-COVID-Lite-updated.py --gzip_output --path_ngg_barabasi_albert ngg_barabasi_albert --path_gemf GEMF --path_coatran_constant coatran_constant --path_seqgen seq-gen --cn_n 5000000 --cn_m 8 --tn_s_to_e_seed 0 --tn_e_to_p1 125.862069 --tn_p1_to_p2 999999999 --tn_p2_to_i1 23.804348 --tn_p2_to_a1 134.891304 --tn_i1_to_i2 62.931034 --tn_i1_to_h 0.000000 --tn_i1_to_r 62.931034 --tn_i2_to_h 45.061728 --tn_i2_to_r 0.000000 --tn_a1_to_a2 9999999999 --tn_a2_to_r 125.862069 --tn_h_to_r 12.166667 --tn_s_to_e_by_e 0 --tn_s_to_e_by_p1 0 --tn_s_to_e_by_p2 3.513125 --tn_s_to_e_by_i1 6.387500 --tn_s_to_e_by_i2 6.387500 --tn_s_to_e_by_a1 0 --tn_s_to_e_by_a2 3.513125 --tn_freq_s 0.99999980 --tn_freq_e 0.00000020 --tn_freq_p1 0 --tn_freq_p2 0 --tn_freq_i1 0 --tn_freq_i2 0 --tn_freq_a1 0 --tn_freq_a2 0 --tn_freq_h 0 --tn_freq_r 0 --tn_end_time 0.273973 --tn_num_seeds 1 --pt_eff_pop_size 1 --pm_mut_rate 0.00092 --o 
```

That is, the code forces simulations to skip the latent phase of the first case.

For simulations where the stable coalescence is in the primary case (~20% of simulations), this compresses coalescence events around the stable coalescence, reducing the likelihood of an early mutation breaking up a basal polytomy.

For simulations where the stable coalescence is not in the primary case (~80% of simulations), this pushes the estimated time of introduction forward, and decreases its variance.

### Minor errors

The `main` function in the script [stableCoalescence_cladeAnalysis.py](https://github.com/sars-cov-2-origins/multi-introduction/blob/78ec9e3b90215267b45ed34be2720566b7398b77/FAVITES-COVID-Lite/scripts/stableCoalescence_cladeAnalysis.py) restores basal lineages if their MRCA is sufficiently close to the time of stable coalescence.
```
# main function
    ...
    coal_timing = coalescent_timing(time_inf_dict, current_inf_dict, total_inf_dict, subtree, args.num_days)

    # prepare for clade analysis; get the subtree with the stable coalescence (MRCA) root
    eps = 1e-8
    stable_coalescence = coal_timing[1][-1]
    subtree_sc_leaves = []
    for n in subtree.distances_from_root():
        if abs(n[1] - stable_coalescence) < eps:
            # print(n[0].label)
            subtree_sc_leaves += [n.label for n in subtree.extract_subtree(n[0]).traverse_leaves()]
    subtree_sc_leaves = set(subtree_sc_leaves)
    subtree_sc = tree.extract_tree_with(subtree_sc_leaves)
```
This can increase the size of basal polytomies, but only in rare cases where coalescence events are compressed closely enough around the stable coalescence (i.e. < 0.2% of simulations).
In one instance (simulation `0823`) this occured in the primary case, where the compression was amplified by the elision of the latent phase.

Clumsy removal of the primary case sample also removes leaves that include the primary case node number in their label string. (The primary case is always sampled to ensure CoaTran roots the phylogeny in the primary case, and is removed for clade analysis.)
```
    # get subtree that excludes index case
    tree = treeswift.read_tree_newick(args.time_tree)
    subtree_leaves = []
    for n in tree.traverse_leaves():
        if index_case not in n.label:
            subtree_leaves.append(n.label)
    subtree = tree.extract_tree_with(subtree_leaves)
    subtree.suppress_unifurcations()
```
This is a very minor error, mentioned here only because it was corrected in subsequent reanalyses.

## Noise

The Bayes factors are described on pages 11 to 14 of the [Supplementary Materials](https://www.science.org/doi/suppl/10.1126/science.abp8337/suppl_file/science.abp8337_sm.v2.pdf). Combined and expanded, their equations can be written as:

$$
BF = 
\frac{0.25 \cdot P(\tau_P|I_1)^2 \left( P(S_A|Y) + P(S_B|Y) + P(S_{C/C}|Y) + P(S_{T/T}|Y) \right)} 
{0.5 \cdot P(\tau_{c1}|I_1) \left( P(S_A|Y) + P(S_B|Y) \right) + 0.5 \cdot P(\tau_{c2}|I_1) \left( P(S_{C/C}|Y) + P(S_{T/T}|Y) \right)}
$$

where:
- $P(\tau_P|I_1)$, $P(\tau_{1C}|I_1$ and $P(\tau_{2C}|I_1)$ are the likelihoods of the different topologies (taken from the simulations and shown in Fig. 2 and Table S5);
- $P(S_A|Y)$, $P(S_B|Y)$, $P(S_{C/C}|Y)$ and $P(S_{T/T}|Y)$ are the posterior probabilities of the different MRCA haplotypes (taken from BEAST and shown in  Table 1); and
- $0.25$ and $0.5$ are the coefficients of the normalized vectors that distribute the topology likelihoods across the posterior probabilities of compatible MRCA haplotypes.

The data includes RNG seeds for the epidemic simulations, but not for the sample times, coalescence, or mutation simulations. Therefore, the simulations cannot be reproduced. However, assuming the published likelihoods are sufficiently accurate, the results of the 1100 simulations can be replicated by sampling appropriate binomial distributions ($\tau_{2C}$ is neglected here because its measured likelihood was zero).
```
import numpy as np
np.random.seed(42)
def sample_likelihoods():
    p_tau_p_given_i1 = 0.475 # from Fig. 2
    p_tau_1c_given_i1 = 0.031 # from Fig. 2
    # sample the number of simulations that have basal polytomies
    n_tau_p = np.random.binomial(1100,p_tau_p_given_i1)
    # sample the number of basal polytomies that have another on a two mutation branch
    p_tau_1c_given_tau_p = p_tau_1c_given_i1/p_tau_p_given_i1 
    n_tau_1c = np.random.binomial(n_tau_p,p_tau_1c_given_tau_p)
    # return the likelihoods
    return n_tau_p/1100, n_tau_1c/1100
```
The Bayes factor equation can be written in terms of the likelihoods and posterior probabilities.
```
def compute_bfs(p_tau_p_given_i1, p_tau_1c_given_i1, posteriors):
     bf = 0.25*p_tau_p_given_i1**2*sum(posteriors)/(0.5*p_tau_1c_given_i1*sum(posteriors[:2]))
     return bf
```
Repeatedly resampling the likelihoods and computing the Bayes factors then provides a distribution of results that would be expected from replicating the analysis.
```
unconstrained_posteriors = np.array([1.68, 80.85, 10.32, 0.92])/100 # from Table 1
results = []
i in range(20000): # sample 20000 times
    p_tau_p_given_i1, p_tau_1c_given_i1 = sample_likelihoods()
    results.append(compute_bfs(p_tau_p_given_i1, p_tau_1c_given_i1, unconstrained_posteriors))
    results.sort()
    print(f'95% CDI of Bayes factors with unconstrained rooting: {results[500]:.1f}, {results[19500]:.1f}')
95% CDI of Bayes factors with unconstrained rooting: 3.1, 6.0
```
The central 95% of the distribution has a range similar to the size of the Bayes factors. The assumption of sufficient accuracy is invalid.

## Remedies

The epidemic simulation script [FAVITES-COVID-Lite_noSeqgen.py](https://github.com/sars-cov-2-origins/multi-introduction/blob/78ec9e3b90215267b45ed34be2720566b7398b77/FAVITES-COVID-Lite/scripts/FAVITES-COVID-Lite_noSeqgen.py) can be corrected to initialise the primary case in the exposed compartment (`E`). However, rather than repeating the epidemic simulations, the published data stored on  [Zenodo](https://zenodo.org/records/6899613) can be reused by sampling times for the latent periods and adding them to the existing transmission and sampling times. [CoaTran](https://github.com/niemasd/CoaTran) can then be rerun to obtain corrected phylogenies.
```
import os
import numpy as np
from gzip import open as gopen
np.random.seed(42)

# go through each of the 1100 simulations
for i in range(1,1101):
    # generate a latent period
    latent_period = np.random.exponential(2.9/365) #expected value is 2.9/365 days (Table S3)
    # set up file paths, open files and write out the corrected data for transmission network
    old_tn_path = f'simulations/{i:04d}/transmission_network.subsampled.txt.gz'
    new_tn_path = f'simulations/{i:04d}/transmission_network.subsampled.corrected.txt'
    with gopen(old_tn_path, 'rt') as old_tn, open(new_tn_path, 'w') as new_tn:
        for line in old_tn:
            u, v, t = line.strip().split('\t')
            new_tn.write(f'{u}\t{v}\t{float(t) + latent_period:.6f}\n')
    # set up file paths, open files and write out the corrected data for sample times
    old_samples_path = f'simulations/{i:04d}/subsample_times.txt.gz'
    new_samples_path = f'simulations/{i:04d}/subsample_times.corrected.txt'
    with gopen(old_samples_path, 'rt') as old_samples, open(new_samples_path, 'w') as new_samples:
        for line in old_samples:
            u, t = line.strip().split('\t')
            new_samples.write(f'{u}\t{float(t) + latent_period}\n')
    # setup and run CoaTRan to generate a corrected time tree
    command = ['coatran_constant', new_tn_path, new_samples_path, '1']
    env = os.environ.copy()
    coatran_rng_seed = np.random.randint(low=0, high=2**31) # ensure it is reproducible
    env["COATRAN_RNG_SEED"] = str(coatran_rng_seed)
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=env)
    tree_path = f'simulations/{i:04d}/tree.time.subsampled.corrected.nwk'
    with open(tree_path, 'w') as file_handle:
        file_handle.write(result.stdout.decode())
```

The script [stableCoalescence_cladeAnalysis.py](https://github.com/sars-cov-2-origins/multi-introduction/blob/78ec9e3b90215267b45ed34be2720566b7398b77/FAVITES-COVID-Lite/scripts/stableCoalescence_cladeAnalysis.py) can be corrected to determine the stable coalescence properly by:
- breaking the loop in the function `coalescent_timing` once 50,000 individuals have been infected, walking back to find the first day when tMRCA of active sampled cases will jump forward by less than one day, and returning the MRCA of that day as the stable coalescence; and

![](https://github.com/nizzaneela/Programming_error_explanation/blob/f08072e335f3d6988b98300728e831590cdc275d/get_Sc.png)

- modifying the `main` function to extract the subtree rooted at the stable coalescence.

![](https://github.com/nizzaneela/Programming_error_explanation/blob/15248e0d2472ea3dff7a9c6da540f7a58a672cb8/stable_coalescence.png)

Correcting the removal of the primary case sample is not really necessary, but simple enough.

![](https://github.com/nizzaneela/Programming_error_explanation/blob/ff3cab7fca370a03961a6536ffe8740146eab014/primary_removel.png)

Complete code and instructions for reproducibly generating corrected time trees is published in [this branch](https://github.com/nizzaneela/multi-introduction/tree/corrected) of the authors' repository. In order to suppress sampling noise, the code also automates resampling of the mutation simulations and re-analysis of the resulting clades, using the correct stable coalescence, 1000 times. 

![](https://github.com/nizzaneela/Programming_error_explanation/blob/ff3cab7fca370a03961a6536ffe8740146eab014/results.png)

The corrections and resampling reduce the Bayes factors by ~15%.
