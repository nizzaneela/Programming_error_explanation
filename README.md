Errors and noise inflate the corrected Bayes factors.

The [Supplementary Materials](https://www.science.org/doi/suppl/10.1126/science.abp8337/suppl_file/science.abp8337_sm.v2.pdf) describe how short-lived basal lineages are excluded from the clade analysis by rooting each phylogeny at a stable coalescence. The code does not implement the described method. 

The effects of this discrepency are difficult to quantify because the stochastic simulations are not reproducible and the results are sensitive to sampling noise. However, in 1000 resamples of the final stochastic phase of the simulations using corrected code the Bayes factors were reduced by, on average, 20%.

Many more simulations are needed to accurately estimate the Bayes factors, and they should be run according to the method described in the text.

# Explanation

The simulated phylogenies are constructed by coalescing lineages sampled from the first 50,000 simulated infections. Partial and delayed sampling is simulated by sampling only the 15% of simulated infections deemed ascertained, and by ignoring samples that precede the first simulated hospitalization, as described on page 8 of the [Supplementary Materials](https://www.science.org/doi/suppl/10.1126/science.abp8337/suppl_file/science.abp8337_sm.v2.pdf).

![Excerpt from page 8 of the Supplementary Materials](https://github.com/nizzaneela/Programming_error_explanation/blob/4b653347fb1b4642c98d82c50fcea29200c4add1/sample.png)

Samples are also ignored if they are on branches from upstream of a stable coalescence, defined on page 10 of the [Supplementary Materials](https://www.science.org/doi/suppl/10.1126/science.abp8337/suppl_file/science.abp8337_sm.v2.pdf).  

![Excerpt from page 10 of the Supplementary Materials](https://github.com/nizzaneela/Programming_error_explanation/blob/b988d5b5b507d88619c9b9fb9fcaceb5349ff771/sctext.png)

Contrary to this definition, the function `coalescent_timing` in the script [stableCoalescence_cladeAnalysis.py](https://github.com/sars-cov-2-origins/multi-introduction/blob/78ec9e3b90215267b45ed34be2720566b7398b77/FAVITES-COVID-Lite/scripts/stableCoalescence_cladeAnalysis.py) does not stop the tMRCA calculations when 50,000 individuals have been infected. Instead, the calculations always continue until the last day of the simulation.

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

This can be confirmed by inspecting the content of `coalData_parameterized.txt` for each simulation in the data stored on Zenodo [here](https://zenodo.org/records/6899613). For example, the first simulation run `0001` reaches 50,000 infections on day 39, when the tMRCA of active sampled infections is 0.000333 years (~3 hours), but the calculations continue until the end of the simulation, after another 61 days and 1.3 million infections, when only 710 sampled infections are still active and their tMRCA is 0.016277 years (~6 days).

```
time	coalescence time	total infected	currently infected	current samples
...
38	0.000333	45444	33079	5731
39	0.000333	53606	38638	6175
...
99	0.016277	1363477	149166	743
100	0.016277	1371985	144107	710
```

It is the tMRCA from the end of the simulation that is used as the time of stable coalescence, so the code removes basal lineages that do not have active sampled infections at the end of simulation (day 100), even if the lineages do have active sampled infections at the end of sampling period (infection 50,000). 

By removing basal lineages that do not have active sampled infections at the end of the simulation period, and retaining those that do, the code filters out basal lineages that did not undergo early growth, so that the MRCA of the retained lineages is more likely to be associated with a superspreading event, and thus more likely to have a basal polytomy.

Additionally, the `main` function in the script [stableCoalescence_cladeAnalysis.py](https://github.com/sars-cov-2-origins/multi-introduction/blob/78ec9e3b90215267b45ed34be2720566b7398b77/FAVITES-COVID-Lite/scripts/stableCoalescence_cladeAnalysis.py) restores basal lineages if their MRCA is sufficiently close to that of the retained lineages. More specifically, if the MRCA of the retained lineages is on a zero-length branch, the code will add back all basal lineages connected via zero-length branches.
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
This increases the size of basal polytomies when the MRCA is on a zero-length branch.

The zero-length branches occur when coalescent events within a single host are compressed into a short period of time. This compression is exacerbated by the model used for coalescing lineages in the program [CoaTran](https://github.com/niemasd/CoaTran/tree/main), where a fixed effective population size of 1 year is used back until the first transmission event, and then a truncated exponential distribution is used to sample waiting times for coalescent events until the remaining lineages are coalesced down to one lineage at the time of infection. 

This compression is exacerbated further by an error in the epidemic simulation script [FAVITES-COVID-Lite_noSeqgen.py](https://github.com/sars-cov-2-origins/multi-introduction/blob/78ec9e3b90215267b45ed34be2720566b7398b77/FAVITES-COVID-Lite/scripts/FAVITES-COVID-Lite_noSeqgen.py) that skips the latent phase of the primary case. Specifically, the primary case is set to start in the infectious compartment (`P1`).
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
...despite the [command](https://github.com/sars-cov-2-origins/multi-introduction/blob/main/FAVITES-COVID-Lite/commands/command.0.28TF_0.15r.txt) indicating that it should start as exposed but non-infectious (`--tn_freq_e 0.00000020`).
```
~/scripts/FAVITES-COVID-Lite-updated.py --gzip_output --path_ngg_barabasi_albert ngg_barabasi_albert --path_gemf GEMF --path_coatran_constant coatran_constant --path_seqgen seq-gen --cn_n 5000000 --cn_m 8 --tn_s_to_e_seed 0 --tn_e_to_p1 125.862069 --tn_p1_to_p2 999999999 --tn_p2_to_i1 23.804348 --tn_p2_to_a1 134.891304 --tn_i1_to_i2 62.931034 --tn_i1_to_h 0.000000 --tn_i1_to_r 62.931034 --tn_i2_to_h 45.061728 --tn_i2_to_r 0.000000 --tn_a1_to_a2 9999999999 --tn_a2_to_r 125.862069 --tn_h_to_r 12.166667 --tn_s_to_e_by_e 0 --tn_s_to_e_by_p1 0 --tn_s_to_e_by_p2 3.513125 --tn_s_to_e_by_i1 6.387500 --tn_s_to_e_by_i2 6.387500 --tn_s_to_e_by_a1 0 --tn_s_to_e_by_a2 3.513125 --tn_freq_s 0.99999980 --tn_freq_e 0.00000020 --tn_freq_p1 0 --tn_freq_p2 0 --tn_freq_i1 0 --tn_freq_i2 0 --tn_freq_a1 0 --tn_freq_a2 0 --tn_freq_h 0 --tn_freq_r 0 --tn_end_time 0.273973 --tn_num_seeds 1 --pt_eff_pop_size 1 --pm_mut_rate 0.00092 --o 
```

Thus, the code:
- removes basal lineages that do not have active sampled infections at the end of the simulation, effectively filtering out those that did not undergo early growth, thereby increasing the likelihood of basal polytomies,
- adds back basal lineages that are connected to the stable coalescence via zero-length branches, thereby increasing the size of some basal polytomies, and
- skips the latent phase of the primary case, compressing the time for coalescing lineages in the primary case, thereby increasing the likelihood of zero-length branches at the stable coalescence.


This behaviour does not agree with the method defined in the [Supplementary Materials](https://www.science.org/doi/suppl/10.1126/science.abp8337/suppl_file/science.abp8337_sm.v2.pdf).

# Explanation

stableCoalescence_cladeAnalysis.py](https://github.com/sars-cov-2-origins/multi-introduction/blob/78ec9e3b90215267b45ed34be2720566b7398b77/FAVITES-COVID-Lite/scripts/stableCoalescence_cladeAnalysis.py) might be corrected by:
- breaking the loop in the function `coalescent_timing` once 50,000 individuals have been infected, e.g.:
```
def coalescent_timing(time_inf_dict, current_inf_dict, total_inf_dict, tree, num_days=100):
    ...
    for index, time in enumerate(times):
        if time > num_days: # sometimes gemf goes past the limit but we don't always know when
            break
        ### added break condition ###
        elif total_inf_dict[time-1]>50000: # stop after the day of 50000 total infections
            break
        ...
```
- then finding the stable coalescence by walking back as long as the tMRCA of active sampled infections is within one day before the final tMRCA, and returning the stable coalescence, e.g.:
```
    # work back day by day from the last day but stop if
    # the tMRCA of the next earlier day is more than a day earlier than the final tMRCA,
    # or the tMRCA of the next earlier day is later than the final tMRCA
    # or if the next earlier has only one sample, i.e. it doesn't have a tMRCA
    end_tMRCA = heights[-1]
    day = -1
    while (0 <= (end_tMRCA - heights[day-1]) <= (1/365)) and current_samples[day-1] != 1:
        day -= 1
    # stable coalescence is the MRCA of that day
    stable_coalescence = tree.mrca(current_labels[day])

    coalescent_timing_results = [used_times, heights, total_inf, current_inf, current_samples]
    return stable_coalescence, coalescent_timing_results
```
- modifying the `main` function to accept the stable coalescence returned from calling `coalescent_timing`, and estracting the substree rooted at the stable coalescence, e.g.:
```
    # stable coalescence
    time_inf_dict, current_inf_dict, total_inf_dict = tn_to_time_inf_dict(args.transmission_network, subtree)
    stable_coalescence, coal_timing = coalescent_timing(time_inf_dict, current_inf_dict, total_inf_dict, subtree, args.num_days)

    # prepare for clade analysis; get the subtree with the stable coalescence (MRCA) root
    subtree_sc = tree.extract_subtree(stable_coalescence)
    subtree_sc.root.edge_length = 0
    subtree_sc.suppress_unifurcations()
```
