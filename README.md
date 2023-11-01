The Bayes factors are further inflated by another error in the code.

The simulated phylogenies are pruned to remove basal lineages that branch upstream of a "stable coalescence". The "stable coalescence" is defined in the [Supplementary Materials](https://www.science.org/doi/suppl/10.1126/science.abp8337/suppl_file/science.abp8337_sm.v2.pdf) as a tMRCA that ignores basal lineages that do not survive the sampling period. However, the implementation in the function `coalescent_timing` in [stableCoalescence_cladeAnalysis.py](https://github.com/sars-cov-2-origins/multi-introduction/blob/78ec9e3b90215267b45ed34be2720566b7398b77/FAVITES-COVID-Lite/scripts/stableCoalescence_cladeAnalysis.py)) ignores basal lineages that do survive the sampling period. The code removes basal lineages that, according to the text, should be retained.

For the primary analysis, correcting this error shifts 13% of the tMRCAs back in time by an average of 5.5 days, and reduces the Bayes factors by ~6%. 

# Explanation

The stable coalescence was introduced in Pekar et al.'s [Timing the SARS-CoV-2 index case in Hubei province](https://www.science.org/doi/10.1126/science.abf8003), where tMRCAs inferred from observations were matched to tMRCAs from simulations. The authors noted that "coalescent processes can prune basal viral lineages before they have the opportunity to be sampled, potentially pushing SARS-CoV-2 tMRCA estimates forward in time". In other words, the tMRCAs inferred from observations might not account for basal lineages that went extinct before they could be sampled. The stable coalescence is a tMRCA that produces a similar but slightly different effect with the simulations - it ignores simulated basal lineages that went extinct before the end of the sampling period. The effect is shown in [Fig. 2](https://www.science.org/cms/10.1126/science.abf8003/asset/7e12255a-8ddf-4d55-bc59-6644bc8de6e6/assets/graphic/372_412_f2.jpeg) of that paper, reproduced below.

![Fig. 2 of "Timing the index case...](https://github.com/nizzaneela/Programming_error_explanation/blob/dae78dd3e2658b59473d68ce5da2a5c9d2284f8b/timing_f2.jpeg)

In the present analysis, the simulations are sampled from the time of the first hospitalization until the fifty-thousandth infection, as described on page 8 of the [Supplementary Materials](https://www.science.org/doi/suppl/10.1126/science.abp8337/suppl_file/science.abp8337_sm.v2.pdf).

![Excerpt from page 8 of the Supplementary Materials](https://github.com/nizzaneela/Programming_error_explanation/blob/4b653347fb1b4642c98d82c50fcea29200c4add1/sample.png)

The stable coalescence is the tMRCA of sampled infections that are active on the day of the fifty-thousandth infection or, should the epidemic fail to grow so far, the end of the simulation, as described on page 10 of the [Supplementary Materials](https://www.science.org/doi/suppl/10.1126/science.abp8337/suppl_file/science.abp8337_sm.v2.pdf).  

![Excerpt from page 10 of the Supplementary Materials](https://github.com/nizzaneela/Programming_error_explanation/blob/b988d5b5b507d88619c9b9fb9fcaceb5349ff771/sctext.png)


Contrary to this definition, the function `coalescent_timing` in the script [stableCoalescence_cladeAnalysis.py](https://github.com/sars-cov-2-origins/multi-introduction/blob/78ec9e3b90215267b45ed34be2720566b7398b77/FAVITES-COVID-Lite/scripts/stableCoalescence_cladeAnalysis.py) does not stop the tMRCA calculations when 50,000 individuals have been infected. Instead, the calculations always continue until the last (hundredth) day of the simulation.

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

This can be confirmed by inspecting the content of `coalData_parameterized.txt` for each simulation in the data stored on Zenodo [here](https://zenodo.org/records/6899613) (and collected [here](https://github.com/nizzaneela/multi-introduction/blob/6c4c02e1a614d3cf482da76a188729f9c6e1933c/notebooks/0.28TF/simulations.zip) for convenience). For example, the first simulation run `0001` reaches 50,000 infections on day 39, when the tMRCA is 0.000333 years (~3 hours), but the calculations continue until the end of the simulation 61 days later, when the tMRCA is 0.016277 years (~6 days).

```
time	coalescence time	total infected	currently infected	current samples
...
38	0.000333	45444	33079	5731
39	0.000333	53606	38638	6175
...
99	0.016277	1363477	149166	743
100	0.016277	1371985	144107	710
```

It is the tMRCA from the end of the simulation that is used as the stable coalescence. In particular, the `main` function of [stableCoalescence_cladeAnalysis.py](https://github.com/sars-cov-2-origins/multi-introduction/blob/78ec9e3b90215267b45ed34be2720566b7398b77/FAVITES-COVID-Lite/scripts/stableCoalescence_cladeAnalysis.py) extracts the subtree rooted at the tMRCA from the end of the simulation, and uses the extracted subtree for the subsequent analysis.
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

Thus, the code ignores basal branches that do not have active sampled infections at the end of simulation (day 100), even if the branches do have active sampled infections at the end of sampling period (infection 50,000). This behaviour is unjustifiable and contrary to the method defined in the [Supplementary Materials](https://www.science.org/doi/suppl/10.1126/science.abp8337/suppl_file/science.abp8337_sm.v2.pdf).

This error might be corrected by breaking the loop in the function `coalescent_timing` once 50,000 individuals have been infected, e.g.:
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
# Verification

The stable coalescents can be extracted from the `coalData_parameterized.txt` files for each simulation collected at [this repository](https://github.com/nizzaneela/multi-introduction/blob/6c4c02e1a614d3cf482da76a188729f9c6e1933c/notebooks/0.28TF/simulations.zip), or from the summary stored `FAVITES_results` [here](FAVITES-COVID-Lite/cumulative_results/FAVITES_results.zip):
```
wget https://github.com/nizzaneela/multi-introduction/blob/6c4c02e1a614d3cf482da76a188729f9c6e1933c/notebooks/0.28TF/simulations.zip
wget FAVITES-COVID-Lite/cumulative_results/FAVITES_results.zip
unzip simulations.zip
unzip FAVITES_results.zip
```
The stable coalescents stored in `FAVITES_results` can be checked against the voorrect and incorrect values derivable from the `coalData_parameterized.txt` files:



The code simulates mutations through the subtree, starting from the stable coalescence. As the randomw number generator used to simulate the mutations did not have a seed, the simulated mutations cannot be replicated. However, the effect of the correction can s 

