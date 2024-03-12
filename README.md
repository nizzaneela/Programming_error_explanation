Errors and noise inflate the corrected Bayes factors.

The [Supplementary Materials](https://www.science.org/doi/suppl/10.1126/science.abp8337/suppl_file/science.abp8337_sm.v2.pdf) describe how short-lived basal lineages are excluded from the clade analysis by rooting each phylogeny at a stable coalescence. The code does not implement the described method. 

The effects of this discrepency are difficult to quantify because the stochastic simulations are not reproducible and the results are sensitive to sampling noise. However, in 1000 resamples of the final stochastic phase of the simulations using corrected code the Bayes factors were reduced by, on average, 20%.

Many more should simulations are needed to accurately estimate the Bayes factors, and they should be run according to the method described in the text.

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

It is the tMRCA from the end of the simulation that is used as the time of stable coalescence, so that the code removes basal lineages that do not have active sampled infections at the end of simulation (day 100), even if the lineages do have active sampled infections at the end of sampling period (infection 50,000). 

This behaviour does not agree with the method defined in the [Supplementary Materials](https://www.science.org/doi/suppl/10.1126/science.abp8337/suppl_file/science.abp8337_sm.v2.pdf). It is also makes no sense because a lineage can have active infections at the end of the simulation but lack active _sampled_ infections merely because the active infections were not amongst the first 50,000. 

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

By removing basal lineages that do not have active sampled infections at the end of the simulation period, and retaining those that do, the code filters out basal lineages that did not undergo rapid growth, so that the MRCA of the retained lineages is more likely to be associated with a superspreading event, and thus more likely to have a basal polytomy.

Additionally, the `main` function in the script [stableCoalescence_cladeAnalysis.py](https://github.com/sars-cov-2-origins/multi-introduction/blob/78ec9e3b90215267b45ed34be2720566b7398b77/FAVITES-COVID-Lite/scripts/stableCoalescence_cladeAnalysis.py) restores basal lineages if their MRCA is sufficiently close to that of the retained lineages. More specifically, if the MRCA of the retained lineages is on a zero-length branch, the code will add all basal lineages connected via zero-length branches.
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
Thus the code increases the size of a basal polytomy when there MRCA is on a zero-length branch/

The zero-length branches occur when coalescent events are compressed into a short period of time within a single host. This is partly a result of the coalescence model used in 



Thus, the code removes basal lineages that did not undergo rapid growth, thereby increasing 


and by walking back as long as the tMRCA of active sampled infections is within one day before the final tMRCA, e.g.:
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

And by walking back from the final day

However, the MRCA from the end of the simulation is not always used as the stable coalescence. 

In particular, the `main` function of [stableCoalescence_cladeAnalysis.py](https://github.com/sars-cov-2-origins/multi-introduction/blob/78ec9e3b90215267b45ed34be2720566b7398b77/FAVITES-COVID-Lite/scripts/stableCoalescence_cladeAnalysis.py) extracts the subtree rooted at the tMRCA from the end of the simulation, and uses the extracted subtree for the subsequent analysis.
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

# Verification

The stable coalescent was introduced in [Timing the SARS-CoV-2 index case in Hubei province](https://www.science.org/doi/10.1126/science.abf8003), where tMRCAs inferred from observations were aligned with tMRCAs from simulations. The tMRCAs inferred from the observations might have been pushed forward in time by the relatively late and sparse observations failing to detect short-lived basal lineages. The stable coalescence produces a similar effect with the tMRCAs from the simulations by ignoring basal lineages that do not live to the end of the sampling period. The effect is shown in [Fig. 2](https://www.science.org/cms/10.1126/science.abf8003/asset/7e12255a-8ddf-4d55-bc59-6644bc8de6e6/assets/graphic/372_412_f2.jpeg), reproduced below.

![Fig. 2 of "Timing the index case...](https://github.com/nizzaneela/Programming_error_explanation/blob/dae78dd3e2658b59473d68ce5da2a5c9d2284f8b/timing_f2.jpeg)




The script `stableCoalescence_cladeAnalysis.py` uses epidemic simulation output from GEMF, a transmission network from FAVITES, and a time tree from CoaTran. For the 1100 simulations of the main analysis, these are published on Zenodo in twenty-two zip files, each around 7GB compressed. They can be downloaded with:
```
for i in {01..22}; do wget https://zenodo.org/records/6899613/files/simulations_"$i".zip; done
```
and then unzipped and collated with:
```
mkdir ./simulations
for i in {01..22}; do unzip simulations_"$i".zip && mv ./simulations_"$i"/* ./simulations && rm -r ./simulations_"$i"; done
```

stable coalescents can be extracted from the `coalData_parameterized.txt` files for each simulation collected at [this repository](https://github.com/nizzaneela/multi-introduction/blob/6c4c02e1a614d3cf482da76a188729f9c6e1933c/notebooks/0.28TF/simulations.zip), or from the summary stored in `FAVITES_results` [here](FAVITES-COVID-Lite/cumulative_results/FAVITES_results.zip):
```
wget https://github.com/nizzaneela/multi-introduction/blob/6c4c02e1a614d3cf482da76a188729f9c6e1933c/notebooks/0.28TF/simulations.zip
wget FAVITES-COVID-Lite/cumulative_results/FAVITES_results.zip
unzip simulations.zip
unzip FAVITES_results.zip
```
The stable coalescents stored in `FAVITES_results` can be checked against the voorrect and incorrect values derivable from the `coalData_parameterized.txt` files:



The code simulates mutations through the subtree, starting from the stable coalescence. As the randomw number generator used to simulate the mutations did not have a seed, the simulated mutations cannot be replicated. However, the effect of the correction can s 

