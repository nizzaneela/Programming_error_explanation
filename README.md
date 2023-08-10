I’m afraid there is another bug in the notebook [cladeAnalysis.ipynb](https://github.com/sars-cov-2-origins/multi-introduction/blob/71ed420fe11ecdbe589568255ec90ca56d6e221c/notebooks/cladeAnalysis.ipynb). 

The function `calculate_bf` takes an array of values as `simulation_results` and divides them by their sum to produce an array `pr_3_topos`, i.e.:
```
def calculate_bf(asr_results, simulation_results):
...
# FAVITES simulation information
# the 3 trees are in the order (t_p, t_1C, t_2C)
pr_3_topos = np.array(simulation_results)/sum(simulation_results)
```
This would be an appropriate way to derive likelihoods for a specific topologies, were the values the number of simulations with each topology and their sum the total number of simulations.

Unfortunately, the values are likelihoods for a subset of the topologies, as can be seen in `clade_analysis_update`:

```
polytomy_result = count_atLeastMinDescendants/1100
…
cc_result = cc_count_30perc_twoPolytomies/1100
ab_result = ab_count_30perc_twoPolytomies/1100
…
simulation_results = [polytomy_result, ab_result, cc_result]
…
bf_unconstrained = calculate_bf(unconstrained_results, simulation_results)
bf_recCA = calculate_bf(recCA_results, simulation_results)
```
In other words, they are the likelihoods for the subset of topologies shown in Figure 2:

![Figure 2 of Pekar et al. 2022](science.abp8337-f2.jpg)

By dividing these likelihoods by their sum, `calculate_bf` increases them according to the proportion of simulations that do not have a topology in this subset. This increase is applied twice in the likelihoods for the multiple introductions, because they are based on the single introduction likelihoods squared. For example, in the main analysis the likelihoods for the single introduction double from [0.475, 0.03, 0] to [0.94, 0.06, 0], while the likelihoods for the multiple introductions quadruple from [0.226, 0.001, 0] to [0.885, 0.004, 0]. The net result is doubling of the Bayes factors.

This is obviously wrong. The correct calculation, as described on page 13 of the supplementary materials, does not include this increase. 

![Excerpt from page 13 of the supplementary materials](p13.png)

This mistake might have been caused by miscommunication between the programmers of the different functions.
