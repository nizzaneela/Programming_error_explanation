The [Supplementary Materials](https://www.science.org/doi/suppl/10.1126/science.abp8337/suppl_file/science.abp8337_sm.v2.pdf) describe the derivation of the Bayes factors over pages 11 to 14. The approach is set out on page 13:

![](approach.png)

That is, the Bayes factors are the ratio of the likelihood of the sequence data $\bf{Y}$ arising from two introductions $I_2$ vs the likelihood of that sequence data $\bf{Y}$ arising from one introduction $I_1$:   

$$
BF = \frac{P(\mathbf{Y}|I_2)}{P(\mathbf{Y} |I_1)}
$$

Page 14 describes how the Bayes factor is taken from the combination of the posterior and prior odds. Though not stated explicitly, this must refer to multiplication of the posterior odds by the inverse of the prior odds (Bayes' theorem), i.e.:

$$
BF = \frac{P(I_2|\mathbf{Y})}{P(I_1|\mathbf{Y})} \cdot \frac{P(I_1)}{P(I_2)}
$$

Page 12 rewrites the posterior probabilities as joint probabilities and breaks them down in terms of the different MRCA haplotypes $S_{MRCA}$ so that the Bayes factor equations can be written:

$$
BF = \frac{\sum_{S_{MRCA}} P(S_{MRCA} | Y)P(I_2, S_{MRCA})}{\sum_{S_{MRCA}} P(S_{MRCA} | Y)P(I_1, S_{MRCA})} \cdot \frac{P(I_1)}{P(I_2)}
$$

Page 13 breaks down the joint probabilities $P(I_2, S_{MRCA})$ and $P(I_1, S_{MRCA})$ over the possible topologies $\tau$ so that the Bayes factor equations can be written:

$$
BF = \frac{\sum_{S_{MRCA}} \left( P(S_{MRCA} | Y)\left( \sum_{\tau} P(S_{MRCA} | \tau)P(\tau | I_2)P(I_2) \right) \right)}
{\sum_{S_{MRCA}}\left(  P(S_{MRCA} | Y)\left( \sum_{\tau} P(S_{MRCA} | \tau)P(\tau | I_1)P(I_1) \right) \right)} 
\cdot \frac{P(I_1)}{P(I_2)}
$$

Page 13 then describes how the likelihood is set for an MRCA under the specific topology (i.e. $P(S_{MRCA}|\tau)$ ):

![](set.png)

This description is misleading in at least two ways.

Firstly, the topology $(\tau_p, \tau_p)$ encompasses the topologies $(\tau_p, \tau_{1c})$, $(\tau_{1c}, \tau_p)$ and $(\tau_{1c}, \tau_{1c})$. Therefore, the topology likelihoods other than $P((\tau_p, \tau_p)|I_2)$ are redundant. Overlooking this may have lead to the error described in #7.

Secodnly, 



Firstly, 



$\left[ \sum_{\tau} P(S_{MRCA} | \tau)P(\tau | I_n)P(I_n) \right]$

The Bayes factors weigh the single introduction likelihoods according to the posterior probabilites of their compatible MRCA haplotypes. Defining weights $\alpha_{1c}$ and $\alpha_{2c}$ as $\frac{2 \cdot (P(S_A|Y) + P(S_B|Y))}{P(S_A|Y) + P(S_B|Y) + P(S_{C/C}|Y) + P(S_{T/T}|Y)}$ and $\frac{2 \cdot (P(S_{C/C}|Y) + P(S_{T/T}|Y)}{P(S_A|Y) + P(S_B|Y) + P(S_{C/C}|Y) + P(S_{T/T}|Y)}$ , respectively, the Bayes factor can be written as:

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
3. the clades are separated by two mutations (with one being basal for $P(\tau_{1c}|I_1)$ and both being derived for $P(\tau_{2c}|I_1)$ ).

The additional conditions 2 and 3 reduce the single introduction likelihoods. This makes the published Bayes factors (4.2 and 4.3) meaningless, since they may be caused by this reduction rather than a higher plausiblity of the two introduction hypothesis.

This can be corrected by applying conditions 2 and 3 to the two introduction likelihood $P((\tau_P,\tau_P)|I_2)$.

This comment examines the effects of condition 2. A subsequent comment will examine condition 3.

Two simulations can be be tested against condition 2 by going through the first 50,000 infections amongst both, and comparing the number of samples from each, e.g.:
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

The function `clade_analysis_updated` in the notebook [stableCoalescence_cladeAnalysis.py](https://github.com/nizzaneela/multi-introduction/blob/corrected/FAVITES-COVID-Lite/scripts/stableCoalescence_cladeAnalysis.py) can be adapted to do this instead of counting the number of simulations with basal polytomies. The `function calculate_bf` can then be adapted to use the resulting likelihood as $P((\tau_p,\tau_p)|I_2)$. $(\tau_p,\tau_p)$ encompasses all topologies compatible with the two introduction hypothesis, so the likelihoods of the other two introduction topologies ( $(\tau_p,\tau_{1c})$, $(\tau_{1c},\tau_p)$, $(\tau_{1c},\tau_{1c})$ ), can be set to zero.
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

Requiring the two introduction hypothesis to produce similarly sized clades, just as the one introduction hypothesis does already, reduces the Bayes factors by a factor of six.

The paper measures the two week doubling times of the simulated epidemics at the 1000th infection, and found a 95% highest density interval (HDI) of 1.35 to 5.44 days. This range of early growth rates suggests that two introductions are unlikely to grow to similar sizes within the first 50,000 infections, and the results confirm this.

Code and instructions for reproducing these results are available at [this branch](https://github.com/nizzaneela/multi-introduction/tree/corrected_with_relative_size_condition) of the authors' repository.

#### Limitations
1. `test_sizes` assumes both introductions are simultaneous. This maximises the likelihood of the two introductions having similar sizes;
2. Each introduction has a different delay before the first samples due to the different times of first hospitalization. Realistically, the earlier time should apply to both. Therefore, this analysis misses some early samples in the simulation with the later time of first hospitalization. In some cases, missing these samples will drop branch counts below the polytomy threshold, or clade sizes below the relative size threshold, reducing the likelihood of satisfying conditions 1 and 2;
3. Each introduction has a stable coalescence based on when it reaches 50,000 infections. The authors do not explain what aspect of reality the stable coalescence is used to reproduce. Assuming the stable coalescence should instead be based on when simulations reach 50,000 combined infections, this analysis removes more basal lineages than it should, increasing the likelihood of basal polytomies.

Overall, the limitations increase the likelihood of two simulations satisfying conditions 1 and 2, so the resulting Bayes factors may be considered a rough upper bound.

It is theoretically possible that applying condition 3 to the two introduction likelihoods could increase the Bayes factors. This would require an increase in the likelihood of two introductions producing topologies compatible with the more probable MRCA haplotypes (i.e. A or B). This is considered unlikely, but will be confirmed either way in the next comment.
