Key results are inflated by an error in the code.

The authors find that the SARS-CoV-2 lineages "A" and "B" are the result of at least two separate cross-species transmission events into humans. This finding is based on Bayes factors calculated in the Jupyter notebook [cladeAnalysis.ipynb](https://github.com/sars-cov-2-origins/multi-introduction/blob/71ed420fe11ecdbe589568255ec90ca56d6e221c/notebooks/cladeAnalysis.ipynb). This notebook has a bug in the function `clade_analysis_updated`. Re-running the main analysis without the bug reduces the Bayes factors by a factor of six. This significantly reduces confidence in the authors' findings.

The sensitivity analyses should be re-run without the bug (this isn't feasible without the data, which have not been published).

There are other deficiencies. In particular, the different hypotheses are tested against different evidence (more specifically, against different conditions derived from the evidence), so that the Bayes factors are heavily biased in favour of the multiple spillover hypothesis. This issue will be addressed in a separate comment.

## Explanation
The Bayes factor calculations require likelihoods of the observed data being produced by the different hypotheses. For the single introduction hypothesis, this is evaluated by generating 1100 simulated viral phylogenies and counting those with topologies that are deemed compatible with the observed data. The compatible topologies are denoted `CC` and `AB` and correspond, respectively, to the center and right topologies in [Figure 2](https://www.science.org/doi/10.1126/science.abp8337#F2 ) of the paper (reproduced below).

![Figure 2 of Pekar et al. 2022](https://pubpeer.com/storage/image-1690798546245.png)

As shown in Figure 2, the `CC` topology has two clades on branches with one mutation from the MRCA ('one-mutation clades'). The `AB` topology has a clade on a branch with two mutations from the MRCA (a "two-mutation" clade). 

For each simulated viral phylogeny XXXX, the details of the one- and two-mutation clades are collated by the script [stableCoalescence_cladeAnalysis.py](https://github.com/sars-cov-2-origins/multi-introduction/blob/71ed420fe11ecdbe589568255ec90ca56d6e221c/FAVITES-COVID-Lite/scripts/stableCoalescence_cladeAnalysis.py) into the files `XXXX_clade_analysis_CC_polytomy.txt` and `XXXX_clade_analysis_AB_polytomy.txt`, respectively. The function `clade_analysis_updated` extracts these details to `clade_analysis_CC_d[XXXX]` and `clade_analysis_AB_d[XXXX]`.  The details list each clade's size (under the label `'clade_sizes'`) and the sizes of its subclades (under the label `'subclade_sizes`). 

Note that the one-mutation clades include all clades on branches with *at least*  one mutation from the MRCA. Likewise, the two-mutation clades include all clades on branches with *at least* two mutations from the MRCA. This means that the two-mutation clades are a subset of the one mutation-clades.  However, the details of the two-mutation clades will not normally be stored at the same place in the lists of `clade_analysis_CC_d[XXXX]` and `clade_analysis_AB_d[XXXX]`.

A simulated viral phylogeny is deemed to have an `AB` topology if it has:
- a specific size - the number of taxa in the two-mutation clade must be between 30% and 70% of the total number of taxa in the whole phylogeny;
- a basal polytomy - at least 100 clades must descend directly from the MRCA; and 
- a polytomy at the two-mutation clade - at least 100 subclades must descend directly from the two-mutation clade root.

These requirements are tested under `# A/B analysis` in the function `clade_analysis_updated` in the notebook  [cladeAnalysis.ipynb](https://github.com/sars-cov-2-origins/multi-introduction/blob/71ed420fe11ecdbe589568255ec90ca56d6e221c/notebooks/cladeAnalysis.ipynb). The relevant code is reproduced below 
```
    # A/B analysis
    ab_count_30perc = 0 # interested in 2 mutations clade that are at least 30% of all taxa
    ab_count_30perc_polytomy = 0 # interested in 2 mutations clade that are at least 30% of all taxa + has a basal polytomy
    ab_count_30perc_twoPolytomies = 0 # interested in 2 mutations clade that are at least 30% of all taxa + has a basal polytomy + polytomy at 2 mutation clade
    lower_constraint = 0.3 # the 2-mutation clade must be at least 30% of all taxa
    upper_constraint = 0.7 # the 2-mutation clade must be at most 70% of all taxa
    
    for run in clade_analyses_AB_d:
        num_leaves = sum(clade_analyses_CC_d[run]['clade_sizes'])
        base_polytomy_size = len(clade_analyses_CC_d[run]['clade_sizes'])
        clade_sizes = clade_analyses_AB_d[run]['clade_sizes']
        subclade_sizes = clade_analyses_CC_d[run]['subclade_sizes'].copy()
        if not clade_sizes: # no 2 mutation clades
            continue
        if max(clade_sizes) >= lower_constraint*num_leaves and max(clade_sizes) <= upper_constraint*num_leaves: # clades match size restrictions
            if len(clade_sizes) == 1:
                ab_count_30perc += 1
                if base_polytomy_size >= min_polytomy_size: # basal polytomy
                    ab_count_30perc_polytomy += 1
                    if len(subclade_sizes[0]) >= min_polytomy_size: # two-mutation clade has polytomy
                        ab_count_30perc_twoPolytomies += 1
            else:
                clade2 = sorted(clade_sizes, reverse=True)[1]
                ab_count_30perc += 1
                if base_polytomy_size >= min_polytomy_size: # basal polytomy
                    ab_count_30perc_polytomy += 1
                    max2mutCladeLoc = clade_sizes.index(max(clade_sizes))
                    if len(subclade_sizes[max2mutCladeLoc]) >= min_polytomy_size: # two mutation clade has polytomy
                        ab_count_30perc_twoPolytomies += 1
```
Although the code is somewhat messy, it can be seen that the size of the polytomy at the two-mutation clade is determined in two different ways:
- if there is only one two-mutation clade (i.e. `len(clade_sizes) == 1`),  the size of the polytomy is taken as `len(subclade_sizes[0])`; and
- otherwise, the index of the largest two-mutation clade is determined (i.e. `max2mutCladeLoc = clade_sizes.index(max(clade_sizes)`) and the size of  the polytomy is taken as `len(subclade_sizes[max2mutCladeLoc])`.

It can also be seen that  `clade_sizes` is from the two-mutation list (i.e. `clade_sizes = clade_analyses_AB_d[run]['clade_sizes']`), but `subclade_sizes` is from the one-mutation list (i.e. `subclade_sizes = clade_analyses_CC_d[run]['subclade_sizes'].copy()`). 

This means that, rather than testing if the two-mutation has a sufficiently large polytomy, the code instead tests a one-mutation clade that is randomly selected, depending on how the clades are ordered in `clade_analysis_CC_d[run]` and where the largest (or only) clade was located in `clade_analysis_AB_d[run]`.

This does not agree with he `AB` topology requirements described in text and the supplementary materials. Moreover, it makes no sense. There can be no valid reason for introducing randomness in this way. Rather, this seems to be the result of a simple copy-paste error, where `subclade_sizes = clade_analyses_CC_d[run][
'subclade_sizes'].copy()` should have been changed to `subclade_sizes = clade_analyses_AB_d[run]['subclade_sizes'].copy()`. 

## Verification
The bug can be reproduced with code and data from the GitHub [repository](https://github.com/sars-cov-2-origins/multi-introduction).
```
wget https://raw.githubusercontent.com/sars-cov-2-origins/multi-introduction/main/notebooks/cladeAnalysis.ipynb
wget https://github.com/sars-cov-2-origins/multi-introduction/blob/71ed420fe11ecdbe589568255ec90ca56d6e221c/FAVITES-COVID-Lite/cumulative_results/clade_analyses_CC.zip
wget https://github.com/sars-cov-2-origins/multi-introduction/blob/71ed420fe11ecdbe589568255ec90ca56d6e221c/FAVITES-COVID-Lite/cumulative_results/clade_analyses_AB.zip
unzip simulations_cumulative_results/clade_analyses_CC.zip
unzip simulations_cumulative_results/clade_analyses_AB.zip
```
After launching the notebook, it may be necessary to comment out any unnecessary imports that cause errors, e.g.
```
#import treeswift
#from utils import *
.
#import seaborn as sns
```
The second and sixth code cells must be updated so that `clade_analyses_CC_dir` and `clade_analyses_AB_dir` point to the newly unpzipped `clade_analyses_CC/` and `clade_analyses_AB/` directories, e.g.
```
clade_analyses_CC_dir = 'clade_analyses_CC/'
clade_analyses_AB_dir = 'clade_analyses_AB/'
```
The seventh and eighth code cells may be deleted, because they require data that has not been published (else they may be ignored when they throw errors).

Re-running the notebook should reproduce the published results of the main analysis.

![Screenshot of main analysis output before bugfix](https://pubpeer.com/storage/image-1690805746422.png)

`subclade_sizes_correct` can then be corrected to use the two-mutation list.
```
    for run in clade_analyses_AB_d:
        num_leaves = sum(clade_analyses_CC_d[run]['clade_sizes'])
        base_polytomy_size = len(clade_analyses_CC_d[run]['clade_sizes'])
        clade_sizes = clade_analyses_AB_d[run]['clade_sizes']
####### CC has been changed to AB in the line below###########################    
        subclade_sizes = clade_analyses_AB_d[run]['subclade_sizes'].copy()

```
Re-running the notebook should now produce the intended results of the main analysis.

![Screenshot of main analysis output before bugfix](https://pubpeer.com/storage/image-1690806023778.png)

Correcting the error increases the likelihood of a single introduction producing the `AB` topology from 0.5% to 3%, and reduces the Bayes factors from ~60 to less than 10.
