# Key results are inflated by an error in the code
The authors find that the SARS-CoV-2 lineages "A" and "B" are the result of at least two separate cross-species transmission events into humans. This finding is based on Bayes factors calculated in the Jupyter notebook [cladeAnalysis.ipynb](https://github.com/sars-cov-2-origins/multi-introduction/blob/71ed420fe11ecdbe589568255ec90ca56d6e221c/notebooks/cladeAnalysis.ipynb). In this notebook, a bug in the function `clade_analysis_updated` causes the main analysis to output:

![Screenshot of main analysis output before bugfix](https://github.com/nizzaneela/Drafting/blob/afeaa3058f8185152451b8751a843924b7f93ef6/main_result_original.png)

After fixing the bug, this output becomes:

![Screenshot of main analysis output after bugfix](https://github.com/nizzaneela/Drafting/blob/abbf139cce753e0cbc5988b1a144761272c4a367/main_result_fixed.png)

Thus, after re-running the main analysis without the bug, the Bayes factors are reduced by a factor of six. This significantly reduces confidence in the authors' findings.

The sensitivity analyses should be re-run without the bug (re-running the sensitivity analyses isn't feasible without the data, which have not been published).

A published correction may be necessary.

There are other deficiencies. In particular, the different hypotheses are tested against different evidence (more specifically, against different conditions derived from the evidence), so that the Bayes factors are heavily biased in favour of the multiple spillover hypothesis. This issue will be addressed in a separate comment.

## Explanation
The Bayes factor calculations require likelihoods of the different hypotheses (single introduction or multiple introductions) producing the observed data. For the single introduction hypothesis, this is evaluated by generating 1100 simulated viral phylogenies and counting those with topologies that are deemed compatible with the observed data. The compatible topologies are denoted `CC` and `AB` and correspond, respectively, to the center and right topologies in Figure 2 (reproduced below).

![Figure 2 of Pekar et al. 2022](https://github.com/nizzaneela/Drafting/blob/a0ca501f04b8937370334bcacde35d906abd5288/science.abp8337-f2.jpg)  

As shown in Figure 2, the `CC` topology has two clades on branches with one mutation from the MRCA ('one-mutation clades'). The `AB` topology has a clade on a branch with two mutations from the MRCA (a "two-mutation" clade). The `AB` topology also requires:
- a specific size - the number of taxa in the two-mutation clade must be between 30% and 70% of the total number of taxa in the whole phylogeny;
- a basal polytomy - at least 100 clades must descend directly from the MRCA; and 
- a polytomy at the two-mutation clade - at least 100 subclades must descend directly from the two-mutation clade root.

For each run, the function `clade_analysis_updated` stores details of the one- and two-mutation clades:
- for each one-mutation clade, `clade_analysis_CC_d[run]` stores that one-mutation clade's size in a list labelled `'clade_sizes'` and stores a sublist with the sizes of that one-mutation clade's subclades in another list labelled `subclade_sizes`; and
- for each two-mutation clade, `clade_analysis_AB_d[run]` stores that two-mutation clade's size in a list labelled `'clade_sizes'` and stores a sublist with the sizes of that one-mutation clade's subclades in another list labelled `subclade_sizes`. 

Then, under the comment `# A/B analysis`, each phylogeny is tested as follows. 

For the size requirement:
```
# get the list with the sizes of the two-mutation clades
clade_sizes = clade_analyses_AB_d[run]['clade_sizes']
.
# get the size of the largest two-mutation clade and test it against the size requirement
if max(clade_sizes) >= lower_constraint*num_leaves and max(clade_sizes) <= upper_constraint*num_leaves: 
```

For the basal polytomy requirement:
```
# get the number of the one-mutation clades from the length of the list with the sizes of the one-mutation clades
base_polytomy_size = len(clade_analyses_CC_d[run]['clade_sizes'])
.
# test the number of one-mutation clades against the minimum polytomy size
if base_polytomy_size >= min_polytomy_size:
```

For the two-mutation polytomy requirement, the relevant code is:
```
# get the list with the subclade sizes of the one-mutation clades
subclade_sizes = clade_analyses_CC_d[run]['subclade_sizes'].copy()
.
# if there is only one two-mutation clade, get a number from the length of the sublist with the subclade sizes for the first one-mutation clade,
# and test the number against minimum polytomy size
    if len(subclade_sizes[0]) >= min_polytomy_size: 
# else get the index of the largest two-mutation clade
    max2mutCladeLoc = clade_sizes.index(max(clade_sizes))
# then get a number from the length of the sublist with the subclade sizes for the corresponding one-mutation clade,
# and test the number against minimum polytomy size
    if len(subclade_sizes[max2mutCladeLoc]) >= min_polytomy_size: 
```

Thus, rather than checking if the two-mutation clade has enough subclades, the code checks to see if a one-mutation clade has enough subclades.

It seems 
```
subclade_sizes = clade_analyses_CC_d[run]['subclade_sizes'].copy()
```
should have been changed to 
```
subclade_sizes = clade_analyses_AB_d[run]['subclade_sizes'].copy()
```
but wasn't.
## Verification
The bug can be reproduced, analysed and fixed as follows.

Get the `cladeAnalysis.ipynb` from the GitHub.
```
wget https://raw.githubusercontent.com/sars-cov-2-origins/multi-introduction/main/notebooks/cladeAnalysis.ipynb
```
Get and unzip the necessary data from Zenodo (checksum is md5:217dfad4e075cdba908268116f43a45e ).
```
wget https://zenodo.org/record/6899613/files/simulations_cumulative_results.zip
unzip simulations_cumulative_results.zip
unzip simulations_cumulative_results/clade_analyses_CC.zip
unzip simulations_cumulative_results/clade_analyses_AB.zip
```
Launch the notebook
```
jupyter-notebook cladeAnalysis.ipynb
```
Comment out (using `#`) any unnecessary imports that cause errors, e.g.
```
#import treeswift
#from utils import *
.
#import seaborn as sns
```
Update the second and sixth code cells so that `clade_analyses_CC_dir` and `clade_analyses_AB_dir` point to the newly unpzipped `clade_analyses_CC/` and `clade_analyses_AB/` directories, e.g.
```
clade_analyses_CC_dir = 'clade_analyses_CC/'
clade_analyses_AB_dir = 'clade_analyses_AB/'
```
Delete the seventh and eighth code cells, because they require data that has not been published (or just ignore them when they throw errors).

Run the notebook to reproduce the published results of the main analysis.
![Screenshot of main analysis output before bugfix](https://github.com/nizzaneela/Drafting/blob/afeaa3058f8185152451b8751a843924b7f93ef6/main_result_original.png)

To test the bug, modify the A/B analysis by adding a list `subclade_sizes_correct` that takes the lists of subclade sizes for the two-mutation clades.
```
for run in clade_analyses_AB_d:
        num_leaves = sum(clade_analyses_CC_d[run]['clade_sizes'])
        base_polytomy_size = len(clade_analyses_CC_d[run]['clade_sizes'])
        clade_sizes = clade_analyses_AB_d[run]['clade_sizes']
        subclade_sizes = clade_analyses_CC_d[run]['subclade_sizes'].copy()
### add the line below
        subclade_sizes_correct = clade_analyses_AB_d[run]['subclade_sizes'].copy()
```
Every time the two-mutation polytomy test fails, re-run the test using `subclade_sizes_correct`. If it then passes, print the size of the two-mutation polytomy and to the count.
```
# single two-mutation clade case
if len(subclade_sizes[0]) >= min_polytomy_size: # two-mutation clade has polytomy
    ab_count_30perc_twoPolytomies += 1
### add the four lines below
elif len(subclade_sizes_correct[0]) >= min_polytomy_size:
    print('Two-mutation polytomy test failed when the two-mutation polytomy was:')
    print(len(subclade_sizes_correct[0]))
    ab_count_30perc_twoPolytomies += 1
.
# multiple two-mutation clade case
if len(subclade_sizes[max2mutCladeLoc]) >= min_polytomy_size: # two mutation clade has polytomy
    ab_count_30perc_twoPolytomies += 1
### add the four lines below
elif len(subclade_sizes_correct[max2mutCladeLoc]) >= min_polytomy_size:
    print('Two-mutation polytomy test failed when the two-mutation polytomy was:')
    print(len(subclade_sizes_correct[max2mutCladeLoc]))
    ab_count_30perc_twoPolytomies += 1
```
Re-run the notebook to print sizes of the two-mutation polytomies that were being wronjgly rejected.

![Output printing sizes of two-mutation polytomies that were rejected](https://github.com/nizzaneela/Drafting/blob/57d6769b4d8e63ec59565706a92694197ac29b88/main_result_with_rejected_polytomies.png)

And to get the corrected Bayes factors, which have been reduced by a factor of six.

![Output including the corrected Bayes factors](https://github.com/nizzaneela/Drafting/blob/f001930392bd187cd5259b14c6ec0d40b1958122/main_result_with_rejected_polytomies_and_results.png)

