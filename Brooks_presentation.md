## Enhancing the coevolutionary signal
----------
### [Travis Hoppe](http://thoppe.github.io/)
NIH/NIDDK/LCP Postdoctoral Fellow

====*

# Outline
MSA
Coevolution
(scoring!)
structure prediction

====

### Sequence to structure prediction
Sequence / Mean sequence alignment 
    DTSGVQGIDVSHWQGSINWSSVKSAGMSFAYIKATEGTNYKDDRFSANYTNAYNAGIIRGAYHFARPNASSGTAQADYFASNGGGWSRDNRTLPGVLDIEHNPSGAMCYGLSTTQMRTWINDFHARYKARTTRDVVIYTTASWWNTCTGSWNGMAAKSPFWVAHWGVSAPTVPSGFPTWTFWQYSATGRVGGVSGDVDRNKFNGSAARLLALANNTA
    
    ----DYGIDVSSSTSQSQWSCLAGKN-QRAIIQVWSGGYGLNSQASSIISAAKSAGFQVDVYAFLCNQCSPSSNVIQQIVNSL---GGQFGT--LWIDVEQCS---GCWG-DVNDNAAFVAEAVQTAAS-LGVTVGVYSSLGEWPQTVGSL-SSLSSYPQWYAHYDGVAASQYGGWDNPEMKQYVGNTNECGV--SVDLDYYG--------------
    ----ELGIDVSSATSQSQWSCLAQKN-QRAIIQVWSGGYGMNNGVVSAIQAAQNAGFQVDLYAFLCNQCSPSSNVIQQIVSKIKQSGVSFGT--LWIDVEQCS---GCWG-STSANAAFVVEAVQTAAS-LGVRVGVYSSSGEWPQTVGTL-TSLSSYPQWYAHYDGVPAGQYGGWNNPEMKQYVGNTNQCGV--SVDLDFYG--------------
    ----TYGVDL------AGFQCLVGKGF-FAIVRCYMSSGGIDPNCASSVSAAWAGGMTVDLYLFPCFSCG----SLVQFAQS---NGVNFGK--IWLDIEGPG---TYWG-DQGANQQFFEGLVQGL--S-GVSVGIYTSESQWSPIMGDY-SGGSNFPLWYANYDGSPN-PFGGWSTPTMKQFDDPSN-CGI--GIDENWIG--------------
    ----GTGIDISSPTSKTQWSCLAKQN-TKAIIQVWSGGYGYNTNIASSVSAAKSAGIQVDLYAFLCSQCSPSSSAIKTLVSNLRSQNVEFGT--LWIDVEQCS---NCWG-STSTNAQFVVEAVQTAQQ-LGVSVGVYSSIGEWSQTVGSL-NSLSSFPLWYAHYDNVPASQFGSWSSPAMKQYAGNTQQCGV--SVDLDFFQ--------------
Contact maps / Structure
!(figures/1jfx_cmap.png) <<height:400px;>>
!(figures/1jfx.png)      <<height:400px;>>

====

## What is coevolution?

====

# Mutual information
(naïve attempt) 

# $MI_{ij} = \sum f_{ij} (A_i, A_j) \ln \left ( \frac{f_{ij}(A_i,A_j)}{f_i(A_i)f_j(A_j)} \right )$

$f_i, f_{ij}$ are observed frequencies and co-frequencies respectively. 

works poorly due to transitivity. e.g. A-B and B-C, this model predicts A-C.

====

### Maximum-entropy model / Markov Random Field
# $P(X=x) = \frac{1}{Z} \exp \left ( \sum_{i=1}^L \left [ v_i(x_i) + \sum_{j>i}^L w_{ij}(x_i,x_j) \right ] \right ) $
$X_i$ is a random variable that represents 
amino acid composition at position $i$ from MSA.

##### Learned parameters
$v_i$ encodes individual propensity of each amino acid at position $i$
$w_{ij}$ statistical coupling of amino acid propensities between positions $i,j$

&& [Assessing the utility of coevolution-based residue–residue contact predictions in a sequence- and structure-rich era](http://www.pnas.org/content/110/39/15674.abstract) Kamisetty, Ovchinnikova, and Baker.

====*
# DCA
====*
# PSICOV
====*

## GREMLIN

Optimize the _pseduolikelihood_ of $v,w$
# $ P(v,w | D) = \sum_{n=1}^N \sum_{i=1}^L \log P (x_i^n | x_{i'}, v, w)$
# $ P (x_i^n | x_{i'}, v, w) = \frac{1}{Z_i} \exp \left ( v_i(x_i^n) + \sum_{j=1,j \neq i}^L w_{ij}(x_i^n,x_j^n) \right )$

Regularization (and priors), prevent overfitting
## $R(v,w) = \lambda_v ||v||^2 + \sum_{i,j} \lambda_w^{i,j} || w_{i,j} || ^2$

====
# Target dataset

150 monomeric proteins $50< N <250$ residues
====
How good does it do?
====
GREMLIN model
!(figures/GREMLIN_only_Acc_Pre.png)

% ~/git-repo/GREMLIN_RF/analysis/plot_stats.py

====
RF model
!(figures/GREMLIN_RF_Acc_Pre.png)
====

# Scoring
For a given protein, alignment, GREMLIN gives $(N,N,21,21)$ tensor.

Reduce GREMLIN's tensor output:
# $S_{N,N,21,21} \rightarrow S_{N,N,20,20} \rightarrow  G_{N,N} \rightarrow G^{\text{APC}}_{N,N} = g$

Drop information about gaps.
Compute the Frobenius norm over each position.
Apply "average product correlation" (necessary!)

!! Show image of g map

====*

## Top score model

Rank sort $g$, take top $\alpha N$ contacts.

Typically $\alpha \in [0.1, 1.5]$.

Accuracy vs. Precision

====

# Hypothesis:
Local structure can enhance contact prediction
!(figures/local_structure_zoom.png) <<width:500px; transparent>>

Secondary structure is local (helices, sheets)

====

Machine learn local pixel maps

### Normalize data
Subtract mean, scale to unit variance

Use Random Forests to predict

====

Example proteins
!(figures/1a3a/1a3a_cartoon.png) <<width:500px; transparent>> [1a3a](http://www.rcsb.org/pdb/explore.do?structureId=1a3a) IIA MANNITOL FROM ESCHERICHIA COLI
!(figures/1a3a/1avs_cartoon.png) <<width:500px; transparent>> [1avs](http://www.rcsb.org/pdb/explore.do?structureId=1avs) CALCIUM-SATURATED N-TERMINAL DOMAIN OF TROPONIN C
====*
Example proteins, GREMLIN APC corrected score
!(figures/1a3a_GREMLIN.png) <<height:500px; transparent>> [1a3a](http://www.rcsb.org/pdb/explore.do?structureId=1a3a) IIA MANNITOL FROM ESCHERICHIA COLI
!(figures/1avs_GREMLIN.png) <<height:500px; transparent>> [1avs](http://www.rcsb.org/pdb/explore.do?structureId=1avs) CALCIUM-SATURATED N-TERMINAL DOMAIN OF TROPONIN C

====
!! Results plot
====
## Contact map vs cutoff length (1a3a)
!(figures/1a3a/animated_1a3a.gif) <<width:1200px; transparent>>
====*
## Contact map vs cutoff length (1avs) 
!(figures/1a3a/animated_1avs.gif) <<width:1200px; transparent>>
====*
### Folding simulations
!(figures/methods_pairplot.png)
====
### Folding 1a3a
!(figures/folding/folding_1a3a.gif)
folding is rapid in coarse grained simulation
====

What is being predicted?
ADD: Gaussian kernels (averaged decision trees)

====

What is being predicted?

!(figures/local_structure_distance.png) <<height:450px; transparent>> 
!(figures/FP_distance.png) <<height:450px; transparent>> 
====


# Improvement in folding
!(figures/folding/Q_avg.png)

====


Additional value 

$g$ can be used as an effective Hamiltonian for evolutionary movement

====

# Future work

Convolutional neural networks
Correct structure prediction (ROSETTA?) (Go?)

Disambiguation of intra/inter predictions

Binding partners

====

## Thank you.

####  Robert Best (NIH/NIDDK)
Wenwei Zheng
Travis Hoppe
*Pengfei Tian*
Jan Domanski
Mathias Bellaiche

#### Jeff Gray (John Hopkins)
Julia Joehler















