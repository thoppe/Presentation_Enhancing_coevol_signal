# Enhancing the coevolutionary signal
----------
### [Travis Hoppe](http://thoppe.github.io/)
NIH/NIDDK/LCP Postdoctoral Fellow

====*

# Outline

### MSA / coevolution
### *Score functions*
### Structure prediction

====

### Sequence to structure prediction
Sequence $\rightarrow$ Mean Sequence Alignment (MSA) $\rightarrow$ ...
    DTSGVQGIDVSHWQGSINWSSVKSAGMSFAYIKATEGTNYKDDRFSANYTNAYNAGIIRGAYHFARPNASSGTAQADYFASNGGGWSRDNRTLPGVLDIEHNPSGAMCYGLSTTQMRTWINDFHARYKARTTRDVVIYTTASWWNTCTGSWNGMAAKSPFWVAHWGVSAPTVPSGFPTWTFWQYSATGRVGGVSGDVDRNKFNGSAARLLALANNTA
    
    ----DYGIDVSSSTSQSQWSCLAGKN-QRAIIQVWSGGYGLNSQASSIISAAKSAGFQVDVYAFLCNQCSPSSNVIQQIVNSL---GGQFGT--LWIDVEQCS---GCWG-DVNDNAAFVAEAVQTAAS-LGVTVGVYSSLGEWPQTVGSL-SSLSSYPQWYAHYDGVAASQYGGWDNPEMKQYVGNTNECGV--SVDLDYYG--------------
    ----ELGIDVSSATSQSQWSCLAQKN-QRAIIQVWSGGYGMNNGVVSAIQAAQNAGFQVDLYAFLCNQCSPSSNVIQQIVSKIKQSGVSFGT--LWIDVEQCS---GCWG-STSANAAFVVEAVQTAAS-LGVRVGVYSSSGEWPQTVGTL-TSLSSYPQWYAHYDGVPAGQYGGWNNPEMKQYVGNTNQCGV--SVDLDFYG--------------
    ----TYGVDL------AGFQCLVGKGF-FAIVRCYMSSGGIDPNCASSVSAAWAGGMTVDLYLFPCFSCG----SLVQFAQS---NGVNFGK--IWLDIEGPG---TYWG-DQGANQQFFEGLVQGL--S-GVSVGIYTSESQWSPIMGDY-SGGSNFPLWYANYDGSPN-PFGGWSTPTMKQFDDPSN-CGI--GIDENWIG--------------
    ----GTGIDISSPTSKTQWSCLAKQN-TKAIIQVWSGGYGYNTNIASSVSAAKSAGIQVDLYAFLCSQCSPSSSAIKTLVSNLRSQNVEFGT--LWIDVEQCS---NCWG-STSTNAQFVVEAVQTAQQ-LGVSVGVYSSIGEWSQTVGSL-NSLSSFPLWYAHYDNVPASQFGSWSSPAMKQYAGNTQQCGV--SVDLDFFQ--------------
... $\rightarrow$ Contact maps $\rightarrow$ Structure
!(figures/1jfx_cmap.png) <<height:400px;transparent>>
!(figures/1jfx.png)      <<height:400px;transparent>>

====*

## What is coevolution?

ADD DETAILS HERE!

====

# Mutual information
(naïve attempt) 

# $MI_{ij} = \sum f_{ij} (A_i, A_j) \ln \left ( \frac{f_{ij}(A_i,A_j)}{f_i(A_i)f_j(A_j)} \right )$

$f_i, f_{ij}$ are observed frequencies and co-frequencies respectively. 

works poorly due to transitivity. e.g. A-B and B-C, this model predicts A-C.

&& [Covariation of mutations in the V3 loop of human immunodeficiency virus: an information theoretic analysis](http://www.pnas.org/content/110/39/15674.abstract) </br> Korber, Farber, Wolpert, & Lapedes

====*

### Maximum-entropy model / Markov Random Field
# $P(X=x) = \frac{1}{Z} \exp \left ( \sum_{i=1}^L \left [ v_i(x_i) + \sum_{j>i}^L w_{ij}(x_i,x_j) \right ] \right ) $
$X_i$ is a random variable that represents 
amino acid composition at position $i$ from MSA.

##### Learned parameters
$v_i$ encodes individual propensity of each amino acid at position $i$
$w_{ij}$ statistical coupling of amino acid propensities between positions $i,j$

&& [Learning generative models for protein fold families](http://www.ncbi.nlm.nih.gov/pubmed/21268112) </br> Balakrishnan, Kamisetty, Carbonell, Lee, and Langmead

====*
# DCA

&& [Direct-coupling analysis of residue coevolution captures native contacts across many protein families](http://www.pnas.org/content/108/49/E1293) </br> Morcosa, Pagnanib, Lunta, Bertolinoc, Marksd, Sandere, Zecchinab,  Onuchica, Hwaa, and Weigt
====*
# PSICOV

&& [PSICOV: precise structural contact prediction using sparse inverse covariance estimation](http://bioinformatics.oxfordjournals.org/content/28/2/184.abstract) </br> Jones, Buchan, Cozzetto, and Pontil


====*

## GREMLIN

Optimize the _pseduolikelihood_ of $v,w$
# $ P(v,w | D) = \sum_{n=1}^N \sum_{i=1}^L \log P (x_i^n | x_{i'}, v, w)$
# $ P (x_i^n | x_{i'}, v, w) = \frac{1}{Z_i} \exp \left ( v_i(x_i^n) + \sum_{j=1,j \neq i}^L w_{ij}(x_i^n,x_j^n) \right )$

Regularization (and priors), prevent overfitting
## $R(v,w) = \lambda_v ||v||^2 + \sum_{i,j} \lambda_w^{i,j} || w_{i,j} || ^2$

&& [Assessing the utility of coevolution-based residue–residue contact predictions in a sequence- and structure-rich era](http://www.pnas.org/content/110/39/15674.abstract) Kamisetty, Ovchinnikova, and Baker.
====
# Target dataset

150 monomeric proteins $50< N <250$ residues
====*
# Data pipeline

Download, parse, and clean PDB.
Build FASTA and original contact map.
Align each FASTA using `HHBLITS`*.
Score alignments with GREMLIN$\dagger$.
Build contact maps from GREMLIN.
(optional) Optimize contact map score with RF.
Fold coarse-grained protein from contact map.

&& * `hhblits -i input.seq -n 4 -diff inf -cov 75 -e 0.0000000001` </br> $\dagger$ Dockerize GREMLIN's MATLAB for maximum performance.
====*
## Performance measurements

_Accuracy_    : Predictions that are correct : $(TP+TN)/(TP+FP+FN+TN)$
_Specificity_ : Non-contacts identified :  $TN / (TN+FP)$
_Precision _: Contacts identified that are true : $TP/(TP+FP)$
_Sensitivity_ : True contacts identified :  $TP / (TP+FN)$

ROC curves measure *Sensitivity* vs *Specificity*.

False positives (FP) are worse than false negatives (FN).

We measure *Precision* vs *Sensitivity*.

====

# Scoring
For a given protein and alignment GREMLIN gives $(N,N,21,21)$ tensor.

Reduce GREMLIN's tensor output:
# $S_{N,N,21,21} \rightarrow S_{N,N,20,20} \rightarrow  G_{N,N} \rightarrow G^{\text{APC}}_{N,N} = g$

Drop information about gaps.
Compute the Frobenius norm over each position.
Apply "average product correlation" (necessary!)

====*

## Top score model
Rank sort top diagonal of $g$, take top $\alpha N$ contacts.
Typically values for $\alpha \in [0.1, 3.0]$.

!(figures/1a3a/1a3a_cartoon.png) <<width:500px; transparent>> [1a3a](http://www.rcsb.org/pdb/explore.do?structureId=1a3a) IIA MANNITOL FROM ESCHERICHIA COLI
!(figures/1a3a/1avs_cartoon.png) <<width:500px; transparent>> [1avs](http://www.rcsb.org/pdb/explore.do?structureId=1avs) CALCIUM-SATURATED N-TERMINAL DOMAIN OF TROPONIN C
====*
Example proteins, GREMLIN APC corrected score
!(figures/1a3a_GREMLIN.png) <<height:500px; transparent>> [1a3a](http://www.rcsb.org/pdb/explore.do?structureId=1a3a) IIA MANNITOL FROM ESCHERICHIA COLI
!(figures/1avs_GREMLIN.png) <<height:500px; transparent>> [1avs](http://www.rcsb.org/pdb/explore.do?structureId=1avs) CALCIUM-SATURATED N-TERMINAL DOMAIN OF TROPONIN C
====*

## GREMLIN Predictions
!(figures/GREMLIN_only_Acc_Pre.png) <<height:750px; transparent>>
% ~/git-repo/GREMLIN_RF/analysis/plot_stats.py

====

Can we do better?

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
## Improved RF model Predictions
!(figures/GREMLIN_RF_Acc_Pre.png)  <<height:750px; transparent>>
====

## Contact map vs cutoff length (1a3a)
!(figures/1a3a/animated_1a3a.gif) <<width:1200px; transparent>>
====*
## Contact map vs cutoff length (1avs) 
!(figures/1a3a/animated_1avs.gif) <<width:1200px; transparent>>

====
## Folding simulations
$C_\alpha$ coarse-grained MD simulation

Unbiased estimate of contact map $\rightarrow$ fold.

No prior knowledge (ROSETTA fragments, SS pred., etc...).

Potential = Backbone + smoothed well with range ~ $8\AA$


====*
### Rapid collapse to contact potential
!(figures/folding/folding_1a3a.gif) <<height:600px; transparent>>
$C_\alpha$ coarse-grained MD simulations
====*
## Folding simulations, $Q_\text{native}$
!(figures/pairplot_Q.png) <<height:750px; transparent>>
====*
## Folding simulations, RMSD
!(figures/pairplot_RMSD.png) <<height:750px; transparent>> RMSD
====*
Decompose Random Forest features
!(figures/SVD_RF.png) <<height:725px; transparent>> SVD of Decision Tree weights
====

### Predicted contacts are closer to true contacts
## $d(x) = \min_{y \in \text{TP}}|\vec{x}-\vec{y}|$

!(figures/local_structure_distance.png) <<height:450px; transparent>>  Distance from predicted contact to true contact 
!(figures/FP_distance.png) <<height:450px; transparent>> Average distance for all proteins
====
## Improvement in folding
!(figures/folding/Q_avg.png) <<height:700px; transparent>>
====


Additional value 

$g$ can be used as an effective Hamiltonian for evolutionary movement

====

# Future work
Convolutional neural networks
!(figures/example_CNN.png) <<transparent>>

Enhanced structure prediction eg. ROSETTA?
Disambiguation of intra/inter predictions
Binding partners

====

## Thank you.

####  *Robert Best* (NIH/NIDDK)
Wenwei Zheng
*Travis Hoppe*
*Pengfei Tian*
Jan Domanski
Mathias Bellaiche

#### Jeff Gray (John Hopkins, Chemical Engineering)
Julia Joehler















