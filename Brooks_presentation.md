# Enhancing the coevolutionary signal
----------
### [Travis Hoppe](http://thoppe.github.io/)
NIH/NIDDK/LCP Postdoctoral Fellow

====*

# Outline

### Alignment / coevolution
### *Score functions*
### Structure prediction

====

### Sequence to structure prediction
Sequence $\rightarrow$ Multiple Sequence Alignment $\rightarrow$ ...
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

Observation: Homologous proteins impose strong constraints 
on their sequence variability.

Assume: If two residues form a contact, a destabilizing substitution at one position is expected to be compensated by a substitution of the other position over the evolutionary timescale, in order for the residue pair to maintain attractive interaction.

====

# Mutual information
(naïve attempt) 

# $MI_{ij} = \sum f_{ij} \ln \left ( \frac{f_{ij}}{f_i f_j} \right )$

$f_i, f_{ij}$ are observed frequencies and co-frequencies respectively. 

works poorly due to transitivity. e.g. A-B and B-C, this model predicts A-C.

&& [Covariation of mutations in the V3 loop of human immunodeficiency virus: an information theoretic analysis](http://www.pnas.org/content/110/39/15674.abstract) </br> Korber, Farber, Wolpert, & Lapedes

====*

### Maximum-entropy model / Markov Random Field
# $P(X=x) = \frac{1}{Z} \exp \left ( \sum_{i=1}^L \left [ v_i(x_i) + \sum_{j>i}^L w_{ij}(x_i,x_j) \right ] \right ) $
Least-constraint model that matches marginal distributions of $f_i$ and $f_{ij}$.
%$X_i$ is a random variable that represents 
%amino acid composition at position $i$ from MSA.

Brute force computational complexity of pairwise $Z$ is $O(N^2 21^N)$.

##### Learned parameters
$v_i$ encodes individual propensity of each amino acid at position $i$
$w_{ij}$ statistical coupling of amino acid propensities between positions $i,j$

&& [Learning generative models for protein fold families](http://www.ncbi.nlm.nih.gov/pubmed/21268112) </br> Balakrishnan, Kamisetty, Carbonell, Lee, and Langmead

====*
## DCA (direct coupling analysis)*
Focus on high MI pairs, use reduced two-residue systems.

## PSICOV $\dagger$
Compute pairwise covariance $\text{Cov}(X,Y)=E(XY)-E(X)E(Y)$ over all pairs of sites for all residues from MSA. Invert the matrix with tricks to avoid singular matrices (assume sparsity, most entries are zero in inversion).


&& *[Identification of direct residue contacts in protein–protein interaction by message passing] (http://www.pnas.org/content/106/1/67.full) </br> Weigt, White, Szurmant, Hoch, and Hwa </br> $\dagger$ [PSICOV: precise structural contact prediction using sparse inverse covariance estimation](http://bioinformatics.oxfordjournals.org/content/28/2/184.abstract) </br> Jones, Buchan, Cozzetto, and Pontil

% && [Direct-coupling analysis of residue coevolution captures native contacts across many protein families](http://www.pnas.org/content/108/49/E1293) </br> Morcos, Pagnani, Lunta, Bertolinoc, Marksd, Sandere, Zecchinab,  Onuchica, Hwaa, and Weigt



====*

## GREMLIN
Optimize the _pseudolikelihood_ of $v,w$
# $ P(v,w | D) = \sum_{n=1}^N \sum_{i=1}^L \log P (x_i^n | x_{i'}^n, v, w)$
# $ P (x_i^n | x_{i'}^n, v, w) = \frac{1}{Z_i} \exp \left ( v_i(x_i^n) + \sum_{j=1,j \neq i}^L w_{ij}(x_i^n,x_j^n) \right )$
Models conditional distribution of the original joint distribution 
instead of the joint distribution itself. Can add regularization 
to prevent overfitting and prior knowledge.

$v_i$ encodes individual propensity of each amino acid at position $i$
$w_{ij}$ statistical coupling of amino acid propensities between positions $i,j$

%## $R(v,w) = \lambda_v ||v||^2 + \sum_{i,j} \lambda_w^{i,j} || w_{i,j} || ^2$
&& [Assessing the utility of coevolution-based residue–residue contact predictions in a sequence- and structure-rich era](http://www.pnas.org/content/110/39/15674.abstract) Kamisetty, Ovchinnikova, and Baker.
====
# Target dataset

Pfam families with $\ge 1000$ sequences with high resolution $\le 1.9 \AA$.
150 monomeric proteins $50< N <275$ residues; [diverse set](http://bioinformatics.oxfordjournals.org/content/suppl/2011/11/29/btr638.DC1/target.txt).

    PDB-ID  Pfam-ID Nseq    Length  Description
    ========================================================================================
    1GUUA   PF00249 10393   50      Myb-like DNA-binding domain
    1BRFA   PF00301 1430    53      Rubredoxin
    1AAPA   PF00014 2256    56      Kunitz/Bovine pancreatic trypsin inhibitor domain
    1JO8A   PF00018 6287    58      SH3 domain
    1KU3A   PF04545 8439    61      Sigma-70, region 4
    1M8AA   PF00048 1062    61      Small cytokines (intecrine/chemokine), interleukin-8 like
    1C9OA   PF00313 6807    66      'Cold-shock' DNA-binding domain
    1VFYA   PF01363 1645    67      FYVE zinc finger
    1CTFA   PF00542 2390    68      Ribosomal protein L7/L12 C-terminal domain
    1KW4A   PF07647 1192    70      SAM domain (Sterile alpha motif)
    1CC8A   PF00403 9383    72      Heavy-metal-associated domain
    1ATZA   PF00092 7567    75      von Willebrand factor type A domain
    1TIFA   PF05198 1947    76      Translation initiation factor IF-3, N-terminal domain
    1H98A   PF00037 10421   77      4Fe-4S binding domain
    1T8KA   PF00550 20685   77      Phosphopantetheine attachment site
    1BDOA   PF00364 11826   80      Biotin-requiring enzyme
    1AVSA   PF00036 13234   81      EF hand
    1CXYA   PF00173 3200    81      Cytochrome b5-like Heme/Steroid binding domain
    1I71A   PF00051 1082    83      Kringle domain
    1ABAA   PF00462 5749    87      Glutaredoxin
    1DSXA   PF02214 1372    87      K+ channel tetramerisation domain
    1SMXA   PF10150 2203    87      Ribonuclease E/G family
    1NPSA   PF00030 1153    88      Beta/Gamma crystallin
    1PCHA   PF00381 3344    88      PTS HPr component phosphorylation site
    1VJKA   PF02597 3283    88      ThiS family
    1FNAA   PF00041 17137   91      Fibronectin type III domain
    1G9OA   PF00595 14944   91      PDZ domain (Also known as DHR or GLGF)
    1FK5A   PF00234 3346    93      Protease inhibitor/seed storage/LTP family

====*
# Data pipeline

Download, parse, and clean PDB.
Build FASTA and reference contact map.
Align each FASTA using `HHBLITS`*.
Score alignments with GREMLIN$\dagger$.
Build contact maps from GREMLIN.
(optional) Optimize contact map score with RF.
Fold coarse-grained protein from contact map.

&& * `hhblits -i input.seq -n 4 -diff inf -cov 75 -e 0.0000000001` </br> $\dagger$ Dockerize GREMLIN's MATLAB for maximum performance.
====

# Scoring
For a given protein and alignment GREMLIN gives $(N,N,21,21)$ tensor.

Reduce GREMLIN's tensor output:
# $S_{N,N,21,21} \rightarrow S_{N,N,20,20} \rightarrow  G_{N,N} \rightarrow G^{\text{APC}}_{N,N} = g$
Drop information about gaps.
Compute the Frobenius norm over each position.
Subtract _average product correlation_*, structural vs. shared ancestry
#### $APC(a,b) = MI(a,\overline{x}) MI(b,\overline{x}) / \overline{\text{MI}}$

&& *[Mutual information without the influence of phylogeny or entropy dramatically improves residue contact prediction](http://bioinformatics.oxfordjournals.org/content/24/3/333.abstract) <br/> Dunn, Wahl, and Gloor

====*

## Top score model
Rank sort top diagonal of $g$, take top $L N$ contacts.
Typically values for $L \in [0.1, 3.0]$.

!(figures/1a3a/1a3a_cartoon.png) <<width:500px; transparent>> [1a3a](http://www.rcsb.org/pdb/explore.do?structureId=1a3a) IIA MANNITOL FROM ESCHERICHIA COLI
!(figures/1a3a/1avs_cartoon.png) <<width:500px; transparent>> [1avs](http://www.rcsb.org/pdb/explore.do?structureId=1avs) CALCIUM-SATURATED N-TERMINAL DOMAIN OF TROPONIN C
====*
Example proteins, GREMLIN APC corrected score
!(figures/1a3a_GREMLIN.png) <<height:500px; transparent>> [1a3a](http://www.rcsb.org/pdb/explore.do?structureId=1a3a) IIA MANNITOL FROM ESCHERICHIA COLI
!(figures/1avs_GREMLIN.png) <<height:500px; transparent>> [1avs](http://www.rcsb.org/pdb/explore.do?structureId=1avs) CALCIUM-SATURATED N-TERMINAL DOMAIN OF TROPONIN C
====*
## Performance measurements

_Accuracy_    : Predictions that are correct : $(TP+TN)/(TP+FP+FN+TN)$
_Specificity_ : Non-contacts identified :  $TN / (TN+FP)$
_Precision _: Contacts identified that are true : $TP/(TP+FP)$
_Sensitivity_ : True contacts identified :  $TP / (TP+FN)$

ROC curves measure *Sensitivity* vs *Specificity*.

False positives (FP) are worse than false negatives (FN).

We measure *Precision* vs *Sensitivity*.
====*

## GREMLIN Predictions
!(figures/GREMLIN_only_Acc_Pre.png) <<height:750px; transparent>>
% ~/git-repo/GREMLIN_RF/analysis/plot_stats.py

====
# Hypothesis:
Local structure can enhance contact prediction.
!(figures/local_structure_zoom.png) <<width:500px; transparent>>

Secondary structure is local (helices, sheets, turns).
====*
## Random forest (RF) score model
Machine learn local pixel maps for contact/non-contact.

Normalize data: subtract mean, scale to unit variance.

Train with extremely random forests*, variant with dropout.

(e)RF's were more robust than traditional shallow learning like SVM.

&& *RF parameters `kernel_window=2`, `n_trees=200`, `ratio_TP_to_TN=20`
====*
## What are Random Forests?
Decision trees are good for simple data, but tend to overfit.
!(figures/example_Dtree.png) <<transparent>>

Random forests are multiple decision trees with 1] "random splits", 
2] selective subsets, (each tree only gets to see a subset of the data). 
This increases individual bias but the average corrects for overfitting.


====*
## Improved RF model Predictions
!(figures/GREMLIN_RF_Acc_Pre.png)  <<height:750px; transparent>>
====

## Contact map vs cutoff length (1a3a)
!(figures/1a3a/animated_1a3a.gif) <<width:1200px; transparent>>
====*
## Contact map vs cutoff length (1avs) 
!(figures/1a3a/animated_1avs.gif) <<width:1200px; transparent>>

====
# Folding simulations
$C_\alpha$ coarse-grained MD simulation

Unbiased estimate of contact map $\rightarrow$ fold.

No prior knowledge (ROSETTA fragments, SS pred., etc...).

Potential = Backbone + smoothed well with range ~ $8\AA$.

%Poor RMSD due to loose interations.
====*
### Rapid collapse to contact potential
!(figures/folding/folding_1a3a.gif) <<height:600px; transparent>>
$C_\alpha$ coarse-grained MD simulations
====*
## Folding simulations, $Q_\text{native}$
!(figures/pairplot_Q.png) <<height:750px; transparent>>
%====*
%## Folding simulations, RMSD
%!(figures/pairplot_RMSD.png) <<height:750px; transparent>> RMSD
====
# Features of the RF model
scientific insight beyond predictive capability
====*
### Predicted contacts are closer to true contacts
## $d(x) = \min_{y \in \text{TP}}|\vec{x}-\vec{y}|$

!(figures/local_structure_distance.png) <<height:450px; transparent>>  Distance from predicted contact to true contact 
!(figures/FP_distance.png) <<height:450px; transparent>> Average distance for all proteins
====*
## Improvement in folding
In $C_\alpha$ potential, more contacts $\approx$ better RF fold.
!(figures/folding/Q_avg.png) <<height:675px; transparent>>
====*
Random Forest features (central difference most important)
!(figures/SVD_RF.png) <<height:725px; transparent>> SVD of Decision Tree weights
====
## Future work & Extensions

Convolutional neural networks to improve prediction:
!(figures/example_CNN.png) <<height:300px; transparent>>


$g$ can be used as an effective Hamiltonian for evolutionary movement.

====*
## Future work & Extensions

Enhanced structure prediction, ROSETTA et. al.


Disambiguation of intra/inter predictions.


Estimation of binding partners and hetrodimers from $g$.
====

# Thanks, you.

####  *Robert Best* (NIH/NIDDK)
Wenwei Zheng
*Travis Hoppe*
*Pengfei Tian*
Jan Domanski
Mathias Bellaiche

#### Jeff Gray (John Hopkins, Chemical Engineering)
Julia Joehler















