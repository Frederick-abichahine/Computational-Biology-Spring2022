# Computational-Biology-Final

---

## Goal
Improving the results of multi-network diffusion (mND) by adjusting scores using geneSurrounder (GS).

---

## Methods & Materials

### Required Libraries

1. Multi-Network Diffusion (mND): https://doi.org/10.1093/bioinformatics/btz652 
>> mND uses multi-layer Network Diffusion to find gene networks that contain high scoring genes. It considers two or more "layers" of genome-wide scores and an interactome.

```
## Install using:
install.packages(c("devtools", "igraph", "parallel"))
library(devtools)
install_github("emosca-cnr/mND", build_vignettes = TRUE)
```


2. GeneSurrounder (GS): https://doi.org/10.1186/s12859-019-2829-y 
>> Gene Surrounder is an analysis method that integrates expression data and network information in a novel procedure to detect genes that are sources of dysregulation on the network.

```
get from github https://github.com/sahildshah1/gene-surrounder
```

 3. Dependencies - other Libraries:
 
```
require(isma)
require(mND)
require(limma)
require(igraph)
require(pcaPP)
require(SciViews)
require(tidyverse)
```


---

### Strategy

#### 1. Pre-Adjustment of Scores: 
a. Run mND alone - without GS ==> Save Results mND_scores   
b. <Future improvements?>   

#### 2. Adjusting Scores:
a. Run GS & obtain p-decay (ideally, we wish to obtain p-fisher).   
b. Run mND on data with adjusted scores by p-decay   

```
Adjusted_scores = -log_10(p-decay)*score_original
```

#### 3. Analyzing & Comparing:
a. Compare Results & Validate any Improvement
b. Visualization & Reporting

---

### Future Plans:
> Apply heuristics (Genetic Algorithm?)

---



---



