# Study of the effect of EZH2 inhibition in hESC
## **Gene regulatory network (GRN) inference:**

### GRN inference and _in silico_ perturbation simulation by [CellOracle](https://morris-lab.github.io/CellOracle.documentation/) (v0.16.0)

Codes: 
_utils/CellOracle_scripts/*.py_

${\color{red}Important:}$

For CellOracle dependencies it's worth making a [pyenv](https://github.com/pyenv/pyenv) virtual environment with Python version 3.8 (or above). CellOracle installation is documented thouroughly [here](https://morris-lab.github.io/CellOracle.documentation/installation/index.html).

**Necessary inputs for hESC EZH2i CellOracle analysis:** 
- scRNA-Seq

In our case scRNA-Seq was initially processed by [Seurat](https://satijalab.org/seurat/articles/get_started.html) (see _utils/sc_rna_workflow.R_). As CellOracle needs scRNA-Seq data in h5ad format ([AnnData](https://anndata.readthedocs.io/en/latest/)), at first we must convert Seurat data object into h5ad. See: _utils/create_h5ad_from_seurat.R_ script.
Additionally, we must export the raw read count tables of the Seurat objects (e.g. seurat_obj@assays$RNA@counts), as CellOracle object needs the raw count table. 

Seurat version: Seurat_5.0.2

- scATAC-Seq

For CellOracle's base GRN inference, we must quantify the co-accessibility sites by [cicero](https://cole-trapnell-lab.github.io/cicero-release/). The cicero outputs (connections and the peaks used by cicero) are mandatory inputs for our _celloracle_atac_prep_ Python file. The cicero analysis (_utils/cicero.R_) needs computational resources (worth running on Uppmax). Also, base GRN configuration script needs a reference genome in fasta format (in our case, it was a [UCSC hg38](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/) fasta). For the cicero analysis, scATAC-Seq data were processed by standard scATAC Seurat workflow (e.g. [here](https://satijalab.org/seurat/archive/v3.1/atacseq_integration_vignette)) and afterwards only the most variable scATAC-Seq peaks were considered in the co-accessibility estimation. 

cicero version: cicero_1.3.9, with monocle3_1.3.5 + _celloracle_atac_prep.py_ needs available BEDTools in the path  

- pseudotime analysis

Ideally, it is estimated by scanpy (see here: [Pseudotime calculation](https://morris-lab.github.io/CellOracle.documentation/tutorials/pseudotime.html)). Pseudotime calculation (_celloracle_pseudotime.py_) needs the processed scRNA-Seq object in h5ad format. 

Optional input:

The GRN construction might be supported by manually adding TF-target gene pairs (e.g. statistically significant TF-target gene predictions of [Pando](https://github.com/quadbio/Pando))

**CellOracle run:**

If the inputs are properly generated, we must run the individual Python scripts sequentially: 

1. Pseudotime calculation
2. scATAC-Seq preprocessing ( = cicero output processing)
3. base GRN construction
4. _in silico_ perturbations





