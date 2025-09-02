# __[The Alzheimer's Cell Atlas](https://taca.lerner.ccf.org)__

This is the official code for manuscript 'A Scalable Human Brain Cell Atlas Generative AI Model Dissects the Cellular Immune Codes and Cognitive Resilience'.<br>

Here, we present the largest human brain single-nucleus RNA-seq atlas to date (~14 million nuclei across 31 regions from 2,239 deeply phenotyped postmortem samples, including cognitively normal (CN), AD (clinico-pathologically confirmed), cognitively resilient (CR), and primary age-related tauopathy (PART)).Integrating deep clinical and genomic data under AI-powered single-cell model, this comprehensive brain cell atlas deciphers brain-immune interactions underlying resilience and offers an open-access platform to accelerate therapies for neurodegeneration.<br>

### __Installation__

```bash
git clone https://github.com/ChengF-Lab/TACA
cd TACA
conda create -n TACA python=3.11
conda activate TACA
pip install .
```
### __Integration & Analyses__

- TACA
  - utils
  - pyproject.toml
  - QC
    - CxG_qc.ipynb
    - GSE_qc.ipynb
    - CxG_metadata.ipynb
  - Preintegration
    - preintegration.ipynb
  - Integration_Evaluation
    - integration_unsupervised.ipynb
    - integration_semisupervised.ipynb
    - evaluation.ipynb

#### __TACA - Quality Control and Preprocessing__
&#128073;QC/CxG_qc.ipynb, GSE_qc.ipynb<br>
- Read data into AnnData object;<br>
- Convert original gene representation (e.g., symbol, NCBI_ID) into ensemblID;<br>
- Remove ambient RNA if necessary;<br>
- Universal quality control for all datasets;<br>
- Adaptive quality control for each dataset;<br>
- Remove doublet if any using [Scrublet](https://github.com/swolock/scrublet);<br>
- Add sampleID and original cell type annotations if any into AnnData object.<br>

&#128073;QC/CxG_metadata.ipynb<br>
- For each dataset, we created two types of metadata information, i.e., 'standard metadata' and 'original metadata';<br>
- 'standard metadata' are the same for each dataset. It includes information including: 'SampleID', 'DonorID', 'Organism', 'Age', 'Sex', 'Race', 'Ethnicity', 'Tissue', 'Region', 'Disease Type', and 'Assay'.
- 'original metadata' are adaptive and represent any additional meta information provided for each dataset.
- For multiple datasets with original data from [Gene Expression Omnibus](https://www.ncbi.nlm.nih.gov/gds), their meta information can be retrieved from [TACA](https://taca.lerner.ccf.org).

#### __TACA - Preintegration: Merge datasets into one AnnData object__
&#128073;Preintegration/preintegration.ipynb<br>
- Import references cell type annotations if any;<br>
- Import 'SampleID' from each dataset into AnnData object;<br>
- Import 'BatchID' from each dataset into AnnData object;<br>
- Import 'Dataset' from each dataset into AnnData object;<br>
- Import additional meta information (user defined) from each dataset into AnnData object;<br>
- Delete all pre-existing cell embeddings if any;<br>
- Merge all datasets;<br>
- Save the merged AnnData object.<br>

#### __TACA - Integration and Evaluation: Integration with either unsupervised ([scVI](https://www.nature.com/articles/s41592-018-0229-2)) or semi-supervised ([scANVI](https://www.embopress.org/doi/full/10.15252/msb.20209620)) model__

&#128073;Integration_Evaluation/integration_unsupervised.ipynb & integration_semisupervised.ipynb<br>
- Compute highly variable genes (HVGs) and differential expressed genes (DEGs, if using semi-supervised model);<br>
- Split adata into adata_train, adata_val, adata_test using train_ratio and test_ratio predefined;<br>
- Tune hyperparameter parameters ([OPTUNA](https://optuna.org)) based on evidence lower bound (ELBO) value, and generate best hyperparameters;<br>
- Generate entire nuclei embedding (and cell type predictions if using semi-supervised model) using best hyperparameters.<br>

&#128073;Integration_Evaluation/evaluation.ipynb<br>
- Evaluate cell type annotations by comparing consistencies between predictions and reference labels (for semi-supervised model only)

### __Dataset Access__

We integrated 21 Alzheimer's Disease Related Dementia (ADRD) datasets from the following four sources:

| [CELLXGENE](https://cellxgene.cziscience.com)   | [Gene Expression Omnibus (GEO)](https://www.ncbi.nlm.nih.gov/geo/index.cgi?qb=dat) | [SYNAPSE](https://www.synapse.org)        | Cheng's Lab   |
|:------:|:---:|:-----------:|:-----------:|
| [30747918](https://pubmed.ncbi.nlm.nih.gov/30747918/)  |  [GSE148822](https://pubmed.ncbi.nlm.nih.gov/33609158/) | [Syn52293417](https://www.cell.com/cell/fulltext/S0092-8674(23)00973-X)    | ADRD_multiome    |
| [34051145](https://pubmed.ncbi.nlm.nih.gov/34051145/)    |  [GSE157827](https://pubmed.ncbi.nlm.nih.gov/32989152/) | [Syn51123521](https://www.nature.com/articles/s41586-024-07871-6) | ALS_multiome |
| [35513515](https://pubmed.ncbi.nlm.nih.gov/35513515/)  |  [GSE160936](https://pubmed.ncbi.nlm.nih.gov/34767070/) | [Syn51105515](https://www.cell.com/cell/fulltext/S0092-8674(24)00234-4)    |     |
| [36007006](https://pubmed.ncbi.nlm.nih.gov/36007006/)    | [GSE174367](https://pubmed.ncbi.nlm.nih.gov/34239132/) |  |  |
| [37217978](https://pubmed.ncbi.nlm.nih.gov/37217978/)  |  [GSE243639](https://pubmed.ncbi.nlm.nih.gov/38245794/) |     |     |
| [37824655](https://pubmed.ncbi.nlm.nih.gov/37824655/)    |  [GSE254569](https://pubmed.ncbi.nlm.nih.gov/40053590/) |  |  |
| [37824663](https://pubmed.ncbi.nlm.nih.gov/37824663/)    |   |  |  |
| [GSE129308](https://www.sciencedirect.com/science/article/pii/S0896627322006006?via%3Dihub)    |   |  |  |
| [SEA-AD_DLPFC](https://www.nature.com/articles/s41593-024-01774-5)    |   |  |  |
| [SEA-AD_MTG](https://www.nature.com/articles/s41593-024-01774-5)    |   |  |  |


### __Citing TACA__
- [Zhou et. al., 2022, The Alzheimer's Cell Atlas (TACA): A single-cell molecular map for translational therapeutics accelerator in Alzheimer's disease](https://pubmed.ncbi.nlm.nih.gov/36254161/)<br>

- [Xu et. al., 2025, A Scalable Human Brain Cell Atlas Generative AI Model Dissects the Cellular Immune Codes and Cognitive Resilience](https://sciety-discovery.elifesciences.org/articles/by?article_doi=10.21203/rs.3.rs-6986171/v1)<br>
