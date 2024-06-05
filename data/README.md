# Data provenance

Last upadted 2023/10/01

## Dimensions (r / c)

| File                                          | all_genes  | protein_coding |
| --------------------------------------------- | ---------- | -------------- |
| `DCP_rhythm.rds`                              |            |                |
| `DR_rhythm_params.csv`                        | UNUSED     | UNUSED         |
| `DR_rhythm_params_rounded.csv`                | 2724 / 31  | 2519 / 31      |
| `genemap.csv`                                 | 23631 /  8 | 16754 / 8      |
| `LONG__vs__SHORT-DESeq2-results-all-data.csv` | 24317 / 14 | 16864 / 14     |
| `metadata.csv`                                | 181 / 6    | 181 / 6        |
| `normalized_counts_all.LONG.csv.gz`           | 24317 / 63 | 16864  / 63    |
| `normalized_counts_all.SHORT.csv.gz`          | 24317 / 58 | 16864  / 58    |

## DiffCircaPipeline-derived data

`DCP_rhythm.rds` copied from `results/122-DiffCirca-protein-only/results/`

`DR_rhythm_params.csv` copied from `results/122-DiffCirca-protein-only/results/`

`DR_rhythm_params_rounded.csv` copied from `results/122-DiffCirca-protein-only/results/`


## Genemap
```R
genemap.prot <- genemap %>% filter(gene_biotype == 'protein_coding')
write.csv(genemap.prot, "data/protein_coding/genemap.prot.csv", row.names = FALSE)
```


## DESeq2-derived data

`normalized_counts_all.LONG.csv` copied from `results/070-DE-Analysis.R` or `results/073-DE_analysis-protein-only` 

 `normalized_counts_all.SHORT.csv` copied from `results/070-DE-Analysis.R` or `results/073-DE_analysis-protein-only` 

