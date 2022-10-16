# RFS-SNP-T2D
 Repository for the project "Recommended Food Score (RFS) ‑ SNP interaction on Type 2 Diabetes Survival Time"

## Objectives
- Identify SNPs, genes, and pathways having significant interaction effect with (modified)RFS on T2D survival time.
- Significant interaction implies that RFS level can affects mechanism of the pathway.
- Unlike previous studies, focused on “survival time” rather than incidence and “interaction effect” between RFS and SNP rather than main effect
- Performed pathway analysis and identified two KEGG pathways related with type 2 diabetes with consideration of RFS using PLINK and MAGMA

## Methods
1. Perform Cox regression by plinkCox in R with the model 
$$\lambda(t|X) = \lambda_0(t) \times \exp{(\beta_0 + X_{cov}^\top\beta_{cov} + \beta_{SNP}X_{SNP} + \beta_{RFS}X_{RFS} + \beta_{inter}X_{SNP}X_{RFS})}$$
2. Using p-value of $\beta_{inter}$ computed by above Cox regression, perform pathway analysis by MAGMA and GSA-SNP2
3. Identify significant pathways and provide their biological interpretation



## :file_folder: code

### :page_facing_up: preprocessing.R
- Preprocessing Korea Association REsource (KARE) Consortium data

### :page_facing_up: RFS_calculation.R
- Recommended Food Score calculation for KARE participants

### :page_facing_up: Coxmodel.R
- Fit Cox proportional hazards model on type 2 diabetes survival time without SNP data

### :page_facing_up: Coxmodel_SNP.R
- Fit Cox proportional hazards model on type 2 diabetes survival time with SNP and RFS*SNP interaction

### :page_facing_up: pathway.R
- Find pathways having p-values smaller than FDR thresholds after performing MAGMA and GSA-SNP2

## :file_folder: pathway
- Pathway analysis results based on various p-values thresholds
