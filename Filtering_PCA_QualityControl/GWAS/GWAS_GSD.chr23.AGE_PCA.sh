 #!/bin/bash
#MSUB -A b1042
#MSUB -q genomics
#MSUB -l walltime=24:00:00
#MSUB -l nodes=1:ppn=1

/projects/b1049/bernabe/GWAS_GSD/snptest_v2.5.2/snptest_v2.5.2 -data /projects/b1049/bernabe/GWAS_GSD/impute_format/chile.X_nonPAR.gen.gz /projects/b1049/bernabe/GWAS_GSD/new_sample_GWAS_wPCA.txt -o /projects/b1049/bernabe/GWAS_GSD/snptest/age_PCA/chile.chrX.age_PCA.snptest -frequentist 1 -method newml -pheno phenotype -assume_chromosome X -hwe -cov_names AGE PC1 PC2 PC3 PC4 PC5 
