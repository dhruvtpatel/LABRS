 #!/bin/bash
#MSUB -A b1042
#MSUB -q genomics
#MSUB -l walltime=24:00:00
#MSUB -l nodes=1:ppn=1

cd 

/projects/b1049/bernabe/GWAS_GSD/snptest_v2.5.2/snptest_v2.5.2 -data /projects/b1049/bernabe/GWAS_GSD/impute_format/chile.21.gen.gz /projects/b1049/bernabe/GWAS_GSD/new_sample_GWAS_wPCA.txt -o /projects/b1049/bernabe/GWAS_GSD/snptest/age_PCA/chile.chr21.age_PCA.snptest -frequentist 1 -method score -pheno phenotype -assume_chromosome 21 -hwe -cov_names AGE PC1 PC2 PC3 PC4 PC5 
