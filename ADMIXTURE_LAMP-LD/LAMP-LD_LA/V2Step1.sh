#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 20:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=10
#SBATCH --mem=20000

for K in {1..20}
do
echo "#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 48:00:00
#SBATCH -N 5
#SBATCH --ntasks-per-node=28
#SBATCH --mail-user=dpatel4@imsa.edu
#SBATCH --mail-type=END
#SBATCH --mem=20000

/projects/b1049/genetics_programs/admixture_linux-1.3.0/admixture --cv myMergedAndPrunedData.bed $K -j28 | tee log${K}.out" > admixture_5cores.$K.cv.sh
done