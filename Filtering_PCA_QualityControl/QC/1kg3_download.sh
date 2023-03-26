cd $refdir

pgen=https://www.dropbox.com/s/j72j6uciq5zuzii/all_hg38.pgen.zst?dl=1
pvar=https://www.dropbox.com/s/vx09262b4k1kszy/all_hg38.pvar.zst?dl=1
sample=https://www.dropbox.com/s/2e87z6nc4qexjjm/hg38_corrected.psam?dl=1

wget $pgen
mv 'all_hg38.pgen.zst?dl=1' all_phase3.pgen.zst
plink2 --zst-decompress all_phase3.pgen.zst > all_phase3.pgen

wget $pvar
mv 'all_hg38.pvar.zst?dl=1' all_phase3.pvar.zst

wget $sample
mv 'hg38_corrected.psam?dl=1' all_phase3.psam
