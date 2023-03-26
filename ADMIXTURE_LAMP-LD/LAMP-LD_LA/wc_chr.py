#!/bin/python

import fileinput
from contextlib import closing
import subprocess

with closing(fileinput.FileInput("../scripts/bash_scripts/wc_ppmi.sh", inplace=True)) as file:
        for line in file:
                num=[0]
                for i in num:
                        if i<1:
				print(line.replace("ppmi_gwas_preQC.vcf.gz", "ppmi_gwas_preQC_chr"+repr(i)+".vcf.gz")
                                subprocess.call(["sbatch", "../scripts/bash_scripts/wc_ppmi.sh"])
                                i+=1
			else i<23:
                                print(line.replace("ppmi_gwas_preQC_chr"+repr(i-1)+".vcf.gz", "ppmi_gwas_preQC_chr"+repr(i)+".vcf.gz"))
                                subprocess.call(["sbatch", "../scripts/bash_scripts/wc_ppmi.sh"])
                                i+=1
