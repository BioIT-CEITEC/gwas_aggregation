########################################
# wrapper for rule: merge_gwas
########################################
import re
import os
import subprocess
from snakemake.shell import shell
shell.executable("/bin/bash")
log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: plink preprocessing \n##\n")
f.close()

command = "Rscript "+os.path.abspath(os.path.dirname(__file__))+"/merge_gwas.R "+\
           str(snakemake.input.assoc) + " " + snakemake.output.combined + " " + snakemake.output.report + " >> " + log_filename + " 2>&1 "

f = open(log_filename, 'a+')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
