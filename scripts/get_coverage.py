################################################################################
#
# Get coverage from low coverage bam files used in phase I
#
################################################################################


import ftplib
import re
from subprocess import call
import time
from os.path import expanduser

ftp = ftplib.FTP("ftp.1000genomes.ebi.ac.uk", "anonymous", "anonymous")

# directory where population files are
popdir = expanduser("~") + "/hla_tools/1kg/data/pops/"

# region encompassing exons 2 and 3 of HLA-A, -B, -C, -DRB1, -DQB1 plus 500 on
# each side
region = "6:29910034-32633344"
positions = expanduser("~") + "/hla_tools/1kg/data/ARS_exons.bed"

for pop in ["CEU", "CHB", "JPT", "YRI", "LWK", "ASW", "MXL", "TSI", "GBR", "FIN", "CHS","PUR", "CLM"]:
    # read each line of file with IDs of samples from this pop
    count = 1
    with open(popdir + "overlap_kg_pag_" + pop + ".txt") as popfile:
        for i in popfile:
        # i = individuals from this pop
            i = i.rstrip() # remove \n
            # get bam file for this individual
            files = []
            ftp.dir("vol1/ftp/phase1/data/" + i + "/alignment/", files.append)
            for f in files:
                m = re.search(i + "\.mapped\.\w+\.\w+\." + pop +
                              "\.low_coverage\.\d{8}\.bam$", f)
                if m:
                    bamfile = m.group()
                    print str(count) + " - Found bam file for individual " + i
                    command = ["samtools view -b ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/data/"+
                               i + "/alignment/" + bamfile + " " + region +
                               " | genomeCoverageBed -ibam 'stdin' -bg " +
                               " | intersectBed -a 'stdin' -b " + positions +
                               " > " + expanduser("~") +
                               "/hla_tools/1kg/data/coverage/coverage_" +
                               i + ".bg"]
                    call(command, shell = True)
                    time.sleep(3)
            count = count + 1
ftp.quit()


