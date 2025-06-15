**TASK 1** 
count sequecnes in fastq (many other methods are valid)
```bash
zgrep -c "^+$" ./data/16S_biochar_run2_10perc_sampled/*.gz
```
check seq length:
- open file page by page without unzipping and find the sequence header
- grep zipped file getting rid of header and +lines, get sequences length with awk
```bash
zless ./data/16S_biochar_run2_10perc_sampled/*.gz #find the header

zgrep -v "^@A00618:" ./data/16S_biochar_run2_10perc_sampled/Bch-16S-V3V4-001-2_S1_L002_R1_001.fastq_10perc.fastq.gz | zgrep -v "^+$" |awk '{print length($1)}'
```


**TASK 2**

- 1.
- Per tile variation are normal, and this run is good quality
-  Per base seq. and GC content cannot be stable and normally distributed (It's a short quite homogenous amplicon, except for hypervariable regions, not a full genome!), maybe worth trying to check shoulder peak from R2.
- Seq. length distribution is totallly fine
- Duplication and overrepresented sequecneare normal for amplicons (and the redundancy is also amplified by PCR cycles). Also overrepresented sequecnes are showing up where the primers were designed, which need to be not variable region of the locus (otherwise primers would not work).
- If curious to know what peak shoulder in R2 could be, unzip some R2 samples in a new folder, put the script /solutions/select_fastq_reads_by_GC_cont.py in the same folder and run it, BLAST the ouptut file (or part of it).
- 2. Not as relevant as would be for older Illumina machines, such as the still commonly used Miseq. These are Novaseq reads, sometimes R1 is slighlty worse than R2, both of them slightly worse at the end, R2 in a liltte more visible way.
- 3. Check out sample 25 with a portion of the reads a bit problematic, the rest is really high quality.
- 4. They clearly start with primer sites, so nothing strange for amplion sequencing (if unsure BLAST some)
- 5. Illumina adapterd seem absent from FastQC report, but primer are clearly there (see overrepresented sequences, or try yourself a grep)
- 6. Reads from Illumina Novaseq6000 are analysed with a 2 color channels chemistry, so all no signal reads are read as G 

**TASK 3**
forward primer search and count with regex
```bash
zgrep --color=auto -B1 -E "CCTACGGG[A|C|G|T][C|G|T|]GCA[C|G]CAG" ./data/16S_biochar_run2_10perc_sampled/*.fastq.gz
```

reverse primer search and count with regex
```bash
zgrep --color=auto -B1 -E "GACTAC[A|C|G|T][A|C|G]GGGTATCTAATCC" ./data/16S_biochar_run2_10perc_sampled/*.fastq.gz
```
if you need to count them drop -B and add -c arguments. 
Some R1 have the reverse primer and viceversa! (very few, it should be close to zero)

if you want only the reads that do not match the pattern of forward primer in R1, so you can check what is "wrong" with them
```bash
 zgrep -A 1 "^@A00618:" ./data/16S_biochar_run2_10perc_sampled/Bch-16S-V3V4-001-2_S1_L002_*R1*.gz | grep -E "^[A|C|G|T|N]" | grep -v -E "CCTACGGG[A|C|G|T][C|G|T|]GCA[C|G]CAG"
```
 Try the same for R2 (but using a regex of the reverse primer) and see how most of them have some mismatch in the primer sequences (other can be polyG polyA... maybe verify this adding a further pipe like ... | grep "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGG")

**TASK 4**
Check the 16S primer used!
Other software are able to use degenerate primer sequences to detect primers. Fastp can detect primer but we would need to create all the possible combinations from the degenerate sequence.

**TASK 5**
```bash 
cd ./results/qiime_artifacts
unzip stats.qza
nano 5940d30a-9bf8-4538-9749-ec95b2140d62/data/stats.tsv
cd -
```

**TASK 6**
Resampling for rarefaction cannot be higher than the total frequency (8000 would work here) of the highest frequency in samples. Also too low would not show well enough how all the samples behave at increased sampling effort. 

**TASK 7**
Chloroplasts and mithochondria ancestors were free living prokaryotes and they have their own 16S rDNA with is captured with universal bacterial 16S primers (also, they can be abundant in marine water samples)
