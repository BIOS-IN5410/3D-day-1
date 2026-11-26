# Calling A/B compartments from Hi-C data
Now, you will reuse the Hi-C contact maps you made with HiC-Pro to find the two large-scale compartments of the genome:
- **A compartment**: open, gene-rich and mostly active chromatin
- **B compartment**: closed, gene-poor and mostly inactive chromatin

Regions in the same compartment interact with each other more than expected, even when they are far apart on the chromosome. This gives the Hi-C map a "checkerboard" (plaid) pattern. To find the compartments you will:
- Remove the effect of genomic distance from the contact map (observed/expected)
- Correlate the contact profiles of all pairs of bins
- Use the first principal component (PC1) of the correlation matrix to split the chromosome in two
- Decide which part is A and which is B, using the genes in the UCSC genome browser

**You need to have finished the previous TAD-exercise first**, because we use the `hicpro_results` folder you made there.

As before, stop at the orange lines
```diff
! Like this
```
and note down your answers.

**1. Log in to Fox (unless you are already logged in), go to your 3D-day-1 folder and allocate resources**
Log in to Fox (if you are not logged in already):
```bash
ssh [ec-username]@fox.educloud.no
```
Then navigate to the exercise folder and allocate resources:
```bash
cd ~/3D-day-1
salloc --ntasks=1 --mem-per-cpu=4G --time=00:30:00 --account=ec34
```

**2. Find the 1 Mb contact map**

In 3D-day-1 you asked HiC-Pro for two bin sizes (`BIN_SIZE = 50000 1000000`). So far we have only used the 50 kb map. Now we will use the 1 Mb map:
```bash
ls hicpro_results/hic_results/matrix/chr18/iced/1000000/
ls hicpro_results/hic_results/matrix/chr18/raw/1000000/
head hicpro_results/hic_results/matrix/chr18/iced/1000000/chr18_1000000_iced.matrix
grep -m 3 chr18 hicpro_results/hic_results/matrix/chr18/raw/1000000/chr18_1000000_abs.bed
```
```diff
! What do the three columns in the .matrix file contain?
! The first chr18 bin does not have ID 1. Why not? (hint: head chrom_hg19.sizes)
! TADs are typically a few hundred kb. Why do you think we use 1 Mb bins for compartments?
```

**3. Look at the compartment.py script**

The script has seven numbered steps, written as comments (lines starting with `#`). Press `q` to leave `less`.
```diff
! In step 4, each diagonal of the matrix is divided by its own mean. 
! What does a diagonal of a Hi-C matrix represent? Why do we need to remove this effect before step 5?
! In step 5, what does a high correlation between two bins mean?
```

**4. Run the script**

The script needs Python with `numpy` and `matplotlib`. We use the same container as the "bonus" in the previous exercise:
```bash
apptainer exec /projects/ec34/biosin5410/HiC/python-healpy.sif python compartments.py \
  hicpro_results/hic_results/matrix/chr18/iced/1000000/chr18_1000000_iced.matrix \
  hicpro_results/hic_results/matrix/chr18/raw/1000000/chr18_1000000_abs.bed \
  chr18 hic/chr18_compartments
```
This makes two files: `hic/chr18_compartments.bedGraph` (the PC1 value for each 1 Mb bin) and `hic/chr18_compartments.png` (a figure).
```bash
head hic/chr18_compartments.bedGraph
```
```diff
! chr18 is about 78 Mb long, which gives 79 bins of 1 Mb. How many bins are in the bedGraph? 
! Which bins are missing, and why? (hint: look at the start coordinates, and think about what lies in the middle of a chromosome)
```

**5. Look at the figure**

Download `hic/chr18_compartments.png` to your own computer and open it. The top panel is the correlation matrix (red = positive, blue = negative correlation). The bottom panel is PC1.
```diff
! Can you see the checkerboard pattern in the correlation matrix?
! Compare the pattern with PC1. What do bins with the same PC1 sign have in common in the matrix?
```

**6. Decide which compartment is A in the UCSC genome browser**

The sign of PC1 is arbitrary: the maths can return the same answer multiplied by -1. So positive does not mean A. You must decide this yourself by comparing with the genes.

- Download `hic/chr18_compartments.bedGraph` to your own computer
- Go to http://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=hg19 and upload the file as a custom track, the same way you uploaded the TADs in the previous exercise (also upload `chr18.consensus.bed` again if it is gone)
- View the whole chromosome by typing `chr18` in the position box
- Make sure a gene track (e.g. NCBI RefSeq or GENCODE) is shown

Positive PC1 values are drawn as red bars and negative values as blue bars.
```diff
! Which colour (sign) lines up with the gene-dense regions? This is the A compartment.
! Roughly what fraction of chr18 is in the A compartment? chr18 is known as a gene-poor chromosome. Does that fit?
```

**7. Compare compartments and TADs**

Zoom in to a region where PC1 changes sign, so you can see both the compartment track and your TAD track.
```diff
! Are compartments larger or smaller than TADs?
! Do compartment switches (where PC1 changes sign) fall at TAD borders?
! Is every TAD border also a compartment switch? What does this tell you about how TADs and compartments relate?
```
