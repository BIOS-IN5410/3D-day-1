# Call A/B compartments from a HiC-Pro (ICE-normalized) matrix
# Usage: python compartments.py <matrix> <abs.bed> <chrom> <output prefix>
# Writes <prefix>.bedGraph (PC1 per bin) and <prefix>.png (correlation matrix + PC1)
from __future__ import print_function
import sys
import numpy as np

matrix_file, bed_file, chrom, prefix = sys.argv[1:5]

# 1. Read the bins of this chromosome (HiC-Pro bin IDs are genome-wide)
bins = [l.split() for l in open(bed_file)]
bins = [(b[0], int(b[1]), int(b[2]), int(b[3])) for b in bins if b[0] == chrom]
first_id = bins[0][3]
n = len(bins)

# 2. Build a dense, symmetric contact matrix
M = np.zeros((n, n))
for line in open(matrix_file):
    i, j, v = line.split()
    i, j = int(i) - first_id, int(j) - first_id
    if 0 <= i < n and 0 <= j < n:
        M[i, j] = M[j, i] = float(v)

# 3. Remove empty bins (e.g. centromere, unmappable regions)
keep = M.sum(axis=0) > 0
M = M[np.ix_(keep, keep)]
m = M.shape[0]

# 4. Observed/expected: divide each diagonal by its mean (removes distance decay)
OE = np.zeros((m, m))
for d in range(m):
    i = np.arange(m - d)
    diag = M[i, i + d]
    if diag.mean() > 0:
        OE[i, i + d] = OE[i + d, i] = diag / diag.mean()

# 5. Pearson correlation matrix, then its first eigenvector (PC1)
C = np.nan_to_num(np.corrcoef(OE))
eigenvalues, eigenvectors = np.linalg.eigh(C)
pc1 = eigenvectors[:, -1]          # eigh sorts eigenvalues ascending: last = largest

# 6. Write PC1 as a bedGraph (positive and negative values = the two compartments)
kept_bins = [b for b, k in zip(bins, keep) if k]
with open(prefix + ".bedGraph", "w") as out:
    print('track type=bedGraph name="PC1 compartments" visibility=full '
          'graphType=bar color=200,0,0 altColor=0,0,200', file=out)
    for (c, start, end, _), value in zip(kept_bins, pc1):
        print("%s\t%i\t%i\t%.4f" % (c, start, end, value), file=out)
print("Wrote %i bins to %s.bedGraph" % (len(kept_bins), prefix))

# 7. Plot the correlation matrix with PC1 underneath
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 7.5), sharex=True,
                               gridspec_kw={"height_ratios": [4, 1]})
ax1.imshow(C, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
ax1.set_title("%s: correlation of observed/expected" % chrom)
x = np.arange(m)
ax2.bar(x, pc1, width=1, color=["#c80000" if v > 0 else "#0000c8" for v in pc1])
ax2.axhline(0, color="black", linewidth=0.5)
ax2.set_xlim(-0.5, m - 0.5)
ax2.set_ylabel("PC1")
ax2.set_xlabel("bin (empty bins removed)")
plt.tight_layout()
plt.savefig(prefix + ".png", dpi=120)
print("Wrote %s.png" % prefix)
