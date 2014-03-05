
from matplotlib.font_manager import FontProperties
import matplotlib.pyplot as plt
from pylab import *
import matplotlib.mlab as mlab
import scipy.stats as spstats
from scipy.stats import norm
import numpy as np
import os
import csv


missing_data = []
full_align_lengths = []
per_bases_trim = []
avg_den_trim_bases = []
density = "0.50"

ifile  = open('sum-stats-0.50.csv', "rb")
reader = csv.reader(ifile)

for r in reader:
    missing_data.append(r[11])
    full_align_lengths.append(r[6])
    per_bases_trim.append(r[9])
    avg_den_trim_bases.append(r[10])

missing_data.pop(0)
missing_data = [float(i) for i in missing_data]
missing_data = [x * 100 for x in missing_data]

full_align_lengths.pop(0)
full_align_lengths = [float(i) for i in full_align_lengths]

per_bases_trim.pop(0)
per_bases_trim = [float(i) for i in per_bases_trim]
per_bases_trim = [x * 100 for x in per_bases_trim]

avg_den_trim_bases.pop(0)
avg_den_trim_bases = [float(i) for i in avg_den_trim_bases]
avg_den_trim_bases = [x * 100 for x in avg_den_trim_bases]

## create scatter plot of align.lenth vs % missing_data
scatter_fig = plt.figure(figsize=(19.2, 10.8))
# Create an Axes object.
ax = scatter_fig.add_subplot(1,1,1) # one row, one column, first plot
# Plot the data.
(m,b) = polyfit(missing_data, full_align_lengths,1)
yp = polyval([m,b],missing_data)
plot(missing_data, yp, alpha=0.5)
ax.scatter(missing_data, full_align_lengths, color="blue", marker="o", s=60, alpha=0.5)

# Add a title.
ax.set_title("Align Length vs Percentage Missing Data")
# Add some axis labels.
ax.set_ylabel("Align Length")
ax.set_xlabel("Percentage Missing Data")
# Produce an image.
scatter_fig.savefig("".join(["scatterplot-",density,".png"]))



# create a histogram of alignment lengths
histo_fig1 = plt.figure(figsize=(19.2, 10.8))
# best fit of data
(mu, sigma) = norm.fit(full_align_lengths)
weights = np.ones_like(full_align_lengths)/len(full_align_lengths)
n, bins, patches = plt.hist(full_align_lengths, weights=weights, normed=True, bins=len(avg_den_trim_bases))
# add a 'best fit' line
y = mlab.normpdf( bins, mu, sigma)
l = plt.plot(bins, y, 'r--', linewidth=3)
# Produce an image.
plt.grid(True)
plt.xlabel('Length (Base Pairs)')
plt.title('Alignment Length Distribution')
plt.ylabel('Frequency')
plt.xlim(0, 5000)
histo_fig1.savefig("".join(["histogram-alignment-lengths-",density,".png"]))

# create a histogram of %bases trimmed
histo_fig2 = plt.figure(figsize=(19.2, 10.8))
# best fit of data
(mu, sigma) = norm.fit(per_bases_trim)
# the histogram of the data
weights = np.ones_like(per_bases_trim)/len(per_bases_trim)
n, bins, patches = plt.hist(per_bases_trim, normed=True, bins=len(avg_den_trim_bases))
# add a 'best fit' line
y = mlab.normpdf( bins, mu, sigma)
l = plt.plot(bins, y, 'r--', linewidth=3)
# Produce an image.
plt.grid(True)
plt.title('Bases Trimmed from Each Alignment')
plt.xlabel('Bases Trimmed (Percent)')
plt.ylabel('Frequency')
plt.ylim(0, 0.5)
histo_fig2.savefig("".join(["histogram-percent-bases-trimmed-",density,".png"]))

# create a histogram of avg. density of trimmed bases
histo_fig3 = plt.figure(figsize=(19.2, 10.8))
# best fit of data
(mu, sigma) = norm.fit(avg_den_trim_bases)
# the histogram of the data
weights = np.ones_like(avg_den_trim_bases)/len(avg_den_trim_bases)
n, bins, patches = plt.hist(avg_den_trim_bases, weights=weights, normed=True, bins=len(avg_den_trim_bases))
# add a 'best fit' line
y = mlab.normpdf( bins, mu, sigma)
l = plt.plot(bins, y, 'r--', linewidth=3)

# Produce an image.
plt.grid(True)
plt.xlabel('Avg Density of Trimmed Bases (Percent)')
plt.ylabel('Frequency')
plt.title('Avg Density of Trimmed Bases')
plt.xlim(50, 100)
histo_fig3.savefig("".join(["histogram-avg-density-trimmed-bases-",density,".png"]))


