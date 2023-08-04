from __future__ import print_function
import sys, os

# Thanks to ChatGPT for fixing indentation issues

'''
This script takes in two folders containing bins from two binning methods, and their complementary 
contamination and completion scores derived from CheckM, then matches corresponding bins in the two
sets based on a minimum of 80% overlap (by length) of the bins, and finally, decides which of the two
bin versions is best in that particular bin.

Then it makes a new folder into which it puts the best version of each bin (changing the naming in
the process), and also makes a new .stats file which is consistent with the new bin folder.

Usage:
./script bin_folder_1 bin_folder_2 stats_file_1 stats_file_2 output_folder min_completion max_contamination
'''

c = float(sys.argv[6])
x = float(sys.argv[7])

# load a list of good bins (>70% complete, <10% contaminated) to save time (won't look at bad bins later on).
print("Loading list of good bins (comp>" + str(c) + "%, cont<" + str(x) + "%)")
good_bins_1 = {}
good_bins_2 = {}
for line in open(sys.argv[3]):
    if "completeness" in line:
        continue
    cut = line.strip().split('\t')
    if float(cut[1]) > c and float(cut[2]) < x:
        good_bins_1[cut[0] + '.fa'] = None
for line in open(sys.argv[4]):
    if "completeness" in line:
        continue
    cut = line.strip().split('\t')
    if float(cut[1]) > c and float(cut[2]) < x:
        good_bins_2[cut[0] + '.fa'] = None

# these are the dictionaries storing the contig names and lengths of each bin
bins_1 = {}
bins_2 = {}

print("load in the info about the contigs in each bin...")
for bin_file in good_bins_1:
    bins_1[bin_file] = {}
    contig_len = 0
    contig_name = ""
    for line in open(sys.argv[1] + '/' + bin_file):
        if not line.startswith('>'):
            contig_len += len(line.strip())
        else:
            if contig_name != "":
                bins_1[bin_file][contig_name] = contig_len
                contig_len = 0
            contig_name = line[1:-1]
    bins_1[bin_file][contig_name] = contig_len

for bin_file in good_bins_2:
    contig_len = 0
    contig_name = ""
    bins_2[bin_file] = {}
    for line in open(sys.argv[2] + '/' + bin_file):
        if not line.startswith('>'):
            contig_len += len(line.strip())
        else:
            if contig_name != "":
                bins_2[bin_file][contig_name] = contig_len
                contig_len = 0
            contig_name = line[1:-1]
    bins_2[bin_file][contig_name] = contig_len

