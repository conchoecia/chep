#!/usr/bin/env python
"""
This script makes several plots that are used to in calculating the heterozygosity of an animal.

To acquire the input data for this script, start with a bam file and the reference fasta file.
samtools mpileup -f assembly.fasta shot_to_assem.sorted.bam | \
"""
import argparse
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import matplotlib.ticker as ticker
from matplotlib.ticker import StrMethodFormatter, NullFormatter
from matplotlib import rc
import numpy as np
import os
import sys

# Set matplotlib style
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Helvetica']
matplotlib.rcParams['axes.facecolor'] = 'white'
matplotlib.rcParams['axes.edgecolor'] = '.95'
matplotlib.rcParams['grid.color'] = '.95'

# Preserve the vertical order of embedded images:
matplotlib.rcParams['image.composite_image'] = False
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

#rc('text', usetex=True)
#plt.rcParams['text.latex.preamble'] = [
#        r"\usepackage{siunitx}",    # micro symbols
#        r"\sisetup{detect-all}",   # ...this to force siunitx to actually use your fonts
#        r"\usepackage{tgheros}",    # helvetica font
#        r"\usepackage{sansmath}",   # math-font matching  helvetica
#        "\sansmath"
#        r"\sisetup{detect-all}",    # force siunitx to use the fonts
#        ]


def argparser():
    """parses the arguments. Returns them as a dictionary."""
    parser = argparse.ArgumentParser(description='clustergenerator')
    parser.add_argument('-f', '--filename',
                        type=str,
                        required=True,
                        help="""the filename to the histogram. Format is:
                          Col1 - read depth
                          Col2 - number of ref alleles
                          Col3 - count""")
    parser.add_argument('-g', '--genome',
                        type=str,
                        required=True,
                        help="""The genome fasta file.""")
    parser.add_argument('-o', '--outprefix',
                        type=str,
                        required=True,
                        help="""The prefix of the output files.""")
    parser.add_argument('-x', '--minX',
                        type=int,
                        help="""only plot the data starting at this X value""")
    parser.add_argument('-X', '--maxX',
                        type = int,
                        help = """only plot the window up to this X value""")
    parser.add_argument('-d', '--dark',
                        action='store_true',
                        help = "make the plot over a dark background")
    args = parser.parse_args()
    args = vars(args)
    return args

def plot_simple_figure(fname, outprefix, xmin, xmax, scale, dark=False):
    maxval = len(scale)-1
    if dark:
        plt.style.use('dark_background')

    df = pd.read_csv(fname, header=None, delim_whitespace=True)
    plt.figure(figsize=(4,4))
    panel1=plt.axes([0.14,0.12,0.8,0.8])

    #this actually produces the plot
    for i in range(0,len(df)):
        x=df[0][i]
        y=df[1][i]
        if x >= xmin and x <= xmax:
            thisfc = [1 - scale[min(df[2][i], maxval)]] *3
            if dark:
                thisfc = [scale[min(df[2][i], maxval)]] *3
            rectangle1=mplpatches.Rectangle((x,y),1,1,
                        linewidth=0,\
                        facecolor=thisfc)
            panel1.add_patch(rectangle1)
    panel1.set_xlim([xmin, xmax])
    panel1.set_ylim([0, xmax*1.1])
    panel1.set_xlabel("read depth")
    panel1.set_ylabel("# of reference bases")
    panel1.tick_params(axis='both', which='both', labelsize=6)

    # now make the color label
    #xpos, ypos, width, height
    panel2=plt.axes([0.16,0.58,0.05,0.3])
    panel2.set_ylim(1,maxval)
    #panel2.set_yscale('log', basey=2)
    panel2.set_xlim(0,1)
    panel2.set_xticklabels([])
    panel2.yaxis.set_label_position("right")
    panel2.yaxis.tick_right()
    plt.tick_params(
            axis='both',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            bottom=False,      # ticks along the bottom edge are off
            top=False,         # ticks along the top edge are off
            labelbottom=False,
            labelsize=6)
    height = 1
    #print("height is :", height)
    #print("maxcol is :", maxcol)
    for i in range(len(scale)):
        thisfc = [1 - scale[i] ] * 3
        if dark:
            thisfc = [scale[i]] * 3
        rectangle1=mplpatches.Rectangle((0,i),1,height,
                    linewidth=0,\
                    facecolor=thisfc)
        panel2.add_patch(rectangle1)

    if dark:
        plt.savefig("{}_simple_het_plot_dark.pdf".format(outprefix))
        plt.savefig("{}_simple_het_plot_dark.png".format(outprefix), dpi=300)
    else:
        plt.savefig("{}_simple_het_plot.pdf".format(outprefix))
        plt.savefig("{}_simple_het_plot.png".format(outprefix), dpi=300)

def figure_with_marginal_histogram(fname, outprefix, xmin, xmax,
                                   scale, dark = False):
    maxval = len(scale)-1
    if dark:
        plt.style.use('dark_background')

    df = pd.read_csv(fname, header=None, delim_whitespace=True)
    plt.figure(figsize=(4,4))
    panel1=plt.axes([0.14,0.12,0.8,0.53])
    panel2=plt.axes([0.14,0.7,0.8,0.2])

    for i in range(0,len(df)):
        x=df[0][i]
        y=df[1][i]
        if x >= xmin and x <= xmax:
            thisfc = [1 - scale[min(df[2][i], maxval)]] * 3
            if dark:
                thisfc = [scale[min(df[2][i], maxval)]] * 3
            rectangle1=mplpatches.Rectangle((x,y),1,1,
                        linewidth=0,\
                        facecolor=thisfc)
            panel1.add_patch(rectangle1)
    panel1.set_xlim([xmin, xmax])
    panel1.set_ylim([0, xmax*1.1])
    panel1.set_xlabel("read depth")
    panel1.set_ylabel("# of reference bases")
    panel1.tick_params(axis='both', which='both', labelsize=6)

    df.rename(columns={0: "depth", 1: "ref", 2: "count"}, inplace=True)
    df2 = df.groupby("depth")["count"].sum()

    panel2.set_xlim([xmin, xmax])
    #panel2.set_ylim(0, max(df2["depth"]))
    panel2.tick_params(
            axis='both',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            bottom=False,      # ticks along the bottom edge are off
            top=False,         # ticks along the top edge are off
            labelbottom=False,
            labelsize=6)
    thismfc = 'black'
    if dark:
        thismfc='white'
    panel2.plot(df2.index.values, df2, mfc =thismfc, mew=0,
                marker='o', linewidth=0, markersize=1)

    if dark:
        plt.savefig("{}_marginal_het_plot_dark.pdf".format(outprefix))
        plt.savefig("{}_marginal_het_plot_dark.png".format(outprefix), dpi=300)
    else:
        plt.savefig("{}_marginal_het_plot.pdf".format(outprefix))
        plt.savefig("{}_marginal_het_plot.png".format(outprefix), dpi=300)

def fig_mhist_hetero(fname, outprefix, xmin, xmax, scale, gsize, dark=False):
    maxval = len(scale)-1
    if dark:
        plt.style.use('dark_background')

    # this df contains the 2D histogram.
    #  - col 0 is the read coverage
    #  - col 1 is the number of reads in that position that had the reference base
    #  - col 2 is the number of sites with that read coverage with that many reads having the reference base
    df = pd.read_csv(fname, header=None, delim_whitespace=True)
    plt.figure(figsize=(4,4))
    panel1=plt.axes([0.14,0.12,0.8, 0.38])
    panel2=plt.axes([0.14,0.74 ,0.8,0.18])
    panel3=plt.axes([0.14,0.53,0.8, 0.18])

    #This plots the histogram below
    for i in range(0,len(df)):
        x=df[0][i]
        y=df[1][i]
        if x >= xmin and x <= xmax:
            thisfc = [1 - scale[min(df[2][i], maxval)]] *3
            if dark:
                thisfc = [scale[min(df[2][i], maxval)]] *3
            rectangle1=mplpatches.Rectangle((x,y),1,1,
                        linewidth=0,\
                        facecolor=thisfc)
            panel1.add_patch(rectangle1)
    panel1.set_xlim([xmin, xmax])
    panel1.set_ylim([0, xmax*1.1])
    panel1.set_xlabel("read depth")
    panel1.set_ylabel("# of reference bases")
    panel1.tick_params(axis='both', which='both', labelsize=6)

    # now plot the depth coverage histogram
    df.rename(columns={0: "depth", 1: "ref", 2: "count"}, inplace=True)
    df2 = df.groupby("depth")["count"].sum()
    panel2.set_xlim([xmin, xmax])
    panel2.set_ylabel("# of bases")
    panel2.tick_params(
            axis='both',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            bottom=False,      # ticks along the bottom edge are off
            top=False,         # ticks along the top edge are off
            labelbottom=False,
            labelsize=6)
    thismfc='black'
    if dark:
        thismfc='white'
    panel2.plot(df2.index.values, df2, mfc =thismfc, mew=0,
                marker='o', linestyle='dashed',
                linewidth=0, markersize=1)

    #make the line dividing the 1x from het sites
    xs = [x for x in range(xmin, int(xmax*1.1))]
    ys = [x*0.75 for x in xs]
    panel1.plot(xs, ys, color="red", ls="--", alpha=0.4, lw = 0.5)

    #make a plot of het sites
    ys = [x*0.5 for x in xs]
    panel1.plot(xs, ys, color="red", ls=":", alpha=0.4, lw = 0.5)

    #make a plot of lower cutoff
    ys = [x*0.25 for x in xs]
    panel1.plot(xs, ys, color="red", ls="--", alpha=0.4, lw = 0.5)

    het_dict = {}
    # now we calculate heterozygosity
    df2.to_csv("{}_depth_num_bases.tsv".format(outprefix), sep = '\t', index=True)
    for i in range(xmin, xmax +1):
        if i in df2.index:
            onex = int(df.query("depth == {} and ref >= {}".format(i, i*0.75))["count"].sum())
            half = int(df.query("depth == {} and ref < {} and ref >= {}".format(
                            i, i*0.75, i*0.25))["count"].sum())
            err  = int(df.query("depth == {} and ref < {}".format(
                            i, i*0.25))["count"].sum())
            het_dict[i] = {"depth": i,
                           "full_cov": onex,
                           "half_cov": half,
                           "err_cov": err,
                           "totdepth": df2.loc[i],
                           "pergenom": (df2.loc[i]/gsize)*100,
                           "het": (half/(onex+half))*100}
    for i in range(xmin, xmax +1):
        if i in het_dict:
            for j in [5,10]:
                full_cov = 0
                half_cov = 0
                err_cov = 0
                totdepth = 0
                for k in range(i-j, i+j+1):
                    if k in het_dict:
                        full_cov   += het_dict[k]["full_cov"]
                        half_cov   += het_dict[k]["half_cov"]
                        err_cov    += het_dict[k]["err_cov"]
                        totdepth   += het_dict[k]["totdepth"]
                het = (half_cov/(full_cov+half_cov))*100
                pergenom = (totdepth/gsize)*100
                het_dict[i]["{}flank_full_cov".format(j)] = full_cov
                het_dict[i]["{}flank_half_cov".format(j)] = half_cov
                het_dict[i]["{}flank_err_cov".format(j)]  = err_cov
                het_dict[i]["{}flank_totdepth".format(j)] = totdepth
                het_dict[i]["{}flank_pergenom".format(j)] = pergenom
                het_dict[i]["{}flank_het".format(j)]      = het
    panel3.set_xlim([xmin, xmax])
    panel3.set_ylabel("% Het")
    panel3.tick_params(
            axis='both',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            bottom=False,      # ticks along the bottom edge are off
            top=False,         # ticks along the top edge are off
            labelbottom=False,
            labelsize=6)

    #unwrap the heterozygosity info into some plottable info
    xs=[key for key in sorted(het_dict)]
    ys=[het_dict[key]["het"] for key in sorted(het_dict)]
    panel3.plot(xs, ys, mfc =thismfc, mew=0,
                marker='o',
                linewidth=0, markersize=1)

    cols = ["depth", "totdepth", "pergenom", "full_cov",
            "half_cov", "err_cov", "het"]
    for j in [5,10]:
        for thing in ["totdepth", "pergenom", "full_cov",
                      "half_cov", "err_cov", "het"]:
            cols.append("{}flank_{}".format(j, thing))

    hdf = pd.DataFrame.from_dict(het_dict, orient="index",
            columns=cols)
    hdf.to_csv("{}_het_by_position.tsv".format(outprefix),
               sep = '\t', index=False)

    if dark:
        plt.savefig("{}_marginal_het_het_plot_dark.pdf".format(outprefix))
        plt.savefig("{}_marginal_het_het_plot_dark.png".format(outprefix), dpi=300)
    else:
        plt.savefig("{}_marginal_het_het_plot.pdf".format(outprefix))
        plt.savefig("{}_marginal_het_het_plot.png".format(outprefix), dpi=300)

def determine_color_scheme(fname, xmin, xmax, dark=False):
    """
    This looks in the middle 50% of the plotting window,
     finds the highest count in a cell at 0.5x coverage,
     and sets that value to maxcolor. Everything from 0 to
     that value is scaled, and everything above that value
     is maxcolor.
    """
    dist = xmax-xmin
    dd=0.25*dist
    df = pd.read_csv(args["filename"], header=None, delim_whitespace=True,
                     names=["depth","ref","count"])
    df2 = df.query("depth <= {} and depth >= {}".format(xmax-dd, xmin+dd))
    df3 = df2.loc[df2["ref"] < df2["depth"]*0.55]
    df4 = df3.loc[df3["ref"] >= df3["depth"]*0.45]
    scalemax = int(max(df4["count"]))
    scale = np.linspace(0,1,scalemax+1)
    return scale

def get_genome_size(genome_file):
    """
    This gets the genome assembly size to calculate the percent of the
    genome in each size bin.
    """
    if not os.path.exists(genome_file):
        raise IOError("The genome assembly file doesn't exist")
    ending = os.path.splitext(genome_file)[-1]
    if ending not in [".fasta", ".fa", ".fna"]:
        raise IOError("The genome assembly file must end in .fa or .fasta")
    gsize = 0
    with open(genome_file, "r") as f:
        for line in f:
            line=line.strip()
            if line[0] != ">":
                gsize += len(line)
    return gsize

def main(args):
    genome_size = get_genome_size(args["genome"])
    #print("genome size: ", genome_size)
    scale = determine_color_scheme(args["filename"], args["minX"], args["maxX"], args["dark"])
    plot_simple_figure(args["filename"], args["outprefix"], args["minX"],
                       args["maxX"], scale, args["dark"])
    figure_with_marginal_histogram(args["filename"], args["outprefix"],
                                   args["minX"], args["maxX"],
                                   scale, args["dark"])
    fig_mhist_hetero(args["filename"], args["outprefix"],
                     args["minX"], args["maxX"],
                     scale, genome_size, args["dark"])

if __name__ == "__main__":
    args = argparser()
    sys.exit(main(args))
