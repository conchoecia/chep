#!/usr/bin/env python
"""
This script makes several plots that are used to in calculating the heterozygosity of an animal.

To acquire the input data for this script, start with a bam file and the reference fasta file.
samtools mpileup -f assembly.fasta shot_to_assem.sorted.bam | \
"""
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
from matplotlib import rc
import numpy as np
import sys

rc('text', usetex=True)
plt.rcParams['text.latex.preamble'] = [
        r'\usepackage{siunitx}',    # micro symbols
        r'\sisetup{detect-all}',   # ...this to force siunitx to actually use your fonts
        r'\usepackage{tgheros}',    # helvetica font
        r'\usepackage{sansmath}',   # math-font matching  helvetica
        r'\sansmath'                # actually tell tex to use it!
        r'\sisetup{detect-all}',    # force siunitx to use the fonts
        ]


def argparser():
    """parses the arguments. Returns them as a dictionary."""
    parser = argparse.ArgumentParser(description='clustergenerator')
    parser.add_argument('-f', '--filename',
                        type=str,
                        help="""the filename to the histogram. Format is:
                          Col1 - read depth
                          Col2 - number of ref alleles
                          Col3 - count""")
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

def plot_simple_figure(fname, xmin, xmax, scale, dark=False):
    maxval = len(scale)-1
    if args["dark"]:
        plt.style.use('dark_background')

    df = pd.read_csv(args["filename"], header=None, delim_whitespace=True)
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
    panel1.set_ylabel("\# of reference bases")
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
            axis='x',          # changes apply to the x-axis
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
        plt.savefig("simple_het_plot_dark.pdf")
        plt.savefig("simple_het_plot_dark.png", dpi=300)
    else:
        plt.savefig("simple_het_plot.pdf")
        plt.savefig("simple_het_plot.png", dpi=300)

def figure_with_marginal_histogram(fname, xmin, xmax, scale, dark = False):
    maxval = len(scale)-1
    if args["dark"]:
        plt.style.use('dark_background')

    df = pd.read_csv(args["filename"], header=None, delim_whitespace=True)
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
    panel1.set_ylabel("\# of reference bases")
    panel1.tick_params(axis='both', which='both', labelsize=6)


    df.rename(columns={0: "depth", 1: "ref", 2: "count"}, inplace=True)
    df2 = df.groupby("depth")["count"].sum()

    panel2.set_xlim([xmin, xmax])
    #panel2.set_ylim(0, max(df2["depth"]))
    panel2.tick_params(
            axis='x',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            bottom=False,      # ticks along the bottom edge are off
            top=False,         # ticks along the top edge are off
            labelbottom=False,
            labelsize=True)
    thismfc = 'black'
    if dark:
        thismfc='white'
    panel2.plot(df2.index.values, df2, mfc =thismfc, mew=0,
                marker='o', linewidth=0, markersize=1)

    if dark:
        plt.savefig("marginal_het_plot_dark.pdf")
        plt.savefig("marginal_het_plot_dark.png", dpi=300)
    else:
        plt.savefig("marginal_het_plot.pdf")
        plt.savefig("marginal_het_plot.png", dpi=300)

def fig_mhist_hetero(fname, xmin, xmax, scale, dark=False):
    maxval = len(scale)-1
    if args["dark"]:
        plt.style.use('dark_background')

    df = pd.read_csv(args["filename"], header=None, delim_whitespace=True)
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
    panel1.set_ylabel("\# of reference bases")
    panel1.tick_params(axis='both', which='both', labelsize=6)

    # now plot the depth coverage histogram
    df.rename(columns={0: "depth", 1: "ref", 2: "count"}, inplace=True)
    df2 = df.groupby("depth")["count"].sum()
    panel2.set_xlim([xmin, xmax])
    panel2.set_ylabel("\# of bases")
    panel2.tick_params(
            axis='x',          # changes apply to the x-axis
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
    df2.to_csv("df2.tsv", sep = '\t', index=True)
    for i in range(xmin, xmax +1):
        onex = int(df.query("depth == {} and ref >= {}".format(i, i*0.75))["count"].sum())
        half = int(df.query("depth == {} and ref < {} and ref >= {}".format(
                        i, i*0.75, i*0.25))["count"].sum())
        err  = int(df.query("depth == {} and ref < {}".format(
                        i, i*0.25))["count"].sum())
        het_dict[i] = {"depth": i,
                       "one": onex,
                       "half": half,
                       "err": err,
                       "totdepth": df2.loc[i],
                       "het": (half/(onex+half))*100}
    panel3.set_xlim([xmin, xmax])
    panel3.set_ylabel("\% Het")
    panel3.tick_params(
            axis='x',          # changes apply to the x-axis
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
    hdf = pd.DataFrame.from_dict(het_dict, orient="index",
            columns=["depth", "totdepth", "one", "half", "err", "het"])
    hdf.to_csv("het_by_position.tsv", sep = '\t', index=False)

    if dark:
        plt.savefig("marginal_het_het_plot_dark.pdf")
        plt.savefig("marginal_het_het_plot_dark.png", dpi=300)
    else:
        plt.savefig("marginal_het_het_plot.pdf")
        plt.savefig("marginal_het_het_plot.png", dpi=300)

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

def main(args):
    scale = determine_color_scheme(args["filename"], args["minX"], args["maxX"], args["dark"])
    plot_simple_figure(args["filename"], args["minX"],
                       args["maxX"], scale, args["dark"])
    figure_with_marginal_histogram(args["filename"], args["minX"],
                                   args["maxX"], scale, args["dark"])
    fig_mhist_hetero(args["filename"], args["minX"],
                     args["maxX"], scale, args["dark"])

if __name__ == "__main__":
    args = argparser()
    sys.exit(main(args))
