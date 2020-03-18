#!/usr/bin/env python

"""
this program looks at every gene model in a transcriptome
 and outputs a bed file of the union of all the introns for each gene

Example

X = exon
- = intron
& = intron

 g1.t1 XXXXXXXX----------XXXXXXXXXXX------------XXXXXXXXX
 g1.t2    XXXXXXXXX-------XXXXXXXXX---------------XXXXXXXXX
 g1.t3  XXXX-------------XXXXXXXXXX-------------XXXXXXXXXX
   bed             &&&&&&           &&&&&&&&&&&&
                   interval1         interval2

The output is:
  - a transcripts.bed file. This is just a bed file of the transcript coords.
  - an intronic.bed file. This is the coordinates of the intronic regions.
"""

import os
import sys
import pandas as pd
import numpy as np

def parse_intervals(TG_list):
    intervals = []
    start = -1
    stop = -1
    in_intron = False
    for i in range(len(TG_list)):
        if TG_list[i] == 1:
            if not in_intron:
                # first time encountering this intronic region
                start = i
                in_intron = True
            else:
                # we're already in the intronic region
                stop = i
        elif TG_list[i] == 0:
            # we may have just stepped out of the intron
            if in_intron:
                if stop != -1:
                    intervals.append([start,stop])
                start = -1
                stop = -1
                in_intron = False
            else:
                # we don't need to do anything
                pass
    if start != -1 and stop != -1:
        # the last intergenic region buts up against the
        #  end of the transcript. This probably won't happen
        intervals.append([start,stop])
    return intervals

def intervals_to_lists(TG_chr, thisgene, TG_start, intervals, strand):
    """
    adds the important information to the df
    """
    return_these = []
    for i in range(len(intervals)):
        entry = intervals[i]
        start = entry[0]
        stop = entry[1]
        return_these.append([TG_chr,
                    TG_start + start,
                    TG_start + stop,
                    "{}.intronic{}".format(thisgene, i+1),
                    ".",
                    strand])
        #add_this = {}
        #add_this["chrom"]      = TG_chr
        #add_this["chromStart"] = TG_start + start
        #add_this["chromEnd"]   = TG_start + stop
        #add_this["name"]       = "{}.intronic{}".format(thisgene, i+1)
        #add_this["score"]      = "."
        #add_this["strand"]     = strand
        #df = df.append(add_this, ignore_index=True)
    return return_these

def gff_to_introns(gff_file):
    add_to_df = []
    with open(gff_file, "r") as f:
        thisgene = ""
        #TG = thisgene
        TG_chr = ""
        TG_start = -1
        TG_stop = -1
        TG_strand = '.'
        TG_list = []
        for line in f:
            if line.strip() and (line.strip()[0] != "#"):
                splitd = line.split()
                # we're now parsing a gff file.
                # fields
                # 0=chromosome
                # 1=source
                # 2=type
                # 3=start
                # 4=stop
                # 5=score
                # 6=strand
                # 7=dunno
                # 8=gene
                if splitd[2] == "gene":
                    #we have found a new gene. This should trigger outputting
                    # the previous gene's introns' bed first.
                    ID = splitd[8].split(";")[0].replace("ID=", "")
                    if thisgene == "":
                        #first thing in the file
                        pass
                    else:
                        # TODO output all the previous gene's info here.
                        intervals = parse_intervals(TG_list)
                        add_to_df += intervals_to_lists(TG_chr, thisgene,
                                                 TG_start, intervals, TG_strand)
                    # No matter what we did, store thisgene,
                    #  the start, the stop, the chromosome
                    TG_chr = splitd[0].strip()
                    TG_start = int(splitd[3])
                    TG_stop  = int(splitd[4])
                    TG_strand = splitd[6]
                    thisgene = ID
                    # turn off the 1s if there is exon there.
                    TG_list = [1]*(TG_stop - TG_start + 1)
                elif splitd[2] == "exon":
                    # we only care about where there are exons
                    sta = int(splitd[3]) - TG_start
                    sto  = int(splitd[4]) - TG_start + 1
                    for i in range( sta, sto):
                        # There is an exon here, so turn off the intron bit
                        TG_list[i] = 0
        intervals = parse_intervals(TG_list)
        add_to_df += intervals_to_lists(TG_chr, thisgene,
                                 TG_start, intervals, TG_strand)
    df = pd.DataFrame(add_to_df,
                      columns=['chrom', 'chromStart', 'chromEnd',
                               'name', 'score', 'strand'])
    return df

def df_to_995p(df):
    """
    this removes the introns in the 0.5th largest percentile
    percentile   length
    99.1          13396
    99.2          14470
    99.3          16174
    99.4          18512
    99.5          22109
    99.6          26908
    99.7          36173
    99.8          50502
    99.9          68541
    """
    df["length"] = df["chromEnd"] - df["chromStart"] + 1
    filter_cutoff = int(np.percentile(df["length"], 99.5))
    return df.loc[df["length"] <= filter_cutoff, ]

def print_df_to_bed(df, handle):
    """
    prints out the intronic regions for each gene in bed format
    """
    for index, row in df.iterrows():
        print("{}\t{}\t{}\t{}\t{}\t{}".format(
              row['chrom'],
              row['chromStart'],
              row['chromEnd'],
              row['name'],
              row['score'],
              row['strand']), file=handle)

def gff_to_spliced_transcripts_df(gff_file):
    """
    just makes a dataframe of the gene starts and stops
    """
    add_to_df = []
    most_recent_gene = []
    thisgene = ""
    exon_counter = 0
    with open(gff_file, "r") as f:
        for line in f:
            if line.strip() and line.strip()[0] != '#':
                splitd = line.split()
                # we're now parsing a gff file.
                # fields
                # 0=chromosome
                # 1=source
                # 2=type
                # 3=start
                # 4=stop
                # 5=score
                # 6=strand
                # 7=dunno
                # 8=gene
                if splitd[2] == "transcript":
                    #we have found a new gene. This should trigger outputting
                    # the previous gene's introns' bed first.
                    if exon_counter >= 2:
                        add_to_df.append(most_recent_gene)
                    exon_counter = 0
                    thisgene = splitd[8].split(";")[0].replace("ID=", "")
                    most_recent_gene = [splitd[0].strip(), #chrom
                                        int(splitd[3]), # start
                                        int(splitd[4]), # stop
                                        thisgene, # name
                                        '.', # score
                                        splitd[6]] #strand
                elif splitd[2] == "exon":
                    exon_counter += 1
    df = pd.DataFrame(add_to_df,
                      columns=['chrom', 'chromStart', 'chromEnd',
                               'name', 'score', 'strand'])
    return df

def transcript_95per_start_stop(df):
    """
    prints out the intronic regions for each gene in bed format
    """
    df["length"] = df["chromEnd"] - df["chromStart"] + 1
    df["chromStart_new"] = 1
    df["chromEnd_new"] = 1
    for index, row in df.iterrows():
        twopfive = int(row["length"] *0.05)
        df.loc[index, "chromStart_new"] = row["chromStart"] + twopfive - 1
        df.loc[index, "chromEnd_new"]   = row["chromEnd"] - twopfive + 1
    return df

def find_antisense_spliced_in_intron(intron_df, tx_df):
    """
    takes a df of introns and of transcripts, then finds
     the transcripts that are antisense, spliced, and in introns
    """
    gene_in_intron_pairs = []
    for index, row in tx_df.iterrows():
        chrom = intron_df.loc[intron_df["chrom"] == row["chrom"],]
        start = chrom.loc[chrom["chromStart"] <= row["chromStart_new"],]
        stop = start.loc[start["chromEnd"] >= row["chromEnd_new"],]
        if row["strand"] == "+":
            dfdir = stop.loc[stop["strand"] == '-',]
        elif row["strand"] == "-":
            dfdir = stop.loc[stop["strand"] == '+',]
        for i2, r2 in dfdir.iterrows():
            gene_in_intron_pairs.append([row["name"], r2["name"]])
            #print(gene_in_intron_pairs[-1][0], gene_in_intron_pairs[-1][1])
    return gene_in_intron_pairs

def main():
    # This part generates the intron regions
    df  = gff_to_introns(sys.argv[1])
    df2 = df_to_995p(df)
    with open("intronic.bed", "w") as f:
        print_df_to_bed(df2, f)

    # This part looks in the GFF and the introns, and finds the ones that
    # have 95% within an intron, and antisense.
    # first get the genes
    transcript_df = gff_to_spliced_transcripts_df(sys.argv[1])
    with open("transcripts.bed", "w") as f:
        print_df_to_bed(transcript_df, f)
    tx_df_ss = transcript_95per_start_stop(transcript_df)

    # now get the pairs of transcripts
    pairs = find_antisense_spliced_in_intron(df2, tx_df_ss)
    with open("antisense_spliced_in_intron_pairs.txt", "w") as f:
        for entry in pairs:
            print("{}\t{}".format(
                entry[0],
                entry[1]), file = f)

if __name__== "__main__":
    main()
