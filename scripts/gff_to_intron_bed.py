#!/usr/bin/env python

"""
This program looks at every gene model in a transcriptome
 and outputs a bed file of the different genomic regions.
 It also outputs a list of genes that are contained within
 introns.

Usage:
  chep_gff2regions outprefix ref gff

Regions:
  - Whole genome .bed file
  - Genic Regions
    - Includes introns and exons
  - Exonic regions
    - All exons, overlaps merged
  - Intronic regions
    - All intronic regions merged, No overlap whatsoever with exons.
  - Intergenic regions
    - All regions that do not contain genic regions
  - Noncoding
    - Introns and intergenic

Explanation of intronic regions

X = exon
- = intron
& = intron

 g1.t1 XXXXXXXX----------XXXXXXXXXXX------------XXXXXXXXX
 g1.t2    XXXXXXXXX-------XXXXXXXXX---------------XXXXXXXXX
 g1.t3  XXXX-------------XXXXXXXXXX-------------XXXXXXXXXX
   bed             &&&&&&           &&&&&&&&&&&&
                   interval1         interval2

Parameters:
  chep_gff2intron out_prefix annotation.gff

Required Software:
  - samtools
  - bedtools
  - awk

The output is:
  # files related to getting bed files of different genomic regions
  - prefix_whole_genome.bed
    - just a bed file of the whole genome
  - prefix_whole_genome.genfile
    - the bedtools genFile format. Like a bed file, but no start.
  - prefix_exonic.bed
    - The exonic regions
  - prefix_genic.bed
    - The transcribed regions
  - prefix_intergenic.bed
    - The regions not transcribed
  - prefix_intronic.bed
    - The intronic regions
  - prefix_noncoding.bed
    - The intronic and intergenic regions
  - prefix_genome_stats.txt
    - Stats of the the percent of the genome in each category
  # files related to finding transcripts contained in introns
  - prefix_custom_introns.bed
    - introns with special cutoffs to remove the longest of introns
  - prefix_sense_spliced_in_intron_pairs.txt
    - genes that are contained in introns that are the same sense as the intron
    - 98% of the gene must be in a single intron
  - prefix_antisense_spliced_in_intron_pairs.txt
    - genes that are contained in introns that are antisense to the intron
    - 98% of the gene must be in a single intron
"""

import os
import sys
import numpy as np
import pandas as pd
import subprocess

def runner(run_this):
    """Just runs a process on the system"""
    process = subprocess.Popen(run_this,
                     stdout=subprocess.PIPE,
                     stderr=subprocess.PIPE,
                     shell = True)
    stdout, stderr = process.communicate()
    for print_this in [stdout, stderr]:
        print_this = print_this.decode('utf-8').strip()
        if print_this:
            print(print_this)

def runner_w_output(run_this):
    """Runs a process on the system and prints the output"""
    process = subprocess.Popen(run_this,
                     stdout=subprocess.PIPE,
                     stderr=subprocess.PIPE,
                     shell = True)
    stdout, stderr = process.communicate()
    return stdout.decode('utf-8').strip()

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

    functions to remove false positives, or genes that are obviously too long.
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

def transcript_98per_start_stop(df):
    """
    prints out the intronic regions for each gene in bed format
    """
    df["length"] = df["chromEnd"] - df["chromStart"] + 1
    df["chromStart_new"] = 1
    df["chromEnd_new"] = 1
    for index, row in df.iterrows():
        twopfive = int(row["length"] *0.02)
        df.loc[index, "chromStart_new"] = row["chromStart"] + twopfive - 1
        df.loc[index, "chromEnd_new"]   = row["chromEnd"] - twopfive + 1
    return df

def find_antisense_spliced_in_intron(intron_df, tx_df):
    """
    takes a df of introns and of transcripts, then finds
     the transcripts that are antisense, spliced, and in introns
    """
    gene_in_intron_pairs = []
    gene_in_intron_same_direction = []
    for index, row in tx_df.iterrows():
        chrom = intron_df.loc[intron_df["chrom"] == row["chrom"],]
        start = chrom.loc[chrom["chromStart"] <= row["chromStart_new"],]
        stop = start.loc[start["chromEnd"] >= row["chromEnd_new"],]

        # this gets the pairs of antisense transcripts
        if row["strand"] == "+":
            dfdir = stop.loc[stop["strand"] == '-',]
        elif row["strand"] == "-":
            dfdir = stop.loc[stop["strand"] == '+',]
        for i2, r2 in dfdir.iterrows():
            gene_in_intron_pairs.append([row["name"], r2["name"]])
            #print(gene_in_intron_pairs[-1][0], gene_in_intron_pairs[-1][1])

        # this gets the pairs of sense transcipts
        if row["strand"] == "+":
            dfdir = stop.loc[stop["strand"] == '+',]
        elif row["strand"] == "-":
            dfdir = stop.loc[stop["strand"] == '-',]
        for i2, r2 in dfdir.iterrows():
            gene_in_intron_same_direction.append([row["name"], r2["name"]])

    return gene_in_intron_pairs, gene_in_intron_same_direction

def main():
    # make sure the program has the right commands
    if len(sys.argv) != 4:
        raise IOError("Run the program like: chep_gff2regions outprefix ref gff")
    # now make sure that the files are correct
    out_prefix = sys.argv[1]
    reference = sys.argv[2]
    if not os.path.exists(reference):
        raise OSError("  - The reference genome does not exist %s" % reference)
    if os.path.splitext(reference)[1] not in [".fasta", ".fa"]:
        raise IOError("  - The reference genome must be a .fasta or .fa file. It cannot be gzipped")
    gff_file = sys.argv[3]
    if not os.path.exists(gff_file):
        raise OSError("  - The gff file does not exist %s" % gff_file)
    if os.path.splitext(gff_file)[1] not in [".gff"]:
        raise IOError("  - The gff must have the ending '.gff'. It cannot be gzipped.")
    # make sure that all the files in the path exist
    split_paths = out_prefix.split('/')[:-1]
    for i in range(len(split_paths)):
        thisdir = "/".join(split_paths[0:i+1])
        print("trying to make directory {}".format(thisdir))
        if not os.path.exists(thisdir):
            try:
                os.mkdir(thisdir)
            except OSError:
                print ("  - Creation of the directory %s failed" % thisdir)
            else:
                print ("  - Successfully created the directory %s " % thisdir)

    # first generate the bed file of the whole genome
    whole_genome_bed = "{}_whole_genome.bed".format(out_prefix)
    whole_genome_coords = "{}_whole_genome.genfile".format(out_prefix)
    print("Running samtools faidx")
    run_this = "samtools faidx {}".format(reference)
    runner(run_this)
    print("generating whole-genome bed")
    fai = "{}.fai".format(reference)
    write_here = open(whole_genome_bed, "w")
    with open(fai, "r") as f:
        for line in f:
            line = line.strip()
            if line:
                splitd = line.split()
                print("{}\t0\t{}".format(
                       splitd[0], splitd[1]),
                       file = write_here)
    write_here.close()
    run_this = """sort -k1,1 -k2,2n {} | cut -f1,3 > {}""".format(
            whole_genome_bed, whole_genome_coords)
    runner(run_this)

    assert os.path.exists(whole_genome_bed)
    assert os.path.exists(whole_genome_coords)

    # now generate a bed file of the transcript regions
    print("generating a bed file of genic regions")
    genic_bed = "{}_genic.bed".format(out_prefix)
    temp_bed = "{}_temp.bed".format(out_prefix)
    write_here = open(temp_bed, "w")
    with open(gff_file, "r") as f:
        for line in f:
            line = line.strip()
            if line:
                splitd = line.split()
                if splitd[2] == "transcript":
                    print("{}\t{}\t{}".format(
                         splitd[0], splitd[3], splitd[4]),
                         file = write_here)
    write_here.close()
    run_this = """sort -k1,1 -k2,2n {} | \
               bedtools merge | bedtools sort > {}""".format(temp_bed, genic_bed)
    runner(run_this)
    os.remove(temp_bed)
    assert os.path.exists(genic_bed)

    # now generate a bed file of the intergenic regions
    print("generating a bed file of intergenic regions")
    intergenic_bed = "{}_intergenic.bed".format(out_prefix)
    run_this = "bedtools complement -i {} -g {} > {}".format(
        genic_bed, whole_genome_coords, intergenic_bed)
    runner(run_this)
    assert os.path.exists(intergenic_bed)

    # now merge the intergenic and genic, make sure it matches the genFile exactly
    print("verifying that the genic and intergenic regions are complete")
    temp_file = "{}_intergenic_genic_merge.temp".format(out_prefix)
    run_this = """cat {} {} | bedtools sort | bedtools merge | \
                sort -k1,1 -k2,2n | cut -f1,3 > {}""".format(
                genic_bed, intergenic_bed, temp_file)
    runner(run_this)
    assert os.path.exists(temp_file)
    run_this = """diff {} {}""".format(whole_genome_coords, temp_file)
    process = subprocess.Popen(run_this,
                     stdout=subprocess.PIPE,
                     stderr=subprocess.PIPE,
                     shell = True)
    stdout, stderr = process.communicate()
    assert stdout.decode('utf-8').strip() == ""
    os.remove(temp_file)

    # now generate a bed file of the exonic regions
    print("generating a bed file of exonic regions")
    exonic_bed = "{}_exonic.bed".format(out_prefix)
    temp_bed = "{}_temp.bed".format(out_prefix)
    write_here = open(temp_bed, "w")
    with open(gff_file, "r") as f:
        for line in f:
            line = line.strip()
            if line:
                splitd = line.split()
                if splitd[2] == "exon":
                    print("{}\t{}\t{}".format(
                         splitd[0], splitd[3], splitd[4]),
                         file = write_here)
    write_here.close()
    run_this = """sort -k1,1 -k2,2n {} | \
               bedtools merge | bedtools sort > {}""".format(temp_bed, exonic_bed)
    runner(run_this)
    os.remove(temp_bed)
    assert os.path.exists(exonic_bed)

    # now generate a bed file of the intronic regions
    print("generating a bed file of intronic regions")
    intronic_bed = "{}_intronic.bed".format(out_prefix)
    run_this = """cat {} {} | bedtools sort | bedtools merge | \
                sort -k1,1 -k2,2n | bedtools complement -i - -g {} > {}""".format(
        exonic_bed, intergenic_bed, whole_genome_coords, intronic_bed)
    runner(run_this)
    assert os.path.exists(intronic_bed)

    # now merge the exonic and intronic, make sure it matches genic
    print("verifying that the exonic and intronic regions are complete")
    temp_file = "{}_exonic_intronic_merge.temp".format(out_prefix)
    run_this = """cat {} {} | bedtools sort | bedtools merge | \
                sort -k1,1 -k2,2n > {}""".format(
                exonic_bed, intronic_bed, temp_file)
    runner(run_this)
    assert os.path.exists(temp_file)
    run_this = """diff {} {}""".format(genic_bed, temp_file)
    process = subprocess.Popen(run_this,
                     stdout=subprocess.PIPE,
                     stderr=subprocess.PIPE,
                     shell = True)
    stdout, stderr = process.communicate()
    assert stdout.decode('utf-8').strip() == ""
    os.remove(temp_file)

    # now generate a bed file of the noncoding regions
    print("generating a bed file of noncoding regions")
    noncoding_bed = "{}_noncoding.bed".format(out_prefix)
    run_this = """cat {} {} | bedtools sort | bedtools merge | \
                sort -k1,1 -k2,2n > {}""".format(
        intronic_bed, intergenic_bed, noncoding_bed)
    runner(run_this)
    assert os.path.exists(noncoding_bed)

    #calculate % of the genome stats
    run_this = """awk '{{sum = sum + $2}} END{{print(sum)}}' {}""".format(whole_genome_coords)
    whole_genome_size = int(runner_w_output(run_this))
    run_this = """awk '{{sum = sum + $3 - $2 }} END{{print(sum)}}' {}""".format(exonic_bed)
    exonic_size = int(runner_w_output(run_this))
    run_this = """awk '{{sum = sum + $3 - $2 }} END{{print(sum)}}' {}""".format(genic_bed)
    genic_size = int(runner_w_output(run_this))
    run_this = """awk '{{sum = sum + $3 - $2 }} END{{print(sum)}}' {}""".format(intergenic_bed)
    intergenic_size = int(runner_w_output(run_this))
    run_this = """awk '{{sum = sum + $3 - $2 }} END{{print(sum)}}' {}""".format(intronic_bed)
    intronic_size = int(runner_w_output(run_this))
    run_this = """awk '{{sum = sum + $3 - $2 }} END{{print(sum)}}' {}""".format(noncoding_bed)
    noncoding_size = int(runner_w_output(run_this))

    genome_stats = "{}_genome_stats.txt".format(out_prefix)
    outfile = open(genome_stats, "w")
    for writehere in [sys.stdout, outfile]:
        print("# whole_genome_size: {}".format(whole_genome_size),
              file=writehere)
        print("region\tnum_bases\tpercent_of_total", file = writehere)
        print("exonic\t{}\t{:.4f}".format(
              exonic_size,
              (exonic_size/whole_genome_size)*100), file = writehere)
        print("genic\t{}\t{:.4f}".format(
              genic_size,
              (genic_size/whole_genome_size)*100), file = writehere)
        print("intergenic\t{}\t{:.4f}".format(
              intergenic_size,
              (intergenic_size/whole_genome_size)*100), file = writehere)
        print("intronic\t{}\t{:.4f}".format(
              intronic_size,
              (intronic_size/whole_genome_size)*100), file = writehere)
        print("noncoding\t{}\t{:.4f}".format(
              noncoding_size,
              (noncoding_size/whole_genome_size)*100), file = writehere)
        intergenic_genic = intergenic_size + genic_size
        print("intergenic_genic\t{}\t{:.4f}".format(
              intergenic_genic,
              (intergenic_genic/whole_genome_size)*100), file = writehere)
        intergenic_exonic_intronic = intergenic_size + exonic_size + intronic_size
        print("intergenic_exonic_intronic\t{}\t{:.4f}".format(
              intergenic_exonic_intronic,
              (intergenic_exonic_intronic/whole_genome_size)*100), file = writehere)

    outfile.close()

    print("now finding transcripts located within the introns of other transcripts")
    # THIS PART FINDS TRANSCRIPTS that are contained within transcripts
    # This part generates the intron regions
    df  = gff_to_introns(gff_file)
    df2 = df_to_995p(df)
    with open("{}_custom_introns.bed".format(out_prefix), "w") as f:
        print_df_to_bed(df2, f)

    # This part looks in the GFF and the introns, and finds the ones that
    # have 95% within an intron, and antisense.
    # first get the genes
    transcript_df = gff_to_spliced_transcripts_df(gff_file)
    tx_df_ss = transcript_98per_start_stop(transcript_df)

    # now get the pairs of transcripts
    antisense_pairs, sense_pairs = find_antisense_spliced_in_intron(df2, tx_df_ss)
    with open("{}_antisense_spliced_in_intron_pairs.txt".format(out_prefix), "w") as f:
        for entry in antisense_pairs:
            print("{}\t{}".format(
                entry[0],
                entry[1]), file = f)

    with open("{}_sense_spliced_in_intron_pairs.txt".format(out_prefix), "w") as f:
        for entry in sense_pairs:
            print("{}\t{}".format(
                entry[0],
                entry[1]), file = f)

if __name__== "__main__":
    main()
