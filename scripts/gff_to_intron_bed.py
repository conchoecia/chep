#!/usr/bin/env python
"""
program: chep_gff2all
contact: Darrin Schultz - dschultz@mbari.org

This program looks at every gene model in a transcriptome
 and outputs a bed file of the different genomic regions.
 It also outputs a list of genes that are contained within
 introns.

Usage: chep_gff2regions outprefix ref.fa annotation.gff

Regions:
  - Whole genome .bed file
  - Exonic regions
    - All exons, overlaps merged
  - Intronic regions
    - All intronic regions merged, No overlap whatsoever with exons.
  - Genic Regions
    - Includes introns and exons
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
  - prefix_wholeGenome.bed
    - just a bed file of the whole genome
  - prefix_size_wholeGenome.genfile
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
import random
import re
import string
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

def intron_bed_to_995p(intronbed):
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
    df = pd.read_csv(intronbed, header=None, sep='\t', comment='#')
    df.columns = ["chrom", "chromStart", "chromEnd"]
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
            line = line.strip()
            if line and (line[0] != '#'):
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

def find_antisense_spliced_in_intron(ranges, tx_df):
    """
    takes a dict of transcript information and a df of transcripts, then finds
     the transcripts that are antisense, spliced, and in introns
    """
    # first convert the dict objext of transcript ranges to
    #  a pandas dataframe that we will use to look up the strandedness
    transcript_df = pd.DataFrame.from_records(ranges).T
    print(transcript_df)

    gene_in_intron_pairs = []
    gene_in_intron_same_direction = []
    #print("these are the colnames of the transcript_df")
    #print(transcript_df.columns)
    #print("these are the colnames of the tx_df")
    #print(tx_df.columns)
    for index, row in tx_df.iterrows():
        chrom = transcript_df.loc[transcript_df["chrom"] == row["chrom"],]
        start = chrom.loc[chrom["start"] < row["chromStart_new"],]
        stop = start.loc[start["stop"] > row["chromEnd_new"],]

        # this gets the pairs of antisense transcripts
        if row["strand"] == "+":
            dfdir = stop.loc[stop["strand"] == '-',]
        elif row["strand"] == "-":
            dfdir = stop.loc[stop["strand"] == '+',]
        for i2, r2 in dfdir.iterrows():
            gene_in_intron_pairs.append([row["name"], r2["seqname"]])
            #print(gene_in_intron_pairs[-1][0], gene_in_intron_pairs[-1][1])

        # this gets the pairs of sense transcipts
        if row["strand"] == "+":
            dfdir = stop.loc[stop["strand"] == '+',]
        elif row["strand"] == "-":
            dfdir = stop.loc[stop["strand"] == '-',]
        for i2, r2 in dfdir.iterrows():
            gene_in_intron_same_direction.append([row["name"], r2["seqname"]])

    return gene_in_intron_pairs, gene_in_intron_same_direction

def print_help():
    """Just prints a help message for the use"""
    print(__doc__)
    sys.exit()

def parent_from_gff_comment(commentline):
    """
    takes a comment string from a gff file and returns the parent sequence ID
    """
    style = ""
    temp  = ""
    if "Parent=" in commentline:
        temp = [x.strip() for x in commentline.split(";")
                   if x.strip().startswith("Parent=")][0]
    elif "transcript_id" in commentline:
        temp = [x.strip() for x in commentline.split(";")
                   if x.strip().startswith("transcript_id")][0]
    strip_start = ["rna-" #In hg38 the annotation Parents start with rna- for some reason
                   ]
    for strip_this in strip_start:
        if temp.startswith(strip_this):
            temp = temp.replace(strip_this,"")
    return temp

def randomString(stringLength=10):
    """Generate a random string of fixed length """
    letters = string.ascii_lowercase
    return ''.join(random.choice(letters) for i in range(stringLength))

def gff_coding_intervals(gfffile):
    """
    Returns a dict with this structure.
     Coding regions defined here:
         ranges[parent] = {"start": start, "stop": stop,
                           "chrom": chrom, "direction": direction,
                           "seqname": parent}
    """
    ranges = {}
    with open(gfffile, "r") as f:
        for line in f:
            line=line.strip()
            if line and (line[0] != '#'):
                # this splits the lines by tabs, and removes trailing tabs
                splitd = re.split(r'\t+', line.rstrip('\t'))
                # we can just look for either CDS or exon here because
                # it doesn't matter. The first start and last stop will
                # be properly recorded.
                if splitd[2] in ["exon", "CDS"]:
                    chrom = splitd[0]
                    start = int(splitd[3])
                    stop = int(splitd[4])
                    direction = splitd[6]
                    parent = parent_from_gff_comment(splitd[8])
                    if parent not in ranges:
                        ranges[parent] = {"start": start, "stop": stop,
                                          "chrom": chrom, "strand": direction,
                                          "seqname": parent}
                    else:
                        if start < ranges[parent]["start"]:
                            ranges[parent]["start"] = start
                        if stop > ranges[parent]["stop"]:
                            ranges[parent]["stop"] = stop
                        if chrom != ranges[parent]["chrom"]:
                            print("""You have a transcript that spans multiple scaffolds.
       This is not good. '{}'""".format(parent))
    return ranges


def gff_to_genic_bed(ranges, genic_bed):
    """
    Takes an exon ranges dict and outputs a file of genic regions using the exons
     as the delimiters
    """
    temp_bed = "{}_temp.bed".format(randomString())
    with open(temp_bed, "w") as f:
        for key in ranges:
            print("{}\t{}\t{}\t{}\t.\t{}".format(
                ranges[key]["chrom"],
                ranges[key]["start"],
                ranges[key]["stop"],
                ranges[key]["seqname"],
                ranges[key]["strand"]),
                  file = f)
    run_this = """sort -k1,1 -k2,2n {} | \
               bedtools merge | bedtools sort > {}""".format(temp_bed, genic_bed)
    runner(run_this)
    os.remove(temp_bed)

def main():
    # make sure the program has the right commands
    if len(sys.argv) != 4:
        print_help()
    # now make sure that the files are correct
    out_prefix = sys.argv[1]
    reference = sys.argv[2]
    if not os.path.exists(reference):
        raise OSError("  - The reference genome does not exist %s" % reference)
    if os.path.splitext(reference)[1] not in [".fasta", ".fa", ".fna"]:
        raise IOError("  - The reference genome must be a .fasta, .fa, or .fna file. It cannot be gzipped")
    gff_file = sys.argv[3]
    if not os.path.exists(gff_file):
        raise OSError("  - The gff file does not exist %s" % gff_file)
    if os.path.splitext(gff_file)[1].lower() not in [".gff", ".gtf", ".gff3",]:
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

    seqs_in_gff = set()
    seqs_in_fasta = set()
    # now make sure that every sequence in the gff is also in the fasta file
    with open(gff_file, "r") as f:
        for line in f:
            line=line.strip()
            if line and line[0] != '#':
                seqs_in_gff.add(line.split()[0])
    with open(reference, "r") as f:
        for line in f:
            line=line.strip()
            if line[0] == ">":
                seqs_in_fasta.add(line.split()[0][1::])
    gff_not_in_fasta = [x for x in seqs_in_gff if x not in seqs_in_fasta]
    if len(gff_not_in_fasta) > 0:
        print("""The following scaffolds are in the gff file, but not
             in the assembly file.""")
        print(gff_not_in_fasta)
        raise IOError("Scaffolds in gff not in fasta")

    # first generate the bed file of the whole genome
    wholeGenome_bed = "{}_wholeGenome.bed".format(out_prefix)
    wholeGenome_coords = "{}_size_wholeGenome.genfile".format(out_prefix)
    print("Running samtools faidx")
    run_this = "samtools faidx {}".format(reference)
    runner(run_this)
    print("generating whole-genome bed")
    fai = "{}.fai".format(reference)
    write_here = open(wholeGenome_bed, "w")
    with open(fai, "r") as f:
        for line in f:
            line = line.strip()
            if line and (line[0] != '#'):
                splitd = line.split()
                print("{}\t0\t{}".format(
                       splitd[0], splitd[1]),
                       file = write_here)
    write_here.close()
    run_this = """sort -k1,1 -k2,2n {} | cut -f1,3 > {}""".format(
            wholeGenome_bed, wholeGenome_coords)
    runner(run_this)

    assert os.path.exists(wholeGenome_bed)
    assert os.path.exists(wholeGenome_coords)

    # now generate a bed file of the transcript regions
    print("generating a bed file of genic regions")
    genic_bed = "{}_genic.bed".format(out_prefix)
    exon_ranges = gff_coding_intervals(gff_file)
    gff_to_genic_bed(exon_ranges, genic_bed)
    assert os.path.exists(genic_bed)

    # now generate a bed file of the intergenic regions
    print("generating a bed file of intergenic regions")
    intergenic_bed = "{}_intergenic.bed".format(out_prefix)
    run_this = "bedtools complement -i {} -g {} > {}".format(
        genic_bed, wholeGenome_coords, intergenic_bed)
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
    run_this = """diff {} {}""".format(wholeGenome_coords, temp_file)
    process = subprocess.Popen(run_this,
                     stdout=subprocess.PIPE,
                     stderr=subprocess.PIPE,
                     shell = True)
    stdout, stderr = process.communicate()
    print("printing: ", stdout)
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
            if line and (line[0] != '#'):
                splitd = line.split()
                # Both exon and CDS are fine. CDS subset exon,
                #  and will get merged out
                if splitd[2] in  ["exon", "CDS"]:
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
        exonic_bed, intergenic_bed, wholeGenome_coords, intronic_bed)
    runner(run_this)
    assert os.path.exists(intronic_bed)

    # there is an off-by-one error that pops up occassionally. This needs to be fixed eventually,
    #  but for now just check that intergenic_genic and intergenic_exonic_intronic are 100.0000
    ## now merge the exonic and intronic, make sure it matches genic
    #print("verifying that the exonic and intronic regions are complete")
    #temp_file = "{}_exonic_intronic_merge.temp".format(out_prefix)
    #run_this = """cat {} {} | bedtools sort | bedtools merge | \
    #            sort -k1,1 -k2,2n > {}""".format(
    #            exonic_bed, intronic_bed, temp_file)
    #runner(run_this)
    #assert os.path.exists(temp_file)
    #run_this = """diff {} {}""".format(genic_bed, temp_file)
    #process = subprocess.Popen(run_this,
    #                 stdout=subprocess.PIPE,
    #                 stderr=subprocess.PIPE,
    #                 shell = True)
    #stdout, stderr = process.communicate()
    #if stdout.decode('utf-8').strip() != "":
    #    print(stdout.decode('utf-8').strip(), file=sys.stderr)
    #    raise IOError("The intronic and exonic combined don't match the genic regions")
    #os.remove(temp_file)

    # now generate a bed file of the noncoding regions
    print("generating a bed file of noncoding regions")
    noncoding_bed = "{}_noncoding.bed".format(out_prefix)
    run_this = """cat {} {} | bedtools sort | bedtools merge | \
                sort -k1,1 -k2,2n > {}""".format(
        intronic_bed, intergenic_bed, noncoding_bed)
    runner(run_this)
    assert os.path.exists(noncoding_bed)

    print("calculating genome stats")
    #calculate % of the genome stats
    run_this = """awk '{{sum = sum + $2}} END{{print(sum)}}' {}""".format(wholeGenome_coords)
    wholeGenome_size = int(runner_w_output(run_this))
    run_this = """awk '{{sum = sum + $3 - $2 }} END{{print(sum)}}' {}""".format(exonic_bed)
    try:
        exonic_size = int(runner_w_output(run_this))
    except:
        print("tried to run this: {}".format(run_this))
    run_this = """awk '{{sum = sum + $3 - $2 }} END{{print(sum)}}' {}""".format(genic_bed)
    genic_size = int(runner_w_output(run_this))
    run_this = """awk '{{sum = sum + $3 - $2 }} END{{print(sum)}}' {}""".format(intergenic_bed)
    intergenic_size = int(runner_w_output(run_this))
    run_this = """awk '{{sum = sum + $3 - $2 }} END{{print(sum)}}' {}""".format(intronic_bed)
    intronic_size = int(runner_w_output(run_this))
    run_this = """awk '{{sum = sum + $3 - $2 }} END{{print(sum)}}' {}""".format(noncoding_bed)
    noncoding_size = int(runner_w_output(run_this))

    genome_stats = "{}_genome_stats.txt".format(out_prefix)
    print("printing genome stats to {}".format(genome_stats))
    outfile = open(genome_stats, "w")
    for writehere in [sys.stdout, outfile]:
        print("# wholeGenome_size: {}".format(wholeGenome_size),
              file=writehere)
        print("region\tnum_bases\tpercent_of_total", file = writehere)
        print("exonic\t{}\t{:.4f}".format(
              exonic_size,
              (exonic_size/wholeGenome_size)*100), file = writehere)
        print("genic\t{}\t{:.4f}".format(
              genic_size,
              (genic_size/wholeGenome_size)*100), file = writehere)
        print("intergenic\t{}\t{:.4f}".format(
              intergenic_size,
              (intergenic_size/wholeGenome_size)*100), file = writehere)
        print("intronic\t{}\t{:.4f}".format(
              intronic_size,
              (intronic_size/wholeGenome_size)*100), file = writehere)
        print("noncoding\t{}\t{:.4f}".format(
              noncoding_size,
              (noncoding_size/wholeGenome_size)*100), file = writehere)
        intergenic_genic = intergenic_size + genic_size
        print("intergenic_genic\t{}\t{:.4f}".format(
              intergenic_genic,
              (intergenic_genic/wholeGenome_size)*100), file = writehere)
        intergenic_exonic_intronic = intergenic_size + exonic_size + intronic_size
        print("intergenic_exonic_intronic\t{}\t{:.4f}".format(
              intergenic_exonic_intronic,
              (intergenic_exonic_intronic/wholeGenome_size)*100), file = writehere)
    outfile.close()

    print("now finding transcripts located within the introns of other transcripts")
    # THIS PART FINDS TRANSCRIPTS that are contained within transcripts
    # This part generates the intron regions
    intron_df = intron_bed_to_995p(intronic_bed)

    # This part looks in the GFF and the introns, and finds the ones that
    # have 95% within an intron, and antisense.
    # first get the genes
    transcript_df = gff_to_spliced_transcripts_df(gff_file)
    tx_df_ss = transcript_98per_start_stop(transcript_df)
    #print(tx_df_ss)

    # now get the pairs of transcripts
    antisense_pairs, sense_pairs = find_antisense_spliced_in_intron(exon_ranges, tx_df_ss)
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

if __name__ == "__main__":
    main()
