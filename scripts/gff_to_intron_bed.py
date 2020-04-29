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
"""
import gzip
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

def print_help():
    """Just prints a help message for the use"""
    print(__doc__)
    sys.exit()

def df_to_999p(df):
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
    df["length"] = df["stop"] - df["start"] + 1
    filter_cutoff = int(np.percentile(df["length"], 99.9))
    return df.loc[df["length"] <= filter_cutoff, ]

def transcript_95per_start_stop(df):
    """
    This trims off 2.5% of the start and stop
    """
    df["length"] = df["stop"] - df["start"] + 1
    df["start_new"] = df["start"] + (df["length"] * 0.025) - 1
    df["stop_new"] = df["stop"]   - (df["length"] * 0.025) + 1
    return df

def find_gene_in_intron(intron_df, tx_df):
    """
    takes a df of introns and of transcripts, then finds
     the transcripts that are in introns
    """
    gene_in_intron_entries = []
    # the goal of this is to find introns that encapsulate the transcript
    for index, row in tx_df.iterrows():
        intron_sub = intron_df.loc[intron_df["chrom"] == row["chrom"],]
        intron_sub = intron_sub.loc[row["start_new"] >= intron_sub["start"],]
        intron_sub = intron_sub.loc[row["stop_new"] <= intron_sub["stop"],]
        # if we weren't able to filter down to one
        if intron_sub.shape[0] >= 1 and isinstance(intron_sub, pd.DataFrame):
            for i2, r2 in intron_sub.iterrows():
                if row["strand"] != r2["strand"]:
                    sense_anti = "antisense"
                elif row["strand"] == r2["strand"]:
                    sense_anti = "sense"
                gene_in_intron_entries.append(
                    {"chrom" : row["chrom"],
                     "gene" : row["seqname"],
                     "gene_start" : row["start"],
                     "gene_stop"  : row["stop"],
                     "surrounding_gene" : r2["seqname"],
                     "surrounding_start" : r2["start"],
                     "surrounding_stop"  : r2["stop"],
                     "sense_anti" : sense_anti,
                     "num_exons": row["num_exons"]})

    df_gene_in_introns = pd.DataFrame.from_records(gene_in_intron_entries)
    return df_gene_in_introns

def parent_from_gff_comment(commentline):
    """
    takes a comment string from a gff file and returns the parent sequence ID
    """
    style = ""
    temp  = ""
    if "Parent=" in commentline:
        temp = [x.strip() for x in commentline.split(";")
                   if x.strip().startswith("Parent=")][0].replace("Parent=","")
    elif "transcript_id" in commentline:
        temp = [x.strip() for x in commentline.split(";")
                   if x.strip().startswith("transcript_id")][0]
    # this is the penultimate last-ditch effort to try to get the
    elif "ID=" in commentline:
        temp = [x.strip() for x in commentline.split(";")
                   if x.strip().startswith("ID=")][0]
    # More conditions might need to be added later on.
    strip_start = ["rna-", #In hg38 the annotation Parents start with rna-
                   "braker1_"]
    for strip_this in strip_start:
        if temp.startswith(strip_this):
            temp = temp.replace(strip_this,"")
    if temp == "":
        raise IOError("This gene's parent was returned as \"\". : {}".format(commentline))
    return temp

def randomString(stringLength=10):
    """Generate a random string of fixed length """
    letters = string.ascii_lowercase
    return ''.join(random.choice(letters) for i in range(stringLength))

def num_exon_bases_in_genelist(gfffile, gene_list):
    """
    takes a gff file, and a list of genes.
    This returns the number of bases in exons that those genes occupy.
    """
    gzipped = gfffile.endswith(".gz")
    # these are things that span multiple scaffolds and should be ignored

    if len(gene_list) == 0:
        return 0
    else:
        random_seed = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(25))
        temp_bed = "{}.temp.bed".format(random_seed)
        temp_bed_handle = open(temp_bed, "w")
        if gzipped:
            f = gzip.open(gfffile, "rt")
        else: f = open(gfffile, "r")
        for line in f:
            line=line.strip()
            if line and (line[0] != '#'):
                # this splits the lines by tabs, and removes trailing tabs
                splitd = re.split(r'\t+', line.rstrip('\t'))
                # we can just look for either CDS or exon here because
                # it doesn't matter. The first start and last stop will
                # be properly recorded.
                if splitd[2] in ["exon", "CDS"]:
                    skip_this = False
                    if "tRNAscan" in splitd[8]:
                        skip_this = True
                    if not skip_this:
                        chrom = splitd[0]
                        start = int(splitd[3])
                        stop = int(splitd[4])
                        parent = parent_from_gff_comment(splitd[8])
                        if parent in gene_list:
                            print("{}\t{}\t{}".format(chrom, start, stop), file = temp_bed_handle)
        f.close()
        temp_bed_handle.close()

        # now merge all the exons from this genelist
        temp_bed2 = "{}.temp.2.bed".format(random_seed)
        run_this = """sort -k1,1 -k2,2n {} | \
               bedtools merge | bedtools sort | \
               sort -k1,1 -k2,2n > {}""".format(temp_bed, temp_bed2)
        runner(run_this)

        # now sum up all the lengths
        exon_df = pd.read_csv(temp_bed2, header = 0, sep = "\t")
        exon_df.columns = ["chrom", "start", "stop"]
        exon_bp = int(np.sum(exon_df["stop"] - exon_df["start"] + 1))
        os.remove(temp_bed)
        os.remove(temp_bed2)
        return exon_bp

def gff_coding_intervals(gfffile):
    """
    Returns a dict with this structure.
     Coding regions defined here:
         ranges[parent] = {"start": start, "stop": stop,
                           "chrom": chrom, "direction": direction,
                           "seqname": parent}
    """
    ranges = {}
    exons = {}
    gzipped = gfffile.endswith(".gz")
    # these are things that span multiple scaffolds and should be ignored
    purge_set = set()
    if gzipped:
        f = gzip.open(gfffile, "rt")
    else: f = open(gfffile, "r")
    for line in f:
        line=line.strip()
        if line and (line[0] != '#'):
            # this splits the lines by tabs, and removes trailing tabs
            splitd = re.split(r'\t+', line.rstrip('\t'))
            # we can just look for either CDS or exon here because
            # it doesn't matter. The first start and last stop will
            # be properly recorded.
            if splitd[2] in ["exon", "CDS"]:
                skip_this = False
                if "tRNAscan" in splitd[8]:
                    skip_this = True
                if not skip_this:
                    chrom = splitd[0]
                    start = int(splitd[3])
                    stop = int(splitd[4])
                    direction = splitd[6]
                    parent = parent_from_gff_comment(splitd[8])
                    if parent not in ranges:
                        ranges[parent] = {"start": start, "stop": stop,
                                          "chrom": chrom, "strand": direction,
                                          "seqname": parent, "num_exons": 1}
                        exons[parent] = [ [start, stop] ]
                    else:
                        exons[parent].append([start, stop])
                        if start < ranges[parent]["start"]:
                            ranges[parent]["start"] = start
                            ranges[parent]["num_exons"]+=1
                        if stop > ranges[parent]["stop"]:
                            ranges[parent]["stop"] = stop
                            ranges[parent]["num_exons"]+=1
                        if chrom != ranges[parent]["chrom"]:
                            purge_set.add(parent)
                            print("""You have a transcript that spans multiple scaffolds.
       This is not good. '{}'""".format(parent))
    f.close()

    for get_rid_of_this in purge_set:
        del ranges[get_rid_of_this]
        del exons[get_rid_of_this]

    # for each parent, use the exons to find introns
    intron_list = []
    for parent_name in ranges:
        # make sure it isn't single-exon
        if len(exons[parent_name]) > 1:
            parent_start = ranges[parent_name]["start"]
            parent_stop  = ranges[parent_name]["stop"]
            exon_ranges = {x:0 for x in range(parent_start, parent_stop + 1)}
            for exon_start, exon_stop in list(exons[parent_name]):
                for e in range(exon_start, exon_stop+1):
                    exon_ranges[e] = 1

            in_intron = False
            intron_start = -1
            for x in range(parent_start, parent_stop + 1):
                if exon_ranges[x] == 0:
                    # if exon_ranges == 0, it must be an intronic region
                    if in_intron:
                        # if we're already in an intron, don't do anything.
                        pass
                    else:
                        # we just steppd into an intron
                        intron_start = x
                        in_intron = True
                elif exon_ranges[x] == 1:
                    # this is an exon
                    if in_intron:
                        # we have just stepped out of an intron
                        intron_list.append({"start":   intron_start,
                                            "stop":    x-1,
                                            "chrom":   ranges[parent_name]["chrom"],
                                            "strand":  ranges[parent_name]["strand"],
                                            "seqname": parent_name})
                        in_intron = False
                    else:
                        # we were already in an exon. don't do anything
                        pass
            #just stepped out of the range. if we're still in intron, raise error
            # transcript ranges are capped by exons
            if in_intron:
                raise IOError("This transcript is capped by an intron. {}".format(parent_name))
    introns_df = pd.DataFrame.from_records(intron_list)
    return ranges, introns_df

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
               bedtools merge | bedtools sort | \
               sort -k1,1 -k2,2n > {}""".format(temp_bed, genic_bed)
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
    reference_ok = False
    for chance in [".fasta", ".fa", ".fna",
                   ".fasta.gz", ".fa.gz", ".fna.gz"]:
        if chance in reference.split("/")[-1]:
            reference_ok = True
    if not reference_ok:
        raise IOError("  - The reference genome must be a .fasta, .fa, .fna file, or a gzipped variant.")
    ref_gzipped = False
    if reference.endswith(".gz"):
        reference_seed = ''.join(random.choice(string.ascii_uppercase + string.digits + string.ascii_lowercase) for _ in range(50))
        run_this = """zcat {} > {}.temp.fa""".format(
            reference, reference_seed)
        runner(run_this)
        assert os.path.exists("{}.temp.fa".format(reference_seed))
        reference = "{}.temp.fa".format(reference_seed)
        ref_gzipped = True

    #make sure the gff file is ok
    gff_file = sys.argv[3]
    if not os.path.exists(gff_file):
        raise OSError("  - The gff file does not exist %s" % gff_file)
    gff_ok=False
    for chance in [".gff", ".gtf", ".gff3",
                   ".gff.gz", ".gtf.gz", ".gff3.gz"]:
        if chance in gff_file.split("/")[-1]:
            gff_ok=True
    if not gff_ok:
        raise IOError("  - The gff must have the ending '.gff', '.gtf', '.gff3', or be a gzipped variant.")

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
    gzipped = gff_file.endswith(".gz")
    if gzipped:
        f = gzip.open(gff_file, "rt")
    else:
        f = open(gff_file, "r")
    for line in f:
        line=line.strip()
        if line and line[0] != '#':
            seqs_in_gff.add(line.split()[0])
    f.close()

    with open(reference, "r") as f:
        for line in f:
            line=line.strip()
            if line[0] == ">":
                seqs_in_fasta.add(line.split()[0][1::])
    gff_not_in_fasta = [x for x in seqs_in_gff if x not in seqs_in_fasta]
    if len(gff_not_in_fasta) > 0:
        print("""The following scaffolds are in the gff file, but not
             in the assembly file.""", file=sys.stderr)
        print(gff_not_in_fasta)
        raise IOError("Scaffolds in gff not in fasta")

    # generate a fai file
    fai = "{}.fai".format(reference)
    if os.path.exists(fai):
        print("  - .fai index of reference exists", file=sys.stderr)
    else:
        print("  - Running samtools faidx", file=sys.stderr)
        run_this = "samtools faidx {}".format(reference)
        runner(run_this)
    assert os.path.exists(fai)

    # first generate the bed file of the whole genome
    print("  - Generating whole-genome bed.", file=sys.stderr)
    wholeGenome_bed = "{}_wholeGenome.bed".format(out_prefix)
    if os.path.exists(reference):
        print("    - Whole-genome bed exists.")
    else:
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
    assert os.path.exists(wholeGenome_bed)

    # now generate a coordinate file. This is used for bedtools, maybe
    print("  - Generating a whole-genome coordinates file", file=sys.stderr)
    wholeGenome_coords = "{}_size_wholeGenome.genfile".format(out_prefix)
    if os.path.exists(wholeGenome_coords):
        print("    - file already exists", file=sys.stderr)
    else:
        run_this = """sort -k1,1 -k2,2n {} | cut -f1,3 > {}""".format(
                wholeGenome_bed, wholeGenome_coords)
        runner(run_this)
    assert os.path.exists(wholeGenome_coords)

    # now generate a bed file of the transcript regions
    print("  - Generating a bed file of genic regions", file=sys.stderr)
    genic_bed = "{}_genic.bed".format(out_prefix)
    if os.path.exists(genic_bed):
        print("    - file already exists", file=sys.stderr)
    else:
        exon_ranges, introns_df = gff_coding_intervals(gff_file)
        gff_to_genic_bed(exon_ranges, genic_bed)
    assert os.path.exists(genic_bed)

    print(" - Making a file of exon/cds groups to verify parsing is correct.",
          file=sys.stderr)
    parentlist = "{}_exoncds_group_list.txt".format(out_prefix)
    if os.path.exists(parentlist):
        print("    - file already exists", file=sys.stderr)
    else:
        with open(parentlist, "w") as f:
            parent_list = sorted([x for x in exon_ranges])
            for parent in parent_list:
                print(parent, file=f)
    assert os.path.exists(parentlist)

    # now generate a bed file of the intergenic regions
    print("  - Generating a bed file of intergenic regions", file=sys.stderr)
    intergenic_bed = "{}_intergenic.bed".format(out_prefix)
    if os.path.exists(intergenic_bed):
        print("    - file already exists", file=sys.stderr)
    else:
        run_this = "bedtools complement -i {} -g {} > {}".format(
            genic_bed, wholeGenome_coords, intergenic_bed)
        runner(run_this)
    assert os.path.exists(intergenic_bed)

    # now merge the intergenic and genic, make sure it matches the genFile exactly
    print("  - Verifying that the genic and intergenic regions are complete", file=sys.stderr)
    temp_file = "{}_intergenic_genic_merge.temp".format(out_prefix)
    if os.path.exists(temp_file):
        print("    - file exists", file=sys.stderr)
    else:
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
    print("    - printing: ", stdout, file=sys.stderr)
    assert stdout.decode('utf-8').strip() == ""
    os.remove(temp_file)

    # now generate a bed file of the exonic regions
    print("  - Generating a bed file of exonic regions", file =sys.stderr)
    exonic_bed = "{}_exonic.bed".format(out_prefix)
    if os.path.exists(exonic_bed):
        print("    - file exists", file=sys.stderr)
    else:
        temp_bed = "{}_temp.bed".format(out_prefix)
        write_here = open(temp_bed, "w")

        gzipped = gff_file.endswith(".gz")
        if gzipped:
            f = gzip.open(gff_file, "rt")
        else:
            f = open(gff_file, "r")
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
        f.close()
        write_here.close()
        run_this = """sort -k1,1 -k2,2n {} | \
                   bedtools merge | bedtools sort > {}""".format(temp_bed, exonic_bed)
        runner(run_this)
        os.remove(temp_bed)
    assert os.path.exists(exonic_bed)

    # now generate a bed file of the intronic regions
    print("  - Generating a bed file of intronic regions", file=sys.stderr)
    intronic_bed = "{}_intronic.bed".format(out_prefix)
    if os.path.exists(intronic_bed):
        print("    - file exists.", file=sys.stderr)
    else:
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
    print("  - Generating a bed file of noncoding regions", file=sys.stderr)
    noncoding_bed = "{}_noncoding.bed".format(out_prefix)
    if os.path.exists(noncoding_bed):
        print("    - file exists", file=sys.stderr)
    else:
        run_this = """cat {} {} | bedtools sort | bedtools merge | \
                    sort -k1,1 -k2,2n > {}""".format(
            intronic_bed, intergenic_bed, noncoding_bed)
        runner(run_this)
    assert os.path.exists(noncoding_bed)


    # THIS WHOLE SECTION IS:
    #  for finding the transcripts that are inside the introns of other sequences
    print("  - Converting exon ranges to df.", file = sys.stderr)
    raw_tx_file = "{}_transcripts.txt".format(out_prefix)
    introns_file = "{}_introns.txt".format(out_prefix)
    if os.path.exists(raw_tx_file):
        print("    - file exists: {}".format(raw_tx_file), file=sys.stderr)
        raw_tx_df = pd.read_csv(raw_tx_file, sep="\t", header=0)
        if not os.path.exists(introns_file):
            raise IOError("introns file doesn't exist, but should")
        else:
            introns_df = pd.read_csv(introns_file, sep = "\t", header = 0)
    else:
        print("    - Couldn't find the raw transcript file: {}".format(raw_tx_file), file=sys.stderr)
        exon_ranges_exists = 'exon_ranges' in locals() or 'exon_ranges' in globals()
        if exon_ranges_exists:
            pass
        else:
            print("    - No information on exon ranges. Generating from gff.", file=sys.stderr)
            exon_ranges, introns_df = gff_coding_intervals(gff_file)
            print("    - Converting exon ranges to a dataframe.", file=sys.stderr)
        raw_tx_df = pd.DataFrame.from_records(exon_ranges).T
        raw_tx_df["start"] = raw_tx_df["start"].astype(int)
        raw_tx_df["stop"]  = raw_tx_df["stop"].astype(int)
        raw_tx_df.reset_index(inplace=True, drop=True)
        print("    - saving the transcripts file", file = sys.stderr)
        raw_tx_df.to_csv(raw_tx_file, index=False, sep="\t")
        print("    - saving the introns file", file = sys.stderr)
        introns_df.to_csv(introns_file, index=False, sep="\t")
    assert os.path.exists(raw_tx_file)
    assert os.path.exists(introns_file)

    print("  - Clipping off a few percent of the ends of exons for tolerance", file = sys.stderr)
    tx_df_ss = transcript_95per_start_stop(raw_tx_df)

    for intron_filter in ["noIFilt", "iFilt"]:
        if intron_filter == "iFilt":
            print("  - Getting rid of the 0.1% largest introns. These are often splice-lead variants, misannotations.", file =sys.stderr)
            introns_df = df_to_999p(introns_df)
        else:
            print("  - Running without filtering introns", file =sys.stderr)

        # now get the pairs of transcripts
        print("    - Finding genes in transcripts", file = sys.stderr)
        df_genes_in_introns = find_gene_in_intron(introns_df, tx_df_ss)
        genes_in_introns_file = "{}_genes_in_introns_{}.txt".format(out_prefix, intron_filter)
        df_genes_in_introns.to_csv(genes_in_introns_file, index=False, sep="\t")
        # End of section about finding transcripts inside the introns of other transcripts
        #######################################################################3

        print("    - Calculating genome stats", file=sys.stderr)
        #calculate % of the genome stats
        run_this = """awk '{{sum = sum + $2}} END{{print(sum)}}' {}""".format(wholeGenome_coords)
        wholeGenome_size = int(runner_w_output(run_this))
        run_this = """awk '{{sum = sum + $3 - $2 }} END{{print(sum)}}' {}""".format(exonic_bed)
        try:
            exonic_size = int(runner_w_output(run_this))
        except:
            print("    - Tried to run this: {}".format(run_this), file=sys.stderr)
        run_this = """awk '{{sum = sum + $3 - $2 }} END{{print(sum)}}' {}""".format(genic_bed)
        genic_size = int(runner_w_output(run_this))
        run_this = """awk '{{sum = sum + $3 - $2 }} END{{print(sum)}}' {}""".format(intergenic_bed)
        intergenic_size = int(runner_w_output(run_this))
        run_this = """awk '{{sum = sum + $3 - $2 }} END{{print(sum)}}' {}""".format(intronic_bed)
        intronic_size = int(runner_w_output(run_this))
        run_this = """awk '{{sum = sum + $3 - $2 }} END{{print(sum)}}' {}""".format(noncoding_bed)
        noncoding_size = int(runner_w_output(run_this))

        # num bp sense
        if "sense_anti" in df_genes_in_introns.columns:
            sense_gene_list = list(df_genes_in_introns.loc[df_genes_in_introns["sense_anti"] == "sense", "gene"].unique())
            bp_sense_exons = num_exon_bases_in_genelist(gff_file, sense_gene_list)
            # num bp antisense
            antisense_gene_list = list(df_genes_in_introns.loc[df_genes_in_introns["sense_anti"] == "antisense", "gene"].unique())
            bp_antisense_exons = num_exon_bases_in_genelist(gff_file, antisense_gene_list)
            # num bp both
            all_gene_list = list(df_genes_in_introns["gene"].unique())
            bp_all_exons = num_exon_bases_in_genelist(gff_file, all_gene_list)

            # now get in-intron genes that break splicing
            sense_df = df_genes_in_introns.loc[df_genes_in_introns["sense_anti"] == "sense",]
            breaks_splicing_df = sense_df.loc[sense_df["num_exons"] > 1, ]
            breaks_splicing_genelist = list(breaks_splicing_df["gene"].unique())
            bp_breaks_splicing = num_exon_bases_in_genelist(gff_file, breaks_splicing_genelist)

            # now get in-intron genes that don't break splicing
            doesnt_break_splicing_df = df_genes_in_introns[~df_genes_in_introns.index.isin(breaks_splicing_df.index)]
            doesnt_break_splicing_genes = list(doesnt_break_splicing_df["gene"].unique())
            doesnt_break_splicing_genes_bp = num_exon_bases_in_genelist(gff_file, doesnt_break_splicing_genes)
            # now get the num bp of single sense exons
            single_exon_sense_df = sense_df.loc[sense_df["num_exons"] == 1, ]
            single_exon_sense_genelist = list(single_exon_sense_df["gene"].unique())
            single_exon_sense_bp = num_exon_bases_in_genelist(gff_file, single_exon_sense_genelist)
        else:
            bp_sense_exons                 = 0
            bp_antisense_exons             = 0
            bp_all_exons                   = 0
            bp_breaks_splicing             = 0
            doesnt_break_splicing_genes_bp = 0
            single_exon_sense_bp           = 0

        genome_stats = "{}_genome_stats_{}.txt".format(out_prefix, intron_filter)
        print("    - Printing genome stats to {}".format(genome_stats), file=sys.stderr)
        outfile = open(genome_stats, "w")
        for writehere in [sys.stderr, outfile]:
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
            print("num_sense_exon_in_intron_bp\t{}\t{:.4f}\t(This is the percent of the exonic size)".format(
                  bp_sense_exons,
                  (bp_sense_exons/exonic_size)*100), file = writehere)
            print("num_antisense_exon_in_intron_bp\t{}\t{:.4f}\t(This is the percent of the exonic size)".format(
                  bp_antisense_exons,
                  (bp_antisense_exons/exonic_size)*100), file = writehere)
            print("num_antisense_and_sense_exon_in_intron_bp\t{}\t{:.4f}\t(This is the percent of the exonic size)".format(
                  bp_all_exons,
                  (bp_all_exons/exonic_size)*100), file = writehere)
            if bp_all_exons == 0:
                print("num_exonic_bp_breaks_splicing\t{}\tNaN\t(This is the percent of the exonic-inside-intron bp)".format(
                      bp_breaks_splicing), file= writehere)
                print("num_exonic_bp_doesnt_break_splicing\t{}\tNaN\t(This is the percent of the exonic-inside-intron bp)".format(
                      doesnt_break_splicing_genes_bp), file= writehere)
            else:
                print("num_exonic_bp_breaks_splicing\t{}\t{:.4f}\t(This is the percent of the exonic-inside-intron bp)".format(
                      bp_breaks_splicing,
                      (bp_breaks_splicing/bp_all_exons)*100), file= writehere)
                print("num_exonic_bp_doesnt_break_splicing\t{}\t{:.4f}\t(This is the percent of the exonic-inside-intron bp)".format(
                      doesnt_break_splicing_genes_bp,
                      (doesnt_break_splicing_genes_bp/bp_all_exons)*100), file= writehere)
            if bp_sense_exons == 0:
                print("single_exon_sense_bp\t{}\tNaN\t(This is the percent of the total sense exonic bp)".format(
                    single_exon_sense_bp), file= writehere)
            else:
                print("single_exon_sense_bp\t{}\t{:.4f}\t(This is the percent of the total sense exonic bp)".format(
                    single_exon_sense_bp,
                    (single_exon_sense_bp/bp_sense_exons)*100), file= writehere)
        outfile.close()
        # if we had to make a temp fasta, remove it and its index now.
    if ref_gzipped:
        os.remove(reference)
        os.remove("{}.fai".format(reference))

if __name__ == "__main__":
    main()
