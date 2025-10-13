#include <algorithm>
#include <fstream>
#include <iostream>
#include <list>
#include <locale>
#include <numeric>
#include <sstream>
#include <stdio.h>
#include <string>
#include <unistd.h>
#include <math.h>
#include <vector>
#include <map>
#include "windowed_het.hpp"

/* in this array we keep track of the following varibles
   0 - count of As Aa
   1 - count of Cs Cc
   2 - count of Ts Tt
   3 - count of Gs Gg
   4 - matches to ref ,. 44 46
   5 - insertions +
   6 - deletions -
   7 - indel *#  42 35
   8 - junk - everything else
   */
const unsigned char seq_nt4_table[256] = { // translate ACGT to 0123
    0, 1, 2, 3,  8, 8, 8, 8,  8, 8, 8, 8,  8, 8, 8, 8,
    8, 8, 8, 8,  8, 8, 8, 8,  8, 8, 8, 8,  8, 8, 8, 8,
    8, 8, 8, 7,  8, 8, 8, 8,  8, 8, 7, 5,  4, 6, 4, 8,
    8, 8, 8, 8,  8, 8, 8, 8,  8, 8, 8, 8,  8, 8, 8, 8,
    8, 0, 8, 1,  8, 8, 8, 2,  8, 8, 8, 8,  8, 8, 8, 8,
    8, 8, 8, 8,  3, 3, 8, 8,  8, 8, 8, 8,  8, 8, 8, 8,
    8, 0, 8, 1,  8, 8, 8, 2,  8, 8, 8, 8,  8, 8, 8, 8,
    8, 8, 8, 8,  3, 3, 8, 8,  8, 8, 8, 8,  8, 8, 8, 8,
    8, 8, 8, 8,  8, 8, 8, 8,  8, 8, 8, 8,  8, 8, 8, 8,
    8, 8, 8, 8,  8, 8, 8, 8,  8, 8, 8, 8,  8, 8, 8, 8,
    8, 8, 8, 8,  8, 8, 8, 8,  8, 8, 8, 8,  8, 8, 8, 8,
    8, 8, 8, 8,  8, 8, 8, 8,  8, 8, 8, 8,  8, 8, 8, 8,
    8, 8, 8, 8,  8, 8, 8, 8,  8, 8, 8, 8,  8, 8, 8, 8,
    8, 8, 8, 8,  8, 8, 8, 8,  8, 8, 8, 8,  8, 8, 8, 8,
    8, 8, 8, 8,  8, 8, 8, 8,  8, 8, 8, 8,  8, 8, 8, 8,
    8, 8, 8, 8,  8, 8, 8, 8,  8, 8, 8, 8,  8, 8, 8, 8
};

const unsigned char index_to_char[256] = { // translate 0123 to ACGT
  'A', 'C', 'G', 'T'};


void help(void){
  printf(
  " \n"
  "usage (single-threaded): \n"
  "  samtools mpileup -af ref.fa reads_to_ref.bam | chep_windowed_het -f 10 -p 178 -w 50000\n"
  "\n"
  "usage (multi-threaded): \n"
  "  cat genome.bed | parallel -P 23 \"echo {} > {}.temp; samtools mpileup -af ref.fa -l {}.temp reads_to_ref.bam | chep_windowed_het -f 10 -p 178 -w 50000 >> unsorted.het; rm {}.temp\"\n"
  "  For a better explanation, see https://github.com/conchoecia/chep \n"
  "Program: windowed het\n"
  "Author: DTS@ucsc.edu\n\n"

  "Description: \n"
  "            \n\n"

  "Params:\n"
  "  -h        Print this help message\n"
        );
  exit(3);
}

Vars process_cl_args(int argc, char **argv){
  Vars vars;
  int c;
  // args
  if (argc == 1){
    help();
  }
  while( (c=getopt( argc, argv, "dp:f:w:h" )) != -1 ) {
    switch(c) {
    case 'd' :
      vars.debug=1;
      break;
    case 'p' :
      vars.peak=strtoul(optarg, NULL, 0);
      break;
    case 'f' :
      vars.flank=strtoul(optarg, NULL, 0);
      break;
    case 'w' :
      vars.window=strtoul(optarg, NULL, 0);
      break;
    case 'h' :
      help();
    }
  }
  if (vars.peak == 0){
    std::cout << "You must provide a parameter for -p. The value 0 is not allowed.\n";
    help();}
  if (vars.flank == 0){
    std::cout << "You must provide a parameter for -f. The value 0 is not allowed.\n";
    help();}
  if (vars.window == 0){
    std::cout << "You must provide a parameter for -f. The value 0 is not allowed.\n";
    help();}
  return(vars);
}


int test_table(){
  /* prints out the indices of the characters matching the tables.
     Should be
   0 - Aa
   1 - Cc
   2 - Tt
   3 - Gg
   4 - ,.
   5 -  +
   6 -  -
   7 -  *#
  */
  char print_these[] = {'A', 'a', 'C', 'c',
                        'G', 'g', 'T', 't',
                        ',', '.', '+', '-',
                        '*', '#'};
  for (char c : print_these){
    std::cout << c << ": " << (int)seq_nt4_table[c] << "\n";
  }
  return 0;
}


int main(int argc, char **argv) {
  //test_table();
  Vars vars = process_cl_args(argc, argv);

  std::cerr << "peak: " << vars.peak << "\n";
  std::cerr << "flank: " << vars.flank << "\n";
  std::cerr << "window: " << vars.window << "\n";

  uint tempmax = vars.peak;
  uint acceptable_max = ceil(tempmax * 0.75);
  uint acceptable_min = tempmax * 0.25;


  //make data structure to store counts
  std::cout << "# peak=" << vars.peak << "\n";
  std::cout << "# flank=" << vars.flank << "\n";
  std::cout << "# window=" << vars.window << "\n";
  std::cout << "# Heterozygous sites are considered as:\n";
  std::cout << "#  For a read depth at col1(depth), a site with >= col2(min) and <= col(max)3\n";
  std::cout << "#   reads matching the ref is considered a heterozygous site.\n"; 
  std::cout << "# depth min max\n";
  std::map<uint,std::tuple<uint,uint>> table;
  for (uint i = vars.peak-vars.flank ; i <= vars.peak + vars.flank +1 ; i++){
    table[i] = std::make_tuple(static_cast<uint>(i*0.25), static_cast<uint>(ceil(i*0.75)));
    std::cout << "# " << " " << i << " " << std::get<0>(table[i]) << " " << std::get<1>(table[i]) << "\n";
  }

  //parse the mpileup input
  std::string line;
  std::list<uint> depth_counter;
  std::string prev_chrom = "";
  std::string base = "";
  uint depth       = 0;
  uint start       = 1;
  uint stop        = 0;
  uint next_stop   = vars.window;
  uint this_pos    = 0;
  uint num_sites_measured = 0;
  uint num_het_sites      = 0;
  uint num_insertion_sites= 0;
  uint num_deletion_sites= 0;
  double num_insertion_bases = 0;
  double num_deletion_bases = 0;
  uint num_tt_AC = 0;
  uint num_tt_AG = 0;
  uint num_tt_AT = 0;
  uint num_tt_CG = 0;
  uint num_tt_CT = 0;
  uint num_tt_GT = 0;
  uint num_GC      = 0;
  uint num_non_N   = 0;
  double pergc = 0;
  double het = 0;
  //print out the header
  std::cout << "chrom\ttarg_start\t";
  std::cout << "start\tstop\ttarg_stop\t";
  std::cout << "num_sites_mes\thet_sites\t";
  std::cout << "het\tdepth_med\tdepth_mean\t";
  std::cout << "depth_sd\tper_gc\t";
  std::cout << "num_insertion_sites\tnum_deletion_sites\t";
  std::cout << "num_insertion_bases\tnum_deletion_bases\t";
  std::cout << "num_indel_sites\tnum_indel_bases\t";
  std::cout << "per_insertion_sites\tper_insertion_bases\t";
  std::cout << "per_deletion_sites\tper_deletion_bases\t";
  std::cout << "per_indel_sites\tper_indel_bases\t";
  std::cout << "num_tt_AC\tnum_tt_AG\t";
  std::cout << "num_tt_AT\tnum_tt_CG\t";
  std::cout << "num_tt_CT\tnum_tt_GT\t";
  std::cout << "per_tt_AC\tper_tt_AG\t";
  std::cout << "per_tt_AT\tper_tt_CG\t";
  std::cout << "per_tt_CT\tper_tt_GT\t";
  std::cout << "num_transitions\tnum_transversions\t";
  std::cout << "per_transitions\tper_transversions\t";
  std::cout << "transition_div_transversion\n";

  //this iterates through every line
  while(std::getline(std::cin, line)) {     // '\n' is the default delimiter
    std::vector<std::string> tokens;
    std::istringstream iss(line);
    std::string token;
    while(std::getline(iss, token, '\t'))   // but we can specify a different one
      tokens.push_back(token);
    //this_chrom = tokens[0];
    this_pos   = std::stoi(tokens[1]);

    // first check if we need to output the data
    // outformat is:
    // chrom start stop size num_sites_measured het_sites het
    if (this_pos > next_stop || tokens[0].compare(prev_chrom) != 0){
      //we have reached the window size or have switched to a new chromosome
      // all this code just outputs all the stats for this window
      if (prev_chrom.compare("") == 0){
        //we've just stepped into the mpileup output. don't do anything
      } else {
        // calculate the heterozygosity
        het = 100.0*((float)num_het_sites/(float)num_sites_measured);
        if (isnan(het)){
          het = 0;
        }
        //calculate the median read depth
        depth_counter.sort();
        uint median_depth = *std::next(depth_counter.begin(), depth_counter.size()/2);
        //calculate the mean read depth
        double mean_depth = std::accumulate(depth_counter.begin(), depth_counter.end(), 0.0) / depth_counter.size();
        // calculate the sd of the read depth
        double sd = 0;
        for (auto const& v : depth_counter) {
          sd += pow(v - mean_depth, 2);
        }
        sd=sqrt(sd/depth_counter.size());
        // calculate the percent GC
        pergc = 100*((float)num_GC/(float)num_non_N);
        //calculate the target start
        uint targ_start = next_stop - vars.window + 1;
        //other calculations
        uint num_indel_sites = num_insertion_sites + num_deletion_sites;
        uint num_indel_bases = num_insertion_bases + num_deletion_bases;
        double per_insertion_sites = 100*((float)num_insertion_sites/(float)num_sites_measured);
        double per_insertion_bases = 100*((float)num_insertion_bases/(float)num_sites_measured);
        double per_deletion_sites = 100*((float)num_deletion_sites/(float)num_sites_measured);
        double per_deletion_bases = 100*((float)num_deletion_bases/(float)num_sites_measured);
        double per_indel_sites = 100*((float)num_indel_sites/(float)num_sites_measured);
        double per_indel_bases = 100*((float)num_indel_bases/(float)num_sites_measured);
        uint num_transitions = num_tt_AG + num_tt_CT;
        uint num_transversions = num_tt_AC + num_tt_AT + num_tt_CG + num_tt_GT;
        double per_tt_AC = 100*((float)num_tt_AC/((float)num_het_sites + 0.000000001));
        double per_tt_AT = 100*((float)num_tt_AT/((float)num_het_sites + 0.000000001));
        double per_tt_CT = 100*((float)num_tt_CT/((float)num_het_sites + 0.000000001));
        double per_tt_AG = 100*((float)num_tt_AG/((float)num_het_sites + 0.000000001));
        double per_tt_CG = 100*((float)num_tt_CG/((float)num_het_sites + 0.000000001));
        double per_tt_GT = 100*((float)num_tt_GT/((float)num_het_sites + 0.000000001));
        double per_transitions = 100*((float)num_transitions/((float)num_het_sites + 0.000000001));
        double per_transversions = 100*((float)num_transversions/((float)num_het_sites + 0.000000001));
        double transition_div_transversion = (float)num_transitions/((float)num_transversions + 0.000000001);
        std::cout << prev_chrom << "\t" << targ_start << "\t"
                  << start << "\t" << stop << "\t" << next_stop << "\t"
                  << num_sites_measured << "\t" << num_het_sites << "\t"
                  << het << "\t" << median_depth << "\t" << mean_depth << "\t"
                  << sd << "\t" << pergc << "\t"
                  << num_insertion_sites << "\t" << num_deletion_sites << "\t"
                  << num_insertion_bases << "\t" << num_deletion_bases << "\t"
                  << num_indel_sites << "\t" << num_indel_bases << "\t"
                  << per_insertion_sites << "\t" << per_insertion_bases << "\t"
                  << per_deletion_sites << "\t" << per_deletion_bases << "\t"
                  << per_indel_sites << "\t" << per_indel_bases << "\t"
                  << num_tt_AC << "\t" << num_tt_AG << "\t"
                  << num_tt_AT << "\t" << num_tt_CG << "\t"
                  << num_tt_CT << "\t" << num_tt_GT << "\t"
                  << per_tt_AC << "\t" << per_tt_AG << "\t"
                  << per_tt_AT << "\t" << per_tt_CG << "\t"
                  << per_tt_CT << "\t" << per_tt_GT << "\t"
                  << num_transitions << "\t" << num_transversions << "\t"
                  << per_transitions << "\t" << per_transversions << "\t"
                  << transition_div_transversion << "\n";
        if (tokens[0].compare(prev_chrom) != 0){
          next_stop = vars.window;
        } else {
          next_stop += vars.window;
        }
      }
      depth_counter.clear();
      num_sites_measured = 0;
      num_het_sites = 0;
      num_insertion_sites= 0;
      num_deletion_sites= 0;
      num_insertion_bases = 0;
      num_deletion_bases = 0;
      num_tt_AC = 0;
      num_tt_AG = 0;
      num_tt_AT = 0;
      num_tt_CG = 0;
      num_tt_CT = 0;
      num_tt_GT = 0;
      het = 0.0;
      prev_chrom = tokens[0];
      start = this_pos;
      stop = this_pos;
      num_GC=0;
      num_non_N=0;
    }

    // do some things to help calculate %GC
    base = tokens[2];
    std::transform(base.begin(), base.end(), base.begin(), ::toupper);
    if (base.compare("N") != 0){
      num_non_N++;
      if (base.compare("G") == 0 || base.compare("C") == 0){
          num_GC++;
      }
    }
    // this section determines whether the site is heterozygous or not
    depth = std::stoi(tokens[3]);
    depth_counter.push_back(depth);
    if (depth >= vars.peak-vars.flank && depth <= vars.peak+vars.flank){
      //this site is acceptable to measure. add one to num_sites_measured
      num_sites_measured++;
      uint acceptable_min = 0;

      // these variables are used to keep track of what type of site this is
      uint ref_match_count = 0;
      /* in this array we keep track of the following varibles
         0 - count of As Aa
         1 - count of Cs Cc
         2 - count of Ts Tt
         3 - count of Gs Gg
         4 - matches to ref ,. 44 46
         5 - insertions +
         6 - deletions -
         7 - indel *#  42 35
         8 - junk - everything else
         */
      uint type_counter[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0};

      // these variables are used in parsing the indel information
      uint bp_insertions = 0; //just the sum of the total bp of indels. Get mean later
      uint bp_deletions = 0;
      std::string indel_string = "";
      bool in_indel = 0;
      bool insertion_or_deletion = 0; //insertion is 0, deletion is 1

      uint skip_counter = 0;
      //this block parses the string and turns it into raw counts
      for (char const &c: tokens[4]) {
        if (in_indel == 1){
          //std::cout << c << ": " << skip_counter << "\n";
          if (skip_counter != 0){
            //raise an exception. we should never be skipping characters while we're inside an indel
            std::cerr << " - The program attempted to skip lines while parsing the output of mpileup. This is not possible.\n";
            std::cerr << tokens[4] << "\n";
            return 3;
          }
          //this will be a number, or we have finished recording the number and stepped into the indel string
          if ( isdigit(c) == 1){
            //this digit is one of the numbers showing the length of the insertion
            indel_string.push_back(c);
          }
          else{
            //otherwise we just stepped into the actual indel. record the indel size, increase the count, skip the next size(indel)-1 chars
            in_indel = 0;
            // we already recorded the occurence of the indel, just record the length now
            if ( insertion_or_deletion == 0 ){
              bp_insertions += stoi(indel_string);
            } else {
              bp_deletions += stoi(indel_string);
            }
            skip_counter = stoi(indel_string) - 1;
            indel_string = "";
            //std::cout << " skip counter: " << skip_counter << "\n";
          }
        }
        else {
          // we're not parsing an indel
          if (skip_counter != 0){
            //We're just skipping this character for some reason
            skip_counter--;
          }
          else {
            //we're not skipping the character
            if ( c == '^' ){
              skip_counter = 1;
            }
            else if ( c == '+' || c == '-' ) {
              type_counter[seq_nt4_table[c]]++;
              insertion_or_deletion = seq_nt4_table[c] - 5; //uses the lookup table //insertions 0, del 5
              in_indel = 1;
              //std::cout << " found indel: " << insertion_or_deletion << "\n";
            } else{
              //it's just a plain character that we can parse
              type_counter[seq_nt4_table[c]]++;
            }
          }
        }
      }
      // mark if this was an indel-originating site or a het-originating site
      if (type_counter[5] >= std::get<0>(table[depth]) && type_counter[5] <= std::get<1>(table[depth])){
        num_insertion_bases += (float)bp_insertions / (float)type_counter[5];
        num_insertion_sites++;
      }
      if (type_counter[6] >= std::get<0>(table[depth]) && type_counter[6] <= std::get<1>(table[depth])){
          num_deletion_bases += (float)bp_deletions / (float)type_counter[6];
          num_deletion_sites++;
      }

      // not figure out if this is a het site
      if (type_counter[4] >= std::get<0>(table[depth]) && type_counter[4] <= std::get<1>(table[depth])){
        // the site is heterozygous, but we don't know if it is from a base or indel at this point
        uint max_int = 0;
        int i_of_max_int = -1;
        uint indices[] = {0,1,2,3};
        for(uint i: indices){
          if (type_counter[i] > max_int){
            max_int = type_counter[i];
            i_of_max_int = i;
          }
        }

        // print this for some diagnostics
        //std::cout << "[ ";
        //for (uint n : type_counter){
        //  std::cout << n << ", ";
        //}
        //std::cout << "]\n";
        //std::cout << "Maxind: " << i_of_max_int << "\n" << tokens[4] << "\n\n";

        if (i_of_max_int != -1){
          //figure out what type of transition to transversion this is
          std::string mut_type = tokens[2];
          mut_type[0] = std::toupper(mut_type[0]);
          mut_type.push_back(index_to_char[i_of_max_int]);
          std::toupper(mut_type[1]);
          //std::for_each(mut_type.begin(), mut_type.end(), [](char & c){
          //  c = ::toupper(c)}
          std::sort(mut_type.begin(), mut_type.end());
          //can't use switch with std::strings
          if (mut_type == "AC"){      num_tt_AC++; }
          else if (mut_type == "AG"){ num_tt_AG++; }
          else if (mut_type == "AT"){ num_tt_AT++; }
          else if (mut_type == "CG"){ num_tt_CG++; }
          else if (mut_type == "CT"){ num_tt_CT++; }
          else if (mut_type == "GT"){ num_tt_GT++; }
          else {
            std::cout << "found some string we couldnt' parse: " << mut_type << "\n";
          }
          num_het_sites++; // this is some sort of het site, so add it
        }
      }
    }

    //exit modifications. now log where we were.
    stop = this_pos;
  }
  //we're at the end. print out the het for the last window
  het = 100.0*((float)num_het_sites/(float)num_sites_measured);
  if (isnan(het)){
    het = 0;
  }
  //calculate the median read depth
  depth_counter.sort();
  uint median_depth = *std::next(depth_counter.begin(), depth_counter.size()/2);
  //calculate the mean read depth
  double mean_depth = std::accumulate(depth_counter.begin(), depth_counter.end(), 0.0) / depth_counter.size();
  // calculate the sd of the read depth
  double sd = 0;
  for (auto const& v : depth_counter) {
    sd += pow(v - mean_depth, 2);
  }
  sd=sqrt(sd/depth_counter.size());
  // calculate the percent GC
  pergc = 100*((float)num_GC/(float)num_non_N);
  //calculate the target start
  uint targ_start = next_stop - vars.window + 1;
  //other calculations
  uint num_indel_sites = num_insertion_sites + num_deletion_sites;
  uint num_indel_bases = num_insertion_bases + num_deletion_bases;
  double per_insertion_sites = 100*((float)num_insertion_sites/(float)num_sites_measured);
  double per_insertion_bases = 100*((float)num_insertion_bases/(float)num_sites_measured);
  double per_deletion_sites = 100*((float)num_deletion_sites/(float)num_sites_measured);
  double per_deletion_bases = 100*((float)num_deletion_bases/(float)num_sites_measured);
  double per_indel_sites = 100*((float)num_indel_sites/(float)num_sites_measured);
  double per_indel_bases = 100*((float)num_indel_bases/(float)num_sites_measured);
  uint num_transitions = num_tt_AG + num_tt_CG;
  uint num_transversions = num_tt_AC + num_tt_AT + num_tt_CG + num_tt_GT;
  double per_tt_AC = 100*((float)per_tt_AC/((float)num_het_sites + 0.000000001));
  double per_tt_AT = 100*((float)per_tt_AT/((float)num_het_sites + 0.000000001));
  double per_tt_CT = 100*((float)per_tt_CT/((float)num_het_sites + 0.000000001));
  double per_tt_AG = 100*((float)per_tt_AG/((float)num_het_sites + 0.000000001));
  double per_tt_CG = 100*((float)per_tt_CG/((float)num_het_sites + 0.000000001));
  double per_tt_GT = 100*((float)per_tt_GT/((float)num_het_sites + 0.000000001));
  double per_transitions = 100*((float)num_transitions/((float)num_het_sites + 0.000000001));
  double per_transversions = 100*((float)num_transversions/((float)num_het_sites + 0.000000001));
  double transition_div_transversion = (float)num_transitions/((float)num_transversions + 0.000000001);
  std::cout << prev_chrom << "\t" << targ_start << "\t"
            << start << "\t" << stop << "\t" << next_stop << "\t"
            << num_sites_measured << "\t" << num_het_sites << "\t"
            << het << "\t" << median_depth << "\t" << mean_depth << "\t"
            << sd << "\t" << pergc << "\t"
            << num_insertion_sites << "\t" << num_deletion_sites << "\t"
            << num_insertion_bases << "\t" << num_deletion_bases << "\t"
            << num_indel_sites << "\t" << num_indel_bases << "\t"
            << per_insertion_sites << "\t" << per_insertion_bases << "\t"
            << per_deletion_sites << "\t" << per_deletion_bases << "\t"
            << per_indel_sites << "\t" << per_indel_bases << "\t"
            << num_tt_AC << "\t" << num_tt_AG << "\t"
            << num_tt_AT << "\t" << num_tt_CG << "\t"
            << num_tt_CT << "\t" << num_tt_GT << "\t"
            << per_tt_AC << "\t" << per_tt_AG << "\t"
            << per_tt_AT << "\t" << per_tt_CG << "\t"
            << per_tt_CT << "\t" << per_tt_GT << "\t"
            << num_transitions << "\t" << num_transversions << "\t"
            << per_transitions << "\t" << per_transversions << "\t"
            << transition_div_transversion << "\n";

  return 0;
}
