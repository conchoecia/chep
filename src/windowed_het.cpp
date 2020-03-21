#include <fstream>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <string>
#include <unistd.h>
#include <math.h>
#include <vector>
#include <map>
#include "windowed_het.hpp"

void help(void){
  printf(
  " \n"
  "usage (single-threaded): \n"
  "  samtools mpileup -f ref.fa reads_to_ref.bam | chep_windowed_het -f 10 -p 178 -w 50000\n"
  "\n"
  "usage (multi-threaded): \n"
  "  cat genome.bed | parallel -P 23 \"echo {} > {}.temp; samtools mpileup -f ref.fa -l {}.temp reads_to_ref.bam | chep_windowed_het -f 10 -p 178 -w 50000 >> unsorted.het; rm {}.temp\"\n"
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


int main(int argc, char **argv) {
  Vars vars = process_cl_args(argc, argv);

  std::cerr << "peak: " << vars.peak << "\n";
  std::cerr << "flank: " << vars.flank << "\n";
  std::cerr << "window: " << vars.window << "\n";

  uint tempmax = 178;
  uint acceptable_max = ceil(tempmax * 0.75);
  uint acceptable_min = tempmax * 0.25;


  //make data structure to store counts
  std::cout << "# eak=" << vars.peak << "\n";
  std::cout << "# flank=" << vars.flank << "\n";
  std::cout << "# window=" << vars.window << "\n";
  std::cout << "# Heterozygous sites are considered as:\n";
  std::cout << "#  For a read depth at col1(depth), a site with >= col2(min) and <= col(max)3\n";
  std::cout << "#   reads matching the ref is considered a heterozygous site.\n"; 
  std::cout << "# depth min max\n";
  std::map<uint,std::tuple<uint,uint>> table;
  for (uint i = vars.peak-vars.flank ; i <= vars.peak + vars.flank +1 ; i++){
    table[i]= {i*0.25, ceil(i*0.75)};
    std::cout << "# " << " " << i << " " << std::get<0>(table[i]) << " " << std::get<1>(table[i]) << "\n";
  }

  //parse the mpileup input
  std::string line;
  std::string prev_chrom = "";
  uint depth       = 0;
  uint start       = 1;
  uint stop        = 0;
  uint next_stop   = vars.window;
  uint this_pos    = 0;
  uint num_sites_measured = 0;
  uint num_het_sites      = 0;
  double het = 0;
  //print out the header
  std::cout << "chrom\ttarg_start\tstart\tstop\ttarg_stop\tnum_sites_mes\thet_sites\thet\n";

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
      if (prev_chrom.compare("") == 0){
        //we've just stepped into the mpileup output. don't do anything
      } else {
        het = 100.0*((float)num_het_sites/(float)num_sites_measured);
        if (isnan(het)){
          het = 0;
        }
        uint targ_start = next_stop - vars.window + 1;
        std::cout << prev_chrom << "\t" << targ_start << "\t" << start << "\t" << stop << "\t" << next_stop << "\t" << num_sites_measured << "\t" << num_het_sites << "\t" << het <<"\n";
        if (tokens[0].compare(prev_chrom) != 0){
          next_stop = vars.window;
        } else {
          next_stop += vars.window;
        }
      }
      num_sites_measured = 0;
      num_het_sites = 0;
      het = 0.0;
      prev_chrom = tokens[0];
      start = this_pos;
      stop = this_pos;
    }

    // this section determines whether the site is heterozygous or not
    depth = std::stoi(tokens[3]);
    if (depth >= vars.peak-vars.flank && depth <= vars.peak+vars.flank){
      //this site is acceptable to measure. add one to num_sites_measured
      num_sites_measured++;
      uint count = 0;
      uint acceptable_min = 0;
      for (char const &c: tokens[4]) {
        if (c == ',' || c == '.'){
          count += 1;
        }
      }
      if (count >= std::get<0>(table[depth]) && count <= std::get<1>(table[depth])){
        num_het_sites++;
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
  uint targ_start = next_stop - vars.window + 1;
  std::cout << prev_chrom << "\t" << targ_start << "\t" << start << "\t" << stop << "\t" << next_stop << "\t" << num_sites_measured << "\t" << num_het_sites << "\t" << het <<"\n";

  return 0;
}
