#include <fstream>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <string>
#include <unistd.h>
#include <vector>
#include <map>
#include "pileup_to_array.hpp"

void help(void){
  printf(
  "Program: pileup_to_array\n"
  "Author: DTS@ucsc.edu\n\n"

  "Description: \n"
  "             \n\n"

  "Params:\n"
  "  -h        Print this help message\n"
        );
}

Vars process_cl_args(int argc, char **argv){
  Vars vars;
  int c;
  // args
  vars.min=5;
  vars.max=500;
  while( (c=getopt( argc, argv, "dm:M:h" )) != -1 ) {
    switch(c) {
    case 'd' :
      vars.debug=1;
      break;
    case 'm' :
      vars.min=atoi(optarg);
      break;
    case 'M' :
      vars.max=atoi(optarg);
      break;
    case 'h' :
      help();
    }
  }
  return(vars);
}

int main(int argc, char **argv) {
  Vars vars = process_cl_args(argc, argv);
  printf("min = %d\n", vars.min);
  printf("max = %d\n", vars.max);

  //make data structure to store counts
  std::map<uint,std::map<uint,uint>> table;
  for (uint i = vars.min ; i <= vars.max ; i++){
    table[i]= {};
    for (uint j = 0 ; j <= i ; j++){
      table[i][j] = 0;
    }
  }

  //for (std::string line; std::getline(std::cin, line);) {
  //  std::cout << line << std::endl;
  //}

  //parse the mpileup input
  std::string line;
  while(std::getline(std::cin, line)) {     // '\n' is the default delimiter
    std::vector<std::string> tokens;
    std::istringstream iss(line);
    std::string token;
    while(std::getline(iss, token, '\t'))   // but we can specify a different one
      tokens.push_back(token);
    std::cout << tokens[3] << "\t" << tokens[4] << "\n";
  }

  //print out the entire table
  for(auto const &key1 : table) {
    for(auto const &key2 : key1.second) {
      std::cout << key1.first << "\t" << key2.first << "\t" << table[key1.first][key2.first] << "\n";
    }
  }
  return 0;
}
