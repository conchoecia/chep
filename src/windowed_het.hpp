typedef struct Vars {
  uint debug=0;
  uint peak=0;
  uint flank=0;
  uint window=0;
} Vars;

Vars process_cl_args(int argc, char **argv);
