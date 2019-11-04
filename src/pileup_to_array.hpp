typedef struct Vars {
  uint debug=0;
  uint min;
  uint max;
} Vars;

Vars process_cl_args(int argc, char **argv);
