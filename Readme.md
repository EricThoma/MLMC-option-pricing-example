Requires GSL for random number generation.

Compile with 

$ clang main.c -l gsl -O2

Run with

$ GSL_RNG_SEED=4234672 ./a.out

where you replace the random seed with any number. Or to put in the current time:

$ GSL_RNG_SEED="$(date +%s%3)" ./a.out
