(Completed for a mini course project back in 2017 when I was an undergraduate. The program implements option pricing using a multi-level Monte Carlo simulation. The idea was to generate some plots and compare to theory, but admittedly the code itself never really got much attention after this happened.)

Requires GSL for random number generation.

Compile with 

$ clang main.c -l gsl -O2

Run with

$ GSL_RNG_SEED=4234672 ./a.out

where you replace the random seed with any number. Or to put in the current time:

$ GSL_RNG_SEED="$(date +%s%3)" ./a.out
