# Vectorized-Array-Multiplication-Reduction-using-SSE
Introduction

There are many problems in scientific and engineering computing where you want to multiply arrays of numbers together and add up all the multiplies to produce a single sum (Fourier transformation, convolution, autocorrelation, etc.): sum = ΣA[i]*B[i]

This project is to test array multiplication/reduction usinf SIMD and non-SIMD.

For the "control groups" benchmarks, do not use OpenMP parallel for-loops. Just use straight C/C++ for-loops. In this project, we are only using OpenMP for the timing.

Requirements

    Use the supplied SIMD SSE intrinsics code to run an array multiplication/reduction timing experiment. Run the same experiment a second time using your own C/C++ array multiplication/reduction code.

    Use different array sizes from 1K to 8M. The choice of in-between values is up to you, but pick something that will make for a good graph.

    Run each array-size test a certain number of trials. Use the peak value for the performance you record.

    Create a table and a graph showing SSE/Non-SSE speed-up as a function of array size. Speedup in this case will be S = Psse/Pnon-sse = Tnon-sse/Tsse (P = Performance, T = Elapsed Time).

    Note: this is not a multithreading assignment, so you don't need to worry about a NUMT. Don't use any OpenMP-isms except for getting the timing.

    The Y-axis performance units in this case will be "Speed-Up", i.e., dimensionless.

    Parallel Fraction doesn't apply to SIMD parallelism, so don't compute one. 
