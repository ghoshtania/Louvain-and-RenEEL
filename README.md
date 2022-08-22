# Louvain-and-RenEEL
Usage: To use the RenEEL scheme for maximizing modularity (Q) using Louvain algorithm as base algorithm, follow the steps below:

1.Prepare data
1.1 Use an edgelist file including 2 columns separated by space or tab with no header. The network should be unweighted, undirected. Self-loops will be ignored.
(See example input file karate.txt)
1.2 Use bash script (work.sh) to generate the three files required by the program.

Example: sh work.sh graph.txt

2.Compile code
compile main.c, help.c and rg.c with required libraries (math).

Example: g++ main.cpp louvain.cpp community.cpp graph.cpp -lm 

This will generate file a.out.

3. Run a.out with 2 arguments.
argument 1: Positive Integer, maximum and initial ensemble size of partitions used in RenEEL
argument 2: Positive Integer, ensemble size of partitions of the reduced network for iteration part in RenEEL
(seed of random number will be generated using system time)

Example: ./a.out 10 5
