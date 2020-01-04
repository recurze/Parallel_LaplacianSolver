# Parallel Laplacian Solver using Random Walk
This repository will contain the implementation of the algorithm presented in
[this paper][paper]. The paper describes a random walk based method to solve an
important class of Laplacian Systems (Lx = b), called "one-sink" systems,
where exactly one of the coordinates of b is negative. The documentation of the
source code and the results of the expermentation can also be found in this
repo.

## Problem Statement
You are given an undirected positive weighted connected graph G = (V, E, w) with
adjacency matrix A<sub>uv</sub> = w<sub>uv</sub>. You are required to solve the
system of equations: Lx = b where L is the [Laplacian matrix][lapmat].

## Algorithm
The algorithm works by deriving the canonical solution from the stationary
state of the data collection process: Packets are generated at each node as an
independent Bernoulli process, transmitted to neighbors according to [stochastic
matrix][stomat] where P<sub>uv</sub> is directly propotional to w<sub>uv</sub>
and sunk at sink node only. Naturally, it consists of two phases: find parameter
&beta; such that DCP is ergodic and compute the stationary state, compute the
canonical solution by choosing an appropriate constant offset.

### Computing beta and stationary state
We have a lower limit for &beta;\* below which it's ergodic, so we binary
search this lower limit. Whenever it's not ergodic, there's one &eta; (queue
occupancy probability) which reaches 1, so we simply simulate the DCP at each
&beta; and check for this condition.

```
func estimateQueueOccupancyProbability(P, beta, J, T_mix, T_samp):
    We generate and transmit for T_mix seconds, allowing it to reach
    stationarity and then count the number of seconds queue is not empty for
    T_samp time, to estimate occupancy probability
begin
    t := 0
    Q_t(u) := 0 for all nodes
    cnt(u) := 0 for all nodes

    repeat
        generate packets at all nodes according to Bernoulli(beta J_u)
        transmit packets according to stochastic matrix, P
        increment cnt if queue is not empty and t > t_mix for all nodes
        t := t + 1
    until t <= T_mix + T_samp

    report: cnt/T_samp as the estimate
end

func computeStationaryState(P, b, t_hit):
    We start with beta = 1 and reduce by half everytime we find that it's
    non-ergodic
begin
    T_mix := 64t_hit log(1/e1)
    T_samp := 4logn/(k^2 e2^2)
    J := -b/b_sink

    beta := 1
    repeat
        beta := beta/2
        eta := estimateQueueOccupancyProbability(P, beta, J, T_mix, T_samp)
    until max(eta) < 0.75(1 - e1 - e2)

    report: eta as the stationary state
end
```

### Offset for canonical solution
The solution to Lx=b satisfies <<Lx, 1>> = 0 as &lambda;<sub>1</sub> = 0,
so we need to have an offset for stationary state such that this holds.

```
func computeCanonicalSolution(eta, beta, D, b):
    Shifts the stationary state, eta by z*
begin
    z* := -sum of eta(u)/d(u) for all u
    sum_d := sum of d(u) for all u
    x(u) := -b_sink(eta(u)/d(u) + z*d(u)/sum_d)/beta

    report: x as the solution to Lx=b
end
```

* Refer to the paper for a detailed analysis and description of the algorithm
  regarding constants.

## Implemenatation details

* There's a strong correlation between the fraction of packets sunk and the
  closeness to the stationary state. Thus, the stopping condition is based on
  this fraction C.
* More than 1 packet is transmitted at once.
* `inQ[i][j]` holds the incoming packets to node j coming from thread i
* The algorithm flounders in case of sparse graphs. The k-step speed up is a try
  at mitigating this issue.
* Each packet is moved by more than 1 step at a time; the path is noted and
  occupancy updated accordingly.
* `via[i][j][k]` is 1 if packet of thread i has been to node j in step k. This
  marks the occupancy of the queue of node j at that time step.

## IO format
The input will be of the following format

```
The first line in the file is two integers - n and m, denoting the number of
nodes and the number of edges. The following m lines will have three space
seperated real numbers, the 2 vertices and the weight of the corresponding
edge. The final line will have n space seperated real numbers, the b in Lx=b.

A[i][i] = 0
A[i][j] = A[j][i] >= 0
b[n] = sum(b[1..n-1])
```

Real world graphs can be found [here][konect] and [here][snap].

The output should be of the following format:

```
Print one and only line containing n space seperated real numbers, the solution
to Lx=b
```


## Help?

* Run `./gengraph.py -h` for help in generating the graph
* Run `./genb.py -h` for help in generating the RHS, b
* Run `cat g100.inp b100.inp > 100.inp` to merge the 2 segments of the input
* Run `make` from the parent directory to compile
* Run `./main` for help running the main program
* Use `run.sh` script to run automated tests
* `solve.py` uses least square or Jacobi method to compute the solution
* `compare.py` can be used to compare the output file and the actual answer

---

[paper]: https://arxiv.org/abs/1905.04989
[lapmat]: https://en.wikipedia.org/wiki/Laplacian_matrix#Definition
[stomat]: https://en.wikipedia.org/wiki/Stochastic_matrix#Definition_and_properties
[konect]: http://konect.uni-koblenz.de/networks]
[snap]: https://snap.stanford.edu/
