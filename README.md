#Parallel Laplacian Solver using Random Walk
This repository will contain the implementation of the algorithm presented in
[this paper][paper]. The paper describes a random walk based method to solve an
important class of Laplacian Systems (Lx = b), called "one-sink" systems,
where exactly one of the coordinates of **b** is negative.

##Problem Statement
You are given an undirected positive weighted connected graph G = (V, E, w) with
adjacency matrix A\_uv = w(u, v). You are required to solve the system of
equations: Lx = b where L is the [Laplacian matrix][lapmat].

[paper]: https://arxiv.org/abs/1905.04989
[lapmat]: https://en.wikipedia.org/wiki/Laplacian_matrix#Definition
