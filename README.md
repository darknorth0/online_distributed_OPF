# Online Distributed Optimal Power Flow with Equality Constraints
This repository is the official implementation of the paper [Online Distributed Optimal Power Flow with Equality Constraints.](https://arxiv.org/abs/)

## Abstract :

In this article, we investigate the distributed DC-Optimal Power Flow (DC-OPF) problem against a backdrop of online convex optimization with coupled constraints. While the classical OPF problem refers to a time-invariant optimization environment across the power distribution network, the online variant involves time-varying objectives and equality constraints, with the cost functions gradually being disclosed to the generating agents in the network. The agents (generators and loads) in the network are only privy to their own local objectives and constraints. To this end, we address the problem by proposing a distributed online algorithm based on the modified primal-dual approach, which deals with the online equality-constrained issue. It has been theoretically demonstrated that the proposed algorithm exhibits a sublinearly bounded static regret and constraint violation with a suitable choice of descending stepsizes. Finally, we corroborate the results by illustrating them with numerical examples.

## Usage : 

This implementation provides a setup as demonstrated in the paper for the IEEE-14 bus system of American Electric Power Systems (1962), networked over an undirected graph topology. 

## Requirements :

1. Code is good to go with MATLAB v2015b and above.
2. CVX package is required to be pre-installed with the MATLAB.

## Citation : 


```bibtex
@article{chatterjee2022online,
  title={Online Distributed Optimal Power Flow with Equality Constraints},
  author={Chatterjee, Sushobhan and Kalaimani, Rachel Kalpana},
  journal={arXiv preprint arXiv:},
  year={2022}
}
  }
```
