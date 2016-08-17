# tg

A two dimensional Taylor-Green Vortex with chemical reaction for accessing load balancing tools.

## Prerequisites

- tvmet
- chemkinpp `svn co https://collab.sam.pitt.edu/svn/chemkinpp/trunk chemkinpp`
- paragon git submodule
- zoltan
- metis `svn co https://collab.sam.pitt.edu/svn/metis/trunk metis`

## Getting Started

### Building on mpi.sam.pitt.edu

* `. build-mpi.sh`
* `make -j4`

#### Running interactively

* allocation resources: `salloc --nodes=2 --time=01:00:00`
* submit interactive job: `srun --nodes=1 --tasks-per-node=1 ./tg simparam.in`
