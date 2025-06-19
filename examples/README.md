# Examples

This directory contains a notebook [PottsOn.ipynb](./PottsOn.ipynb) in which crossing symmetry equations are solved for a few four-point correlations of the loop $O(n)$ model. This notebook also includes performance comparisons with the [previous version](https://gitlab.com/s.g.ribault/Bootstrap_Virasoro.git) of this code, which was written in Python.

# Julia setup
The directory also contains an example of a `startup.jl` file. I recommend putting a copy of this file in your `$HOME/.julia/config` folder, where `$HOME` is your home directory. This will ensure you have a few general purpose Julia packages installed and available, which makes the experience with Julia more comfortable.

I also recommend adding a line

```shell
export JULIA_NUM_THREADS=auto
```

in your shell configuration (`$HOME/.bashrc`, `$HOME/.zshrc`, or whatever shell you are using as the default shell on your system), which will automatically start all instances of julia with as many threads as available.
