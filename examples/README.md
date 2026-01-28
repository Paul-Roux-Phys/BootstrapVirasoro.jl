# Julia setup
The directory also contains an example of a `startup.jl` file. I recommend putting a copy of this file in your `$HOME/.julia/config` folder, where `$HOME` is your home directory. This will ensure you have a few general purpose Julia packages installed and available, which makes the experience with Julia more comfortable.

I also recommend adding a line

```shell
export JULIA_NUM_THREADS=auto
```

in your shell configuration (`$HOME/.bashrc`, `$HOME/.zshrc`, or whatever shell you are using as the default shell on your system), which will automatically start all instances of julia with as many threads as available.

# Examples

* The file [compute_blocks.jl](./compute_blocks.jl) contains an example of how to use the package to compute Virasoro conformal blocks. It can be ran from a terminal by doing `julia ./compute_blocks.jl`, or inside of a julia session with `include("bootstrap.jl")`.

* The file [bootstrap.jl](./bootstrap.jl) contains an example of bootstrapping a sphere four-point function of 4 degenerate fields. It can be ran from a terminal by doing `julia ./bootstrap.jl`, or inside of a julia session with `include("bootstrap.jl")`.
