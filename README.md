# local-flags

Experimentations with [flag algebras](http://people.cs.uchicago.edu/~razborov/files/flag.pdf) for graphs of maximum degree `Δ`, when both `Δ` and `n/Δ` can be large.

This module contains is a small set of helper function and a set of example based on the
[flag algebra library](https://docs.rs/flag-algebra).

## Requirements

You need to install cargo to compile Rust.
```
sudo apt install cargo
```
You need to have the `csdp` command line installed to solve semi-definite optimization problems.
```
sudo apt install coinor-csdp
```

## Usage

First go to the directory of a copy of the repo.
```
git clone git@github.com:avangogo/local-flags.git
cd local-flags
```
You can check that `cargo` is installed and that you are at the right place  with
```
cargo test
```

To run one of the examples of the `example/` folder, for instance `example/strong_complete.rs`.
```
cargo run --release --example strong_complete
```
The first compilation may be quite long. The first execution can also take time because the library needs to compute lists of graphs and the matrices of some flag operators. These later are stored in files for later computations. Eventually, the bottleneck is the SDP solver.