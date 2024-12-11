# Coyan

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11191730.svg)](https://doi.org/10.5281/zenodo.11191730)

Coyan is a Rust project that transforms Static Fault Trees into CNF equisatisfiables formulas. Then, with the use of a Model Counter,
computes the Top Event Probability of failure at a given time bound. Also, it allows the creation of random Fault Trees with specific configurations.

## Installation

```bash
git clone https://github.com/NazaGara/Coyan.git
cargo build --release 
```

## Usage

The tool currently posses the modules for `Fault Tree Analysis` [coyan_fta](https://github.com/NazaGara/Coyan/tree/main/coyan_fta) and for `Random Fault Tree Generation` [coyan_rft](https://github.com/NazaGara/Coyan/tree/main/coyan_rft). Both of them can be used via the command line interface [coyan-cli](https://github.com/NazaGara/Coyan/tree/main/coyan-cli).

## License

[MIT](https://choosealicense.com/licenses/mit/)
