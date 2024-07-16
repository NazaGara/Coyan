# Coyan

Coyan is a Rust project that transforms Static Fault Trees into CNF equisatisfiables formulas. Then, with the use of a Model Counter,
computes the Top Event Probability of failure at a given time bound $t$.
Also, it allows the creation of random Fault Trees with specific configurations.

## Installation

```bash
git clone https://github.com/NazaGara/ft_project
cargo build --release 
```

## Usage

The tool currently posses the `Fault Tree to CNF` [coyan_fta](https://github.com/NazaGara/Coyan/tree/main/coyan_fta) module and the `Random Fault Tree Generator` [coyan_rft](https://github.com/NazaGara/Coyan/tree/main/coyan_rft) module.

Each one have their own usage detailed on the correspondent subdirectory.

## License

[MIT](https://choosealicense.com/licenses/mit/)