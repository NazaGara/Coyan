# Coyan: Random Fault Tree Generator

Generator of static Random Fault Trees.

## Build

Create the binary by running:
```bash
cargo build --release -p coyan_rft
```
This creates the binary `/target/release/coyan_rft`.

## Usage

```bash
./target/release/coyan_rft --n-nodes <N_NODES> --output <OUTPUT> [OPTIONS]
```

Options:
  - `-n, --n-nodes <N_NODES>`: Total number of nodes.
  - `-o, --output <OUTPUT>`: Output file, writes .dft files
  - `--rate-be <RATE_BE:>`:  Number in `(0,1]`. Rate of Basic Events from the total amount of Nodes. [default: 0.3]
  - `--rate-and <RATE_AND>`: Number in `(0,1]`. Rate of AND gates from the total amount of Nodes. [default: 0.3]
  - `--rate-or <RATE_OR>`:   Number in `(0,1]`. Rate of OR gates from the total amount of Nodes. [default: 0.3]
  - `--rate-vot <RATE_OR>`:   Number in `(0,1]`. Rate of VOT gates from the total amount of Nodes. [default: 0.3]
  - `--prob-multiplier <PROB_MULTIPLIER>`: Number in `(0,1]`, multiplies the random float number of the probability for Basic Events, so it can be a smaller number.[default: 0.0001]
  - `--perc-last <PERC_LAST>`: In which percentage of the last gates start to put the Basic Events if they were not used before [default: 0.6]
  - `-h, --help`: Print help
  - `-V, --version`: Print version

