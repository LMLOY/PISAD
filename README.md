# Peptide Iterative design with Stochastic Algorithm
Peptide Iterative design with Stochastic Algorithm generates and optimizes peptides according to set parameters without structure guidance through three steps: initial sequence generation, sequence mutation, and crossover.

## Table of Contents

- [Overview](#overview)
- [Installation](#installation)
- [Usage](#usage)
  - [Initial Sequence Generation](#initial-sequence-generation)
  - [Mutation Process](#mutation-process)
  - [Crossover Process](#crossover-process)
- [Parameters](#parameters)
- [Contributing](#contributing)

## Overview

The suite includes scripts for generating initial sequences, performing mutations, and creating crossover sequences, all while assessing the sequences against various biological criteria. The primary goals are to predict secondary structures, evaluate stability, and ensure sequences meet specific hydrophilic/hydrophobic residue criteria.

## Installation

Before running the scripts, ensure you have the following prerequisites installed:

- Python 3
- PSIPRED for secondary structure prediction

Clone the repository and navigate to the project directory:

```bash
git clone https://github.com/LMLOY/PISAD.git
cd PISAD
```

## Usage

### Initial Sequence Generation

To generate single chain initial sequences:

```bash
./psipred.sh Initial False 1 50
```

This command generates 50 sequences starting from sequence number 1.

### Mutation Process

To generate mutated  sequences:

```bash
./psipred.sh Mutation False 1 6
```

This command performs 6 mutations per sequence in the specified directory.

### Crossover Process

To generate crossover sequences:

```bash
./psipred.sh Crossover False 1 10
```

This command generates 10 crossover sequences.

## Parameters

Parameters are defined in `singleParameter.py` and can be customized as needed to suit your specific requirements and environment setup:

- `home_dir`: The parent directory where all the codes are located.
- `OutSeq_dir`: The directory where output sequences are saved.
- `psipred_dir`: The directory where PSIPRED is installed.
- `Mut_dir`: The directory containing input sequences for mutation.
- `Cross_dir`: The directory containing input sequences for crossover.
- `seq_length`: The length of the protein sequence.
- `coilPercentage`: The maximum percentage of total coil in a sequence (for filtering based on secondary structure).
- `HHCriteria`: A setting to determine if the sequence should meet specific hydrophobic/hydrophilic residue criteria ("Yes" or "No").
- `minper_philic`: The minimum percentage of hydrophilic residue in the sequence.
- `maxper_philic`: The maximum percentage of hydrophilic residue in the sequence (optional, set to 0 if no maximum is desired).
- `maxper_phobic`: The maximum percentage of hydrophobic residue in the sequence.
- `SSCriteria`: A setting to determine if the sequence should meet specific secondary structure residue criteria ("Yes" or "No").
- `Perc_Helix`: The percentage increase in weight value for amino acids that commonly construct helix.
- `Perc_Beta`: The percentage increase in weight value for amino acids that commonly construct beta sheets.
- `Perc_Turn`: The percentage increase in weight value for amino acids that commonly construct turns.
- `aaWV`: A list of amino acid weight values for sequence generation, where each index corresponds to an amino acid in the order `["T", "S", "P", "G", "D", "K", "Q", "N", "A", "V", "E", "R", "I", "Y", "M", "F", "H", "C", "W", "L"]`.
- `Consecutive`: A parameter to determine the mutation strategy, "Yes" for sequential position mutation, "No" for random, and any other value for randomly choosing between them.
- `mutPosition`: A list of positions where mutations will be applied if `Consecutive` is set to "Yes" or a specific position.
- `mutWay`: The method of mutation, with different values corresponding to different mutation strategies.
- `mutRes`: The number of residues to be mutated in a sequence.

## Contributing

Contributions are welcome! Please fork the repository and submit a pull request with your improvements.
