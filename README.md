# SecondaryStructurePrediction-ChouFasmanAlgorithm
Prediction of secondary structure of a protein, whose primary structure is given as input, by using Chou-Fasman algorithm.

# Requirements
- Python>=3.6
- colorama module: [https://pypi.org/project/colorama/](https://pypi.org/project/colorama/)

# Usage
- python3 chou_fasman.py <sequence_path> <secondary_structure_path>
- An example for protein sequence file is tp53_protein_sequence.txt.
- An example for secondary structure file, which indicates the known secondary structure elements of the protein along with at which position intervals they occur, is tp53_secondary_structure.txt.
- If secondary structure of the protein is specified, the program shows 3x3 and individual confusion matrix computations at the end of the output.
- Specifying secondary structure of the protein is optional. If it is not specified, confusion matrix computations are not shown. So, the program can also be executed as follows:
    > python3 chou_fasman.py <sequence_path>
- When **&#45;&#45;verbose** is added at the end of the command, the program generates a detailed output regarding the steps of the algorithm:
    > python3 chou_fasmen.py <sequence_path> &#45;&#45;verbose
    > python3 chou_fasman.py <sequence_path> <secondary_structure_path> &#45;&#45;verbose
    
# Record
[![asciicast](https://asciinema.org/a/4otZeW06bJi6uapUevriaquRw.svg)](https://asciinema.org/a/4otZeW06bJi6uapUevriaquRw)