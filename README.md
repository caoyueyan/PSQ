# DISP Generator
A Python tool for generating dispensation orders for pyrosequencing and evaluating whether they meet phase conditions.
## File Description
**DISP_generator.py:** Main script for generating sequencing signals, filtering complete signal sequences, and performing phase testing.
**DISP_get.py:** Helper script for randomly generating dispensation orders of specified length and quantity.
## Overview
**Generate Dispensation Orders:** Creates multiple base dispensation sequences based on random rules.
**Simulate Sequencing Signals:** Simulates pyrosequencing signals based on input microhaplotype (MH) sequences.
**Filter Complete Signal Sequences:** Excludes dispensation orders that cannot fully cover the target sequences.
**Phase Testing:** Tests whether each dispensation order can uniquely distinguish different allele combinations to ensure they are phaseable.
## Usage
### Requirements
  Python 3.6+
### Required libraries: 
  pandas, numpy
### Input File
  sequence_list.csv: Contains the list of microhaplotype sequences to be sequenced, in CSV format.
### Run Main Program
bash
python DISP_generator.py
### Parameter Configuration (adjustable in script)
length: Length of dispensation order (default: 50)
Ntest: Number of dispensation orders to generate (default: 1000)
file: Path to input sequence file
### Output Description
Console output:
● List of dispensation orders with complete signals
● Phase test results for each dispensation order
● Number of finally available dispensation orders
### Example Output
‘<text
Complete signal dispensation orders: [['A' 'C' 'G' ... 'T' 'A' 'C']...]
Available dispensation orders: X!
Running time: X.XX Seconds>'
## Notes
Ensure the input sequence file is correctly formatted and the path is accurate.
Adjust length and Ntest to control the complexity and quantity of generated orders.
Current phase testing logic is based on the uniqueness judgment of signal sums.
## License
This project is licensed under the MIT License.
