# bioMAKEFP
bioMAKEFP is a python script designed to build GAMESS MAKEFP input files for amino acids, ligands, and water molecules existing in the solvation shell of interest, using structural data obtained from MD simulations. The script specifically requires `.g96` files representing both the whole structure and the solvation shell structure, as well as a `.itp` file including the atomic charges. For more clarity, see the examples in the `tests` directory. 

## Dependencies
- **Python 3.8 (Anaconda 2020.11)**
- **GAMESS**

## Usage

To run bioMAKEFP, use the following command:

```
python bioMAKEFP.py <input_file_1.g96> <input_file_2.g96> <input_file_2.itp>
```

## Obtaining the Solvation Shell File from GROMACS

To obtain the `.g96` file for the solvation shell from GROMACS, follow these steps:

1. **Create an index file for the solvation shell**:
   ```
   gmx select -f input_file_1.g96 -s input_file_1.g96 -on input_file_2 -select '`...`'
   ```

2. **Extract the solvation shell structure**:
   ```
   gmx editconf -f input_file_1.g96 -n input_file_2.ndx -o input_file_2.g96
   ```

## Author

**Andres S. Urbina**
asurbinab@gmail.com

For any inquiries or feedback, please contact Andres S. Urbina.