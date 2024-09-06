# bioMAKEFP
bioMAKEFP is a python script designed to build GAMESS MAKEFP input files for amino acids, ligands, and water molecules within a specified solvation shell. It employs structural data obtained from MD simulations. The script requires `.g96` files representing for both the entire structure and the solvation shell, along with a `.itp` file that includes atomic charges. For detailed examples and usage, check the `tests` directory. 

## Dependencies
- **Python 3.8 (Anaconda 2020.11)**
- **GAMESS**

## Usage
To run the script, you must have the following files in the working directory:

1. `input_file_1.g96` 
    This structure file corresponds to a snapshot extracted from a curated MD trajectory.

2. `input_file_2.g96`
    This file corresponds to the structure of the solvation shell.  

3. `input_file_3.itp`
    This file contains the atomic charges of the protein and the ligands. The atom ids have to be present and follow the order in the `input_file_1.g96`.

4. `ligands`
    This file contains the residue names of the ligands as they appear in the `input_file_1.g96`. Each ligand should be listed on a new line.

    **Example:**
    ```
    LIG
    ```

5. `taas`
    This file lists the residue names of the terminal amino acids as they appear in the `input_file_1.g96`. Each terminal amino acid should be listed on a new line.

    **Example:**
    ```
    NMET
    CLEU
    ```

6. `settings`
    This file is used to configure key settings. It allows you to enable or disable certain parts of the workflow.

    **Example:**
    ```
    [Settings]
    ligands = yes
    sf = yes
    ```

To run `bioMAKEFP.py`, use the following command:

```
python bioMAKEFP.py <input_file_1.g96> <input_file_2.g96> <input_file_3.itp>
```

## Obtaining the solvation shell using GROMACS

To obtain the `.g96` file for the solvation shell, follow these steps:

1. **Create an index file for the solvation shell**:
   ```
   gmx select -f input_file_1.g96 -s input_file_1.g96 -on input_file_2 -select '`...`'
   ```

2. **Extract the solvation shell structure**:
   ```
   gmx editconf -f input_file_1.g96 -n input_file_2.ndx -o input_file_2.g96
   ```

## Limitations
1. Terminal amino acids present in the solvation shell are not used to generate their associated GAMESS MAKEFP input files.

2. The script currently does not handle sulfur bridges if they are present in the solvation shell.


## Author

**Andres S. Urbina**
asurbinab@gmail.com

For any inquiries or feedback, please contact Andres S. Urbina.