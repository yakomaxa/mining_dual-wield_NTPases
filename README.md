# mining dual-wield NTPase

README for structure-mining against AlphaFold protein structure database using foldcomp to find dual-wield NTPase family.

## contents

find_ploop: contatins the script to find structures containing Walker-A-like sub-structure from foldcomp database. This dumps all the fcz file for the structures contatining Walker-A like sub-structures.

dbmaker: contains shell script to call MMseqs2 to make custom foldcomp database for Walker-A-like-containing structures.

sheet_analysis: contatins the script to analyze the beta-sheet architecture against foldcomp database. In actual analysis, we used custom subset of foldcomp database just owning the structure containing Walker-A-like substructure made by the script in the dbmaker directory. This generates the result of sheet-analysis in csv format.

find_dual: contains an R-script to select structures. We added test data TestData_only_10000line.csv from which you can extract a dwNTPase among example 10000 lines.

md_script: script templates for running MD using Gromacs. Force-field file is excluded.

docked_structure: manually docked structure of representative dwNTPase alphafold prediction and putative ligands 

find_conserved: script to find conserved residue from MSA fasta

The python script in find_ploop (and sheet_analysis as well) should be executed parallely in you need speed. You need to divide index file of foldcomp into chunks.

## license

MIT license.

