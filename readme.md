## Checking and Fetching AlphaFold PDBs

Functions to find check if AlphaFold entries exist for uniprot IDs and retrieve PDBs for IDs with entries.

Three execution options:

- *Uniprot_Full*: Checks accession ID database for an exhaustive search of AlphaFold entries relating to uniprot IDs, for those with matches retrieves the AlphaFold PDBs
  - Hence requires local database AlphaFold accession IDs
- *Uniprot_Quick*: Takes uniprot IDs, *guesses* the name of AlphaFold PDB and attempts to retrieve it
  - Convenient if database check if not necessary

- *AlphaFoldPDBs*: Takes a list of AlphaFold PDB names and attempts to retrieve them
  - Useful when AlphaFold PDB name is known with confidence


All functions capture information on the success of ID matching/PDB fetching



### WalkthroughDemoAndTest

Explains functions within FetchAlphaFoldPDBs, and tests the demo case.

some functions require local database AlphaFold accession IDs (csv available from [AlphaFolds FTP server](http://ftp.ebi.ac.uk/pub/databases/alphafold/)).



### Run

Notebook to execute each of the three main functions.
