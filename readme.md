## Checking and Fetching AlphaFold PDBs

Functions to find check if AlphaFold entries exist for uniprot IDs and retrieve PDBs for IDs with entries.

Three execution options:

- *Uniprot_Full*: Checks accession ID database for an exhaustive search of AlphaFold entries relating to uniprot IDs, for those with matches retrieves the AlphaFold PDBs where possible
  - Hence requires local database AlphaFold accession IDs
  - Also catches where uniprot ID has a '-' present and will search a 'reduced' ID without hyphen, e.g. PXXXXX-Y will attempt to match to PXXXXX, this is captured in output (and ultimately sequences are checked anyway).
- *Uniprot_Guess*: Takes uniprot IDs, *guesses* the name of AlphaFold PDB and attempts to retrieve it
  - Convenient if database check if not necessary
  - Guess is of the form 'AF-{uniprotID}-F1-model_v{latestVersion}.pdb'

- *AlphaFoldPDBs*: Takes a list of AlphaFold PDB names and attempts to retrieve them
  - Useful when AlphaFold PDB name is known with confidence


All functions capture information on the success of ID matching and PDB fetching and saves to a new csv



### WalkthroughDemoAndTest

Explains functions within FetchAlphaFoldPDBs, and tests the demo case.

some functions require local database AlphaFold accession IDs (csv available from [AlphaFolds FTP server](http://ftp.ebi.ac.uk/pub/databases/alphafold/)).



### Run

Notebook to execute each of the three main functions.
