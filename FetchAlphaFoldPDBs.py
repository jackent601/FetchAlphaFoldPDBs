import sqlite3
import pandas as pd
import os
import shutil
import requests

def url_ok(url):
    r = requests.head(url)
    return r.status_code == 200

# Hard Coded For Local DataBase
ACCESSION_DB_PATH = ".\\AlphaFold\\accession_id_db.db"
ACCESSION_ID_TABLE_NAME = "accession_ids"
UNI_PROT_ID_FEATURE_NAME = "UniProtAccessionID"
ALPHAFOLD_DB_ID_FEATURE_NAME = "AlphaFoldDBID"

# ================================================================================================================================
#   RETRIEVING ALPHAFOLD DATABASE IDENTIFIER FROM UNIPROT IDS
# ================================================================================================================================

def getMatchingUniProtIDEntries(cursor, UniProtID, ACCESSION_ID_TABLE_NAME="accession_ids", UNI_PROT_ID_FEATURE_NAME="UniProtAccessionID"):
    """
    Cursor should be a connection to database with accession information, found on AlphaFolds FTP server.
    Uses this database to query and return AF info matching Uniprot ID
    """
    sqlQuery = f'SELECT * FROM {ACCESSION_ID_TABLE_NAME} WHERE {UNI_PROT_ID_FEATURE_NAME} = \'{UniProtID}\''
    return cursor.execute(sqlQuery).fetchall()

def getMatchingManyUniProtIDEntries(cursor, UniProtIDs, ACCESSION_ID_TABLE_NAME="accession_ids", UNI_PROT_ID_FEATURE_NAME="UniProtAccessionID"):
    """
    Cursor should be a connection to database with accession information, found on AlphaFolds FTP server.
    Uses this database to query and return AF info matching Uniprot ID list all at once avoiding slow multiple DB queries (as DB is huge)
    """
    # Generate SQL QuerySet
    INquery = "(\'" + "\', \'".join(UniProtIDs) + "\')"
    sqlQuery = f'SELECT * FROM {ACCESSION_ID_TABLE_NAME} WHERE {UNI_PROT_ID_FEATURE_NAME} IN {INquery}'
    return cursor.execute(sqlQuery).fetchall()

def getAFinfoForListOfUniProtIDs(cursor, UniProtIDs, debug=False, ACCESSION_ID_TABLE_NAME="accession_ids", UNI_PROT_ID_FEATURE_NAME="UniProtAccessionID"):
    """
    Uses getMatchingUniProtIDEntries to return AF info matching Uniprot ID
    Handles what is returned (i.e. if match if found or not) and saves info to a file
        writes to file on the go in case of any error in execution
    Returns info object also
    """
    # Catch hyphened IDs and reduce
    hyphenedIDs = [ID for ID in UniProtIDs if "-" in ID]
    hyphenedReducedIDs = [ID.split("-")[0] for ID in hyphenedIDs]
    extendedUniProtIDs = UniProtIDs + hyphenedReducedIDs

    # Get Query Set for these IDs
    AFentries = getMatchingManyUniProtIDEntries(cursor, extendedUniProtIDs, ACCESSION_ID_TABLE_NAME, UNI_PROT_ID_FEATURE_NAME)
    AF_uniprotIDs = [entry[0] for entry in AFentries]

    AF_information = []
    for ID in UniProtIDs:
        # Catch Reduced Case (search ID used to retain the original hyphenated ID)
        _searchID = ID
        if ID not in AF_uniprotIDs and "-" in ID:
            # Hyphenated ID, check reduced form ID
            IDreduced = ID.split("-")[0]
            if IDreduced in AF_uniprotIDs:
                _searchID = IDreduced
        
        if _searchID in AF_uniprotIDs:
            # Entry Found, Get it
            _index = AF_uniprotIDs.index(_searchID)
            _entry = AFentries[_index]
            # Record AF Info
            AF_info_string = utility_AppendEntryToAFinfoList(AF_information, ID, _entry) 
        else:
            # ID not found in AF DB
            AF_info_string = utility_AppendNullEntryToAFinfoList(AF_information, ID)
        if debug:
            print(AF_info_string.strip()) 

    return pd.DataFrame(AF_information)

def getAndSaveAFinfoForListOfUniProtIDs(cursor, UniProtIDs, outputPath, debug=False, ACCESSION_ID_TABLE_NAME="accession_ids", UNI_PROT_ID_FEATURE_NAME="UniProtAccessionID"):
    """Saves Dataframe result of getAFinfoForListOfUniProtIDs"""
    assert not os.path.exists(outputPath), "output path already exists!"
    AF_info_df = getAFinfoForListOfUniProtIDs(cursor, UniProtIDs, debug, ACCESSION_ID_TABLE_NAME, UNI_PROT_ID_FEATURE_NAME)
    AF_info_df.to_csv(outputPath, index=False)
    return AF_info_df

# ================================================================================================================================
#   RETRIEVING ALPHAFOLD PDBS FROM ALPHAFOLD IDENTIFIERS
# ================================================================================================================================
def fetchPDBsFromAlphaFoldInfoDataFrame(AlphaFoldInfoDataFrame, targetDirectory, localPDBdirectories=None, AlphaFoldFilesHTTPS = 'https://alphafold.ebi.ac.uk/files/', debug=True, outputPath=None):
    """
    Takes AlphaFold Information Dataframe (generated from getAFinfoForListOfUniProtIDs)
    For each entry:
        - checks if PDB present in target diectory i.e. already present, do nothing
        - checks if present in local directories of pdb (provided by localPDBdirectories) and copies pdb to target diectory
    If not found pulls a copy from AlphaFold into target directory
    """
    # Initialise DataFrame that links row to PDB path
    final_paths_to_pdb = []

    # Get pdbs present in target directory
    pdbs_present_in_target_dir = os.listdir(targetDirectory)

    # Get filepaths of all local pdbs (from direcotry list provided)
    if localPDBdirectories is not None:
        local_pdb_paths = [os.path.join(_dir,_name) for _dir in localPDBdirectories for _name in os.listdir(_dir)]
    else:
        local_pdb_paths = ['']
    
    # Loop each entry and conduct pdb checks
    for i, row in AlphaFoldInfoDataFrame.iterrows():
        # Get AlphaFold pdb name
        AF_ID = row['AF_DB_ID']
        uniprotID = row['uniprot_ID_source']

        # Catch where no ID match in AlphaFold
        if not isinstance(AF_ID, str):
            if debug:
                print(f'{uniprotID} has no AF match')
            final_paths_to_pdb.append(None)
            continue

        # Get PDB name
        version = str(int(row['latestVersion']))
        pdb_name = f"{AF_ID}-model_v{version}.pdb"

        # Check if already present in collected dir
        if pdb_name in pdbs_present_in_target_dir:
            if debug:
                print(f'{uniprotID} (AF: {pdb_name}) already present in {targetDirectory}, skipping')
            final_paths_to_pdb.append(os.path.join(targetDirectory, pdb_name))
            continue
        else:
            _TARGET_PATH = os.path.join(targetDirectory, pdb_name)

        # Check all local paths for pdb for match
        AF_FileMatches = [AF_file for AF_file in local_pdb_paths if pdb_name in AF_file]
        if len(AF_FileMatches) >= 1:
            # PDB found! Fetch first PDB match and copy to working space
            if debug:
                print(f'{uniprotID} (AF: {pdb_name}) found locally\n\tcopying {AF_FileMatches[0]} to {_TARGET_PATH}')
            shutil.copy(AF_FileMatches[0], _TARGET_PATH)
            final_paths_to_pdb.append(_TARGET_PATH)
        else:
            # PDB not found in files, pull from AlphaFold to working space
            if debug:
                print(f'{uniprotID} (AF: {pdb_name}) not found locally, pulling from AlphaFold')
            
            # Check URL okay
            pdb_url = f'{AlphaFoldFilesHTTPS}{pdb_name}'
            if url_ok(pdb_url):
                # Pull PDB
                os.system(f'curl {AlphaFoldFilesHTTPS}{pdb_name} -o {_TARGET_PATH}')
                final_paths_to_pdb.append(_TARGET_PATH)
            else:
                # URL not reachable
                final_paths_to_pdb.append(None)
                if debug:
                    print(f'{uniprotID} (AF: {pdb_name}) URL ({pdb_url}) not reachable')

    AlphaFoldInfoDataFrame['PDB_path'] = final_paths_to_pdb
    print('- - - FINISHED - - - ')
    if outputPath is not None:
        AlphaFoldInfoDataFrame.to_csv(outputPath, index=False)
    return AlphaFoldInfoDataFrame

# ================================================================================================================================
#   RUN FUNCTION (End to end)
# ================================================================================================================================

def Uniprot_Full(CONFIG, debug=False):
    """
    See readme - takes list of uniprot IDs, checks for AF entry in accession DB, pulls down PDB for entries present
    """
    # UnPack Config
    RUN_NAME = CONFIG['RUN_NAME']
    OUTPUT_DIR = CONFIG['OUTPUT_DIR']
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    PATH_TO_UNIPROTID_CSV = CONFIG['PATH_TO_UNIPROTID_CSV']
    ID_FEATURE_NAME = CONFIG['ID_FEATURE_NAME']
    LOCAL_ALPHAFOLD_PDB_DIRECTORIES = CONFIG['LOCAL_ALPHAFOLD_PDB_DIRECTORIES']
    
    ACCESSION_DB_PATH = CONFIG['ACCESSION_DB_PATH']
    ACCESSION_ID_TABLE_NAME = CONFIG['ACCESSION_ID_TABLE_NAME']
    TABLE_UNIPROT_ID_FEATURE_NAME= CONFIG['TABLE_UNIPROT_ID_FEATURE_NAME']
    ALPHAFOLD_ENDPOINT = CONFIG['ALPHAFOLD_ENDPOINT']
    
    # Get Unique IDs (shortcut)
    idList = pd.read_csv(PATH_TO_UNIPROTID_CSV)
    uniqueUniProtIDs = list(set(idList[ID_FEATURE_NAME]))
    print(f'Found {len(idList[ID_FEATURE_NAME])} IDs, {len(uniqueUniProtIDs)} unique IDs')
    
    # Set Up DB Connection
    conn = sqlite3.connect(ACCESSION_DB_PATH)
    cur = conn.cursor()
    
    # Determine which have AlphaFold Entries
    print("Checking which have AlphaFold entries")
    AF_INFO_OUTPUT_FILEPATH = os.path.join(OUTPUT_DIR, f'./{RUN_NAME}_AF_info.csv')
    AF_info = getAndSaveAFinfoForListOfUniProtIDs(
        cur, 
        uniqueUniProtIDs, 
        AF_INFO_OUTPUT_FILEPATH,
        ACCESSION_ID_TABLE_NAME=ACCESSION_ID_TABLE_NAME,
        UNI_PROT_ID_FEATURE_NAME=TABLE_UNIPROT_ID_FEATURE_NAME,
        debug=debug)
    
    # Small Debug Info
    match_flags = AF_info['uniprot_ID_match'].values
    numMatches = len([m for m in match_flags if isinstance(m, str)])
    numNoMatches = len([m for m in match_flags if not isinstance(m, str)])
    print(f'{len(match_flags)} IDs ran, {numMatches} matches, {numNoMatches} no match')
    
    # Prep Output Directory
    COLLECTED_PDBS_DIR = os.path.join(OUTPUT_DIR, f'{RUN_NAME}_AF_PDBs')
    os.makedirs(COLLECTED_PDBS_DIR, exist_ok=True)
    FINAL_OUTPUT_PATH = os.path.join(OUTPUT_DIR, f'./{RUN_NAME}_AF_info_with_PDB_paths.csv')
    
    # Fetch PDBs for those with entries
    print("Fetching PDBs for IDs with entries")
    AF_info_with_PDB_paths = fetchPDBsFromAlphaFoldInfoDataFrame(
        AF_info, 
        COLLECTED_PDBS_DIR, 
        LOCAL_ALPHAFOLD_PDB_DIRECTORIES, 
        outputPath=FINAL_OUTPUT_PATH,
        debug = debug,
        AlphaFoldFilesHTTPS=ALPHAFOLD_ENDPOINT)
    
    # Small Debug Info
    pdb_flags = AF_info_with_PDB_paths['PDB_path'].values
    numPDBsFound = len([m for m in pdb_flags if isinstance(m, str)])
    numNoPDBFound = len([m for m in pdb_flags if not isinstance(m, str)])
    print(f'{len(pdb_flags)} IDs ran, {numPDBsFound} PDBs found, {numNoPDBFound} no PDB found')

def Uniprot_Quick(CONFIG, assumedVersion=4, debug=False):
    """
    See readme, takes list of uniprot IDs, guesses AlphaFold PDB name, attempts to pull down PDB
    """
    # TODO - functionalise unpacking
    RUN_NAME = CONFIG['RUN_NAME']
    PATH_TO_UNIPROTID_CSV = CONFIG['PATH_TO_UNIPROTID_CSV']
    ID_FEATURE_NAME = CONFIG['ID_FEATURE_NAME']
    OUTPUT_DIR = CONFIG['OUTPUT_DIR']
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    LOCAL_ALPHAFOLD_PDB_DIRECTORIES = CONFIG['LOCAL_ALPHAFOLD_PDB_DIRECTORIES']
    ALPHAFOLD_ENDPOINT = CONFIG['ALPHAFOLD_ENDPOINT']

    # Read list of IDs (TODO - Functionalise)
    idList = pd.read_csv(PATH_TO_UNIPROTID_CSV)
    uniqueUniProtIDs = list(set(idList[ID_FEATURE_NAME]))
    print(f'Found {len(idList[ID_FEATURE_NAME])} IDs, {len(uniqueUniProtIDs)} unique IDs')

    # Prep Output Directory
    COLLECTED_PDBS_DIR = os.path.join(OUTPUT_DIR, f'{RUN_NAME}_AF_PDBs')
    os.makedirs(COLLECTED_PDBS_DIR, exist_ok=True)
    FINAL_OUTPUT_PATH = os.path.join(OUTPUT_DIR, f'./{RUN_NAME}_AF_info_with_PDB_paths.csv')

    # Assume PDB entry is of form "AF-{uniprotID}-F1-model_v{assumedVersion}.pdb" to construct AF info dataframe
    AF_DB_ID_guesses = [f'AF-{id}-F1' for id in uniqueUniProtIDs]
    _AlphaFoldInfoDataFrame = pd.DataFrame({'uniprot_ID_source':uniqueUniProtIDs, 'AF_DB_ID':AF_DB_ID_guesses, 'latestVersion':assumedVersion})
    AF_info_with_PDB_paths = fetchPDBsFromAlphaFoldInfoDataFrame(
        _AlphaFoldInfoDataFrame, 
        COLLECTED_PDBS_DIR, 
        localPDBdirectories=LOCAL_ALPHAFOLD_PDB_DIRECTORIES, 
        AlphaFoldFilesHTTPS = ALPHAFOLD_ENDPOINT, 
        debug=debug, 
        outputPath=FINAL_OUTPUT_PATH)

    # Rename AF_ID_DB column to include that it was a guess
    AF_info_with_PDB_paths = AF_info_with_PDB_paths.rename(columns={"AF_DB_ID":"AF_DB_ID_guess"})
    AF_info_with_PDB_paths.to_csv(FINAL_OUTPUT_PATH, index=False)
    
    # Small Debug Info (TODO - Functionalise)
    pdb_flags = AF_info_with_PDB_paths['PDB_path'].values
    numPDBsFound = len([m for m in pdb_flags if isinstance(m, str)])
    numNoPDBFound = len([m for m in pdb_flags if not isinstance(m, str)])
    print(f'{len(pdb_flags)} IDs ran, {numPDBsFound} PDBs found, {numNoPDBFound} no PDB found')



# ================================================================================================================================
#   UTILITIES (Small functions either to tidy code for readability, or of limited/specific use)
# ================================================================================================================================

def utility_fetchSingleAlphaFoldPDB(uniprotID, AlphaFoldID, AlphaFoldVersion, targetDirectory, localPDBdirectories=None):
    """
    Small utility for when a single PDB is needed fetched rather than a list (dataframe) if AF info is known
    """
    AlphaFoldInfoDataFrame = pd.DataFrame(
        {
            "uniprot_ID_source": uniprotID,
            "uniprot_ID_match": None,
            "AF_DB_ID": AlphaFoldID,
            "firstResidueIndex": None,
            "lastResidueIndex": None,
            "latestVersion": AlphaFoldVersion
        })
    fetchPDBsFromAlphaFoldInfoDataFrame(AlphaFoldInfoDataFrame, targetDirectory, localPDBdirectories, AlphaFoldFilesHTTPS = 'https://alphafold.ebi.ac.uk/files/')

def utility_AppendNullEntryToAFinfoList(AF_information, uniProtID):
    """Small utility to keep AF info scrapping code a bit cleaner"""
    # Add Null Entry to info list
    AF_information.append(
        {
        "uniprot_ID_source": uniProtID,
        "uniprot_ID_match": None,
        "AF_DB_ID": None,
        "firstResidueIndex": None,
        "lastResidueIndex": None,
        "latestVersion": None
        }
    )
    # Get write string
    return f'{uniProtID}, NA, NA, NA, NA, NA\n'

def utility_AppendEntryToAFinfoList(AF_information, uniProtID, entry):
    """Small utility to keep AF info scrapping code a bit cleaner"""
    # Unpack Entry
    uID, firstResidueIndx, lastResidueIndx, AF_DB_ID, latestVersion = entry
    # Add to info list
    AF_information.append(
        {
        "uniprot_ID_source": uniProtID,
        "uniprot_ID_match": uID,
        "AF_DB_ID": AF_DB_ID,
        "firstResidueIndex": firstResidueIndx,
        "lastResidueIndex": lastResidueIndx,
        "latestVersion": latestVersion
        }
    )
    # Get Write String
    return f'{uniProtID}, {uID}, {AF_DB_ID}, {firstResidueIndx}, {lastResidueIndx}, {latestVersion}\n'
    
# ================================================================================================================================
#   COMMAND LINE UTILITY
# ================================================================================================================================
# Not Needed Currently
# if __name__ == '__main__':
