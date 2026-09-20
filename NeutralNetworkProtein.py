import json
import sys
import requests
import time
import random
import csv 

import Bio
import Bio.PDB
import Bio.SeqRecord
from Bio.PDB import PDBParser, PDBIO
from pdbfixer import PDBFixer
from openmm.app import PDBFile

amino_acids = ["A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
aa_codes = {
        "A": "ALA",
        "R": "ARG",
        "N": "ASN",
        "D": "ASP",
        "C": "CYS",
        "Q": "GLN",
        "E": "GLU",
        "G": "GLY",
        "H": "HIS",
        "I": "ILE",
        "L": "LEU",
        "K": "LYS",
        "M": "MET",
        "F": "PHE",
        "P": "PRO",
        "S": "SER",
        "T": "THR",
        "W": "TRP",
        "Y": "TYR",
        "V": "VAL"
    }

reduced_amino_acids = ['D','G','L','A']

qmean_url = "https://swissmodel.expasy.org/qmean/submit/"
path_1LYZ = "./1LYZ.pdb"

#Command Line Argumesnts:
# 1. mutation_folder: The folder where the mutated PDB files will be saved. Default is "./mutations/"
# 2. original_pdb: The path to the original PDB file to be mutated. Default is "./1LYZ.pdb"

def main():

    # Check command line arguments
    mutation_folder = sys.argv[1] if len(sys.argv) > 1 else "./mutations/"
    original_pdb = sys.argv[2] if len(sys.argv) > 2 else path_1LYZ
    seq = getProteinSequenceFromFile(original_pdb)

    # Current Number of mutations and maximum number of mutation steps
    mutations = 0
    max_steps = 100

    cur_seq_path = original_pdb
    cur_pdb_mut_path = mutation_folder + str(mutations) + ".pdb"

    #get the initial QMean score of the original PDB file
    cur_high_qmean = getQMean(cur_seq_path)
    loops = 0
    results = []

    #While we have not reached the maximum number of mutations and the number of loops is less than 720, we will continue to mutate the protein sequence and evaluate the QMean score of the mutated structure.
    while mutations < max_steps and loops <= 720:
        #try block to catch any errors that may occur during the mutation process, such as network errors or issues with the PDB file.
        try:
            #Mutate our protein sequence by 1 amino acid and get the new sequence, position of mutation, new amino acid, and old amino acid.
            new_seq , position, new_aa, old_aa = mutate(seq)
            #Generate a new PDB file corrosponding to the mutated sequence
            mutate_pdb(cur_seq_path, cur_pdb_mut_path, position, new_aa, chain_id="A")
            #get the Qmean of the mutated PDB file
            cur_qmean = getQMean(cur_pdb_mut_path)
            #If the Qmean score is within .005 if the highest Qmean score we have got so far set it as a permenent step and move on
            if cur_qmean > (cur_high_qmean - .005):
                mutations += 1
                cur_seq_path = cur_pdb_mut_path
                cur_pdb_mut_path = mutation_folder + str(mutations) + ".pdb"
                seq = new_seq
                #add new sequence to results list with the Qmean score and number of loops
                results.append([mutation_folder + str(mutations) + ".pdb", cur_qmean, loops])

                #if the new Qmean score is higher than the highest Qmean score we have got so far, set it as the new highest Qmean score
            if cur_qmean > cur_high_qmean:
                cur_high_qmean = cur_qmean
            loops += 1
            time.sleep(10)

        #catch errors and print them out, then sleep for 10 seconds before continuing the loop
        except Exception as e:
            print(f"Error during mutation {mutations}: {e}")
            time.sleep(10)
            continue
    print(results)

    #Write results to a csv file
    with open(mutation_folder + 'results.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',',
                                quotechar='|', quoting=csv.QUOTE_MINIMAL)
        for row in results:
            writer.writerow(row)

    reduce_sequence(getProteinSequenceFromFile(original_pdb))

# Function to extract protein sequence from PDB file
# pdbPath: The path to the PDB file
# pdbcode: The PDB code of the structure (default is "1LYZ")
def getProteinSequenceFromFile(pdbPath: str,pdbcode = "1LYZ") -> str:
    pdbparser = Bio.PDB.PDBParser(QUIET=True)   # suppress PDBConstructionWarning
    struct = pdbparser.get_structure(pdbcode, pdbPath)

    
    ppb = Bio.PDB.PPBuilder()
    seqrecords = []
    for i, chain in enumerate(struct.get_chains()):
        # extract and store sequences as list of SeqRecord objects
        pps = ppb.build_peptides(chain)    # polypeptides
        seq = Bio.Seq.Seq("".join([str(i.get_sequence()) for i in pps]))
    return str(seq)

# Function to get QMean score from PDB file
# pdbPath: The path to the PDB file
def getQMean(pdbPath: str) -> str:
    poll_interval = 10
    timeout = 600
    #Queries the qmean api with our file

    success = False
    loop = 0

    #queries the qmean api with our file, if there is a network error it will retry up to 300 times with a 10 second delay between each retry
    while success == False and loop < 300:
        try:
            with open(pdbPath, 'rb') as f:
                response = requests.post(url=qmean_url,
                                        data={
                                            "email": "andrewsvester@gmail.com"
                                        },
                                        files={
                                            "structure": f
                                        })

            response.raise_for_status()
            success = True

        except requests.exceptions.RequestException as e:
            loop += 1
            print("Network error:", e)
            print("Retrying in", poll_interval, "seconds...")
            time.sleep(poll_interval)
            continue
        

    results_url = response.json()["results_json"]

    print(f"Submitted {pdbPath}")
    print(f"Checking: {results_url}")

    start_time = time.time()
    
    #loop to wait untill we get the response correctly
    while True:
        try:
            status_response = requests.get(
                results_url,
                timeout=30
            )
            status_response.raise_for_status()

            results = status_response.json()

        except requests.exceptions.RequestException as e:
            print("Network error:", e)
            print("Retrying in", poll_interval, "seconds...")
            time.sleep(poll_interval)
            continue

        status = results.get("status")

        print(f"Status: {status}")
        
        if status == "COMPLETED":
            qmean = results["models"]["model_001"]["scores"]["global_scores"]

            qmean_disco = qmean["avg_local_score"]

            print(f"QMEANDisCo: {qmean_disco:.4f}")

            return qmean_disco
        
        if status in ["FAILED", "ERROR"]:
            raise RuntimeError(f"QMEAN processing failed: {results}")
            
        if time.time() - start_time > 300:
            raise TimeoutError(
                f"QMEAN did not finish within {timeout} seconds"
            )
        
        time.sleep(1)

# Function to mutate a protein sequence by changing one amino acid to another
# seq: The protein sequence to be mutated
def mutate(seq):
    i = random.randrange(len(seq))
    a = random.randrange(20)
    while seq[i] == amino_acids[a]:
        a = random.randrange(20)
    
    old_aa = seq[i]
    new_aa = amino_acids[a]
    seq = list(seq)
    seq[i] = amino_acids[a]
    seq = "".join(seq)
    return seq, i, new_aa, old_aa

# Function to mutate a PDB file by changing one amino acid to another
# input_pdb: The path to the input PDB file
# output_pdb: The path to the output PDB file
# position: The position of the amino acid to be mutated (0-indexed)
# new_aa: The new amino acid to be introduced
# chain_id: The chain ID of the protein (default is "A")
def mutate_pdb(input_pdb, output_pdb, position, new_aa, chain_id="A"):
    

    new_resname = aa_codes[new_aa.upper()]

    # Read structure
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", input_pdb)

    model = structure[0]
    chain = model[chain_id]

    # Get standard amino-acid residues
    residues = [
        residue
        for residue in chain
        if residue.id[0] == " "
    ]

    if position < 0 or position >= len(residues):
        raise IndexError(
            f"Position {position} is outside chain {chain_id} "
            f"(length {len(residues)})"
        )

    residue = residues[position]

    old_resname = residue.resname

    print(
        f"Mutating chain {chain_id}, "
        f"position {position + 1}: "
        f"{old_resname} -> {new_resname}"
    )

    # Change residue identity
    residue.resname = new_resname

    # Save temporary PDB
    temp_pdb = output_pdb + ".temp.pdb"

    io = PDBIO()
    io.set_structure(structure)
    io.save(temp_pdb)

    # Use PDBFixer to add atoms missing from the new residue
    fixer = PDBFixer(filename=temp_pdb)

    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()

    fixer.findMissingAtoms()
    fixer.addMissingAtoms()

    with open(output_pdb, "w") as output_file:
        PDBFile.writeFile(
            fixer.topology,
            fixer.positions,
            output_file
        )

    return output_pdb

#Function to create as similar of a protein structure as possible with only the reduced amino acids D, G, L, and A. This is done by mutating all other amino acids to one of the reduced amino acids and then evaluating the QMean score of each mutated structure. The amino acid that results in the highest QMean score is then used to mutate the original structure.
def reduce_sequence(seq):
    for i in range(len(seq)):
        if seq[i] in reduced_amino_acids:
            continue
        seqa = list(seq)
        seqa[i] = 'A'
        seqa = "".join(seqa)
        mutate_pdb("./1LYZ_Results/1LYZ_Reduced_AA/1LYZ_reduced.pdb", f"./1LYZ_Results/1LYZ_Reduced_AA/1LYZ_reduced_a.pdb", i, 'A', chain_id="A")
        seqd = list(seq)
        seqd[i] = 'D'
        seqd = "".join(seqd)
        mutate_pdb("./1LYZ_Results/1LYZ_Reduced_AA/1LYZ_reduced.pdb", f"./1LYZ_Results/1LYZ_Reduced_AA/1LYZ_reduced_d.pdb", i, 'D', chain_id="A")
        seql = list(seq)
        seql[i] = 'L'
        seql = "".join(seql)
        mutate_pdb("./1LYZ_Results/1LYZ_Reduced_AA/1LYZ_reduced.pdb", f"./1LYZ_Results/1LYZ_Reduced_AA/1LYZ_reduced_l.pdb", i, 'L', chain_id="A")
        seqg = list(seq)
        seqg[i] = 'G'
        seqg = "".join(seqg)
        mutate_pdb("./1LYZ_Results/1LYZ_Reduced_AA/1LYZ_reduced.pdb", f"./1LYZ_Results/1LYZ_Reduced_AA/1LYZ_reduced_g.pdb", i, 'G', chain_id="A")

        qmeans = {}
        qmeans['A'] = getQMean(f"./1LYZ_Results/1LYZ_Reduced_AA/1LYZ_reduced_a.pdb")
        qmeans['D'] = getQMean(f"./1LYZ_Results/1LYZ_Reduced_AA/1LYZ_reduced_d.pdb")
        qmeans['L'] = getQMean(f"./1LYZ_Results/1LYZ_Reduced_AA/1LYZ_reduced_l.pdb")
        qmeans['G'] = getQMean(f"./1LYZ_Results/1LYZ_Reduced_AA/1LYZ_reduced_g.pdb")
        maxqmean = max(qmeans, key=qmeans.get)
        mutate_pdb("./1LYZ_Results/1LYZ_Reduced_AA/1LYZ_reduced.pdb", f"./1LYZ_Results/1LYZ_Reduced_AA/1LYZ_reduced.pdb", i, maxqmean, chain_id="A")





if __name__ == "__main__":
    main()