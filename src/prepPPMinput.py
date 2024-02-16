import shutil
import os
import argparse

def separate_chain(s):
    if len(s) > 1:
        return ','.join(list(s))
    else:
        return s

def write_in_file(protein, destination_directory, ppm_directory, chains):
    ppminput = protein.replace(".pdb", ".in")
    destination_in_file_path = os.path.join(destination_directory, ppminput)
    chains_list = separate_chain (chains)
    
    with open(destination_in_file_path, "w") as f:
        f.write("2\n")
        f.write("no\n")
        f.write(f"{protein}\n")
        f.write("2\n")
        f.write("OPC\n")
        f.write("planar\n")
        f.write("in\n")
        f.write(f"{chains_list}\n")
    
    print(f".in file written to {destination_in_file_path}")
    
    source_res_lib = os.path.join(ppm_directory, "res.lib")
    destination_res_lib = os.path.join(destination_directory, "res.lib")
    
    shutil.copy(source_res_lib, destination_res_lib)
    
    print(f"res.lib copied to {destination_directory}")

def main():
    parser = argparse.ArgumentParser(description="Generate .in file and copy res.lib")
    parser.add_argument("-p", "--pdb", required=True, help="Protein name")
    parser.add_argument("-d", "--dest", required=True, help="Destination directory")
    parser.add_argument("-s", "--ppm", required=True, help="PPM directory")
    parser.add_argument("-c", "--chain", required=True, help="Chains of the protein (especially if the PDB contains multiple chains)")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.dest):
        os.makedirs(args.dest)
        
    write_in_file(args.pdb, args.dest, args.ppm, args.chain)

if __name__ == "__main__":
    main()
