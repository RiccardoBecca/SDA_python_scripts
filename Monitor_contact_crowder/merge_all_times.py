import numpy as np
import argparse
import os
from tqdm import tqdm

def main(args):

    data_merged=np.array([])
    for file in tqdm(os.listdir(args.folder)):
        if args.files_name in file:
            data=np.load(args.folder+"/"+file)
            data_merged=np.hstack((data_merged, data))

    np.save(args.folder+"/"+args.files_name+"_merged", data_merged)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description = "Python script to merge npy bounded_times files",
                                     epilog="Example usage:\n"
                                     "python merge_all_times.py --files_name bounded_times --folder crowder_contacts")
    parser.add_argument("--files_name", help = "string prefix for files to be merged", type=str, required = True)
    parser.add_argument("--folder", help = "folder where all files are saved", type=str, required = True)

    args = parser.parse_args()

    main(args)
