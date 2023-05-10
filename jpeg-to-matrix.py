from pathlib import Path
from PIL import Image 
import numpy as np
import os
import argparse

"""
example usage of code:
python3 jpeg-to-matrix --convert_dir images/custom --save_dir .
"""


"""Convert JPEG images to grayscale matrices and save them in CSV format.

Usage:
python3 jpeg-to-matrix --convert_dir images/custom --save_dir .

Args:
--convert_dir (str): The directory to convert images from.
--save_dir (str): The directory to save the converted image.

Returns:
None
"""
def main():
    """Main function that converts JPEG images to grayscale matrices and saves them in CSV format."""
    # Parse command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--convert_dir', type=str, help='The directory to convert images from.')
    parser.add_argument('--save_dir', type=str, help='directory to save the converted image.')
    args = parser.parse_args()

    # Convert images to matrices and save them in CSV format
    source_dir = Path(args.convert_dir)
    files = source_dir.iterdir()
    for file in files:
        img = Image.open(file).convert('L')  # Convert image to grayscale
        img = np.array(img)  # Convert image to numpy array
        file_name = os.path.basename(file)  # Get file name
        file_name, file_extension = os.path.splitext(file_name)  # Get file name and extension
        np.savetxt(f'{args.save_dir}/{file_name}.csv', img, delimiter=",")  # Save matrix in CSV format


if __name__ == "__main__":
    main()