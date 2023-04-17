from pathlib import Path
import sys
from PIL import Image 
import numpy as np
import os


def main():
    source_dir = Path('images/train/')
    files = source_dir.iterdir()
    for file in files:
        img = Image.open(file).convert('L')
        img = np.array(img)
        
        np.savetxt(f'images/mat_files/{os.path.basename(file)}.csv', img, delimiter=",")
        np.save(f'images/mat_files/{os.path.basename(file)}', img)

if __name__ == "__main__":
    main()