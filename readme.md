# sobel-filter-MPI-implementation

This github repo contains a MPI implemention of the edge-detector [Sobel filter](https://en.wikipedia.org/wiki/Sobel_operator), which was made originally to run on [Dardel](https://www.pdc.kth.se/hpc-services/computing-systems/about-the-dardel-hpc-system-1.1053338), the fastest computer in Sweden, at the time of writing this.

## Convert `.jpg` to `.csv`

The original code was made to run on `.csv` file, representing a grayscale image. Thus to filter an image, we provide python script [jpeg-to-matrix.py](https://github.com/antongorillam/sobel-filter-MPI-implementation/blob/master/jpeg-to-matrix.py) which converts all `.jpg` in a given folder to grayscale `.csv`. To do so, run:

``` bash
mkdir -p images/custom # Creates folder if it doesn't exist
python3 jpeg-to-matrix --convert_dir images/custom --save_dir .
```

This will convert all the files in the folder `images/custom` to the right format.

## Filtering

To run sobel filtering, we will provide instruction how to run it in Dardel, but simple MPI should also work. In Dardel, simply run:

``` bash
srun sobel_filter_mpi.x images/mat_files/8049.csv -05
```

This will filter the image and output it (assuming the input image is called `8049.csv`) to `images/mat_filtered/filtered_8049.csv`.

## Converting filtered `.csv` to viewable image

This python code would do the trick:

``` python
import pandas as pd
from PIL import Image

image_filtered = pd.read_csv('images/mat_filtered/filter_8049', header=None).values

normalized_data = ((image_filtered - image_filtered.min()) * (1/(image_filtered.max() - image_filtered.min()) * 255)).astype('uint8')

# Create an image from the numpy array
img = Image.fromarray(normalized_data)

img.save("output.png")

```
