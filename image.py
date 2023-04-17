import numpy as np
from PIL import Image
from matplotlib import pyplot as plt

def arrays_equal(a, b):
    if a.shape != b.shape:
        return False
    for ai, bi in zip(a.flat, b.flat):
        if ai != bi:
            return False
    return True


csv1 = np.genfromtxt('images/mat_files/12074.jpg.csv', delimiter=",")
csv2 = np.load('images/mat_files/12074.jpg.npy')
# csv = csv.astype(int)
print(arrays_equal(csv1, csv2))
# img = Image.fromarray(csv1, 'L')
# img.show()

plt.imshow(csv1, cmap='gray')
plt.show()