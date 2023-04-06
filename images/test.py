import scipy
import sklearn
from sklearn.feature_extraction import image
from scipy.io import loadmat
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
images = loadmat('IMAGES.mat',variable_names='IMAGES',appendmat=True).get('IMAGES')

