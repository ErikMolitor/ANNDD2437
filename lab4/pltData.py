from util import *
import numpy as np 
import matplotlib.pyplot as plt

image_size = [28, 28]
train_imgs, train_lbls, test_imgs, test_lbls = read_mnist(
        dim=image_size, n_train=20, n_test=10000)

i = 7
# idx = 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16
# rsl = 5 0 4 1 9 2 1 3 1 4 3  5  3  6   1  7  2

lbl = np.where(train_lbls[i,:]==1)[0]
img = train_imgs[i,:]

img = np.reshape(img,[28,28])

plt.imshow(img)
plt.title(lbl)
plt.show()

print('DONE')
