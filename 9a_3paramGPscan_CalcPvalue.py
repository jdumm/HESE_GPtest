import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii

# After having run scrambled trials saved in .out/.txt files, calculate the p-value.

#data = ascii.read('all_4year.out')
data = ascii.read('all_6year.out')
p = data['col1']

# From 4year, 1856 trials so far
#myp = 0.000101817314443 # 4-year result
myp = 0.000578051702569 # 6-year result

len(p[p<=myp])/float(len(p))
# Trial-corrected pvalue for 4-year result:
# 0.034482758620689655
# Trial-corrected pvalue for 6-year result:
#0.1557112068965517

plt.hist(-1*np.log10(p),bins=20)
plt.semilogy()
plt.show()

