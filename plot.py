__author__ = 'kcho'
import matplotlib.pyplot as plt
import numpy as np


#plt.plot([1,2,3,4])
#plt.plot([1,2,3,4], [1,4,9,16])
# plt.plot([1,2,3,4], [1,4,9,16], 'ro')
# plt.axis([0, 6, 0, 20])
# plt.ylabel('some numbers')
# plt.show()




t = np.arange(0., 5., 0.2)
plt.plot(t, t, 'r--', t, t**2, 'bs', t, t**3, 'g^')
plt.axis([0, 20, 0, 20])

plt.show()

import matplotlib.pyplot as plt

t11 = ['00', '01', '02', '03', '04', '05', '10', '11', '12', '13', '14', '15',
       '20', '21', '22', '23', '24', '25', '30', '31', '32', '33', '34', '35',
       '40', '41', '42', '43', '44', '45', '50', '51', '52', '53', '54', '55']

t12 = [173, 135, 141, 148, 140, 149, 152, 178, 135, 96, 109, 164, 137, 152,
       172, 149, 93, 78, 116, 81, 149, 202, 172, 99, 134, 85, 104, 172, 177,
       150, 130, 131, 111, 99, 143, 194]


plt.bar(range(len(t12)), t12, align='center')
plt.xticks(range(len(t12)), t11, size='small')
plt.show()