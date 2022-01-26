import pandas as pd
import matplotlib.pyplot as plt


print(plt.style.available)
data1 = pd.read_csv('gdd25N50.txt', sep='\t', header=None)
data2 = pd.read_csv('gdd25N75.txt', sep='\t', header=None)
data3 = pd.read_csv('gdd25N100.txt', sep='\t', header=None)
#data4 = pd.read_csv('gdd25N125.txt', sep='\t', header=None)

fig, ax = plt.subplots(figsize=(8, 5.5), dpi=100)
fig.subplots_adjust(bottom=0.17, top=0.95)
plt.style.use('tableau-colorblind10')
plt.plot(data1[0], data1[1], '-', label=r'$N_{\rm tot} = 50$', color='black', linewidth=1.5, alpha=0.65)
plt.plot(data2[0], data2[1], '-', label=r'$N_{\rm tot} = 75$', linewidth=1.5, alpha=0.65)
plt.plot(data3[0], data3[1], '-', label=r'$N_{\rm tot} = 100$', linewidth=1.5, alpha=0.65)
#plt.plot(data4[0], data4[1], '-', label=r'$N_{\rm tot} = 125$', linewidth=1.5, alpha=0.65)
plt.xlabel(r'$T\;\left(\frac{\hbar^2}{mL^2 k_B}\right)$', fontsize = 16)
plt.ylabel(r'$W$', fontsize = 16)
plt.legend(loc="upper right", fontsize = 16)
plt.grid(color='gray', linestyle='-.', linewidth=0.5)
plt.tick_params(axis='both', which='major', labelsize=16)
plt.savefig('various_N.png', transparent=True)

# plt.style.use('ggplot')
# fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, figsize=(12, 4), dpi=100)
# fig.subplots_adjust(left=0.06, bottom=0.125, right=0.94, wspace=0.3)

# ax1.plot(data1[0], data1[1], '-',
#          label='N = 50', color='black', linewidth=1.25, alpha=0.65)
# ax1.plot(data2[0], data2[1], '-',
#          label='N = 75', color='orange', linewidth=1.25, alpha=0.65)
# ax1.plot(data3[0], data3[1], '-',
#          label='N = 100', color='blue', linewidth=1.25, alpha=0.65)
# ax1.plot(data4[0], data4[1], '-',
#          label='N = 125', color='green', linewidth=1.25, alpha=0.65)

# ax1.set_title('Width', fontname='Times New Roman')
# ax1.set_ylabel('a', fontname='Times New Roman')
# ax1.set_xlabel('T', fontname='Times New Roman')
# ax1.legend(loc="upper right")


# ax2.plot(data1[0], data1[2], '-',
#          label='N = 50', color='black', linewidth=1.25, alpha=0.65)
# ax2.plot(data2[0], data2[2], '-',
#          label='N = 75', color='orange', linewidth=1.25, alpha=0.65)
# ax2.plot(data3[0], data3[2], '-',
#          label='N = 100', color='blue', linewidth=1.25, alpha=0.65)
# ax2.plot(data4[0], data4[2], '-',
#          label='N = 125', color='green', linewidth=1.25, alpha=0.65)

# ax2.set_title('Width Fluctuations',
#               fontname='Times New Roman')
# ax2.set_ylabel(r'$\sigma^2_{\leftangle a \rightangle}$',
#                fontname='Times New Roman')
# ax2.set_xlabel('T', fontname='Times New Roman')
# # ax2.legend(loc="upper right")

# ax3.plot(data1[0], data1[3], '-',
#          label='N = 50', color='black', linewidth=1.25, alpha=0.65)
# ax3.plot(data2[0], data2[3], '-',
#          label='N = 75', color='orange', linewidth=1.25, alpha=0.65)
# ax3.plot(data3[0], data3[3], '-',
#          label='N = 100', color='blue', linewidth=1.25, alpha=0.65)
# ax3.plot(data4[0], data4[3], '-',
#          label='N = 125', color='green', linewidth=1.25, alpha=0.65)

# ax3.set_title('Helmholtz Free Energy', fontname='Times New Roman')
# ax3.set_ylabel('F(T)', fontname='Times New Roman')
# ax3.set_xlabel('T', fontname='Times New Roman')
# # ax3.legend(loc="upper right")

plt.show()
