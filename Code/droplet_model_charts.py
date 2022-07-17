import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
#plt.rcParams['text.usetex'] = True

# Read data from files
# data1 = pd.read_csv('C:/Users/Mateusz/Downloads/Studnia/DATA_BEC_3D_N1000.txt', sep='\t', header=None)
# data2 = pd.read_csv('C:/Users/Mateusz/Downloads/Studnia/DATA_BEC_3D_N1000.txt', sep='\t', header=None)
# data3 = pd.read_csv('C:/Users/Mateusz/Downloads/Studnia/DATA_BEC_3D_N1000.txt', sep='\t', header=None)
# data4 = pd.read_csv('C:/Users/Mateusz/Downloads/Studnia/DATA_BEC_3D_N1000.txt', sep='\t', header=None)
#data1 = pd.read_csv('../Data/gdd1000N200_well.txt', sep='\t', header=None)
data2 = pd.read_csv('../Data/gdd25N50_new_well.txt', sep='\t', header=None)
data3 = pd.read_csv('../Data/gdd25N50_old_well.txt', sep='\t', header=None)
#data4 = pd.read_csv('../Data/gdd1000N200_msc_new.txt', sep='\t', header=None)

# Set plot
plt.figure(figsize=(8, 5.5), dpi = 100)
plt.subplots_adjust(bottom=0.17, top=0.95)

# Set seaborn color palette
sns.set_palette("mako_r")

# Plot data
#sns.lineplot(data=data1, x=data1[0], y=data1[1], label=r'$g_{\rm dd} = 60$', linewidth=2)
sns.lineplot(data=data2, x=data2[0], y=data2[1], label=r'$g_{\rm dd} = 1000\;N=200\;Studnia$', linewidth=2)
sns.lineplot(data=data3, x=data3[0], y=data3[1], label=r'$g_{\rm dd} = 1000\;N=200\;Maciek$', linewidth=2)
#sns.lineplot(data=data4, x=data4[0], y=data4[1], label=r'$g_{\rm dd} = 120$', linewidth=2)
# Set grid on the plot
plt.grid(color = 'gray', linestyle = '-.', linewidth = 0.5)
# Set x label description
plt.xlabel(r'Temperatura $\left[\frac{\hbar^2}{k_BmL^2}\right]$', fontsize = 16)
# Set y label description
plt.ylabel(r'Szerokość kropli', fontsize = 16)
#plt.yticks(fontname = "Times New Roman")
#plt.xticks(fontname = "Times New Roman")
plt.tick_params(axis = 'both', which = 'major', labelsize = 16)
plt.tick_params(which = 'both', direction='in', top=True, right=True, length = 5, width = 1)
# plt.tick_params(which = 'minor', direction='in', top=True, right=True, length = 3, width = 1)
# Legend
plt.legend(loc = "upper right", fontsize = 16)
# Set fontsize of ticks
# Save and show
plt.savefig('various_N.png', transparent=True)
plt.show()
