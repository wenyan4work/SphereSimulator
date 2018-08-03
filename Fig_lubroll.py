import matplotlib as mpl
import os

mpl.rcParams[r'font.family'] = r'serif'

import numpy as np
import matplotlib.pyplot as plt

# rc params
fig_width_pt = 384  # Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 1.0 / 72.27  # Convert pt to inch
aspectratio = 0.5  # Aesthetic ratio
fig_width = fig_width_pt * inches_per_pt  # width in inches
fig_height = fig_width * aspectratio  # height in inches
fig_size = [fig_width, fig_height]
# define rc params
params = {r'backend': r'Agg',
          r'axes.labelsize': 9,
          r'font.size': 9,
          r'legend.fontsize': 9,
          r'xtick.labelsize': 9,
          r'ytick.labelsize': 9,
          r'text.usetex': True,
          r'figure.figsize': fig_size,
          r'text.latex.preamble': r'\usepackage{amsmath},\usepackage[notextcomp]{stix},'
                                  + r'\usepackage[T1]{fontenc}, \usepackage{bm}'}
plt.rcParams.update(params)

# nd solve data A
# alpha=delta^2/ell^2
transfactor = (6 * np.pi * 0.001 * 1)
rotfactor = (8 * np.pi * 0.001 * 1)


# load data from text file
def getData(filename):
    v0file = 'tempv0' + filename
    v1file = 'tempv1' + filename
    os.system('grep "id 0" ./' + filename + ' > ' + v0file)
    os.system('grep "id 1" ./' + filename + ' > ' + v1file)

    os.system("tr ' ' ',' <"+v0file+" > " + "new"+v0file)
    os.system("tr ' ' ',' <"+v1file+" > " + "new"+v1file)

    v0 = np.loadtxt("new"+v0file, delimiter=",", usecols=[1, 13])
    v1 = np.loadtxt("new"+v1file, delimiter=",", usecols=[1, 13])

    os.system("rm ./"+v0file)
    os.system("rm ./"+v1file)
    os.system("rm ./new"+v0file)
    os.system("rm ./new"+v1file)

    return v0, v1


v0p6, v1p6 = getData("lubrollp6.txt")
v0p8, v1p8 = getData("lubrollp8.txt")
v0p10, v1p10 = getData("lubrollp10.txt")
v0p12, v1p12 = getData("lubrollp12.txt")
v0p16, v1p16 = getData("lubrollp16.txt")
v0p24, v1p24 = getData("lubrollp24.txt")

print(v0p24)

v0true = v0p24

fig = plt.figure()

ax1 = plt.subplot(1, 2, 1)

Lv0p6, = ax1.plot(v0p6[:, 0], np.abs(v0p6[:, 1]-v0true[:, 1])*rotfactor, '-',
                  markeredgecolor='#FF0000', markerfacecolor='none', color='#FF0000')
Lv0p8, = ax1.plot(v0p8[:, 0], np.abs(v0p8[:, 1]-v0true[:, 1])*rotfactor, '-',
                  markeredgecolor='#00DD00', markerfacecolor='none', color='#00DD00')
Lv0p10, = ax1.plot(v0p10[:, 0], np.abs(v0p10[:, 1]-v0true[:, 1])*rotfactor, '-',
                   markeredgecolor='#0000FF', markerfacecolor='none', color='#0000FF')
Lv0p12, = ax1.plot(v0p12[:, 0], np.abs(v0p12[:, 1]-v0true[:, 1])*rotfactor, '-',
                   markeredgecolor='#00DDDD', markerfacecolor='none', color='#00DDDD')
Lv0p16, = ax1.plot(v0p16[:, 0], np.abs(v0p16[:, 1]-v0true[:, 1])*rotfactor, '-',
                   markeredgecolor='#FF00DD', markerfacecolor='none', color='#FF00DD')
Lv0p24, = ax1.plot(v0p24[:, 0], np.abs(v0p24[:, 1]-v0true[:, 1])*rotfactor, '-',
                   markeredgecolor='#000000', markerfacecolor='none', color='#000000')


ax1.set_ylim(1e-14, 1e-2)
ax1.set_yscale('log')
ax1.set_ylabel(r'$\lvert \omega_{y,error} \rvert/(T/8\pi\eta a)$')

ax1.set_xlim(0.01, 5.0)
ax1.set_xscale('log')
ax1.set_xlabel(r'$\epsilon/a$')

# formater = mpl.ticker.FormatStrFormatter('%2.1f')
# ax1.xaxis.set_major_locator(plt.FixedLocator([0.1,0.2,0.4,0.8,1.0]))
# ax1.xaxis.set_major_formatter(formater)
# formatery = mpl.ticker.FormatStrFormatter('%3.2f')
# ax1.yaxis.set_major_locator(plt.FixedLocator([0.02,0.05,0.1,0.2,0.5,1.0]))
# ax1.yaxis.set_major_formatter(formatery)

ax1.legend((Lv0p6, Lv0p8, Lv0p10, Lv0p12, Lv0p16, Lv0p24),
           (r'$p=6$',
            r'$p=8$',
            r'$p=10$',
            r'$p=12$',
            r'$p=16$',
            r'$p=24$',
            ), frameon=True, loc=r'upper right', fontsize=7, numpoints=1, scatterpoints=1, ncol=1)

ax2 = plt.subplot(1, 2, 2)
Lve01p6, = ax2.plot(v0p6[:, 0], np.abs(v0p6[:, 1]-v1p6[:, 1])*rotfactor, '--',
                    markeredgecolor='#FF0000', markerfacecolor='none', color='#FF0000')
Lve01p8, = ax2.plot(v0p8[:, 0], np.abs(v0p8[:, 1]-v1p8[:, 1])*rotfactor, '--',
                    markeredgecolor='#00DD00', markerfacecolor='none', color='#00DD00')
Lve01p10, = ax2.plot(v0p10[:, 0], np.abs(v0p10[:, 1]-v1p10[:, 1])*rotfactor, '--',
                     markeredgecolor='#0000FF', markerfacecolor='none', color='#0000FF')
Lve01p12, = ax2.plot(v0p12[:, 0], np.abs(v0p12[:, 1]-v1p12[:, 1])*rotfactor, '--',
                     markeredgecolor='#00DDDD', markerfacecolor='none', color='#00DDDD')
Lve01p16, = ax2.plot(v0p16[:, 0], np.abs(v0p16[:, 1]-v1p16[:, 1])*rotfactor, '--',
                     markeredgecolor='#FF00DD', markerfacecolor='none', color='#FF00DD')
Lve01p24, = ax2.plot(v0p24[:, 0], np.abs(v0p24[:, 1]-v1p24[:, 1])*rotfactor, '--',
                     markeredgecolor='#000000', markerfacecolor='none', color='#000000')

ax2.set_ylim(1e-14, 1e-2)
ax2.set_yscale('log')
ax2.set_ylabel(r'$\lvert \Delta \omega_y \rvert/(T/8\pi\eta a)$')

ax2.set_xlim(0.01, 5.0)
ax2.set_xscale('log')
ax2.set_xlabel(r'$\epsilon/a$')
plt.setp(ax2.get_yticklabels(), visible=False)

ax2.legend((Lve01p6, Lve01p8, Lve01p10, Lve01p12, Lve01p16, Lve01p24),
           (
               r'$p=6$',
               r'$p=8$',
               r'$p=10$',
               r'$p=12$',
               r'$p=16$',
               r'$p=24$',
), frameon=True, loc=r'upper right', fontsize=7, numpoints=1, scatterpoints=1, ncol=1)


# this is an inset axes over the main axes
axin = plt.axes([.23, .3, .1, .15])

axin.plot(v0p6[:, 0], v0p6[:, 1]*rotfactor, '-',
          markeredgecolor='#FF0000', markerfacecolor='none', color='#FF0000', linewidth=1)
axin.plot(v0p8[:, 0], v0p8[:, 1]*rotfactor, '-',
          markeredgecolor='#00DD00', markerfacecolor='none', color='#00DD00', linewidth=1)
axin.plot(v0p10[:, 0], v0p10[:, 1]*rotfactor, '-',
          markeredgecolor='#0000FF', markerfacecolor='none', color='#0000FF', linewidth=1)
axin.plot(v0p12[:, 0], v0p12[:, 1]*rotfactor, '-',
          markeredgecolor='#00DDDD', markerfacecolor='none', color='#00DDDD', linewidth=1)
axin.plot(v0p16[:, 0], v0p16[:, 1]*rotfactor, '-',
          markeredgecolor='#FF00DD', markerfacecolor='none', color='#FF00DD', linewidth=1)
axin.plot(v0p24[:, 0], v0p24[:, 1]*rotfactor, '-',
          markeredgecolor='#000000', markerfacecolor='none', color='#000000', linewidth=1)

# axin.set_ylim(-1.41, -1.39)
axin.set_ylabel(r'$\omega_y/(T/8\pi\eta a)$')
axin.set_xlim(0.01, 0.2)
axin.set_xscale('log')
axin.set_xlabel(r'$\epsilon/a$')

ax1.tick_params(direction='in',which='both')
ax2.tick_params(direction='in',which='both')
axin.tick_params(direction='in',which='both')

# manual page 121
# The pads are specified in fraction of fontsize.
# Meanwhile, use of pad at least larger than 0.3 is recommended.
plt.tight_layout(pad=0.3, w_pad=0.3, h_pad=0.3)

# export
plt.savefig('Fig_Lubroll.png', bbox_inches=0, dpi=300)
# plt.show()
