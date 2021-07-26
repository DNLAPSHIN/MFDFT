import numpy as np
import subprocess
import random
import itertools
import os
import glob


# input
ilength = 8        # the length of the intrawall pore ilength*2=16
p_length = 5       # the length of the single pore
sec_length = 20    # the length of the pore section between intrawall pores
iwn = 8            # the number of intrawall pores of different size, depends on sigma

# the main channel
l1 = np.full(0, 6)
l2 = np.full(66, 7)
l3 = np.full(368, 8)
l4 = np.full(66, 9)
l5 = np.full(0, 10)
pore = list(l1)+list(l2)+list(l3)+list(l4)+list(l5)
random.shuffle(pore)
pore = np.array(pore)
psd1, psd2 = [], []
for i in pore:
    if (i %2) == 0:
        x = i/2
        y = i-x
        psd1.append(int(x))
        psd2.append(int(y))
    else:
        xx = (i+random.randrange(-1,2,2))/2
        yy = i-xx
        psd1.append(int(xx))
        psd2.append(int(yy))
psd = [psd1, psd2]

# intrawall pore distribution
iwd = int((len(pore) * p_length / sec_length - 1) / iwn)
# the number of intrawall pores with different widths
l6 = np.full(iwd, 2)
l7 = np.full((iwd+1), 3)
l8 = np.full(iwd, 4)
l9 = np.full((iwd+1), 5)
l10 = np.full((4*iwd+2), 6)
iw = list(l6) + list(l7) + list(l8) + list(l9) + list(l10)
random.shuffle(iw)
iwidth = np.array(iw)
iwalls = []
for j in iwidth:
    iwalls.append(list(itertools.repeat(ilength, j)))
# the total number of pore sites
npsites = len(pore) * p_length + sum(iw)
# network
network = []
for ipsd in psd:
    # generate segments of a specific length
    segments = [k for k in ipsd for _ in range(p_length)]
    # split segments into sections of 20 length
    sections = [segments[l:l + sec_length] for l in range(0, len(segments), sec_length)]
    # generate the network
    def netw(sections, iwalls):
        outcome = []
        ls = len(sections)
        lw = len(iwalls)
        for i in range(max(ls, lw)):
            if i < ls:
                outcome.append(sections[i])
            if i < lw:
                outcome.append(iwalls[i])
        return outcome
    network.append(list(itertools.chain.from_iterable(netw(sections, iwalls))))

# update the number of sites in *.dat file
with open("mft_nvt.dat") as file:
    content = file.read()
    lines = content.splitlines()
    elements = lines[1].split()
    elements[0] = str(npsites)
    elements[2] = str(npsites)
    nline = " ".join(elements)
    ncontent = content.replace(lines[1], nline)
with open("mft_nvt.dat", 'w') as ufile:
    ufile.write(ncontent)

# save psd to file
np.savetxt('mft_psd.dat', np.vstack((network[0], network[1])), fmt='%1.0f')

# run MFT
# current_working_dir = os.getcwd()
# subprocess.call(['ifort', '-O2', '-o', 'mft.exe', 'mft_nvt.f'])
# subprocess.call(['./mft.exe'])