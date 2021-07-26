import pathlib
import collections
import subprocess
import sys
import random
import numpy as np
import itertools


# input_data
el = ['C', 'O', 'H', 'N', 'S']
mz = 16
psd, psd_channel, psd_up, psd_down = [], [], [], []
psd_iwall, psd_channel2 = [], []
pores, iwall = [], []

def main():
    # read data
    file = '{0}.dat'.format(str(sys.argv[1]))
    source_data = pathlib.Path(file).read_text()
    out_data = collections.defaultdict(list)
    for line in source_data.splitlines()[2:]:
        element, column1, _, column2, *_ = line.split()
        out_data[column1].append([element, column2])

    # new psd
    for pore in out_data.values():
        pores.append(pore)
        pore_down = [i for i in pore if i[0] in el and int(i[1])<=(mz/2)]
        pore_up = [i for i in pore if i[0] in el and int(i[1])>(mz/2)]
        psd_up.append(len(pore_up))
        psd_down.append(len(pore_down))
        pore_all = [i for i in pore if i[0] in el]
        psd.append(len(pore_all))
        # choose the main channel
        if len(pore) == mz:
            pore_iwall = [i for i in pore if i[0] in el and int(i[1]) in np.arange(5, 13)]
            psd_channel2.append(len(pore_iwall))
        else:
            psd_channel.append(len(pore_all))

    # add 10 gas sites for each side of the pore
    lgas = np.full(10, (mz/2))
    psd_up[:0] = lgas
    psd_up[len(psd_up):] = lgas
    psd_down[:0] = lgas
    psd_down[len(psd_up):] = lgas

    # pore lengths
    pore_len = [len(x) for x in pores]
    len_dist = []
    for value, counter in itertools.groupby(pore_len):
        len_dist.append(len(list(counter)))

    # choose iwalls only
    iwalls2 = []
    inputt = iter(pores)
    pore_dist = [list(itertools.islice(inputt, elem)) for elem in len_dist]
    for iwalls in pore_dist:
        if len(iwalls[0]) == mz:
            iwalls2.append(np.hstack(iwalls))

    # reshape the matrix
    iw_new = []
    for pp in iwalls2:
        for p in pp:
            p_new = [list(item) for item in zip(p[::2], p[1::2])]
            iw_new.append(p_new)

    # choose only iwalls
    for n in iw_new:
        iw = []
        for ii in n:
            if ii[0] in el:
                if int(ii[1]) < 5 or int(ii[1]) > 12:
                    iw.append(ii)
        iwall.append(iw)
    iwall_sorted = [x for x in iwall if x!= []]
    for iw_sorted in iwall_sorted:
        psd_iwall.append(len(iw_sorted))

    # count the number of pores of a particular size
    psd_count = collections.Counter(psd_channel) + collections.Counter(psd_channel2) + collections.Counter(psd_iwall)

    # update psite in mft_nvt.dat
    psite = sum(psd)
    tsite = psite + 2*len(lgas)
    # with open("mft.dat") as file:
    #     content = file.read()
    #     lines = content.splitlines()
    #     elements = lines[3].split()
    #     elements[0] = str(psite)
    #     elements[2] = str(tsite)
    #     nline = " ".join(elements)
    #     ncontent = content.replace(lines[1], nline)
    #     with open("mft.dat", 'w') as file:
    #         file.write(ncontent)

    # save to file
    np.savetxt('mft_psd.dat', np.vstack((psd_up, psd_down)), fmt='%1.0f')
    np.savetxt('{0}_mft_psd.dat'.format(str(sys.argv[1])), np.vstack((psd_up, psd_down)), fmt='%1.0f')

    with open("{0}_mft_psd_count.dat".format(str(sys.argv[1])), 'w') as file:
        for key, val in sorted(psd_count.items()):
            file.writelines(str(key) + '\t' + str(val) + '\n')
        file.writelines(str(psite))
        file.close()
main()

# run mft with new psd
current_working_dir = pathlib.Path.cwd()
subprocess.call(['ifort', '-O2', '-o', 'mft.exe', 'mft2.f'])
subprocess.call(['./mft.exe'])