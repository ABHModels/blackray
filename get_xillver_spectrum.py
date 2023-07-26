from astropy.io import fits
import numpy as np
from sys import argv

# path_to_xillver_file = "/Users/askarabd/Documents/fudan/FITS/xillver-a-Ec5.fits"

gamma = 2
afe = 1.5
logxi = 1.5
ecut = 270
incl = 45

# Input parameters: gamma, afe, logxi, ecut, incl, path_to_xillver_file

gamma = float(argv[1])
afe = float(argv[2])
logxi = float(argv[3])
ecut = float(argv[4])
incl = float(argv[5])
path_to_xillver_file = argv[6]

hdu = fits.open(path_to_xillver_file)

# print(gamma, afe, logxi, ecut, incl)

input_par = np.array([gamma, afe, logxi, ecut, incl])

data = hdu[1].data
pars, pmax = [], []
for item in data:
    # print(item[0], item[-1][:item[-2]])
    pars.append(np.array(item[-1][:item[-2]]))
    pmax.append(item[-2])

data = hdu[2].data
energy_array = np.zeros(len(data))
for i in range(len(data)):
    energy_array[i] = data[i][0]

indices, fracs = [], []
# print(input_par)
for i in range(len(input_par)-1):
    if input_par[i] < pars[i][0]:
        input_par[i] = pars[i][0]
    if input_par[i] > pars[i][-1]:
        input_par[i] = pars[i][-1]
    val = input_par[i]

    inds = np.argwhere(pars[i] > val).flatten()
    if not len(inds):
        inds = len(pars[i]) - 1
    else:
        inds = inds[0]
    indices.append([inds-1, inds])
    fracs.append((val - pars[i][indices[i][0]])/(pars[i][indices[i][1]] - pars[i][indices[i][0]]))    
    # print(val, pars[i][inds-1], pars[i][inds], fracs[i], i)

indices.append(np.arange(len(pars[-1])))
fracs.append(1)

gami, afei, xii, ecti, inci  = indices

sets = []
for g in gami:
    for a in afei:
        for x in xii:
            for e in ecti:
                for i in inci:
                    sets.append([g, a, x, e, i])
                    # print(sets[-1])

data = []
d3 = hdu[3].data
for iset in sets:
    ii, jj, kk, ll, mm = iset
    rownum = (((ii * pmax[1] + jj) * pmax[2] + kk) * pmax[3] + ll) * pmax[4] + mm
    data.append(d3[rownum][1])
    p1 = d3[rownum][0]
    p2 = [pars[0][ii], pars[1][jj], pars[2][kk], pars[3][ll], pars[4][mm]]
    check = np.array(p1 == p2)
    if np.sum(check-1) > 0:
    	print("Error during extracting xillver spectra! ", iset)

data = np.array(data)
sets = np.array(sets)
xill_spec = np.zeros([10, 2999])
for i in range(10):
    inds = 10 * np.arange(16) + i
    length = len(sets[inds])
    sets2 = sets[inds]
#     print(sets[inds])
    inds1 = np.arange(0, length-1, 2)
    inds2 = np.arange(1, length, 2)
#     sets1 = sets[inds][inds1]
#     sets2 = sets[inds][inds2]
#     print(sets1, sets2)
    data2 = data[inds][inds1] + fracs[3] * (data[inds][inds2] - data[inds][inds1])
    # print(data2.shape)    
    inds1 = np.arange(0, int(length/2)-1, 2)
    inds2 = np.arange(1, int(length/2), 2)
    data2 = data2[inds1] + fracs[2] * (data2[inds2] - data2[inds1])
    # print(data2.shape)
    inds1 = np.arange(0, int(length/4)-1, 2)
    inds2 = np.arange(1, int(length/4), 2)
    data2 = data2[inds1] + fracs[1] * (data2[inds2] - data2[inds1])
    # print(data2.shape)
    xill_spec[i] = data2[0] + fracs[0] * (data2[1] - data2[0])
    # print(xill_spec[i].shape)


filename = "data/xill_spec_gam_%.2f_afe_%.2f_xi_%.2f_ecut_%.2f_incl_%.2f.dat" % (gamma, afe, logxi, ecut, incl)
file = open(filename, "w")

for i in range(10):
	for j in range(len(xill_spec[i])):
		line = str(energy_array[j]) + " " + str(xill_spec[i][j]) + "\n"
		file.write(line)

file.close()

