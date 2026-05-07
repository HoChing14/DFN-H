import numpy as np
import math
import matplotlib.pyplot as plt
import subprocess

from scipy.interpolate import interp1d

def ana_result(t, x, son, mot, v, T0, Ti):
    h = t - x / v
    if h <= 0:
        return T0
    final = son * x / mot / np.sqrt(v * h)
    T = T0 + (Ti - T0) * math.erfc(final) * np.heaviside(h, 1)
    return T


def read_foft(foft_path, total_time, n_mesh, T0, t_pts=None, P=None, T=None, t_unit=8.64E4):
    if t_pts is None:
        t_pts, T_total = [0.0], [[T0] for i in range(n_mesh)]
    with open(foft_path) as foft:
        for ct in foft.readlines():
            fs = ct.strip().split(',')
            ti = float(fs[1].strip()) / t_unit
            t_pts.append(ti)
            for i in range(n_mesh):
                start = 5 * i + 2
                temps = float(fs[start + 2].strip())
                T_total[i].append(temps)
    if abs(1 - max(t_pts) * 8.64E4 / total_time) > 1e-4:
        t_pts, T_total = [], [[] for i in range(n_mesh)]
    return t_pts, T_total


def read_toft(toft_path, n_time, model_mesh, skipevery=0):  # temperature distribution at a certain time
    Tt_total = [[] for i in range(n_time)]
    with open(toft_path) as toft:
        fc = toft.readlines()
    for i in range(n_time):
        start = (model_mesh + 2) * i + 2
        data = fc[start: start + model_mesh]
        count = 0
        for tt in data:
            count += 1
            if skipevery > 1 and (count - 1) % skipevery == 0:
                voids = tt.strip().split()
                Tt_total[i].append(float(voids[-1]))
    return Tt_total


xcoords = []
skip_space, counts = 5, 0
with open('MESH') as mesh:  # get xcoords for analytical solution
    for eleme in mesh.readlines()[1:]:
        try:
            info = eleme.strip()[5:].split()[1]
        except IndexError:
            break
        counts += 1
        if skip_space > 1 and (counts - 1) % skip_space == 0:
            xcoords.append(float(info[:10]))
fig = plt.figure(figsize=(8, 6), dpi=100)
ax1, ax2 = fig.add_subplot(121), fig.add_subplot(122)
xlocks, tlocks = [36.25], [300.0]  # space and time for temperature distributions
tcoords = [i * 40 for i in range(11)]  # time points for scatter
rouw, cw, uf = 1000.0, 4200.0, 0.008
roum, cm, lamdam = 2600.0, 920.0, 3.0
T_ini, T_inj, b = 200.0, 30.0, 0.001
erfc_son = lamdam / rouw / cw / b
erfc_mot = np.sqrt(uf * lamdam / roum / cm)
subprocess.call('main_tough_2', shell=True)
t_sin, T_sin_total = read_foft('FOFT', 3.456E7, len(xlocks), T_ini)  # read PT data of multiple points
Tx_sim = read_toft('TOFT', len(tlocks), 80, skip_space)  # manual separation of PT distribution file
T_tt, T_xx = open('temporal.txt', 'w'), open('spatial.txt', 'w')
for T_sin, xlock in zip(T_sin_total, xlocks):  # calculate analytical results and compare with numerical results
    T_obs = [T_ini]
    T_tt.write(f'temperature versus time at x={xlock}m\n')
    T_tt.write(f'0.0,{T_ini},{T_ini}\n')
    sum_error, max_error, max_e_time = 0.0, 0.0, 0.0
    f1 = interp1d(t_sin, T_sin, kind='cubic')
    T_interp = list(f1(tcoords))
    ax1.plot(tcoords, T_interp, label='calculated')
    for cc, cd in zip(tcoords[1:], T_interp[1:]):
        ccs = cc * 8.64E4
        T_real = ana_result(ccs, xlock, erfc_son, erfc_mot, uf, T_ini, T_inj)
        T_obs.append(T_real)
        ab_error = abs(cd - T_real)
        s_e = ab_error ** 2
        if ab_error > max_error:
            max_error = ab_error
            max_e_time = cc
        sum_error += s_e
        T_tt.write(f'{cc},{T_real},{cd}\n')
    RR = np.sqrt(sum_error / len(T_interp))
    ax1.scatter(tcoords, T_obs, label=f'x={xlock}m')
    #print(f'result temperature .vs. time at x={xlock} m:')
    #print(f'simulation error={RR} degC ({RR / (max(T_obs) - min(T_obs)) * 100}%)')
    #print(f'max error={max_error} degC ({max_error / (max(T_obs) - min(T_obs)) * 100}%) at {max_e_time} days\n')
T_tt.close()
for tlock, Ts in zip(tlocks, Tx_sim):  # calculate analytical results and compare with numerical results
    tlocka = tlock * 8.64E4
    T_xx.write(f'temperature versus space at t={tlock}d\n')
    T_obs2 = []
    sum_error, max_error, max_e_time = 0.0, 0.0, 0.0
    ax2.plot(xcoords, Ts, label='calculated')
    for cc, cd in zip(xcoords, Ts):
        T_real2 = ana_result(tlocka, cc, erfc_son, erfc_mot, uf, T_ini, T_inj)
        T_obs2.append(T_real2)
        ab_error = abs(cd - T_real2)
        s_e = ab_error ** 2
        if ab_error > max_error:
            max_error = ab_error
            max_e_time = cc
        sum_error += s_e
        T_xx.write(f'{cc},{T_real2},{cd}\n')
    RR = np.sqrt(sum_error / len(xcoords))
    ax2.scatter(xcoords, T_obs2, label=f't={tlock}d')
    #print(f'result temperature .vs. space at t={tlock} d:')
    #print(f'simulation error={RR} degC ({RR / (max(T_obs) - min(T_obs)) * 100}%)')
    #print(f'max error={max_error} degC ({max_error / (max(T_obs) - min(T_obs)) * 100}%) at {max_e_time} m\n')
T_xx.close()
ax1.legend()
ax2.legend()
plt.savefig('fit_ana.png')
plt.savefig('fit_ana.pdf')

