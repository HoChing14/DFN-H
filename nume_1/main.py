import copy
import random
import subprocess
import os
import shutil

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from scipy.interpolate import interp1d
from time import time
from tough_script import load_uge_mesh, uge2tough2, parse_con_p_meshes
from tough_script import main_script, write_string, read_foft, run_tough


def read_obs(filename, short='0'):
    t_obs, c_obs, c_obs2, sh_t_obs, sh_c_obs = [], [], [], [], []
    with open(filename) as obs:
        for iss in obs.readlines():
            pte = iss.strip().split(',')
            t_obs.append(float(pte[0]))
            c_obs.append(float(pte[1]))
            if len(pte) > 2:
                c_obs2.append(float(pte[2]))
            if float(short) != 0 and float(pte[0]) > float(short):
                continue
            sh_t_obs.append(float(pte[0]))
            sh_c_obs.append(float(pte[1]))
    return t_obs, c_obs, sh_t_obs, sh_c_obs, c_obs2


def compliklihood(t_cal, c_cal, t_obs, c_obs, t_inj=0, t_ini=0, tem=False):
    ab_errors = []
    sum2, max_err1 = 0.0, 0.0
    f1 = interp1d(t_cal, c_cal, kind='cubic')
    ct_cal = list(f1(t_obs))
    div_fac = max(c_obs) - min(c_obs)
    if tem:
        div_fac = abs(t_ini - t_inj)
    for cal, obs in zip(ct_cal, c_obs):
        abs_e = abs(cal - obs)
        ab_errors.append(abs_e)
        e = abs_e ** 2
        sum2 += e
        rel_e = abs_e / div_fac
        if rel_e > max_err1:
            max_err1 = rel_e
    rmse1 = np.sqrt(sum2 / len(c_obs)) / div_fac
    return rmse1, max_err1, ct_cal, ab_errors


def gen_ran_with_acc(minv, maxv, acc, longs):
    rands = random.randint(0, longs)
    get = minv + acc * rands
    rest = maxv - get
    tail = min(acc, rest)
    get += random.uniform(0, tail)
    return get


def write_current_accept(iters, rmse, ref_new, valid_ratio):
    fst, ipt = '', 0
    for st in rmse:
        fst += str(st)
        ipt += 1
        if ipt != 7:
            fst += ','
    fst += '\n'
    os.mkdir(f'{os.getcwd()}/accept/result_{iters}_{rmse[0]:.5e}')
    shutil.copy('cal_con.txt', f'{os.getcwd()}/accept/result_{iters}_{rmse[0]:.5e}/cal_con_{rmse[1]:.5e}.txt')
    shutil.copy('cal_P_T.txt', f'{os.getcwd()}/accept/result_{iters}_{rmse[0]:.5e}/cal_P_T_{rmse[2]:.5e}.txt')
    write_string('accept.txt', 'a', f'{iters}\t{rmse[0]}({rmse[1]},{rmse[2]})\t{valid_ratio}\t{len(ref_new)}\n')
    ipt = 0
    for fs in ref_new:
        ip = 0
        for fg in fs[:-1]:
            fst += str(fg)
            ip += 1
            if ip != 6:
                fst += ','
        ipt += 1
        if ipt != len(ref_new):
            fst += ';'
    fst += f'\n{valid_ratio}\n'
    write_string('current_ref.txt', 'w', fst)
    write_string('fit_ref_log.txt', 'a', f'{iters}\n{fst}\n\n')


def accept_or_reject(iters, pr_new, pr_old=None):
    pr_acc = copy.deepcopy(pr_new)
    if pr_old is None:
        write_current_accept(1, pr_new[0], pr_new[1], pr_new[2])
    else:
        if pr_new[0][0] < pr_old[0][0]:
            print(f'new comprehensive rmse={pr_new[0][0]}')
            write_current_accept(iters, pr_new[0], pr_new[1], pr_new[2])
        else:
            print(f'present rmse={pr_old[0][0]}, rejected rmse={pr_new[0][0]}')
            pr_acc = copy.deepcopy(pr_old)
    return pr_acc


def read_coord(filename):
    coords, vecs, nors = [], [], []  # vectors of edges, normal vectors
    with open(filename) as frac:
        dpd = frac.readlines()[3:][:-1]
    for dp in dpd:
        pots = []
        dpt = dp.strip().split('\t')
        for dta in dpt:
            dga = dta[1:-1].split(',')
            dpa = list(map(float, dga))
            pots.append(dpa)
        v1 = [pots[1][0] - pots[0][0], pots[1][1] - pots[0][1], pots[1][2] - pots[0][2]]
        v2 = [pots[2][0] - pots[1][0], pots[2][1] - pots[1][1], pots[2][2] - pots[1][2]]
        nor1 = np.cross(v1, v2)
        nor = list(nor1 / np.linalg.norm(nor1))
        trant = - np.dot(nor, pots[0])
        nor.append(trant)
        v3 = [pots[3][0] - pots[2][0], pots[3][1] - pots[2][1], pots[3][2] - pots[2][2]]
        v4 = [pots[0][0] - pots[3][0], pots[0][1] - pots[3][1], pots[0][2] - pots[3][2]]
        vecs.append([v1, v2, v3, v4])  # BC vectors of a fracture
        coords.append(pots)
        nors.append(nor)
    return coords, vecs, nors


fig = plt.figure(figsize=(8, 6), dpi=100)
ff, scto = os.getcwd(), time()
td = open('param.txt', 'r')
T_range = td.readline().strip().split(':')[1].split(',')
time_split = td.readline().strip().split(':')[1].split(',')
max_time = td.readline().strip().split(':')[1].split(',')
tracer = ['1.d-10' for i in range(len(max_time))]
inj_param = td.readline().strip().split(':')[1].split(',')
tracer[0], inj_T, T_loss = inj_param[0], float(inj_param[2]), float(inj_param[3])
short_time = str(float(td.readline().strip().split(':')[1]) / 8.64E4)
td.close()
top_T, bot_T, acc_T = float(T_range[0]), float(T_range[1]), float(T_range[2])
h_top, h_bot = float(T_range[3]) / 2, - float(T_range[3]) / 2
ss_g = []  # locations of inj/ext
with open('source_sink.txt') as sst:
    for jis in sst.readlines():
        jct = jis.strip().split(',')
        jcs = list(map(float, jct[:3]))
        jcs.append(jct[3])
        jcs.append(float(jct[4]))
        ss_g.append(jcs)
c_t, c_c, c11, c12, sas = read_obs('obs_con.txt')
T_t, T_c, sh_T_t, sh_T_c, T_c2 = read_obs('obs_T.txt', short_time)
p_t, p_c, p11, p12, sas = read_obs('obs_p.txt')
plt.scatter(T_t, T_c, label='observed')
#ax3.scatter(sh_T_t, sh_T_c, label='observed')
dfn_mesh = load_uge_mesh(f'{ff}/f/full_mesh.uge')
times = {'total': time_split, 'max_step': max_time} 
cells, connections = dfn_mesh[1:(int(dfn_mesh[0][1]) + 1)], dfn_mesh[(int(dfn_mesh[0][1]) + 2):]
b = 0.001
perm = [f'{b ** 2 / 12:.2e}', '1.00e-22']
with open(f'{ff}/f/normal_vectors.dat') as oc:
    novs = oc.readlines()  # normal vectors
with open(f'{ff}/f/translations.dat') as lo:
    locs = lo.readlines()  # fractures locations
sog, post_normal = 0, []
for ov in locs[1:-2]:
    transl = ov.strip().split(' ')
    if transl[-1] != 'R':  # exclude isolated fractures
        g = list(map(float, transl))
        norsl = novs[sog].strip().split(' ')
        nov = list(map(float, norsl))
        dm = - nov[0] * g[0] - nov[1] * g[1] - nov[2] * g[2]
        nov.append(dm)
        post_normal.append(nov)
        sog += 1
# main inversion
if os.path.isfile('next.txt'):
    shutil.copy('next.txt', 'next_bak.txt')
    write_string('next.txt', 'w')
ele, con, inj, ext = uge2tough2(cells, connections, ss_g, b, post_normal)  # mesh generation and conversion
w_c, w_s, w_st = parse_con_p_meshes(ele, con, top_T, bot_T, h_top, h_bot, T_loss, inj_T, ref_d=3500.0, out1=None, out2=None)
write_string(f'{ff}/treact/MESH', 'w', w_c)
write_string(f'{ff}/treact/INCON', 'w', w_s)
write_string(f'{ff}/treact2/MESH', 'w', w_c)
write_string(f'{ff}/treact2/INCON', 'w', w_s)
if os.path.isfile('SAVE'):
    shutil.copy('SAVE', f'{ff}/treact/INCON')
    shutil.copy('SAVE', f'{ff}/treact2/INCON')
write_string(f'{ff}/treact/TQLOS', 'w', w_st)
write_string('next.txt', 'a', f'\naperture={b}\t')
print(f'fracture aperture={b}, permeability={perm[0]}')
valid_rates = f'{15.0 * float(inj_param[1]):.2f}'
print(f'assumed valid injection rates={valid_rates} kg/s')
gen = [valid_rates for i in range(len(tracer))]
p_ini, T_ini = [3.005747509168E1, 3.001043649727E1], [193.5, 193.5]
fl, ch = main_script(f'{ff}/treact', inj, ext, perm, times, gen, tracer)
fl0, ch0 = main_script(f'{ff}/treact2', inj, ext, perm, times, gen, tracer)
sct = time()
tc, cc, ttem, ptem2, tem2 = run_tough(f'{ff}/treact', times['total'], fl, ch, 'd', p_ini, T_ini, len(p_ini))
sib = time() - sct
exit(0)
sct0 = time()
tc0, cc0, ttem0, ptem20, tem20 = run_tough(f'{ff}/treact2', times['total'], fl0, ch0, 'd', p_ini, T_ini, len(p_ini))
sib0 = time() - sct0
print(f'forward model run 1 for {sib} seconds')
print(f'forward model run 2 for {sib0} seconds')
write_string('next.txt', 'a', f'\ttime={sib}')
tem, ptem = tem2[1], ptem2[0]
tem0, ptem0 = tem20[1], ptem20[0]
tcs, ccs, tcs0, ccs0, shttem, shtem = [], [], [], [], [], []
w_con, w_tem, w_tep = 'time,concentration\n', 'time,temperature\n', 'time,pressure\n'
plt.plot(ttem, tem, label='calculated', c='r')
plt.plot(ttem0, tem0, label='no rock matrix')
if tc is not None:
    for ac, ad in zip(tc, cc):
        if ac - max(c_t) > 1e-6:
            break
        if ac not in tcs:
            tcs.append(ac)
            ccs.append(max(0, ad))
    for ac, ad in zip(tc0, cc0):
        if ac - max(c_t) > 1e-6:
            break
        if ac not in tcs0:
            tcs0.append(ac)
            ccs0.append(max(0, ad))
    rr_c, max_c, interp_c, td = compliklihood(tcs, ccs, c_t, c_c)
    rr_c0, max_c0, interp_c0, td0 = compliklihood(tcs0, ccs0, c_t, c_c)
    for ac, ad, ae in zip(c_t, interp_c, interp_c0):
        w_con += f'{ac},{ad},{ae}\n'
    write_string('cal_con.txt', 'w', w_con)
for ae, af in zip(ttem, tem):
    if ae <= float(short_time):
        shttem.append(ae)
        shtem.append(af)
rr_T, max_T, interp_T, abe = compliklihood(ttem, tem, T_t, T_c, min(T_c), max(T_c), True)
rr_T0, max_T0, interp_T2, abe0 = compliklihood(ttem0, tem0, T_t, T_c, min(T_c), max(T_c), True)
ccsc, ccsp, interp_p, tcg = compliklihood(ttem, ptem, p_t, p_c)
ccsc0, ccsp0, interp_p2, tca = compliklihood(ttem0, ptem0, p_t, p_c)
#sh_rr_T, sh_max_T, sh_interp_T, tcgg = compliklihood(shttem, shtem, sh_T_t, sh_T_c, min(sh_T_c), max(sh_T_c), True)
#ax3.scatter(T_t, abe)
for ae, af, ad in zip(T_t, interp_T, interp_T2):
    w_tem += f'{ae},{af},{ad}\n'
for ag, ah, ai in zip(p_t, interp_p, interp_p2):
    w_tep += f'{ag},{ah},{ai}\n'
print(f'real injection pressure={p_c[-1]}\ncalculated injection pressure={ptem[-1]}')
write_string('cal_T.txt', 'w', w_tem)
write_string('cal_P.txt', 'w', w_tep)
print(f'long temperature error={rr_T}')
print(f'maximum error in long temperature error={max_T}')
plt.legend()
plt.ylim(75.0, 200.0)
plt.savefig('fit_result.png')
plt.savefig('fit_result.pdf')

