import numpy as np
import math
import random
import shutil
import subprocess
import pandas as pd
import os
import copy

from scipy.interpolate import griddata


def load_uge_mesh(fileName):
    dataSet = []
    with open(fileName) as fr:
        for line in fr.readlines():
            line.replace(' ', '')
            if line.startswith('#'):
                continue
            curline = line.strip().split()
            try:  # test if all elements in the list are numbers
                fltline = list(map(float, curline))  # if so, convert string numbers to numbers
                fltline[0] = int(fltline[0])  # convert sequence number of cells to integers
                if len(fltline) == 6:
                    fltline[1] = int(fltline[1])
            except:  # if not, append the original list
                dataSet.append(curline)
            else:
                dataSet.append(fltline)
    return dataSet


def dist(t1, t2):
    d, se = 0.0, []
    for i in range(min(len(t1), len(t2))):
        sea = t1[i] - t2[i]
        d += sea ** 2
        se.append(sea)
    return math.sqrt(d), se


def data_convert(mat, flow=4):
    string = '0.0000E+00'
    if mat != 0:
        magnitude, mac = 0, abs(mat)
        while mac < 0.1 or mac >= 1.0:
            if mac < 0.1:
                mac *= 10
                magnitude -= 1
            if mac >= 1.0:
                mac /= 10
                magnitude += 1
        pre = int(round(mac, flow) * 10 ** flow)
        if pre == 10 ** flow:  # solve the bug of mac=9999
            pre = 10 ** (flow - 1)
            magnitude += 1
        string_mag, string_flo = f'E+{str(magnitude).zfill(2)}', f'-.{pre}'
        if magnitude < 0:
            string_mag = f'E-{str(-magnitude).zfill(2)}'
        if mat > 0:
            string_flo = f'0.{pre}'
        string = string_flo + string_mag
    return string


def cell_num_trans(numcells):
    a3, pret = '', int(numcells / 100)
    if pret > 0:
        n = int(math.log(pret, 36))
        for i in range(n):
            a3 += num_tran_letter(int(pret / 36 ** (n - i)))
            pret = int(pret % (36 ** (n - i)))
    a3 += num_tran_letter(pret)
    cell_str = a3.zfill(3) + str(numcells % 100).rjust(2)
    return cell_str


def num_tran_letter(num):
    if num >= 10:
        return str(chr(87 + num).upper())
    return str(num)


def uge2tough2(cell, connection, ss, aperture, anay, scale=1.0):
    cell_number, tough_conn, tough_mesh = [], [], []
    inject, extract, inja, exta, inj_p, ext_p = [], [], [], [], [], []
    inj, ext = np.inf, np.inf
    num_inj, num_ext = 0, 0
    for d in ss:
        if d[4] > 0:
            inject = d
        elif d[4] < 0:
            extract = d
    for t in range(len(cell)):
        cell_name, rock_type = f'{cell_num_trans(cell[t][0])}', ' ' * 4 + '1'
        cell_volume, center, belonging = f'{data_convert(cell[t][4] * aperture * scale)}', cell[t][1:4], []
        cell_number.append(cell[t][0])
        cell_heat_area = data_convert(2 * cell[t][4])
        if len(anay) > 1:  # judge which fracture the point belongs to, and solve the heat exchange area with confining beds
            proj_cos = -1.0
            for ac in anay:
                check = np.dot(ac[:3], center) + ac[3]
                if abs(check) <= 1e-6:
                    if proj_cos < 0:  # no fracture has included this point
                        proj_cos = 1.0
                    else:
                        cell_heat_area = '0.0000E+00'
        tocal = [cell_name, ' ' * 10, rock_type, cell_volume, cell_heat_area, ' ' * 10]
        for ci in center:  # "center" coordinates
            tocal.append(data_convert(ci))
        tough_mesh.append(tocal)
        if inject:
            new_inj, t_1 = dist(center, inject[:3])
            if new_inj < inj:
                inj = new_inj
                inj_p = copy.deepcopy(center)
                num_inj = cell_num_trans(t + 1)
        if extract:
            new_ext, t_2 = dist(center, extract[:3])
            if new_ext < ext:
                ext = new_ext
                ext_p = copy.deepcopy(center)
                num_ext = cell_num_trans(t + 1)
    for i in range(len(connection)):  # read and calculate connection data
        conne_1, conne_2 = connection[i][0], connection[i][1]
        cell_index_1, cell_index_2 = cell_number.index(conne_1), cell_number.index(conne_2)
        pub_pot = connection[i][2:5]
        point_1, point_2 = cell[cell_index_1][1:4], cell[cell_index_2][1:4]
        len_for_conn_1, vec1 = dist(pub_pot, point_1)
        len_for_conn_2, vec2 = dist(pub_pot, point_2)  # connection lengths
        con_conn_1, con_conn_2 = data_convert(len_for_conn_1), data_convert(len_for_conn_2)
        cosine = (point_1[2] - point_2[2]) / dist(point_1, point_2)[0]  # connection direction
        area = connection[i][5] * aperture * scale  # contact area
        trans_co_1, trans_co_2 = cell_num_trans(conne_1), cell_num_trans(conne_2)
        total = f'{trans_co_1}{trans_co_2}' + ' ' * 19  # connect mesh name and blank records
        if abs(cosine) == 1:  # connect in z direction (=3)
            total += f'3{con_conn_1}{con_conn_2}{data_convert(area)}{data_convert(cosine)}'
        elif cosine == 0:  # omit cosine value when it is 0
            if point_2[0] == point_1[0] and point_2[2] == point_1[2]:  # connect in y direction
                total += f'2{con_conn_1}{con_conn_2}{data_convert(area)}'
            else:  # connect in x direction
                total += f'1{con_conn_1}{con_conn_2}{data_convert(area)}'
        else:  # connection not perpendicular or parallel with any of the three axes
            total += f'1{con_conn_1}{con_conn_2}{data_convert(area)}{data_convert(cosine)}'
        tough_conn.append(total)
    write_string('next.txt', 'a', f'source_sinks:inject={num_inj}{inj_p}\textract={num_ext}{ext_p}\n')
    write_string('next.txt', 'a', f'{len(anay)} fractures have eleme={len(cell)} and conne={len(connection)}')
    if inject:
        inja = [num_inj, inject[3], inject[4], inj_p]
    if extract:
        exta = [num_ext, extract[3], extract[4], ext_p]
    return tough_mesh, tough_conn, inja, exta


def write_string(filename, mode, string=''):
    fi = open(filename, mode)
    if not string == '':
        fi.write(string)
    fi.close()


def gen_input(dirt, param, rep_str, formats, orifile, desfile):
    th_cond_const = ''
    for oi, de in zip(orifile, desfile):
        with open(oi) as fch:
            gt = fch.read()
        if oi == f'{dirt}/fl_1.inp':
            for a, b, c in zip(param[:-1], rep_str[:-1], formats):
                gt = gt.replace(b, a.rjust(c))
        elif oi == f'{dirt}/chem.inp':
            gt = gt.replace(rep_str[-1], param[-1])
        write_string(de, 'w', gt)


def parse_con_p_meshes(eleme, conne, T_top, T_bot, top, bot, qloss, inj_T, injs=None, inpt=None, p=1.01325E+05, rho=1000, g=-9.80665, ref_d=0.0, out1=None, out2=None):
    w_init, w_mesh, w_l_t = 'INCON\n', 'ELEME\n', ''
    grad_T = (T_bot - T_top) / (bot - top)
    for single in eleme:
        p0 = p + rho * g * (eval(single[8]) - ref_d)  # reference to static hydraulic pressure
        if out1 is not None:
            cx, cy, cz, cp = out1[0], out1[1], out1[2], out1[3]
            testpt = np.array(list(map(float, single[6:])))
            p0 = griddata((cx, cy, cz), cp, testpt, method='nearest')[0] * 10 ** 6
        T = T_top + grad_T * (eval(single[8]) - top)
        if out2 is not None:
            tcx, tcy, tcz, ct = out2[0], out2[1], out2[2], out2[3]
            testpt = np.array(list(map(float, single[6:])))
            T = griddata((tcx, tcy, tcz), ct, testpt, method='nearest')[0]
        T_start = T * (1 - qloss) + inj_T * qloss
        if injs is not None and inpt is not None and single[0] == injs[0]:
            p0, T = inpt[0], inpt[1]
            single[3] = '0.1000E+40'
        for tc in range(len(single)):
            w_mesh += single[tc]
        w_mesh += '\n'
        w_init += f'{single[0]}\n {data_convert(p0, 13)} {data_convert(T, 13)}\n'
        w_l_t += f'{data_convert(T_start, 13)}\n'
    w_mesh += '\nCONNE\n'
    w_init += '\n'
    w_l_t += '\n'
    for t in conne:
        w_mesh += f'{t}\n'
    w_mesh += '\n'
    return w_mesh, w_init, w_l_t


def main_script(dirt, inj, ext, perm, timet, gener, tracer):
    rep = ['perm1' * 2, 'perm2' * 2, 't' * 10, 'tm' * 5, 'injme', 'extme', 'injra' * 2, 'bon_inj_c']
    with open(f'{dirt}/solu.inp') as solu:
        gb = solu.read()
    gb = gb.replace('extme', ext[0])
    write_string(f'{dirt}/solute.inp', 'w', gb)
    orifilelist, formata = [f'{dirt}/fl_1.inp', f'{dirt}/chem.inp'], [10, 10, 10, 10, 5, 5, 10]
    goal_flo, goal_chem = [], []
    for i, j, l, m in zip(timet['total'], timet['max_step'], gener, tracer):
        tct = timet['total'].index(i) 
        fch_pa = [perm[0], perm[1], i, j, inj[0], ext[0], l, m]
        flo_file, chem_file = f'{dirt}/flow_{tct}.inp', f'{dirt}/chemical_{tct}.inp'
        goal_flo.append(flo_file)
        goal_chem.append(chem_file)
        desfilelist = [flo_file, chem_file]
        gen_input(dirt, fch_pa, rep, formata, orifilelist, desfilelist)
    return goal_flo, goal_chem


def run_tough(dirt, tss, flo_inp, chem_inp, t_unit, inip, iniT, n_mo=1, c=True):
    tg, cg, t_T, pp, tt = None, None, None, None, None
    os.chdir(dirt)
    for i, j, k in zip(flo_inp, chem_inp, tss):
        stat = flo_inp.index(i)
        shutil.copy(i, f'{dirt}/flow.inp')
        shutil.copy(j, f'{dirt}/chemical.inp')
        shutil.copy(f'{dirt}/INCON', f'{dirt}/INCON_{stat}')
        subprocess.call('main_tough_2')  # short period (t1) inpulse tracer injection
        if stat != 2 and c:
            tg, cg = extract_c_t(f'{dirt}/kdd_tim.dat', float(k), t_unit, tg, cg)
            if not tg:  # TOUGH2 simulation
                break
        t_T, pp, tt = read_foft(f'{dirt}/FOFT', float(k), n_mo, iniT, inip, t_T, pp, tt, 8.64E4)
        if not t_T:
            break
        ori_io = [f'{dirt}/kdd_tim.dat', f'{dirt}/iter.dat', f'{dirt}/SAVE', f'{dirt}/FOFT', f'{dirt}/savechem', f'{dirt}/kdd_conc.dat']  # restart
        des_io = [f'{dirt}/kdd_tim_{stat}.dat', f'{dirt}/iter_{stat}.dat', f'{dirt}/INCON', f'{dirt}/FOFT_{stat}', f'{dirt}/inchem', f'{dirt}/kdd_conc_{stat}.dat']
        shutil.copy(f'{dirt}/INCON', f'{dirt}/INCON_{stat}')
        for o, d in zip(ori_io, des_io):
            shutil.copy(o, d)
        shutil.copy(f'{dirt}/GOFT', f'{dirt}/GOFT_{stat}')
        shutil.copy(f'{dirt}/TOFT', f'{dirt}/TOFT_{stat}')
    os.chdir('../')
    return tg, cg, t_T, pp, tt


def extract_c_t(c_t_path, total_time, unit='day', exist_tim=None, exist_con=None):
    read_c_t = pd.read_table(c_t_path, header=4, sep='\s+').T
    tim, con = read_c_t.iloc[1, :], list(read_c_t.iloc[10, :])
    if unit == 'day' or unit == 'd':
        tims = list(tim)
    elif unit == 'hour' or unit == 'h':
        tims = list(tim * 24)
    elif unit == 'year' or unit == 'yr':
        tims = list(tim / 365.24)
    else:
        tims = list(tim * 8.64E4)
    if exist_tim is None:
        exist_tim, exist_con = [], []
    if abs(1 - max(tim) * 8.64E4 / total_time) <= 1e-4:
        exist_tim.extend(tims)
        exist_con.extend(con)
    else:
        exist_tim, exist_con = [], []
    return exist_tim, exist_con


def read_foft(foft_path, total_time, n_mesh, T0, p0, t_pts=None, p_total=None, T_total=None, t_unit=8.64E4):
    if t_pts is None:
        t_pts, T_total, p_total = [0.0], [[T0[i]] for i in range(n_mesh)], [[p0[i]] for i in range(n_mesh)]
    with open(foft_path) as foft:
        for ct in foft.readlines():
            fs = ct.strip().split(',')
            ti = float(fs[1].strip()) / t_unit
            t_pts.append(ti)
            for i in range(n_mesh):
                start = 5 * i + 2
                temps, press = float(fs[start + 2].strip()), float(fs[start + 1].strip()) / 1.0E6
                T_total[i].append(temps)
                p_total[i].append(press)
    if abs(1 - max(t_pts) * 8.64E4 / total_time) > 1e-4:
        t_pts, T_total, p_total = [], [[] for i in range(n_mesh)], [[] for i in range(n_mesh)]
    return t_pts, p_total, T_total


