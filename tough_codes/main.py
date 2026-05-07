import subprocess
import shutil
import os

outs = []
with open('sequence.txt') as t2react_exe:
    cmg_o = 'gfortran -o main_tough '
    used = t2react_exe.readlines()
for cmd in used[1:]:
    f_files = cmd.strip().split('\t')
    cmg_c = 'gfortran -c '
    for f in f_files:
        cmg_c += f'{f} '
        cmg_o += f'{f[:-2]}.o '
        outs.append(f'{f[:-2]}.o')
    subprocess.call(cmg_c, shell=True)
subprocess.call(cmg_o, shell=True)
shutil.copy('main_tough', f'/home/{used[1].strip()}/.local/bin/main_tough_2')
for de in outs:
    if os.path.isfile(de):
        os.remove(de)
print('compile successful!')
