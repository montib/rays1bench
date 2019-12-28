#!/usr/bin/env python3

import os
import glob
import subprocess
import sys


def print_help():
    print("""\
USAGE:
  bench.py <--all | --single-threaded | --latest>
          [--cxx [cl | clang | g++ | ...]] [--save]
          [--compile-only] [--num NUMBER]

  --all               Run all tests from step1 to latest
  --latest            Run latest only
  --single-threaded   Run latest single-threaded version (step12)
  --step NUMBER       Run stepNUMBER
  --cxx COMPILER      Selects compiler, available options:
                          cl: MSVC compiler
                          clang or clang++ or clang++-9 or similar: clang
                          g++ or gcc or g++-9 or similar: g++

  --quick             Quick bench (lower resolution and spp)
  --compile-only      Compile programs but doesn't run
  --num NUMBER        Number of runs (1..32)
  --save              Save the output to a .tga file (src/step*/*.tga)

EXAMPLES:
  bench.py --all
    Compiles all version from step1 to latest with the default compiler.

  bench.py --latest --cxx clang++-9 --save
    Compiles src/latest with clang++-9 and then start it.
    Output saved to src/latest/out_small.tga, out_medium.tga, out_large.tga.

""")

def del_files(pattern):
    files = glob.glob(pattern)
    for f in files:
        try:
            os.remove(f)
        except Exception:
            pass



cpp_compiler = None
save_picture = False
compile_only = False
num_runs = 1
quick_bench_args_cl = ""
quick_bench_args = ""


dirs = None

# parse arguments

for idx, val in enumerate(sys.argv):

    if val == '--compile-only':
        compile_only = True

    if val == '--save':
        save_picture = True

    if val == '--cxx' and idx < len(sys.argv) - 1:
        cpp_compiler = sys.argv[idx+1]

    if val == '--num' and idx < len(sys.argv) - 1:
        n = int(sys.argv[idx+1])
        num_runs = n

    if val == '--latest':
        dirs = ['src/latest']


    if val == '--step' and idx < len(sys.argv) - 1:
        dirs = ['src/step' + sys.argv[idx+1]]

    if val == '--all':
        dirs = ['src/step1',
                'src/step2',
                'src/step3',
                'src/step4',
                'src/step5',
                'src/step6',
                'src/step7',
                'src/step8',
                'src/step9',
                'src/step10',
                'src/step11',
                'src/step12',
                'src/step13',
                'src/latest']

    if val == '--single-threaded':
        dirs = ['src/step12']

    if val == '--quick':
        quick_bench_args_cl = ' /D "QUICKBENCH" '  # msvc
        quick_bench_args = ' -DQUICKBENCH '     # g++/clang++



if dirs is None:
    print('Either --all, --single or --latest must be specified.\n')
    print_help()
    exit(-1)

for i in dirs:
    os.chdir(i)
    print("\n===============\n" + i + "\n===============")

    del_files('out*.*')
    del_files('a.out')
    del_files('*.exe')

    # run python script to generate SoA sources

    py_files = glob.glob('*.py')
    print('py_files = ' + str(py_files))

    for script in py_files:
        try:
            print('running: ' + script)

            if sys.platform.startswith('win'):
                subprocess.call(('python ' + script).split())
            elif sys.platform.startswith('linux'):
                subprocess.call(('python3 ' + script).split())
            elif sys.platform.startswith('darwin'):
                subprocess.call(('python3 ' + script).split())
        except Exception as e:
            print('ERROR: ' + str(e))
            sys.exit(-1)

    # compile

    cpp_files = " ".join(glob.glob('*.cpp'))

    print('cpp_files = ' + cpp_files)

    compile_cmd = None
    executable = None


    if sys.platform.startswith('win'):

        if cpp_compiler == None:
            cpp_compiler = 'cl'

        if 'clang' in cpp_compiler:
            compile_cmd =  (cpp_compiler + quick_bench_args + ' -ffast-math -O3 -g -fno-rtti -fno-exceptions -std=c++17 -march=native -flto -fuse-ld=lld -DNDEBUG -o rayweek1.exe ' + cpp_files).split()
        elif 'g++' in cpp_compiler or 'gcc' in cpp_compiler:
            compile_cmd = (cpp_compiler + quick_bench_args + ' -pthread -ffast-math -O3 -g -fno-rtti -fno-exceptions -std=c++17 -flto -march=native -m64 -DNDEBUG -o rayweek1.exe ' + cpp_files).split()
        elif cpp_compiler == 'cl':
            compile_cmd = 'cl /Ox /fp:fast /D "NDEBUG" /arch:AVX2 /nologo /GL /Qvec-report:1 ' + quick_bench_args_cl + cpp_files + ' /link /LTCG' # !!! no split(), it causes weird bugs
        else:
            compile_cmd = 'compiler not supported on this platform: ' + cpp_compiler

        executable = 'rayweek1.exe'

    elif sys.platform.startswith('linux'):

        if cpp_compiler == None:
            cpp_compiler = 'g++'

        if 'clang' in cpp_compiler:
            compile_cmd =  (cpp_compiler + quick_bench_args + ' -pthread -ffast-math -O3 -g -fno-rtti -fno-exceptions -std=c++17 -flto -march=native -m64 -fuse-ld=lld -DNDEBUG ' + cpp_files).split()
        elif 'g++' in cpp_compiler or 'gcc' in cpp_compiler:
            compile_cmd = (cpp_compiler + quick_bench_args + ' -pthread -ffast-math -O3 -g -fno-rtti -fno-exceptions -std=c++17 -flto -march=native -m64 -DNDEBUG  ' + cpp_files).split()
        else:
            compile_cmd = 'compiler not supported on this platform: ' + cpp_compiler

        executable = './a.out'
    elif sys.platform.startswith('darwin'):
        if cpp_compiler == None:
            cpp_compiler = 'clang++'

        if 'clang' in cpp_compiler:
            compile_cmd = (cpp_compiler + quick_bench_args + ' -ffast-math -O3 -g -fno-rtti -fno-exceptions -std=c++17 -flto -march=native -m64 -DNDEBUG ' + cpp_files).split()
        else:
            compile_cmd = 'compiler not supported on this platform: ' + cpp_compiler

        executable = './a.out'

    print('compile cmd = ' + str(compile_cmd))

    try:
        subprocess.call(compile_cmd)
    except Exception as e:
        print('FAILED to compile')
        print(str(e))
        sys.exit(-1)


    if save_picture:
        executable = executable + ' -w'

    if num_runs > 1:
        executable = executable + ' -n ' +  str(num_runs)

    if not compile_only:
        print('RUN ' + executable + '\n')
        subprocess.call(executable.split())

    os.chdir('../..')
