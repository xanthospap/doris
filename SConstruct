from __future__ import print_function
import os, sys, glob

## Prefic for install(ed) files
prefix="/usr/local"
if not os.path.isdir(prefix):
    print('[ERROR] Cannot find \'prefix\' directory, aka {:}; aborting'.format(prefix), file=sys.stderr)
    sys.exit(1)

## Library version
lib_version="0.1.0"
## Library name
lib_name="doris"
## Include dir (following prefix) if any
inc_dir="doris"
## the rootdir of the project
root_dir=os.path.abspath(os.getcwd())

## get number of CPUs and use for parallel builds
num_cpu = int(os.environ.get('NUM_CPU', 2))
SetOption('num_jobs', num_cpu)
print("running with -j %s" % GetOption('num_jobs'))

## user can specify the --prefix variables (expanded to $PREFIX)
AddOption('--prefix',
          dest='prefix',
          type='string',
          nargs=1,
          action='store',
          metavar='DIR',
          help='installation prefix',
          default='/usr/local')
AddOption('--cxx',
          dest='cxx',
          type='string',
          nargs=1,
          action='store',
          metavar='CXX',
          help='C++ Compiler',
          default=None)
AddOption('--std',
          dest='std',
          type='string',
          nargs=1,
          action='store',
          metavar='STD',
          help='C++ Standard [11/14/17/20]',
          default='17')
#AddOption('--branchless',
#          dest='branchls',
#          action='store_true',
#          help='Trigger built with BRANCHLESS defined',
#          default=False)
AddOption('--ignore-tests',
          nargs=1,
          type='string',
          action='store',
          metavar='IGNORE_TEST_FILES',
          dest='ignore_tests',
          help='Comma-seperated list of test source files to ignore when building')
AddOption('--include-tests',
          nargs=1,
          type='string',
          action='store',
          metavar='INCLUDE_TEST_FILES',
          dest='include_tests',
          help='Comma-seperated list of test source files to ignore when building')
AddOption('--check-costG',
          dest='costg',
          action='store_true',
          metavar='BUILD_COSTG_TESTS',
          help='Build costG tests')

## Source files (for lib)
lib_src_files = glob.glob(r"src/*.cpp")
lib_src_files += glob.glob(r"src/doris/*.cpp")
lib_src_files += glob.glob(r"src/ids/*.cpp")
lib_src_files += glob.glob(r"src/rinex/*.cpp")
lib_src_files += glob.glob(r"src/planets/*.cpp")
lib_src_files += glob.glob(r"src/gravity/*.cpp")
lib_src_files += glob.glob(r"src/satellite/*.cpp")
lib_src_files += glob.glob(r"src/atmosphere/*.cpp")
lib_src_files += glob.glob(r"src/atmosphere/dtm2020/*.cpp")
lib_src_files += glob.glob(r"src/astrodynamics/*.cpp")
lib_src_files += glob.glob(r"src/web/*.cpp")
lib_src_files += glob.glob(r"src/iers/*.cpp")
lib_src_files += glob.glob(r"src/integrators/*.cpp")
lib_src_files += glob.glob(r"src/relativity/*.cpp")
lib_src_files += glob.glob(r"src/satellites/*.cpp")
lib_src_files += glob.glob(r"src/beacon_tbl/*.cpp")
lib_src_files += glob.glob(r"src/var/*.cpp")
lib_src_files += glob.glob(r"src/filters/*.cpp")

## Headers (for lib)
hdr_src_files = glob.glob(r"src/*.hpp")

## Environments ...
denv = Environment(CXXFLAGS='-std=c++17 -g -pg -Wall -Wextra -Werror -pedantic -W -Wshadow -Wno-error=inline -Wno-class-memaccess -Wdisabled-optimization -DDEBUG -DDORIS_EXTRA_CHECKS -DUSE_OWN_ROTATION_COSTG')
## .. for checking the integrator ...
#penv = Environment(CXXFLAGS='-std=c++17 -Wall -Wextra -Werror -pedantic -W -Wshadow -Wno-error=inline -Wno-class-memaccess -O2 -march=native -DRANDOM_RFO -DINTEGRATOR_CHECK')
penv = Environment(CXXFLAGS='-std=c++17 -Wall -Wextra -Werror -pedantic -W -Wshadow -Wno-error=inline -Wno-class-memaccess -O2 -march=native -DUSE_OWN_ROTATION_COSTG')

## Command line arguments ...
debug = ARGUMENTS.get('debug', 0)
eigen = ARGUMENTS.get('eigen', 0)
branchless = ARGUMENTS.get('branchless', 1)

## Construct the build enviroment
env = denv.Clone() if int(debug) else penv.Clone()

## What compiler should we be using ?
if GetOption('cxx') is not None:
    env['CXX'] = GetOption('cxx')
    if env['CXX'] == "clang++" and '-Wno-class-memaccess' in env['CXXFLAGS']:
        ## LLVM Compiler
        print('Note: Removing flag \'-Wno-class-memaccess\' cause compiler is LLVM')
        env['CXXFLAGS'] = env['CXXFLAGS'].replace('-Wno-class-memaccess', '')
        env['CXXFLAGS'] = ' '.join([env['CXXFLAGS'], '-Wno-error=c++20-attribute-extensions'])

## Set the C++ standard
cxxstd = GetOption('std')
env.Append(CXXFLAGS=' --std=c++{}'.format(cxxstd))

## Various other compilation symobols, for debug builds ...
for key, value in ARGLIST:
    if key == 'count_kepler_iterations':
        env.Append(CXXFLAGS=' -DCOUNT_KEPLER_ITERATIONS')

## !! Warning !!
## This affects build of ggeodesy (header files)
if eigen:
    math_lib = ''
    env.Append(CXXFLAGS=' -DUSE_EIGEN')
#else:
#    math_lib = 'matvec'
#if branchless: env.Append(CXXFLAGS=' -DBRANCHLESS')

## (shared) library ...
vlib = env.SharedLibrary(source=lib_src_files, target=lib_name, CPPPATH=['src/'], SHLIBVERSION=lib_version)

## Build ....
env.Alias(target='install', source=env.Install(dir=os.path.join(prefix, 'include', inc_dir), source=hdr_src_files))
env.Alias(target='install', source=env.InstallVersionedLib(dir=os.path.join(prefix, 'lib'), source=vlib))

## Tests ...
test_sources = glob.glob(r"test/*.cpp")
ignore_test_list = [] if GetOption('ignore_tests') is None else GetOption('ignore_tests').split(',')
if len(ignore_test_list):
    for testscr in ignore_test_list:
        if testscr == "all":
            test_sources = [] 
        elif os.path.isfile(testscr):
            if testscr in test_sources: test_sources.remove(testscr)
include_test_list = [] if GetOption('include_tests') is None else GetOption('include_tests').split(',')
if len(include_test_list):
    test_sources = []
    for testscr in include_test_list:
        if os.path.isfile(testscr):
            test_sources.append(testscr)

env.Append(RPATH=root_dir)
for tsource in test_sources:
    ttarget = tsource.replace('_', '-').replace('.cpp', '.out')
    env.Program(target=ttarget, source=tsource, CPPPATH='src/', LIBS=vlib+['sp3', 'sinex', 'iers2010', 'geodesy', 'datetime', 'yaml-cpp', 'cspice.a', 'csupport', 'curl', 'sofa_c'], LIBPATH='.')

## Unit Tests (only build if user selected)
if ARGUMENTS.get('make-check', 0):
    print('>> Note: Building Unit Tests ...')
    tests_sources = glob.glob(r"unit_test/*.cpp")
    if 'RPATH' not in env or root_dir not in env['RPATH']:
      env.Append(RPATH=root_dir)
    for tsource in tests_sources:
        pth = os.path.dirname(tsource)
        bsn = os.path.basename(tsource)
        ttarget = os.path.join(pth, bsn.replace('_', '-').replace('.cpp', '.out'))
        env.Program(target=ttarget, source=tsource, CPPPATH='src/', LIBS=vlib+['sp3', 'sinex', 'iers2010', 'geodesy', 'datetime', 'yaml-cpp', 'cspice.a', 'csupport', 'curl'], LIBPATH='.')

if GetOption('costg'):
    print(">> Building cost-G test programs ...")
    tests_sources = [r"costg/costg_parsers.cpp"]
    tests_targets = glob.glob(r"costg/*.cpp")
    tests_targets = [ t for t in tests_targets if t not in tests_sources ]
    if 'RPATH' not in env or root_dir not in env['RPATH']:
      env.Append(RPATH=root_dir)
    # env['CXXFLAGS'] = ' '.join([env['CXXFLAGS'], '-DUSE_OWN_ROTATION_COSTG'])
    for prog in tests_targets:
        pth = os.path.dirname(prog)
        bsn = os.path.basename(prog)
        print('>> Building source {:}'.format(prog))
        ttarget = os.path.join(pth, bsn.replace('_', '-').replace('.cpp', '.out'))
        env.Program(target=ttarget, source=[prog]+tests_sources, CPPPATH='src/', LIBS=vlib+['sp3', 'sinex', 'iers2010', 'geodesy', 'datetime', 'yaml-cpp', 'cspice.a', 'csupport', 'curl', 'sofa_c'], LIBPATH='.')
