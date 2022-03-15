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
AddOption('--branchless',
          dest='branchls',
          action='store_true',
          help='Trigger built with BRANCHLESS defined',
          default=False)
AddOption('--ignore-tests',
          nargs=1,
          type='string',
          action='store',
          metavar='IGNORE_TEST_FILES',
          dest='ignore_tests',
          help='Comma-seperated list of test source files to ignore when building')

## Source files (for lib)
lib_src_files = glob.glob(r"src/*.cpp")
lib_src_files += glob.glob(r"src/doris/*.cpp")
lib_src_files += glob.glob(r"src/rinex/*.cpp")
lib_src_files += glob.glob(r"src/planets/*.cpp")
lib_src_files += glob.glob(r"src/gravity/*.cpp")
lib_src_files += glob.glob(r"src/satellite/*.cpp")
lib_src_files += glob.glob(r"src/atmosphere/*.cpp")
lib_src_files += glob.glob(r"src/web/*.cpp")
lib_src_files += glob.glob(r"src/iers/*.cpp")

## Headers (for lib)
hdr_src_files = glob.glob(r"src/*.hpp")

## Environments ...
denv = Environment(CXXFLAGS='-std=c++17 -g -pg -Wall -Wextra -Werror -pedantic -W -Wshadow -Winline -Wdisabled-optimization -DDEBUG')
penv = Environment(CXXFLAGS='-std=c++17 -Wall -Wextra -Werror -pedantic -W -Wshadow -Winline -O2 -march=native')

## Command line arguments ...
debug = ARGUMENTS.get('debug', 0)

## Construct the build enviroment
env = denv.Clone() if int(debug) else penv.Clone()

## What compiler should we be using ?
if GetOption('cxx') is not None: env['CXX'] = GetOption('cxx')

## Set the C++ standard
cxxstd = GetOption('std')
env.Append(CXXFLAGS=' --std=c++{}'.format(cxxstd))

## Get options from command line ...
if GetOption('branchls'): env.Append(CXXFLAGS=' -DBRANCHLESS')

## (shared) library ...
vlib = env.SharedLibrary(source=lib_src_files, target=lib_name, CPPPATH=['src/'], SHLIBVERSION=lib_version)

## Build ....
env.Alias(target='install', source=env.Install(dir=os.path.join(prefix, 'include', inc_dir), source=hdr_src_files))
env.Alias(target='install', source=env.InstallVersionedLib(dir=os.path.join(prefix, 'lib'), source=vlib))

## Tests ...
ignore_test_list = [] if GetOption('ignore_tests') is None else GetOption('ignore_tests').split(',')
print('Note: Ignore Test list: {:}'.format(ignore_test_list))
tests_sources = glob.glob(r"test/*.cpp")
env.Append(RPATH=root_dir)
for tsource in tests_sources:
  if tsource not in ignore_test_list:
    ttarget = tsource.replace('_', '-').replace('.cpp', '.out')
    env.Program(target=ttarget, source=tsource, CPPPATH='src/', LIBS=vlib+['sp3', 'sinex', 'iers2010', 'geodesy', 'datetime', 'cspice.a', 'csupport', 'curl'], LIBPATH='.')
