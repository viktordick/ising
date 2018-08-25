import os
from socket import gethostname

srcdir = 'src'
objdir = '.obj'
bindir = 'bin'
CacheDir('.cache')

env = Environment()
env.VariantDir('.obj', 'src', duplicate=0)

for i in ["PATH", "TERM"]:
    env["ENV"][i] = os.environ[i]


if ARGUMENTS.get('color', 1) == 1:
    try:
        from colorizer import colorizer
        colorizer().colorize(env)
    except ImportError:
        pass

env.Append(CCFLAGS=Split("""
    -march=native
    -Wall
    -Wno-sign-compare
    """) + (["-g","-O0"] if 'debug' in ARGUMENTS else ["-O3"]))


COMPILER = ARGUMENTS.get("cxx", "gcc") #gcc is default

if COMPILER in ["llvm", "clang"]:
    env.Replace(CXX="clang++")
elif COMPILER == "gcc":
    env.Append(CCFLAGS=Split("""
    -fcx-fortran-rules
    -funsafe-math-optimizations
    -finline-limit=800
    -funswitch-loops
    -funsafe-loop-optimizations
    -floop-interchange
    -floop-strip-mine
    -floop-block
    -ftree-loop-distribution
    -ftree-loop-im
    -ftree-loop-linear
    -ftree-loop-ivcanon
    -ftree-vrp
    -fgcse-after-reload
    -fgcse-las
    -fgcse-sm
    -fprefetch-loop-arrays
    """))
else:
    print("cxx argument not recognized.")
    Exit()

if gethostname() == 'gdev1':
    env.Append(CPPDEFINES=['GDEV'])
else:
    env.Append(CXXFLAGS=['-std=c++11'])
try:
    pmin = ARGUMENTS['pmin']
    pmax = ARGUMENTS['pmax']
    bits = ARGUMENTS['bits']
    L = ARGUMENTS['L']
    exe = '{}-{}-{}'.format(pmin, pmax, bits)
    env.Append(CPPDEFINES={'L': L})

    inc=".h/{0}".format(exe)
    h = env.Command(inc+"/random.h", "sc/random", "sc/bitrange {} {} {} | sc/random > $TARGET".format(pmin, pmax, bits))
    o = env.Object("{0}/ising/{1}/{2}.o".format(objdir, L, exe),
        objdir+"/ising.cpp", CPPPATH=inc)
    env.Depends(o,h)
    env.Program("{0}/ising/{1}/{2}".format(bindir,L,exe),o)
except KeyError:
    pass

env.Program(
    bindir+"/analyze", 
    objdir+"/analyze.cpp",
    CPPPATH="/usr/lib/gcc/x86_64-unknown-linux-gnu/5.2.0/include",
    LIBPATH=["/usr/lib/x86_64-linux-gnu"], 
    LIBS=["boost_filesystem","boost_system"]
    )
nc = env.Clone()
nc.CacheDir(None)
nc.Decider('MD5-timestamp')
results = Glob('result/*/*', strings=True)
if os.path.exists('data'):
    for datafile in Glob('data/*/*', strings=True):
        results.append(nc.Command('result'+datafile[4:], 
            [datafile, 'bin/analyze'], 
            'bin/analyze 100 0.01 $SOURCE > $TARGET'))

env["ENV"]["PATH"]+= ":~/bin"
env.Command(["plot/mag.pdf", 'plot/chi.pdf'], ["sc/plot"]+results, "sc/plot")

