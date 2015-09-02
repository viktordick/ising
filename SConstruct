import os

srcdir = 'src'
objdir = '.obj'
bindir = 'bin'
CacheDir('.cache')

env = Environment()
env.VariantDir('.obj', 'src', duplicate=0)

for i in ["PATH", "TERM"]:
    env["ENV"][i] = os.environ[i]


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

if os.path.exists('extent') and os.path.exists('.bits'):
    extent = int(open('extent').read())
    env.Append(CPPDEFINES={'L': extent})

    for line in open('.bits').readlines():
        (beta,bits) = line.split()
        env.Command(".h/{0}/random.h".format(bits), "random", "./random {0} > $TARGET".format(bits))
        o = env.Object(
            "{0}/ising/{1}_{2}.o".format(objdir,extent,bits), 
            objdir+'/ising.cpp', 
            CPPPATH='.h/{0}'.format(bits))
        env.Program("{0}/ising/{1}/{2}".format(bindir,extent,beta), o)

env.Program(
    bindir+"/analyze", 
    objdir+"/analyze.cpp",
    CPPPATH="/usr/lib/gcc/x86_64-unknown-linux-gnu/5.2.0/include",
    LIBPATH=["/usr/lib/x86_64-linux-gnu"], 
    LIBS=["boost_filesystem","boost_system"]
    )
if os.path.exists('data'):
    results = []
    for datafile in Glob('data/*/*', strings=True):
        result = 'result'+datafile[4:]
        env.Command(result, [datafile, 'bin/analyze'], 'bin/analyze 100 1000 $SOURCE > $TARGET')
        results.append(result)

    env["ENV"]["PATH"]+= ":~/bin"
    env.Command(["M.pdf", 'chi.pdf'], [results, "./plot"], "./plot")

