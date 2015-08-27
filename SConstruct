import os

srcdir = 'src'
objdir = '.obj'
bindir = 'bin'
CacheDir('.cache')

env = Environment()
env.VariantDir('.obj', 'src', duplicate=0)


try:
    from colorizer import colorizer
    colorizer().colorize(env)
except ImportError:
    pass

env["ENV"]["TERM"] = os.environ["TERM"] #to get clang to display colored output
env.Append(CCFLAGS=Split("""
    -march=native
    -Wall
    -std=c++11
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

extent = int(open('extent').read())
with open("src/extent.h", "w") as f:
    f.write("const unsigned L = {};\n".format(extent))

for line in open('.bits').readlines():
    (beta,bits) = line.split()
    Command(".h/{}/random.h".format(bits), "random", "./random {} > $TARGET".format(bits))
    o = env.Object(
        "{}/ising/{}_{}.o".format(objdir,extent,bits), 
        objdir+'/ising.cpp', 
        CPPPATH='.h/{}'.format(bits))
    env.Program("{}/ising/{}/{}".format(bindir,extent,beta), o)

env.Program(
    bindir+"/analyze", 
    objdir+"/analyze.cpp",
    CPPPATH="/usr/lib/gcc/x86_64-unknown-linux-gnu/5.2.0/include",
    LIBPATH=["/usr/lib/x86_64-linux-gnu"], 
    LIBS=["boost_filesystem","boost_system"]
    )
results = []
for datafile in Glob('data/*/*', strings=True):
    result = 'result'+datafile[4:]
    Command(result, [datafile, 'bin/analyze'], 'bin/analyze 100 1000 $SOURCE > $TARGET')
    results.append(result)

env["ENV"]["PATH"]+= ":~/bin"
env.Command(["M.pdf", 'chi.pdf'], [results, "./plot"], "./plot")

