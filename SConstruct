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
env.Append(CCFLAGS=( ["-g","-O0"] if 'debug' in ARGUMENTS else ["-O3"]))
cpppath=["/usr/lib/gcc/x86_64-unknown-linux-gnu/5.2.0/include"]
env.Append(
    CCFLAGS=["-march=native", "-Wall", "-std=c++14"],
    LIBPATH=["/usr/lib/x86_64-linux-gnu"], 
    LIBS=["boost_filesystem","boost_system"],
    )

COMPILER = ARGUMENTS.get("cxx", "gcc") #gcc is default

if COMPILER in ["llvm", "clang"]:
    env.Replace(CXX="clang++")
    env.Append(CCFLAGS=Split("""
    -Wno-sign-compare
    """))
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
    -Wno-sign-compare
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
        "{}/ising_{}_{}.o".format(objdir,extent,bits), 
        objdir+'/ising.cpp', 
        CPPPATH=cpppath+['.h/{}'.format(bits)])
    env.Program("{}/ising-{}-{}".format(bindir,extent,beta), o)

env.Program(bindir+"/analyse-ising", objdir+"/analyse-ising.cpp")
# for s in env.Glob(objdir+"/*.cpp", strings=True):
#     env.Program(bindir+s[len(objdir):-4], s)

if os.path.isdir('data'):
    env["ENV"]["PATH"]+= ":~/bin"
    results = []
    for L in os.listdir('data'):
        results.append(Command("result/"+L+".txt", Glob("data/"+L+"/*")+["bin/analyse-ising"], "bin/analyse-ising "+L+" 8000 40"))
    env.Command(["M.pdf", 'chi.pdf'], [results, "./plot"], "./plot")

