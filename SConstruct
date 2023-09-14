import os

srcdir = 'src'
objdir = '.obj'
bindir = 'bin'
CacheDir('.cache')

env = Environment(
    CXX="clang++",
    CCFLAGS="-O3 -Wall",
)
env.VariantDir(objdir, srcdir, duplicate=0)

for i in ["PATH", "TERM"]:
    env["ENV"][i] = os.environ[i]

try:
    from colorizer import colorizer
    colorizer().colorize(env)
except ImportError:
    pass

env.Program(
    bindir+"/analyze", 
    objdir+"/analyze.cpp",
    LIBS=["boost_filesystem","boost_system"],
)

nc = env.Clone()
nc.CacheDir(None)
nc.Decider('MD5-timestamp')
results = Glob('result/*/*', strings=True)
if os.path.exists('data'):
    for datafile in Glob('data/*/*', strings=True):
        results.append(
            nc.Command(
                'result'+datafile[4:], 
                [datafile, 'bin/analyze'], 
                'bin/analyze 100 0.1 $SOURCE > $TARGET',
            )
        )

#env["ENV"]["PATH"]+= ":~/bin"
#env.Command(["plot/mag.pdf", 'plot/chi.pdf'], ["sc/plot"]+results, "sc/plot")
