import os.path
import Options

def get_hg_version() :
    from mercurial import ui, hg
    r = hg.repository(ui.ui(),".")
    ctx = r.parents()[0]
    return "0.1-r" + str(ctx.rev())


APPNAME = "FGMM"
try :
    VERSION = get_hg_version()
except :
    VERSION = "0.1a"

srcdir = '.'
blddir = '_build_'


def set_options(opt):
    # command-line options provided by a waf tool
    opt.tool_options('compiler_cc')
    opt.tool_options('compiler_cxx')
    opt.add_option('--gprof', dest = "use_gprof", action='store_true', 
                   default=False, help='Link with gprof for profiling')

    opt.add_option('--old',dest="oldGMR",action='store_true',
                   default=False, help='Build benchmarks using old GMR lib for comparison')

    opt.add_option('--debug',dest="debug",action="store_true",
                   default=False, help='Build debug symbols')

    opt.add_option('--no-test',dest="test",action="store_false",
                   default=True,help="disable building tests (type ./waf test to run them)")

    opt.add_option('--python',dest="python",action="store_true",
                   default=False,help="build python bindings")

    opt.add_option('--matlab',dest="matlab",action="store_true",
                   default=False,help="build matlab mex files")

    opt.add_option('--octave',dest="octave",action="store_true",
                   default=False,help="build octave mex files")

#     opt.add_option('--check',dest="check",action='store_true',
#                    default=False, help='run the tests')

def configure(conf) :
    #print "configuring"
    conf.check_tool('compiler_cc')  
    if Options.options.debug :
        conf.env.CCFLAGS = conf.env.CCFLAGS_DEBUG
        if conf.env.CC_NAME == "gcc" :
            conf.env.CCFLAGS.append("-Wall")
    else :
        if conf.env.CC_NAME ==  "gcc" :
            conf.env.CCFLAGS = ['-Wall','-O2','--fast-math']
            conf.env.CPPFLAGS = ['-DNDEBUG'] # no more asserts.. 
    if Options.options.use_gprof :
        conf.env.CCFLAGS.append('-pg')
        conf.env.LINKFLAGS = ['-pg']
        
    if Options.options.test :
        conf.env.build_test = True 
        conf.recurse("tests")

    if Options.options.python :
        conf.env.python = True
        conf.recurse("python")
    
    if Options.options.matlab or Options.options.octave :
        conf.recurse("matlab")

    conf.check_cc(lib='m', uselib_store='M') # looking for libm (no need with msvc) 

def build(bld) :
    bld.recurse("src")

    if bld.env.build_test :
        bld.recurse("tests")

    if bld.env.python :
        bld.recurse("python")

    if bld.env.matlab or bld.env.octave :
        bld.recurse("matlab")

def doc(ctx) :
    subprocess.call(["doxygen","Doxfile"])
    subprocess.call(["sphinx-build","-b","html","docsrc", "doc/sphinx/"])

import Scripting, Build, Environment,Utils
import os,sys,subprocess

def test(ctx) :
    # getting the project context 
    try:
        proj=Environment.Environment(Options.lockfile)
    except IOError:
        raise Utils.WafError("Project not configured (run 'waf configure' first)")

    bld = Build.BuildContext()
    bld.load_dirs(proj['srcdir'], proj['blddir']) 
    bld.load_envs()
    
    pybin = sys.executable  # hey you should have a python exe to run this .. 
    for v in bld.lst_variants :
        trunner = subprocess.Popen([pybin,os.path.join(proj['srcdir'],"tests",'run_test.py')],
                                   cwd = os.path.join(bld.bdir,v,"tests"))
        if(trunner.wait() != 0) :
            raise Utils.WafError("test failed " )
        
