import os.path
import Options

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

#     opt.add_option('--check',dest="check",action='store_true',
#                    default=False, help='run the tests')

def configure(conf) :
    #print "configuring"
   
    conf.check_tool('compiler_cc')  
    conf.check_tool('compiler_cxx')  
    if Options.options.debug :
        conf.env.CCFLAGS = ['-Wall','-g'] 
    else :
        conf.env.CCFLAGS = ['-Wall','-O3','--fast-math']
        conf.env.CPPFLAGS = ['-DNDEBUG'] # no more asserts.. 

    if Options.options.use_gprof :
        conf.env.CCFLAGS.append('-pg')
        conf.env.LINKFLAGS = ['-pg']
        
    if Options.options.test :
        conf.env.build_test = True

    if Options.options.oldGMR :
        conf.env.CCXFlAGS = conf.env.CCFLAGS
        conf.env['build_old_gmr'] = True
        conf.env['LIBPATH_MATRIX'] = '/home/fdhalluin/code/MathLib/lib/'
        conf.env['LIB_MATRIX'] = 'Matrix'
        conf.env['LIBPATH_GMR'] = '/home/fdhalluin/code/gmr/lib/'
        conf.env['LIB_GMR'] = 'GMR'
        conf.env['CPPPATH_GMR'] = '/home/fdhalluin/code/gmr/include/'
        conf.env['CPPPATH_MATRIX'] = '/home/fdhalluin/code/MathLib/include/'
        

   
def build(bld) :
    obj_src = bld.path.ant_glob("src/*.c")

    bld(features='cc cstaticlib',
        target = 'fgmm',
        source = obj_src,
        export_incdirs="src/",
        install_path='${PREFIX}/lib') # <-- UGLY NON PORTABLE STUFF 

    bld.install_files('${PREFIX}/include/',  # <-- again ! 
                      ["src/fgmm.h","src/fgmm++.hpp"])

    if bld.env.build_test :
        test_src = bld.path.ant_glob("tests/test_*.c")
        for src in test_src.split():

            targ_name = os.path.basename(src)
            targ_name = os.path.splitext(targ_name)[0]

            bld(features = 'cc cprogram',
                source = src,
                target = targ_name,
                includes = '.',
                lib = ['m'],
                uselib_local="fgmm",
                install_path=False)
            #add_objects=obj)
        bld(target="run_test.py",
            source = "tests/run_test.py",
            rule= "cp ${SRC} ${TGT}")

    if bld.env.build_old_gmr :
        cpgmr = bld(features = 'cxx cprogram',
                    source = 'tests/oldGMR.cpp',
                    target = 'oldGMR',
                    uselib = ['GMR','MATRIX'],
                    install=False)    
        cpgmr = bld(features = 'cxx cprogram',
                    source = 'tests/bench.cpp',
                    target = 'bench',
                    uselib = ['GMR','MATRIX'],
                    uselib_local="fgmm",
                    install=False)  
  
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
        trunner = subprocess.Popen([pybin,os.path.join(bld.bdir,v,'run_test.py')],
                                   cwd = os.path.join(bld.bdir,v))
        trunner.wait()
        
