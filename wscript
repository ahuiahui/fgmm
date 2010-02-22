import os.path
srcdir = '.'
blddir = '_build_'



def set_options(opt):
    # command-line options provided by a waf tool
    opt.tool_options('compiler_cc')

def configure(conf) :
    #print "configuring"
    conf.check_tool('compiler_cc')  
    conf.env.CCFLAGS = ['-Wall','-O2','-g','--fast-math']
   
def build(bld) :
    print "compiling"
    obj_src = bld.path.ant_glob("*.c")
    test_src = bld.path.ant_glob("tests/test_*.c")
    obj = []
    for src in obj_src.split():
        print src
        targ_name = os.path.splitext(src)[0]+'.o'
        bld(features='cc',
            target = targ_name,
            source=src)
        obj.append(targ_name)

    for src in test_src.split():

        targ_name = os.path.basename(src)
        targ_name = os.path.splitext(targ_name)[0]

        bld(features = 'cc cprogram',
            source = src,
            target = targ_name,
            includes = '.',
            lib = ['m'],
            add_objects=obj)
