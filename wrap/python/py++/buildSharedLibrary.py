"""The script should be run with scons"""
import os, sys
import SCons.Util
from distutils.sysconfig import get_python_inc, get_python_lib

def getSMACFiles(src = '.'):
    source_files = filter( lambda s: s.endswith( '.cpp' ), os.listdir(src) )
    return map( lambda fname: os.path.join(src, fname ), source_files )

env = Environment()

# Take the environmental variables into account
if os.environ.has_key('CXX'):
    env['CXX'] = os.environ['CXX']
if os.environ.has_key('CXXFLAGS'):
    env['CXXFLAGS'] += Scons.Util.CLVar(os.environ["CXXFLAGS"])
if os.environ.has_key('LDFLAGS'):
    env['LINKFLAGS'] += Scons.Util.CLVar(os.environ['LDFLAGS'])

# Include python into the library
env.Append(FRAMEWORK_SEARCH_PATH=['/opt/local/*','/System/'])
env.Append(CPPPATH=[get_python_inc(plat_specific = True), '/opt/local/include'])
#env['FRAMEWORKS'] += ['Python']
env['LINKFLAGS'] += ['-lpython2.6']

# Get the argument from SCONS
libs = [ARGUMENTS.get('BOOSTLIB', 'boost_python'), 'smac']
print libs
env.Append(LIBS = libs)

# Build the shared library libsmac
cpp_sources = getSMACFiles()
smac = env.SharedLibrary('smac', cpp_sources, LIBPATH=['../../../lib/','/opt/local/lib/'],SHLIBSUFFIX='.so')
env.Alias('lib', smac)

# Setup the install into the python trunk
install_dir = get_python_lib(plat_specific=True)

# Install it!
env.Install(install_dir, [smac])
modules = ['__init__', 'formats', 'match', 'matchlib']
env.Install(os.path.join(install_dir, 'pysmac'),
            ['pysmac/%s.py' % module for module in modules])
env.Alias('install', install_dir)

#SharedLibrary(target = '',
#              source = [getSMACFiles()],
#              LIBS = ['boost_python'],
#              LIBPATH = [BOOST_LIB_PATH, '/Library/Frameworks/Python.framework/Versions/2.6/Headers/', '/opt/local/lib/', '.'],
#              CPPPATH = ['~/Downloads/boost_1_45_0', '.', '../../../src/', '/Library/Frameworks/Python.framework/Versions/2.6/Headers/'],
#              CCFLAGS = [],
#              SHLIBPREFIX = 'pysmac',
#              SHLIBSUFFIX = 'so'
#              )
