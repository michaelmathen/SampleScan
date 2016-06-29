import os
# Let's define a common build environment first...
common_env = Environment()
common_env.Append(CPPDEFINES={'VERSION': 1})
common_env.Append(CXXFLAGS=['-std=c++11', '-fopenmp'])
common_env.Append(LIBPATH=['/usr/local/lib/', '/home/mmath/lib'])
common_env.Append(CPPPATH=["/usr/local/include", '/home/mmath/include'])
common_env.Append(CPPPATH=["/usr/include/python2.7"])
#common_env.Append( F90FLAGS  = ' -openmp' )
common_env.Append( LINKFLAGS = ' -fopenmp' )

parent_dir = os.getcwd()

release_env = common_env.Clone()
release_env.Append(CPPDEFINES=['RELEASE'])
release_env.Append(CXXFLAGS=['-O3', '-g'])
release_env.VariantDir('build/release', 'src')
release_env.Append(CPPPATH=[os.getcwd() + '/build/release'])

debug_env = common_env.Clone()
debug_env.Append(CPPDEFINES=['DEBUG'])
debug_env.VariantDir('build/debug', 'src')
debug_env.Append(CPPPATH=[os.getcwd() + '/build/debug'])
debug_env.Append(CCFLAGS=['-g','-O2', '-pg'])

# Now that all build environment have been defined, let's iterate over
# them and invoke the lower level SConscript files.
for mode, env in dict(release=release_env, debug=debug_env).iteritems():
    env.SConscript('build/%s/SConscript' % mode, {'env': env})
