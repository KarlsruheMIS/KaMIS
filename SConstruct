import os
import platform
import sys

# Get the current plattform.
SYSTEM = platform.uname()[0]
HOST = platform.uname()[1]

# Get shortcut to $HOME.
HOME = os.environ['HOME']

def GetEnvironment():
    """ Get environment variables from command line and environment

    Returns
        Environment with the configuration from the command line.
    """

    opts = Variables()
    opts.Add('program', 'program or interface to compile', 'redumis, graph_checker, library')
    opts.Add('variant', 'the variant to build, optimized or optimized with output', 'optimized, optimized_output, debug')

    env = Environment(options=opts, ENV=os.environ)

    if not env['variant'] in ['optimized', 'optimized_output', 'debug']:
        print 'Illegal value for variant: %s' % env['variant']
        sys.exit(1)

    if not env['program'] in ['redumis', 'graph_checker', 'library']:
        print 'Illegal value for program: %s' % env['program']
        sys.exit(1)

    if platform.architecture()[0] == '64bit':
        env.Append(CPPFLAGS=['-DPOINTER64=1'])

    return env

# Get the common environment
env = GetEnvironment()

# Setup library and cpp paths
env.Append(CPPPATH = [ './lib' ])
env.Append(CPPPATH = [ './lib/mis' ])
env.Append(CPPPATH = [ './lib/mis/initial_mis' ])
env.Append(CPPPATH = [ './lib/mis/evolutionary' ])
env.Append(CPPPATH = [ './lib/mis/evolutionary/combine' ])
env.Append(CPPPATH = [ './lib/mis/hopcroft' ])
env.Append(CPPPATH = [ './lib/mis/kernel' ])
env.Append(CPPPATH = [ './lib/mis/ils' ])

env.Append(CPPPATH = [ './app' ])
env.Append(CPPPATH = [ './lib/io' ])
env.Append(CPPPATH = [ './interface' ])
env.Append(CPPPATH = [ './lib/algorithms' ])
env.Append(CPPPATH = [ './lib/data_structure' ])
env.Append(CPPPATH = [ './lib/data_structure/priority_queues' ])
env.Append(CPPPATH = [ './lib/data_structure/matrix' ])
env.Append(CPPPATH = [ './lib/partition' ])
env.Append(CPPPATH = [ './lib/partition/uncoarsening/refinement/quotient_graph_refinement/flow_refinement/' ])
env.Append(CPPPATH = [ './lib/parallel_mh' ])
env.Append(CPPPATH = [ './lib/parallel_mh/galinier_combine' ])
env.Append(CPPPATH = [ './lib/tools' ])

env.Append(CPPPATH = [ '../app' ])
env.Append(CPPPATH = [ '../lib' ])
env.Append(CPPPATH = [ '../lib/io' ])
env.Append(CPPPATH = [ '../lib/algorithms' ])
env.Append(CPPPATH = [ '../lib/data_structure' ])
env.Append(CPPPATH = [ '../lib/data_structure/matrix' ])
env.Append(CPPPATH = [ '../lib/data_structure/priority_queues' ])
env.Append(CPPPATH = [ '../lib/partition' ])
env.Append(CPPPATH = [ '../lib/partition/uncoarsening/refinement/quotient_graph_refinement/flow_refinement/' ])
env.Append(CPPPATH = [ '../lib/parallel_mh' ])
env.Append(CPPPATH = [ '../lib/parallel_mh/galinier_combine' ])
env.Append(CPPPATH = [ '../lib/tools' ])

env.Append(CPPPATH = [ './extern/KaHIP' ])
env.Append(CPPPATH = [ './extern/KaHIP/app' ])
env.Append(CPPPATH = [ './extern/KaHIP/interface' ])
env.Append(CPPPATH = [ './extern/KaHIP/lib' ])
env.Append(CPPPATH = [ './extern/KaHIP/lib/algorithms' ])
env.Append(CPPPATH = [ './extern/KaHIP/lib/data_structure' ])
env.Append(CPPPATH = [ './extern/KaHIP/lib/data_structure/priority_queues' ])
env.Append(CPPPATH = [ './extern/KaHIP/lib/data_structure/matrix' ])
env.Append(CPPPATH = [ './extern/KaHIP/lib/io' ])
env.Append(CPPPATH = [ './extern/KaHIP/lib/parallel_mh' ])
env.Append(CPPPATH = [ './extern/KaHIP/lib/parallel_mh/galinier_combine' ])
env.Append(CPPPATH = [ './extern/KaHIP/lib/partition' ])
env.Append(CPPPATH = [ './extern/KaHIP/lib/tools' ])
env.Append(CPPPATH = [ './extern/KaHIP/lib/partition/uncoarsening/refinement/quotient_graph_refinement/flow_refinement/' ])

env.Append(LIBPATH = [ '../extern/KaHIP' ])
env.Append(LIBPATH = [ '../extern/KaHIP' ])
env.Append(CPPPATH = [ './extern/argtable-2.10/include' ])
env.Append(CPPPATH = [ '../extern/argtable-2.10/include' ])

if not SYSTEM == 'Darwin':
    env.Append(LIBPATH = [ './extern/argtable-2.10/lib' ])
    env.Append(LIBPATH = [ '../extern/argtable-2.10/lib' ])
    env.Append(LIBPATH = [ '../../extern/argtable-2.10/lib' ])

conf = Configure(env)

# Set compiler flags
env.Append(CXXFLAGS = '-fopenmp')
env.Append(CXXFLAGS = '-DNDEBUG -Wall -funroll-loops -fno-stack-limit -O3 -std=c++0x')
env.Append(CCFLAGS  = '-DNDEBUG -Wall -funroll-loops -fno-stack-limit -O3 -std=c++0x')

# Execute the SConscript.
SConscript('SConscript', exports=['env'],variant_dir=env['variant'], duplicate=False)
