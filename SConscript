import os

Import('debug', 'compiler')

# Initialize environment
env = Environment(CC = compiler, 
                  CFLAGS = '-std=c99 -fopenmp -D_GNU_SOURCE ',
                  CPPPATH = [os.getcwd(), ARGUMENTS.get('commons-path', os.getcwd() + '/../common-libs/'), '#', '.' ],
                  LIBPATH = ['/usr/lib' ])
env.Decider('MD5-timestamp')

if debug == 1:
    env['CFLAGS'] += ' -O0 -g'
else:
    env['CFLAGS'] += ' -O3'

env['objects'] = []

# Targets
SConscript(['bioformats/SConscript',], exports = ['env'])
SConscript(['aligners/SConscript',], exports = ['env'])

env.Library('bioinfo', env['objects'])

