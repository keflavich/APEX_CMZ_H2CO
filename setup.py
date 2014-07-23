#!/usr/bin/env python

import sys
if 'build_sphinx' in sys.argv or 'develop' in sys.argv:
    from setuptools import setup, Command
else:
    from distutils.core import setup, Command

#with open('README.rst') as file:
#    long_description = file.read()
#
#with open('CHANGES') as file:
#    long_description += file.read()

# no versions yet from agpy import __version__ as version

class PyTest(Command):
    user_options = []
    def initialize_options(self):
        pass
    def finalize_options(self):
        pass
    def run(self):
        import sys,subprocess
        errno = subprocess.call([sys.executable, 'runtests.py'])
        raise SystemExit(errno)

#execfile('image_registration/version.py')

setup(name='apex_shfi_pipeline',
      version='0.0',
      description='SHFI OTF pipeline',
      long_description='SHFI OTF pipeline',
      author='Adam Ginsburg',
      author_email='adam.g.ginsburg@gmail.com',
      url='https://github.com/keflavich/APEX_CMZ_H2CO',
      packages=['shfi_otf_pipeline',], 
      cmdclass = {'test': PyTest},
     )

