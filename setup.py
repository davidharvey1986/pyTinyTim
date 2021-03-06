import sys,os,string,glob,subprocess

from setuptools import setup,Extension
from setuptools.command.build_ext import build_ext
from setuptools.command.install import install

import numpy

long_description = """\
This suite of scripts uses tinytim to produce a pixellised version of 
the ACS PSF
"""
#python setup.py register -r pypi
#sudo python setup.py sdist upload -r pypi

version='0.0.1'
packages = ['pyTinyTim','IDL','cCode']
package_dir =  {'pyTinyTim':'src',
               'IDL':'./lib/IDL',
               'cCode':'./lib/cCode'}
package_data=  {'IDL':'lib/IDL/*',
               'cCode':'lib/cCode/*'}
    
    
INCDIRS=['.']

setup   (       name            = "pyTinyTim",
                version         = version,
                author          = "David Harvey",
                author_email    = "harvey@lorentz.leidenuniv.nl",
                description     = "pyTinyTim module",
                packages        = packages,
                package_dir     = package_dir,
                package_data     = package_data,
                license         = 'MIT',
                url = 'https://github.com/davidharvey1986/pyTinyTim', # use the URL to the github repo
                download_url = 'https://github.com/davidharvey1986/pyTinyTim/archive/'+version+'.tar.gz',
                install_requires=['pyRRG']
        )


