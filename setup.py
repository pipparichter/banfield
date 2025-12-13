import setuptools
import os


setuptools.setup(
    name='banfield',
    version='0.1',    
    description='N/A',
    author='Philippa Richter',
    author_email='prichter@berkeley.edu',
    packages=['src', 'src.files'], 
    include_package_data=True, 
    package_data={'src':['data/*']})


# TODO: What exactly is an entry point?
# https://python-packaging.readthedocs.io/en/latest/command-line-scripts.html
#  