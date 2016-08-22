# -*- coding: utf-8 -*-

from setuptools import setup, find_packages


with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='aesop',
    version='0.0.2',
    description='Module for analyzing electrostatics with protein structures',
    long_description=readme,
    author='Reed Harrison, Rohith Mohan',
    author_email='reed.harrison@email.ucr.edu, rohith.mohan@email.ucr.edu',
    url='https://github.com/rohithmohan/aesop-python',
    license=license,
    packages=find_packages(exclude=('tests', 'docs')),
    install_requires=['prody', 'gridDataFormats', 'scipy', 'numpy', 'python-dateutil']#,
#    zip_safe=False	
)

