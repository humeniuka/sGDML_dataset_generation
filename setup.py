#!/usr/bin/env python
import os
import re
from io import open
from setuptools import setup, find_packages


def get_property(property, package):
    result = re.search(
        r'{}\s*=\s*[\'"]([^\'"]*)[\'"]'.format(property),
        open(package + '/__init__.py').read(),
    )
    return result.group(1)


from os import path

this_dir = path.abspath(path.dirname(__file__))
with open(path.join(this_dir, 'README.rst'), encoding='utf8') as f:
    long_description = f.read()

# Scripts
scripts = []
for dirname, dirnames, filenames in os.walk('scripts'):
    for filename in filenames:
        if filename.endswith('.py') or filename.endswith('.sh'):
            scripts.append(os.path.join(dirname, filename))

setup(
    name='sgdml_dataset_generation',
    version=get_property('__version__', 'sgdml_dataset_generation'),
    description='Parallelize computation of forces using QChem on SLURM queue',
    long_description=long_description,
    long_description_content_type='text/markdown',
    classifiers=[
        'Development Status :: 1 - Beta',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Education',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3.7',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Physics',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ],
    url='https://github.com/humeniuka/sGDML_dataset_generation',
    author='Alexander Humeniuk',
    author_email='alexander.humeniuk@gmail.com',
    license='LICENSE.txt',
    packages=find_packages(),
    install_requires=['numpy', 'ase', 'torchani'],
    scripts=scripts,
    include_package_data=True,
    zip_safe=False,
)
