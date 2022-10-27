from setuptools import setup
import os

this_directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(this_directory, 'README.md')) as f:
    long_description = f.read()



console_scripts = '''
    [console_scripts]
    cinful = cinful.cinful.cinful:main
    cinful_init = cinful.cinful.cinful_initialize:main
'''

setup(
    name= 'cinful',
    version= '1.2.0',
    author= 'Aaron Feller',
    author_email = 'aaronleefeller@gmail.com',
    include_package_data = True,
    url='https://github.com/wilkelab/cinful',
    description=('cinful: A fully automated pipeline to identify microcins'
                 'with associated immunity proteins and export machinery'),
    long_description=long_description,
    long_description_content_type='text/markdown',
    license='GPL-3',
    py_modules=['cinful', 'cinful_init'],
    entry_points=console_scripts,
    packages=['cinful'],
)
