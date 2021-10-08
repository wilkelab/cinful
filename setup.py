from setuptools import setup
import os

this_directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(this_directory, 'README.md')) as f:
    long_description = f.read()

setup(
    name='cinful',
	include_package_data = True,
	version='0.1.0',
	author='T. Jeffrey Cole',
    author_email = "jffcole7@gmail.com",
    url='https://github.com/tijeco/cinful',
	description=('cinful: A fully automated pipeline to identify microcins'
        'along with their associated immunity proteins and export machinery'),
	long_description=long_description,
	long_description_content_type='text/markdown',
	license='GPL-3',
    py_modules=['cinful'],
    entry_points='''
        [console_scripts]
        cinful=cinful.cinful:main
    ''',
    packages=['cinful'],
)