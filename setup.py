from setuptools import setup

setup(
    name='araproc',
    url='https://github.com/clark2688/araproc',
    author='Brian Clark',
    author_email='baclark@umd.edu',
    packages=['araproc'],
    package_data = {'araproc': ['analysis/data/*']},
    include_package_data=True,
    python_requires='>= 3.7',
    install_requires=['numpy', 'scipy', 'matplotlib', 'yaml'],
    version='0.1',
    license='GPLv3',
    description='ARA data analysis framework',
)