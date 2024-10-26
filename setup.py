from setuptools import setup

setup(
    name="araproc",
    url="https://github.com/clark2688/araproc",
    author="Brian Clark",
    author_email="baclark@umd.edu",
    packages=["araproc"],
    package_data = {"araproc": ["analysis/data/*", "framework/config_files/*"]},
    include_package_data=True,
    python_requires=">= 3.9",
    install_requires=["numpy < 2", "scipy >= 1.13.0", 
                      "matplotlib", "pyyaml", "ruff",
                      "scikit-learn", "pandas",
                      ],
    version="0.1",
    license="GPLv3",
    description="ARA data analysis framework",
)