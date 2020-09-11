import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

# version = {}
# with open("seaplan/version.py") as fp:
#     exec(fp.read(), version)

setuptools.setup(
    name="crawtools",
    version='0.0',
    # version=version['__version__'],
    author="Wayne Crawford",
    author_email="crawford@ipgp.fr",
    description="Various python modules",
    long_description=long_description,
    long_description_content_type="text/x-rst; charset=UTF-8",
    # url="https://github.com/WayneCrawford/seaplan",
    packages=setuptools.find_packages(),
    include_package_data=True,
    install_requires=[
          'cartopy==0.18',
          'numpy>=1.17',
          'matplotlib>=3.0',
          'obspy>=1.0'
      ],
    # entry_points={
    #      'console_scripts': [
    #          'seaplan=seaplan.sea_plan:main',
    #          'seaplan-validate=seaplan.validate_json:_console_script'
    #      ]
    # },
    python_requires='>=3.6',
    classifiers=(
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Physics"
    ),
    keywords='software'
)
