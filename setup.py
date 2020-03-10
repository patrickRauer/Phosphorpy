from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='Phosphorpy',
    version='0.7.3',
    python_requires='>=3.6',
    packages=['Phosphorpy', 'Phosphorpy.data',
              'Phosphorpy.data.sub',
              'Phosphorpy.data.sub.plots', 'Phosphorpy.data.sub.interactive_plotting',
              'Phosphorpy.config',
              'Phosphorpy.fitting', 'Phosphorpy.external',
              'Phosphorpy.core', 'Phosphorpy.data.sub.tables',
              'Phosphorpy.local'],
    # package_dir={'': 'Phosphorpy'},
    include_package_data=True,
    # package_data={
    #       'config': ['Phosphorpy/local/survey.conf']
    #   },
    install_requires=['hypothesis', 'seaborn', 'numpy', 'astropy', 'pandas', 'astroquery',
                      'numba', 'scikit-learn', 'armapy', 'requests', 'scipy', 'pyarrow'],
    url='https://github.com/patrickRauer/Phosphorpy',
    license='GPL',
    author='Patrick Rauer',
    author_email='j.p.rauer@sron.nl',
    description='''
        Phosphorpy is python package to mine large photometric sky surveys. 
        It is designed to allow to do common and regular task, which are done if
        large photometric data sets are used, in just few lines. 
        The aim is to provide a simple interface to make such data sets and 
        specially data set combination more accessible to the community.
        
        These task are for example, cross-match different catalogs (photometric and 
        astrometric), make basic selection, plotting results, simple image and/or
        light curves downloads.
    ''',
    long_description=long_description,
    long_description_content_type="text/markdown",

    zip_safe=False,
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Topic :: Scientific/Engineering :: Astronomy",
    ],
)
