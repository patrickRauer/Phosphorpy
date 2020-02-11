from setuptools import setup

setup(
    name='Phosphorpy',
    version='0.6.14',
    python_requires='>=3.6',
    packages=['Phosphorpy', 'Phosphorpy.data',
              'Phosphorpy.data.sub', 'Phosphorpy.data.sub.plots',
              'Phosphorpy.config',
              'Phosphorpy.fitting', 'Phosphorpy.external',
              'Phosphorpy.core', 'Phosphorpy.data.sub.tables', 'Phosphorpy.local'],
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
    description='',
    zip_safe=False
)
