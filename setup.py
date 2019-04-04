from setuptools import setup

setup(
    name='Phosphorpy',
    version='0.5.4.3',
    packages=['Phosphorpy', 'Phosphorpy.data',
              'Phosphorpy.data.sub', 'Phosphorpy.data.sub.plots',
              'Phosphorpy.config', 'Phosphorpy.report',
              'Phosphorpy.fitting', 'Phosphorpy.external'],
    # package_dir={'': 'Phosphorpy'},
    include_package_data=True,
    package_data={
          '': ['local/survey.conf']
      },
    url='https://github.com/patrickRauer/Phosphorpy',
    license='GPL',
    author='Patrick Rauer',
    author_email='j.p.rauer@sron.nl',
    description='',
    zip_safe=False,
    install_requires=['seaborn', 'numpy', 'astropy', 'pandas', 'astroquery', 'numba', 'scikit-learn', 'armapy',
                      'requests', 'extinction', 'scipy']
)