from setuptools import setup, find_packages

setup(name='snappy',
      version='2.2',
      description='Y-chromosome haplogroup inference',
      url='http://github.com/chrisgene/snappy',
      author='Alissa Severson, Jonathan Shortt, Chris Gignoux',
      author_email='jonathan.shortt@cuanschutz.edu',
      license='GPLv3.0',
      packages=['snappy', 'snappy/bin'],
      install_requires=[ #numpy is only module that is not included in standard distributions of python
            'numpy>=1.13.3'
      ],
      entry_points = { 'console_scripts': [
      		'snappy=snappy.main:run_snappy', 
            'snappy-clean=snappy.main:clean_isogg_table',
            'snappy-qc=snappy.main:do_isogg_qc',
            'snappy-build=snappy.main:make_ref_files'
        ],
      },
      zip_safe=False)
