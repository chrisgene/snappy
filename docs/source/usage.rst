Usage
=====

SNAPPY was recently moved to a package installation rather than a standalone script. Some details recorded below may not be applicable to the package installation.

Quick Start:
------------

Download or clone `SNAPPY <https://www.github.com/chrisgene/snappy/>`_, then navigate in the terminal to the folder where SNAPPY is saved. The folder should also contain the files parse_ref_files.py, parse_plink_files.py, and a directory called ‘ref_files’ that contains four additional files. Run SNAPPY with the following command:
::

   snappy --infile plink_library

where ``plink_library`` is the prefix name of the genotypes to be analyzed. Note that SNAPPY uses plink (v1.9) in a preprocessing step so a plink (v1.9) executable must be listed in the user’s path as `plink` for preprocessing steps. If plink (v1.9) is not available in the path, the input can be created manually using plink (v1.9) with the `--recodeAD` option. plink is available for all major operating systems and can be downloaded `here <https://www.cog-genomics.org/plink/1.9/>`_. 

.. _installation:

Installation:
-------------

Download or clone `SNAPPY <https://www.github.com/chrisgene/snappy/>`_, then navigate in the roo of the the local SNAPPY repository in terminal. Run the following commands:
::

   python setup.py sdist
   pip install dist/snappy-X.X.tar.gz
   
SNAPPY is now installed as a command-line tool and can be accessed from any directory. To test the installation, navigate to the root of the local SNAPPY repository, then run:
::

   snappy --infile test_data/TGP_chrY --out test_SNAPPY

Test SNAPPY's output against the output distributed with the software:
::

   diff test_SNAPPY.all test_data/TGP_SNAPPY_results.all > test_SNAPPY_diffs.txt

A working installation of SNAPPY should yield an empty file ``test_SNAPPY_diffs.txt``. 

Citation:
=========

If you use SNAPPY, please cite our `preprint on bioRXiv<https://www.biorxiv.org/content/10.1101/454736v2>`_.

Terms of Use:
-------------

SNAPPY is published under a GPL-3.0 License. More information about the license is available `here <https://opensource.org/licenses/GPL-3.0>`_.
