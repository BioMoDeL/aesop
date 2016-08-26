.. _installation::

Installation
============

Anaconda installation
"""""""""""""""""""""

We recommend installing the Anaconda Python Distribution as it comes with several packages used by AESOP pre-installed.
Anaconda can be downloaded from `here <https://www.continuum.io/downloads>`_.

.. caution::

	The 64-bit version is recommended for Linux and Mac OS while 32-bit is recommended for Windows.

Install PDB2PQR, APBS and Coulomb
"""""""""""""""""""""""""""""""""

Please download and install the appropriate version of `PDB2PQR <https://sourceforge.net/projects/pdb2pqr>`_ and `APBS <https://sourceforge.net/projects/apbs>`_ (Coulomb installs alongside APBS) for your operating system. Ensure that the executables are added to your path by running the executable in terminal/command prompt::

	apbs

	pdb2pqr

Install other dependencies
""""""""""""""""""""""""""

ProDy::

	pip install prody

Modeller::

	conda config --add channels salilab
	conda install modeller

GridDataFormats::

	pip install GridDataFormats

Multiprocessing (optional)::

	pip install multiprocessing

.. note::

	Multiprocessing is only needed if you are planning to utilize multiple cores in your analysis.

Installing AESOP
""""""""""""""""

You can install AESOP from `PyPi <https://pypi.python.org/>`_ using::
	
	pip install aesop

If you already have it installed, you can upgrade to the latest version using::
	
	pip install --upgrade aesop

If you are having issues installing through PyPi, you may try to `Install from source`_.

Install from source
"""""""""""""""""""

To install from source, you can use the following commands to clone the GitHub repository and install manually::

	git clone git@github.com:rohithmohan/aesop-python.git
	cd aesop-python
	pip install setup.py

.. note::

	This may require administrative privileges. 
