.. _installation::

Installation
============

Anaconda Installation
===========
We recommend installing the Anaconda Python Distribution as it comes with several packages used by AESOP pre-installed.
Anaconda can be downloaded from `here <https://www.continuum.io/downloads>`_.

.. caution::

	The 64-bit version is recommended for Linux and Mac OS while 32-bit is recommended for Windows.

Installing AESOP
===========
You can install AESOP from `PyPi <https://pypi.python.org/>`_ using::
	
	pip install aesop

If you already have it installed, you can upgrade to the latest version using::
	
	pip install --upgrade aesop

Install PDB2PQR, APBS and Coulomb
===========
Please download and install the appropriate version of `PDB2PQR <https://sourceforge.net/projects/pdb2pqr>`_ and `APBS <https://sourceforge.net/projects/apbs>`_ (Coulomb installs alongside APBS) for your operating system. Ensure that the executables are added to your path by running the executable in terminal/command prompt::

	apbs

	pdb2pqr

