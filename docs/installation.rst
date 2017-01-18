.. _installation::

Installation
============

AESOP is a Python 2.7 library, and is not compatible with Python 3.

Anaconda installation
"""""""""""""""""""""

We recommend installing the Anaconda Python Distribution as it comes with several packages used by AESOP pre-installed.
Anaconda can be downloaded from `here <https://www.continuum.io/downloads>`_. If you choose this route, install Anaconda2. Anaconda3 is not compatible with AESOP.

.. caution::

	The 64-bit version is recommended for Linux and Mac OS while 32-bit is recommended for Windows.

Install PDB2PQR, APBS and Coulomb
"""""""""""""""""""""""""""""""""

Please download and install the appropriate version of `PDB2PQR <https://sourceforge.net/projects/pdb2pqr>`_ and `APBS <https://sourceforge.net/projects/apbs>`_ (Coulomb installs alongside APBS) for your operating system. Ensure that the executables are added to your path by running the executable in terminal/command prompt::

	apbs

	pdb2pqr

.. note::

	Mac users: If you are using the APBS binary/executable to install (as opposed to compiling from source), you may need to add the following line to your path::

		export PATH=$PATH:`dirname '/Applications/APBS.app/Contents/MacOS/apbs_term'`

Install other dependencies
""""""""""""""""""""""""""

ProDy::

	pip install prody
	
.. note::

	If pip is unable to install prody (which is likely to happen if you are running a Windows operating system), 
	then you may download installation materials appropriate for your OS directly from the 
	`PyPI website <https://pypi.python.org/pypi/ProDy>`_. For Windows, it is typically easier to use the appropriate 
	installation executable appropriate for your Python 2.7 installation (choose amd64 for 64-bit Python 2.7 and choose 
	win32 for 32-bit Python 2.7). As a reminder, AESOP is only compatible with Python 2.7.

Modeller::

	conda config --add channels salilab
	conda install modeller

.. note::

	Modeller will require users to have a license key. 
	Registration is located at the `Sali Lab Website <https://salilab.org/modeller/>`_.
	The Modeller license key will be used in the following manner:
	
	Edit //anaconda/lib/modeller-9.17/modlib/modeller/config.py
	and replace XXXX with your Modeller license key 
	(or set the KEY_MODELLER environment variable before running 'conda install').

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

	git clone git@github.com:rohithmohan/aesop.git
	cd aesop-python
	pip install setup.py

.. note::

	This may require administrative privileges. 
