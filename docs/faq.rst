===========================
FAQ / Troubleshooting Guide
===========================

Installation
------------

Q: How can I avoid polluting my global python installation when installing the SNP pipeline?

A: You can either use a python virtual environment or install into your user area.

To use a python virtual environment, see :ref:`contributing-label`

To install into your user area instead of installing into your global site packages::

	pip install --user snp-pipeline


Running the Pipeline
--------------------

Q: Nothing works?

A: Make sure you have the proper dependencies on your path.  Modify your path if necessary to include bowtie, samtools, and varscan.  Also make sure you install Biopython.  If you are using virtual environments, make sure you issue the command "toggleglobalsitepackages" so Biopython can be found.  See also: the :ref:`installation-label` section of this documentation.

Q: How can I verify the pipeline is installed and working properly?

A: The SNP Pipeline includes a small set of test data with result files.  You can run the pipeline against the test data to verify correct results.  The test data is in the test/testLambdaVirusClean directory.  You can run the test with::

	runLamdaVirusTest.sh   # TODO verify this script works

Q: My results for the included test data do not match the expected results. What is the cause?

A: Different versions of the executable tools can generate different results.  The test data was generated with these versions:
	
	* bowtie2 2.2.2
	* samtools 0.1.18


Developer Questions
-------------------

Q: What causes "ImportError: No module named sphinx_rtd_theme" when building the documentation?

A: The documentation uses the *Read The Docs* theme.  Install it like this::

	pip install sphinx_rtd_theme

Q: I installed sphinx_rtd_theme, but I still get error "ImportError: No module named sphinx_rtd_theme"?

A: Try running sphinx like this::

	python /usr/bin/sphinx-build -b html  .  ./_build