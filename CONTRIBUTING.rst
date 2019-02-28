.. _contributing-label:

============
Contributing
============

.. highlight:: bash

Contributions are welcome, and they are greatly appreciated! Every
little bit helps, and credit will always be given. 

You can contribute in many ways:

Types of Contributions
----------------------

Report Bugs
~~~~~~~~~~~

Report bugs at https://github.com/CFSAN-Biostatistics/snp-pipeline/issues.

If you are reporting a bug, please include:

* Your operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.

Fix Bugs
~~~~~~~~

Look through the GitHub issues for bugs. Anything tagged with "bug"
is open to whoever wants to implement it.

Implement Features
~~~~~~~~~~~~~~~~~~

Look through the GitHub issues for features. Anything tagged with "feature"
is open to whoever wants to implement it.

Write Documentation
~~~~~~~~~~~~~~~~~~~

SNP Pipeline could always use more documentation, whether as part of the 
official SNP Pipeline docs, in docstrings, or even on the web in blog posts,
articles, and such.

Submit Feedback
~~~~~~~~~~~~~~~

The best way to send feedback is to file an issue at https://github.com/CFSAN-Biostatistics/snp-pipeline/issues.

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.


.. _get-started-label:

Get Started!
------------

Ready to contribute? Here's how to set up `snp-pipeline` for local development.

#. Fork the `snp-pipeline` repo on GitHub.
#. Clone your fork locally::

    $ git clone git@github.com:your_name_here/snp-pipeline.git

#. Install your local copy into a virtualenv. Assuming you have virtualenvwrapper installed, this is how you set up your fork for local development::

    $ mkvirtualenv snppipeline
    $ cd snppipeline/
    $ python setup.py develop
    $ pip install sphinx_rtd_theme    # the documentation uses the ReadTheDocs theme

#. Run the unit tests on the supplied data set to verify your installation is working::

    $ python setup.py test

#. Create a branch for local development::

    $ git checkout -b name-of-your-bugfix-or-feature
   
   Now you can make your changes locally.

#. When you're done making changes, check that your changes pass the tests, including testing other Python versions::

    $ python setup.py test
    $ tox

   To get tox, just pip install it into your virtualenv. 

#. Run the regression tests::

    $ test/regression_tests.sh

   To get shunit2, install from https://code.google.com/p/shunit2/

#. Update the documentation and review the changes locally with sphinx::

    $ cd docs
    $ sphinx-build -b html . ./_build
    $ xdg-open _build/index.html

#. Commit your changes and push your branch to GitHub::

    $ git add .
    $ git commit -m "Your detailed description of your changes."
    $ git push origin name-of-your-bugfix-or-feature

#. Submit a pull request through the GitHub website.

Pull Request Guidelines
-----------------------

Before you submit a pull request, check that it meets these guidelines:

#. The pull request should include tests.
#. If the pull request adds functionality, the docs should be updated. Put
   your new functionality into a function with a docstring, and add the
   feature to the list in README.rst.
#. The pull request should work for Python 2.7, 3.4, 3.5, 3.6 and 3.7, and for PyPI.

Tips
----

To run a subset of tests::
  
    $ python -m unittest test.test_snppipeline
    $ python -m unittest test.test_utils
