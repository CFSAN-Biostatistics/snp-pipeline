#!/bin/bash

# This script must be sourced (executed) with command ". test.sh &> test.log"
# because it uses aliases defined in the main shell.

usage()
{
  echo usage: $0 [-h] [-c] [-d]
  echo
  echo 'Test the snp-pipeline against multiple python versions, '
  echo 'each in a separate virtual environment.'
  echo
  echo 'Options:'
  echo '  -h               : Show this help message and exit'
  echo '  -c               : Create temporary virtual environments just for this test'
  echo '  -d               : Developer install, this is implied by option -c'
  echo
  echo 'Example:'
  echo '    . run_tests.sh -c &> test.log'
}


#--------
# Options
#--------

# Need to set these variable because they survive in the shell after execution
# because this scipt is sourced.  Subsequent runs will be affected by these variables.
OPTIND=1
opt_c_set=
opt_d_set=
while getopts ":hcd" option; do
  if [ "$option" = "h" ]; then
    usage
    return 0
  elif [ "$option" = "?" ]; then
    echo "Invalid option -- '$OPTARG'" 1>&2
    usage
    return 1
  elif [ "$option" = ":" ]; then
    echo "Missing argument for option -- '$OPTARG'" 1>&2
    usage
    return 2
  else
    declare opt_"$option"_set="1"
    if [ "$OPTARG" != "" ]; then
      declare opt_"$option"_arg="$OPTARG"
    fi
  fi
done

deactivate

py_ver=`python -c 'import sys; print ".".join(str(x) for x in sys.version_info[:2])'`

for ver in 2.6 2.7 3.3 3.4 3.5
do
    # Create fresh virtual env, return on error
    if [[ "$opt_c_set" == "1" ]]; then
      echo -------------------------------------------------------
      echo  python $ver CREATE venv
      echo -------------------------------------------------------
      rmvirtualenv temp
      if [[ "$ver" == "$py_ver" ]]; then
        mkvirtualenv -p /usr/bin/python temp
      else
        mkvirtualenv -p /usr/bin/python$ver temp
      fi
      status=$?
      if [ "$status" -ne 0 ]; then return $status; fi
      venv="temp"
    else
      workon snp-pipeline-$ver
      venv="snp-pipeline-$ver"
    fi

    # Install the pipeline into the venv, return on error
    if [[ "$opt_c_set" == "1" || "$opt_d_set" == "1" ]]; then
      echo -------------------------------------------------------
      echo  python $ver DEVELOP using venv $venv
      echo -------------------------------------------------------
      python setup.py develop
      status=$?
      if [ "$status" -ne 0 ]; then return $status; fi
    fi

    # Run unit tests, return on error
    echo -------------------------------------------------------
    echo  python $ver TEST using venv $venv
    echo -------------------------------------------------------
    python setup.py test
    status=$?
    deactivate
    if [ "$status" -ne 0 ]; then return $status; fi
done

# Need to set these variable because they survive in the shell after execution
# because this scipt is sourced.  Subsequent runs will be affected by these variables.
OPTIND=1
opt_c_set=
opt_d_set=
