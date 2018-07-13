#!/bin/bash -e

# Set up an environment to run tests under Travis CI (see ../.travis.yml)

if [ $# -ne 2 ]; then
  echo "Usage: $0 conda_dir python_version"
  exit 1
fi

pmi_dir=$(pwd)
conda_dir=$1
python_version=$2
temp_dir=$(mktemp -d)

cd ${temp_dir}

# Use miniconda Python rather than the Travis environment

# Clean up after a potential previous install failure
rm -rf ${conda_dir}
# Save on some downloading if the version is the same
if [ "${python_version}" == "2.7" ]; then
  wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh -O miniconda.sh
else
  wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
fi
bash miniconda.sh -b -p ${conda_dir}
export PATH=${conda_dir}/bin:$PATH
conda update --yes -q conda
conda create --yes -q -n python${python_version} -c salilab python=${python_version} pip scipy matplotlib nose imp-nightly
source activate python${python_version}
pip install coverage

# Replace PMI1 in IMP with that from git
IMP_PATH=$(echo "import IMP, sys, os; sys.stdout.write(os.path.dirname(IMP.__file__))" | python)
cd ${IMP_PATH}
mv pmi1 pmi1.orig
cp -sr ${pmi_dir}/pyext/src pmi1
cp pmi1.orig/__init__.py pmi1.orig/_version_check.py pmi1/

# Also replace PMI1 examples, since some tests use data from them
EXAMPLE_PATH=$(echo "import IMP.pmi1, sys; sys.stdout.write(IMP.pmi1.get_example_path('..'))" | python)
cd ${EXAMPLE_PATH}
mv pmi1 pmi1.orig
cp -sr ${pmi_dir}/examples pmi1

# IMP tests use sys.argv[0] to determine their location, which won't work if
# we use nosetests, so add a workaround
ln -sf $(which nosetests) ${pmi_dir}/test/

rm -rf ${temp_dir}
