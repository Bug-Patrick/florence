before_install:
  FLORENCEPATH=$(pwd)
  echo $FLORENCEPATH

  sudo apt-get update --fix-missing
  sudo apt-get install curl
  # Install BLAS libraries
  echo -ne '\n' | sudo apt-get install gfortran libblas-dev liblapack-dev libatlas-base-dev libopenblas-* libsuitesparse-dev swig

  # Setup python venv
  cd $FLORENCEPATH
  python3 -m venv venvpy
  source venvpy/bin/activate

  # pypi packages - install before installing PostMesh
  pip install -r "requirements.txt"


  # Install Eigen
  EIGENVERSION=3.3.8
  cd ~
  curl https://gitlab.com/libeigen/eigen/-/archive/$EIGENVERSION/eigen-$EIGENVERSION.tar.bz2 -o eigen.tar.bz2
  tar -xvf eigen.tar.bz2
  sudo mv eigen-$EIGENVERSION /usr/local/include/eigen/
  # Install OpenCascade
  echo -ne '\n' | sudo apt-get install libocct-*-dev

  # PostMesh
  pip install PostMeshPy

  # Fastor
  git clone https://github.com/romeric/Fastor
  sudo mv Fastor /usr/local/include/Fastor/

  # update the cache
  sudo ldconfig


install:
  cd $FLORENCEPATH
  echo -ne '\n' | python setup.py build

script:
  cd ~
  export PYTHONPATH=$FLORENCEPATH:$PYTHONPATH
