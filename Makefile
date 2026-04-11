# Disable check for file with target's name / target doesn't produce a file
.PHONY: hello help apt_refresh apt_dependencies before_install dependencies python install

EIGENVERSION:=3.3.8

FLORENCEPATH:=$(CURDIR)

help:
	@echo 'Usage: make [target]'
	@echo ''
	@echo 'Prerequisite:'
	@echo '	Run from directory of cloned repository'
	@echo ''
	@echo 'Targets:'
	@echo '	apt_refresh		Update and Upgrade packages - Unnecessary, use with caution'
	@echo '	before_install		Install dependencies from apt, external sources and Poyas projects'
	@echo '	apt_dependencies	Install apt listed packages'
	@echo '	outer_dependencies	Install Eigen'
	@echo '	custom_dependencies	Install PostMesh and Fastor from Roman Poya'
	@echo '	PostMesh		Install PostMesh from Roman Poya ~ via own Fork'
	@echo '	Fastor			Install Fastor from Roman Poya'
	@echo '	python			DOESNT WORK - COPY TO CMD SHELL - create Virtual Environment with installed packages from requirements'
	@echo '	install			first time cython execution'
	@echo '	check_install		first time cython execution'
	@echo '	run			recurring cython execution'
	@echo '	env			Provide Florence in PythonPath, Usage: eval $$(make env)'

apt_refresh:
	sudo apt update --fix-missing
	apt list --upgradable
	sudo apt upgrade

before_install: apt_dependencies outer_dependencies custom_dependencies
# update the cache
	sudo ldconfig

apt_dependencies:
	sudo apt install curl
# Install BLAS libraries
	sudo apt install -y gfortran libblas-dev liblapack-dev libatlas-base-dev libopenblas-* libsuitesparse-dev swig
# Install OpenCascade - dependency for PostMesh
	sudo apt install -y libocct-*-dev

outer_dependencies:
# Install Eigen
# cd doesn't work in Makefile, as each line runs in an own separate shell!	cd ~
	curl https://gitlab.com/libeigen/eigen/-/archive/$(EIGENVERSION)/eigen-$(EIGENVERSION).tar.bz2 -o eigen.tar.bz2
	tar -xvf eigen.tar.bz2
	sudo mv eigen-$(EIGENVERSION) /usr/local/include/eigen/

# update the cache
	sudo ldconfig

custom_dependencies: PostMesh Fastor


PostMesh:
# pip install PostMeshPy
# better clone updated PostMesh fork
	git clone git@github.com:Bug-Patrick/PostMesh.git ~/repos/PostMesh
# provide own Makefile and Instructions, Is PostMesh in Path? or how can Florence link to it?

Fastor:
	cd ~/repos
	git clone https://github.com/romeric/Fastor ~/repos/Fastor
	sudo mv Fastor /usr/local/include/Fastor/

python:
# Setup python venv
# cd $(FLORENCEPATH) doesn't work
	python3 -m venv venv
# limitation of Makefile: activate only for 1 line, as it is a process run in its own environment! -> .sh better
	. venv/bin/activate
# pypi packages - install before installing PostMesh
	pip install -r "requirements.txt"

install: before_install python
	python setup.py build

check_install:
#check dependencies
	@[ -d "/usr/local/include/eigen/" ] && echo "Eigen exists."
	@[ -d "/usr/local/include/Fastor/" ] && echo "Fastor exists."
	@[ -d "$$HOME/repos/PostMesh/" ] && echo "PostMesh exists."
	@pip check
	@pip install --dry-run -r requirements.txt 2>&1 | grep -q "Would install" \
		&& { echo "Missing packages:"; pip install --dry-run -r requirements.txt; exit 1; } \
		|| echo "All Python packages are installed."

run:
	. venv/bin/activate
	python setup.py build

env:
# limitation of Makefile, as it is a process run in its own environment! -> .sh better
# if-else to prevent a trailing colon :
	@if [ -z "$$PYTHONPATH" ]; then \
		echo "export PYTHONPATH=$(FLORENCEPATH)"; \
	else \
		echo "export PYTHONPATH=$(FLORENCEPATH):$$PYTHONPATH"; \
	fi