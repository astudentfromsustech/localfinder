UV     := uv
PIP    := pip
PYTHON := $(UV) run python


.PHONY: \
	help \
	build install publish  \
	mypy ruff


help:  ## Show this help message
	@awk -F ':|##' '/^[^\t].+?:.*?##/ {printf "\033[36m%7s\033[0m:%s\n", $$1, $$NF}' $(MAKEFILE_LIST)

build:  ## Build package using uv
	$(UV) build

install:  ## Install package in development mode
	# conda create --name lf python=3.9
	# conda install -c bioconda -c conda-forge samtools bedtools ucsc-bigwigtobedgraph
	# conda install -c bioconda -c conda-forge igvtools=2.16.2
	$(PIP) install -e .

publish:  ## Publish package to PyPI using uv
	$(UV) publish

mypy:  ## Run static type checking with mypy
	# uv tool install mypy --with "pip,matplotlib"
	# uv tool run mypy localfinder --install-types
	$(UV) tool run mypy localfinder

ruff:  ## Run code formatting with ruff
	# uv tool install ruff
	$(UV) tool run ruff format
	$(UV) tool run ruff check --fix
