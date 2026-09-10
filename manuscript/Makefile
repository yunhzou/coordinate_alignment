PYTHON ?= ../../.venv/bin/python
.PHONY: all figures animations verify
all:
	bash scripts/build.sh
figures:
	PYTHONPATH=.build-deps $(PYTHON) scripts/build_figures.py
animations:
	PYTHONPATH=.build-deps $(PYTHON) scripts/build_animations.py
verify:
	PYTHONPATH=.build-deps $(PYTHON) scripts/validate_outputs.py
