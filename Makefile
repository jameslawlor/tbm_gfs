.PHONY: all
# Clean temporary and generated files
all: ruff tests clean

.PHONY: clean
# Clean temporary and generated files
clean:
	find . \( -type f -name '*.pyc' -or -type d -name '__pycache__*' \) -delete
	find . \( -type d -name '.eggs' -or -type d -name '*.egg-info' -or -type d -name '.pytest_cache' \) | xargs rm -rf


.PHONY: ruff
# Clean temporary and generated files
ruff:
	ruff format .
	ruff check .
	
.PHONY: tests
tests:
	pytest -vs