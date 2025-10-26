# Simple Makefile for coralsnake

PY=uv run --python 3.13

.PHONY: help format lint test build clean

help:
	@echo "Targets: format lint test build clean"

format:
	$(PY) ruff check --fix .

lint:
	$(PY) ruff check .

# Minimal smoke test: small paired-end run
# Adjust paths if needed
TEST_OUT=/tmp/coralsnake_test.bam

test:
	$(PY) coralsnake map -r test/ref.fa -1 test/test1.fq -2 test/test2.fq -o $(TEST_OUT) -m 0 -t 2
	@echo "Test output: $(TEST_OUT)"

build:
	$(PY) python -m pip install -e .

clean:
	rm -rf dist build *.egg-info __pycache__ **/__pycache__
