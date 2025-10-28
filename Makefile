## Makefile for coralsnake

.PHONY: help install install-dev test test-cov build clean lint format check all

help:  ## Show help
	@echo "Available commands:"
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-15s\033[0m %s\n", $$1, $$2}'

install:  ## Install the package and dependencies
	uv sync
	uv pip install -e . --no-deps

install-dev:  ## Install with development dependencies
	uv sync --extra dev
	uv pip install -e . --no-deps

test:  ## Run tests
	uv run pytest tests/ -v

test-cov:  ## Run tests with coverage
	uv run pytest tests/ --cov=coralsnake --cov-report=html --cov-report=term

build:  ## Build the package (wheel and sdist)
	uv build

clean:  ## Clean build/test artifacts
	rm -rf build/ dist/ *.egg-info/ .pytest_cache/ htmlcov/ .coverage
	find . -type d -name __pycache__ -exec rm -rf {} +
	find . -type f -name "*.pyc" -delete

lint:  ## Run linting
	uv run ruff check coralsnake/ tests/
	uv run ruff format --check coralsnake/ tests/

format:  ## Format code
	uv run ruff format coralsnake/ tests/
	uv run ruff check --fix coralsnake/ tests/

check: lint test  ## Run lint + tests

all: clean install-dev test build  ## Full pipeline
