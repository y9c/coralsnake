PY=uv run --python 3.13

.PHONY: format lint test build clean

format:
	$(PY) ruff check --fix .

lint:
	$(PY) ruff check .

test:
	$(PY) coralsnake map -r test/ref.fa -1 test/test1.fq -2 test/test2.fq -o /tmp/out.bam -m 0 -t 2

build:
	$(PY) python -m build_ext

clean:
	rm -rf dist build *.egg-info

