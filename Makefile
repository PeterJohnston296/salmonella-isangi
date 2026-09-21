check:
	python3 scripts/audit.py --burden
	python3 -m unittest discover -s tests -v

figures:
	Rscript scripts/reproduce.R

setup:
	Rscript environment/install.R

rcheck:
	Rscript tests/check_R.R

archive:
	python3 release/package.py
