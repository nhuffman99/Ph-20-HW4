# USAGE: bash run_HW3.sh
# to produce plots

HW3.pdf: HW3.py
	python HW3.py 0 1

	python HW3.py 1 0

	python HW3.py 1 1

	python HW3.py 2 1

	python HW3.py long 1 0

	python HW3.py long 2 1

	pdflatex HW3.tex