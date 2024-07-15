paths:
	mkdir -p ../results
	mkdir -p ../results/adsbc21/plots
	mkdir -p ../results/dbc19/plots
	mkdir -p ../results/adbc23/plots
	mkdir -p ../data/dbc19_corrected
	mkdir -p ../results/dbc19_corrected/plots

collect_leo13b_surp:
	python 01_collect_surprisal.py -m leo13b

collect_secretgpt2_surp:
	python 01_collect_surprisal.py -m secret-gpt-2

eval_surp:
	python 02_eval_surprisal.py

rERP:
	julia 03_do_rERP.jl

plot_rERP:
	Rscript 04_do_plot_rERP.r

analysis: paths eval_surp rERP plot_rERP

clean:
	rm -rf ../results
	rm -rf ../data/*/*surp*