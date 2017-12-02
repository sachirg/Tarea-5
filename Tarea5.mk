Resultados_hw5.pdf : Resultados_hw5.tex  curvaRotacion.png
	pdflatex $< && rm *.aux *.log 

curvaRotacion.png : plots.py 
	python $<

plots.py: daticos.txt
	./a.x

daticos.txt : a.x
	./a.x

a.x : CurvaRotacion.c
	cc CurvaRotacion.c -lm -o a.x

clean :
	rm *.dat *.png *.pdf *.tex *.log *.aux *.gz *.txt *.dat *.x
