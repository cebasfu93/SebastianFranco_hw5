params.pdf : datos.dat grafica.py
	python grafica.py

datos.dat : a.out
	./a.out > datos.dat

a.out : sismo.c
	gcc -lm sismo.c

clean:
	rm -f a.out datos.dat M.txt params.pdf
