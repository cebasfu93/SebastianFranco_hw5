params.pdf : datos.dat grafica.py
	python grafica.py

datos.dat : a.out
	./a.out > datos.dat

a.out : planetas.c
	gcc planetas.c -lm

clean:
	rm -f a.out datos.dat M.txt params.pdf
