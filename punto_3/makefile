params.pdf : grafica.py datos.dat
	python grafica.py

datos.dat : a.out
	./a.out > datos.dat

a.out : poblaciones.c
	gcc -lm poblaciones.c

clean:
	rm -f a.out datos.dat params.pdf
