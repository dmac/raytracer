.DEFAULT: raytracer
.PHONY: clean run

raytracer: raytracer.c
	cc -o $@ $^ -std=c99 -lm -g -Wall -Werror -pedantic \
		-Wno-missing-braces \
		-Wno-unused-variable \

clean:
	rm -f ./raytracer

run: raytracer
	./raytracer
