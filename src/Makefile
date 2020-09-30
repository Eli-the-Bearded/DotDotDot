
pbmtobraille: constants.c grid.c unicode.c renderer.c pbmtobraille.c
	gcc -g -o pbmtobraille constants.c grid.c unicode.c renderer.c pbmtobraille.c -Wall -Werror -lncursesw -lm -lpbm

sampleprog: examples.c constants.c grid.c unicode.c renderer.c main.c
	gcc -g -o dots examples.c constants.c grid.c unicode.c renderer.c main.c -Wall -Werror -lncursesw -lm

clean:
	rm -f dots pbmtobraille
