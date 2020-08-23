#include <math.h>
#include "grid.h"
#include "renderer.h"
#include <stdio.h>

// print all 256 eight bit patterns
void try_them_all()
{
    int width = 64;
    int height = 64;

    grid *g = grid_new(width, height);
    renderer_new(g);

    renderer_update(g);

    for (int s = 0; s < 256; s ++) {

        for (int b = 0; b < 8; b ++) {
	    int i = 2 * (s % 32) + (b / 4);
	    int j = 4 * (s / 32) + (b % 4);
	    if ( s & (1 << b) ) {
                grid_set_pixel(g, i, j);
	    } 
        }
        renderer_update(g);
    } 

    // Free allocations
    renderer_free();
    grid_free(g);
}

void example_filling_bar()
{
    int width = 100;
    int height = 24;

    grid *g = grid_new(width, height);
    renderer_new(g);

    // Fill grid from left to right (simple animation)
    renderer_update(g);
    for (int i = 0; i < width; i++)
    {
        for (int j = 0; j < height; j++)
        {
            grid_set_pixel(g, i, j);
        }
        renderer_update(g);
    }

    // Free allocations
    renderer_free();
    grid_free(g);
}

void example_build_block()
{
    int width = 100;
    int height = 40;

    grid *g = grid_new(width, height);
    renderer_new(g);

    for (int x = 0; x < width; x++)
    {
        for (int y = 0; y < height; y++)
        {
            grid_set_pixel(g, x, y);
            renderer_update(g);
        }
    }

    // Free allocations
    renderer_free();
    grid_free(g);
}

void example_sine_tracking()
{
    int width = 124;
    int height = 40;

    grid *g = grid_new(width, height);
    renderer_new(g);

    double shift = 0;

    while (1)
    {
        grid_clear(g);

        // Draw line
        grid_draw_line(g, 0, height / 2, width - 1, (height + sin(shift) * height) / 2);

        // Draw curve
        for (int j = 0; j < width; j++)
        {
            grid_set_pixel(g, j, (height / 2 * sin(0.05 * j + shift) + height / 2));
        }

        // Move curve
        shift += 0.05;

        renderer_update(g);
    }

    // Free allocations
    renderer_free();
    grid_free(g);
}

void example_spiral_effect()
{
    int width = 60;
    int height = 32;

    grid *g = grid_new(width, height);
    renderer_new(g);

    // Start with an empty grid
    grid_clear(g);

    int m = width, n = height;
    int sr = 0, sc = 0, er = m - 1, ec = n - 1;
    while (sr <= er && sc <= ec)
    {
        for (int i = sc; i <= ec; ++i)
        {
            grid_set_pixel(g, sr, i);
            renderer_update(g);
        }
        for (int i = sr + 1; i <= er; ++i)
        {
            grid_set_pixel(g, i, ec);
            renderer_update(g);
        }
        for (int i = ec - 1; sr != er && i >= sc; --i)
        {
            grid_set_pixel(g, er, i);
            renderer_update(g);
        }
        for (int i = er - 1; sc != ec && i > sr; --i)
        {
            grid_set_pixel(g, i, sc);
            renderer_update(g);
        }
        sr++, sc++;
        er--, ec--;
    }

    // Free allocations
    renderer_free();
    grid_free(g);
}

void braille_graph( grid *g, double *f, int n )
{
    int    row=0,nrows=24, col=0,ncols=78;
    double bigf=0.0,yval;

    grid_clear(g);

    for ( col=0; col<n; col++ )
       if ( fabs(f[col]) > bigf ) bigf = fabs(f[col]);

    for ( row=0; row<nrows; row++ ) {
        yval = bigf*((double)(nrows/2-row))/((double)(nrows/2));

        for ( col=0; col<ncols; col++ )
	   if (yval*f[(col*(n-1))/(ncols-1)]>=yval*yval)
	       grid_set_pixel(g, col, row);
    }
    
    renderer_update(g);
}

void johnforkosh_function()
{
    int width = 80;
    int height = 48;

    double f[999], pi=3.14159;
    double t=0.0, dt=0.05, w1=16.,w2=3.;  int Nt=50;
    int    i=0, N=511;                    /* f[] index */

    grid *g = grid_new(width, height);
    renderer_new(g);
    renderer_update(g);
    
    while ( --Nt > 0 ) {
	for ( i=0; i<N; i++ ) {
	    double x = 2.0*pi*((double)i)/((double)(N-1));
	    f[i] = .75*sin(2.*x+pi/3.+w1*t) + 1.*sin(1.*x+pi/2.+w2*t);
	}

	braille_graph(g, f, N);
	t += dt;
    }
    renderer_free();
    grid_free(g);
}

