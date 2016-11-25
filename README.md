# HOWTO

Fill a float distance array with the values 0 and INFINITY and then call:

```C
df(squared_distances, nx, ny, closest_points);
```

# Full example

```C
#define DF_IMPLEMENTATION
#include "df.h"

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <math.h>

int main(){
    int i;
    int nx = 1920;
    int ny = 1080;
    FILE *fp;

    /* Pixels for image of distance field. */
    uint8_t *pixels = (uint8_t*)calloc(nx*ny, 3);
    /* Colors for voronoi diagram. */
    uint8_t *colors = (uint8_t*)calloc(nx*ny, 3);

    float *squared_distances = (float*)malloc(nx*ny*sizeof(*squared_distances));
    df_point *closest_points = (df_point*)calloc(nx*ny, sizeof(*closest_points));

    /* Initialize distances to infinity. */
    for (i = 0; i < nx*ny; i++){
        squared_distances[i] = DF_INFINITY;
    }

    /* Except some random points. */
    for (i = 0; i < 100; i++){
        int x = rand() % nx;
        int y = rand() % ny;

        int j = x + y*nx;

        /* Other distance values than 0 and INFINITY work, too! */
        squared_distances[j] = 0.0f;

        colors[j*3 + 0] = rand();
        colors[j*3 + 1] = rand();
        colors[j*3 + 2] = rand();
    }

    /* Calculate squared distance field. */
    /* Last parameter may be NULL if not needed. */
    df(squared_distances, nx, ny, closest_points);

    /* Draw voronoi diagram. */
    for (i = 0; i < nx*ny; i++){
        df_point p = closest_points[i];
        int j = p.x + p.y*nx;
        pixels[i*3 + 0] = colors[j*3 + 0];
        pixels[i*3 + 1] = colors[j*3 + 1];
        pixels[i*3 + 2] = colors[j*3 + 2];
    }

    /* Write ppm image file of voronoi diagram. */
    fp = fopen("voronoi.ppm", "wb");
    fprintf(fp, "P6\n%d %d\n255\n", nx, ny);
    fwrite(pixels, 1, 3*nx*ny, fp);
    fclose(fp);

    /* Draw distance field. */
    for (i = 0; i < nx*ny; i++){
        float d = sqrt(squared_distances[i]);
        pixels[i] = d > 255 ? 255 : d;
    }

    /* Write pgm image file of distance field. */
    fp = fopen("df.pgm", "wb");
    fprintf(fp, "P5\n%d %d\n255\n", nx, ny);
    fwrite(pixels, 1, nx*ny, fp);
    fclose(fp);

    /* Cleanup. */
    free(pixels);
    free(squared_distances);
    free(closest_points);
    free(colors);
    return 0;
}
```

# Distance field

![Image of distance field](https://raw.githubusercontent.com/983/sdf/master/df.jpg)

# Voronoi diagram

![Image of Voronoi diagram](https://raw.githubusercontent.com/983/sdf/master/voronoi.png)
