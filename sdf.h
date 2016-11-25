#ifndef SDF_H
#define SDF_H

typedef struct sdf_point {
    int x;
    int y;
} sdf_point;

#define SDF_INFINITY (1.0f/0.0f)

/* Calculate squared distance field given a partial field of distances.     */
/* Also calculates the closest point with a finite distance.                */
/* Fields for which the squared distance should be calculated should be     */
/* initialized to SDF_INFINITY.                                                 */
/* If you want to calculate the sdf of a binary mask, set the distance to   */
/* 0 if the mask is filled and to SDF_INFINITY if the mask is not filled.       */
/* Values besides 0 and SDF_INFINITY work just as well.                         */
/* If you do not need the closest point to a point, you can pass in NULL.   */
/* 2D array indices are calculated as "i = x + y*nx", i.e. values are       */
/* stored in row-major order.                                               */
void sdf(float *distances, int nx, int ny, sdf_point *closest_points);

#endif /* SDF_H */

#ifdef SDF_IMPLEMENTATION

#include <stdlib.h>

static float sdf_parabola_y(float x0, float dy0, float x){
    float dx0 = x - x0;
    return dx0*dx0 + dy0*dy0;
}

void sdf(float *distances, int nx, int ny, sdf_point *closest_points){
    int x1, y;

    /* There are at most n parabolas and n + 1 intersections between */
    /* consecutive parabolas. */
    int *parabola_x = (int*)malloc((nx + 1)*sizeof(*parabola_x));
    float *x_intersections = (float*)malloc(nx*sizeof(*x_intersections));
    float *temp_row_distances = (float*)malloc(nx*sizeof(*temp_row_distances));
    sdf_point *temp_closest = (sdf_point*)malloc(nx*sizeof(*temp_closest));

    if (closest_points){
        for (y = 0; y < ny; y++){
            for (x1 = 0; x1 < nx; x1++){
                int i = x1 + y*nx;
                closest_points[i].x = x1;
                closest_points[i].y = y;
            }
        }
    }

    /* Find minimum distance in columns. */
    for (x1 = 0; x1 < nx; x1++){
        for (y = 1; y < ny; y++){
            int i = x1 + y*nx;
            int j = i - nx;
            if (distances[i] > distances[j] + 1.0f){
                distances[i] = distances[j] + 1.0f;
                if (closest_points) closest_points[i] = closest_points[j];
            }
        }
        for (y = ny - 2; y >= 0; y--){
            int i = x1 + y*nx;
            int j = i + nx;
            if (distances[i] > distances[j] + 1.0f){
                distances[i] = distances[j] + 1.0f;
                if (closest_points) closest_points[i] = closest_points[j];
            }
        }
    }

    /* Find minimum distance in rows by finding lower envelope of parabolas. */
    for (y = 0; y < ny; y++){
        int n = 0;
        parabola_x[0] = 0;

        /* Find first parabola with finite distance. */
        for (x1 = 0; x1 < nx; x1++){
            if (distances[x1 + y*nx] < SDF_INFINITY){
                parabola_x[0] = x1;
                break;
            }
        }

        for (x1 = x1 + 1; x1 < nx; x1++){
            /* Load new parabola (x1, dy1). */
            float dy1 = distances[x1 + y*nx];

            /* Infinite parabolas are not part of the lower envelope. */
            if (dy1 == SDF_INFINITY) continue;

            while (1){
                /* Load old parabola (x0, dy0). */
                int x0 = parabola_x[n];
                float dy0 = distances[x0 + y*nx];

                /* If the old parabola (x0, dy0) is above the new parabola */
                /* (x1, dy1) at the point of the last intersection */
                if (n > 0 &&
                    sdf_parabola_y(x0, dy0, x_intersections[n-1]) >
                    sdf_parabola_y(x1, dy1, x_intersections[n-1])){
                    /* it will not be in the lower envelope and is discarded. */
                    n--;
                }else{
                    /* Otherwise, add new parabola and its intersection. */
                    x_intersections[n] = 0.5f/(x1 - x0)*
                        (x1*x1 + dy1*dy1 - x0*x0 - dy0*dy0);
                    n++;
                    parabola_x[n] = x1;
                    break;
                }
            }
        }

        /* Fill in distance values based on lower envelopes. */
        for (x1 = nx - 1; x1 >= 0; x1--){
            int i, x0;
            float dy0;
            /* Go to next parabola. */
            while (n > 0 && x1 < x_intersections[n - 1]) n--;

            x0 = parabola_x[n];
            i = x0 + y*nx;
            dy0 = distances[i];
            /* Can not write directly because we also read from same array. */
            temp_row_distances[x1] = sdf_parabola_y(x0, dy0, x1);
            if (closest_points) temp_closest[x1] = closest_points[i];
        }

        /* Copy back temporary values. */
        for (x1 = 0; x1 < nx; x1++){
            int i = x1 + y*nx;
            distances[i] = temp_row_distances[x1];
            if (closest_points) closest_points[i] = temp_closest[x1];
        }
    }

    free(parabola_x);
    free(x_intersections);
    free(temp_row_distances);
    free(temp_closest);
}

#endif /* SDF_IMPLEMENTATION */
