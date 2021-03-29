/*
Copyright (C) 2021 Ellert van der Velden
All rights reserved.

This software is free software, distributed under the BSD-3 License.
You may redistribute and/or modify it without any restrictions, as long as the conditions specified in the terms of the BSD-3 license (included) are met.
*/

// Header file
#include "../include/cube.h"

// Standard libraries
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

// Macros
#define min(A, B) (A < B ? A : B)
#define max(A, B) (A > B ? A : B)

// Declare struct for holding Y and Z coordinate data
struct coord_data{
    // Data value
    double v;
    // Starting index of lower coordinate
    unsigned long i;
};

// Declare struct for holding geometry data
struct rubiks_cube{
    // Dynamic array of cube center X-coordinates
    double *x;
    // Dynamic array of cube center Y-coordinate data
    struct coord_data *y;
    // Dynamic array of cube center Z-coordinate data
    struct coord_data *z;
    // Distance from cube center to cube face in X-direction
    double dx;
    // Distance from cube center to cube face in Y-direction
    double dy;
    // Distance from cube center to cube face in Z-direction
    double dz;
    // Dynamic array of cube densities
    double *density;
    // Dynamic array of cube medium/rock IDs
    int *rock_id;
    // Index of last cube that was checked
    unsigned long cur_idx[3];
    // Index of estimated next cube
    unsigned long next_idx[3];
    // Total number of layers (length of z)
    size_t n_layers;
    // Total number of stacks (length of y)
    size_t n_stacks;
    // Total number of cubes (length of x)
    size_t n_cubes;
};

// Private function for interpreting a line from a CSV file
static int interpret_CSV_line(char *line, double *x, double *y, double *z,
                              double *density, int *rock_id){
    // Initialize variables
    int ntok = 5;
    const char sep[2] = ",";
    char *field;

    // Obtain the x field
    if ((field = strtok(line, sep)) == NULL) return ntok;
    else ntok--;
    *x = strtod(field, NULL);

    // Obtain the y field
    if ((field = strtok(NULL, sep)) == NULL) return ntok;
    else ntok--;
    *y = strtod(field, NULL);

    // Obtain the z field
    if ((field = strtok(NULL, sep)) == NULL) return ntok;
    else ntok--;
    *z = strtod(field, NULL);

    // Obtain the density field
    if ((field = strtok(NULL, sep)) == NULL) return ntok;
    else ntok--;
    *density = strtod(field, NULL);

    // Obtain the rock_id field
    if ((field = strtok(NULL, sep)) == NULL) return ntok;
    else ntok--;
    *rock_id = (int)strtod(field, NULL);

    return ntok;
}

// Create and/or read the meta data file required for the rubiks_cube struct
static enum rubiks_cube_return read_meta_data(const char *file_path, size_t *n_cubes_ptr,
                                              size_t *n_stacks_ptr, size_t *n_layers_ptr){
    // Check if the provided data file exists
    if (access(file_path, R_OK) == -1) {
        // If no existing file was passed, raise error and return
        fprintf(stderr, "ERROR: No valid file path was passed (%s)!\n", file_path);
        return(RUBIKS_CUBE_RETURN_FILE_NOT_FOUND);
    }

    // Obtain name of meta data file
    char file_path_meta[80];
    char line[80];
    sprintf(file_path_meta, "%s.meta", file_path);

    // Check if this meta data file exists
    if (access(file_path_meta, R_OK) == -1) {
        // If not, it must be created
        // Open the data file
        FILE *file = fopen(file_path, "r");

        // Declare some variables required for determining the number of cubes/stacks/layers in the file
        int rock_id_temp, ntok;
        double x_temp, yi, yj, zi, zj, density_temp;
        size_t n_cubes, n_stacks, n_layers;
        char * ret;

        // Initialize numbers to 1 (as the first line is not inspected)
        n_cubes = n_stacks = n_layers = 1;

        // Read in first line of file
        ret = fgets(line, sizeof line, file);
        if (ret == NULL) {
            // If this is the EOF, the file is empty
            fclose(file);
            fprintf(stderr, "ERROR: File is empty!\n");
            return(RUBIKS_CUBE_RETURN_EMPTY_FILE);
        }

        // Check if first character is part of a number
        if (!(line[0] >= '0' && line[0] <= '9')) {
            // If not, the file has a header
            fgets(line, sizeof line, file);
        }

        // Interpret the first line
        ntok = interpret_CSV_line(line, &x_temp, &yi, &zi, &density_temp, &rock_id_temp);
        if (ntok) goto format_error;

        // Read new line
        ret = fgets(line, sizeof line, file);

        // Read all remaining lines in the file and determine number of cubes, stacks and layers
        while (ret != NULL) {
            // Swap yi and zi to yj and zj
            yj = yi;
            zj = zi;

            // Increase n_cubes by 1 
            n_cubes++;

            // Interpret line
            ntok = interpret_CSV_line(line, &x_temp, &yi, &zi, &density_temp, &rock_id_temp);
            if (ntok) goto format_error;

            // Check if a new layer was read in
            if (zi != zj) {
                n_layers++;
                n_stacks++;
            }

            // Check if a new stack was read in
            else if (yi != yj) {
                n_stacks++;
            }

            // Read new line
            ret = fgets(line, sizeof line, file);
        }

        // Close the data file
        fclose(file);
        goto write_metadata;

format_error:
        fclose(file);
        fprintf(stderr, "ERROR: Invalid file format (%s)!\n", file_path);
        return RUBIKS_CUBE_RETURN_FORMAT_ERROR;

write_metadata:
        // Write the meta data to file
        file = fopen(file_path_meta, "w");
        fprintf(file, "%zu,%zu,%zu", n_cubes, n_stacks, n_layers);
        fclose(file);

        // Assign values to proper pointers
        *n_cubes_ptr = n_cubes;
        *n_stacks_ptr = n_stacks;
        *n_layers_ptr = n_layers;
    }
    else {
        // If so, read it
        FILE *file_meta = fopen(file_path_meta, "r");
        fgets(line, sizeof line, file_meta);
        fclose(file_meta);

        // Initialize variables
        const char sep[2] = ",";
        char *field;

        // Obtain the x field
        field = strtok(line, sep);
        sscanf(field, "%zu", n_cubes_ptr);

        // Obtain the y field
        field = strtok(NULL, sep);
        sscanf(field, "%zu", n_stacks_ptr);

        // Obtain the z field
        field = strtok(NULL, sep);
        sscanf(field, "%zu", n_layers_ptr);
    }

    // Return that everything went okay
    return(RUBIKS_CUBE_RETURN_SUCCESS);
}

// Initialize rubiks_cube struct from file
enum rubiks_cube_return rubiks_cube_create(struct rubiks_cube **cube_ptr, 
                                           const char *file_path){
    // Initialize size variables
    size_t n_cubes, n_stacks, n_layers;

    // Obtain values for these variables
    enum rubiks_cube_return rc = read_meta_data(file_path, &n_cubes, &n_stacks, &n_layers);
    if (rc != RUBIKS_CUBE_RETURN_SUCCESS) {
        // If that failed, return the error code
        return(rc);
    }

    // Allocate memory for rubiks_cube struct
    struct rubiks_cube *cube = (struct rubiks_cube *)malloc(sizeof(struct rubiks_cube));

    // Check if memory allocation was done correctly
    if (cube == NULL) {
        // If not, raise error and return
        fprintf(stderr, "ERROR: Failed to allocate memory for rubiks_cube object!\n");
        return(RUBIKS_CUBE_RETURN_MALLOC_ERROR);
    }

    // Allocate memory for member arrays in the rubiks_cube struct
    cube->x = (double *)malloc(sizeof(double)*n_cubes);
    cube->y = (struct coord_data *)malloc(sizeof(struct coord_data)*n_stacks);
    cube->z = (struct coord_data *)malloc(sizeof(struct coord_data)*n_layers);
    cube->density = (double *)malloc(sizeof(double)*n_cubes);
    cube->rock_id = (int *)malloc(sizeof(int)*n_cubes);

    // Check if memory allocation was done correctly
    if (cube->x == NULL || cube->y == NULL || cube->z == NULL ||
        cube->density == NULL || cube->rock_id == NULL) {
        // If not, raise error and return
        fprintf(stderr, "ERROR: Failed to allocate memory for rubiks_cube object members!\n");
        return(RUBIKS_CUBE_RETURN_MALLOC_ERROR);
    }

    // Assign created rubiks_cube object to provided pointer
    *cube_ptr = cube;

    // Initialize some variables required for looping through the file
    unsigned long i;
    double yi, yj, zi, zj;
    char line[80];

    // Open the data file
    FILE *file = fopen(file_path, "r");

    // Read in first line of file
    fgets(line, sizeof line, file);

    // Check if first character is part of a number
    if (!(line[0] >= '0' && line[0] <= '9')) {
        // If not, the file has a header, so read in the actual first line
        fgets(line, sizeof line, file);
    }

    // Interpret the first line of the file outside of the following loop
    interpret_CSV_line(line, &cube->x[0], &yi, &zi,
                       &cube->density[0], &cube->rock_id[0]);

    // Assign values of yi and zi to y and z
    cube->y[0].v = yi;
    cube->y[0].i = 0;
    cube->z[0].v = zi;
    cube->z[0].i = 0;

    // Set the values of dx, dy and dz for checking later
    cube->dx = cube->dy = cube->dz = 0;

    // Loop over all remaining lines in the file and read in their contents
    unsigned long ny = 1;
    unsigned long nz = 1;
    for (i=1; i<n_cubes; i++) {
        // Swap yi and zi to yj and zj
        yj = yi;
        zj = zi;

        // Read line and assign
        fgets(line, sizeof line, file);
        interpret_CSV_line(line, &cube->x[i], &yi, &zi,
                           &cube->density[i], &cube->rock_id[i]);

        // Attempt to determine dx if not done before
        if (cube->dx == 0 && cube->x[i] != cube->x[i-1]) {
            cube->dx = fabs(cube->x[i]-cube->x[i-1])/2;
        }

        // Check if the current z value is different from the previous
        if (zi != zj) {
            // If so, store it in z and add 1 to nz
            cube->z[nz].v = zi;
            cube->z[nz].i = ny;
            nz++;

            // Store new y as well
            cube->y[ny].v = yi;
            cube->y[ny].i = i;
            ny++;

            // Attempt to determine dz if not done before
            if (cube->dz == 0) {
                cube->dz = fabs(zi-zj)/2;
            }
        }

        // Check if the current y value is different from the previous
        else if (yi != yj) {
            // If so, store it in y and add 1 to ny
            cube->y[ny].v = yi;
            cube->y[ny].i = i;
            ny++;

            // Attempt to determine dy if not done before
            if (cube->dy == 0) {
                cube->dy = fabs(yi-yj)/2;
            }
        }
    }
    fclose(file);

    // If dy or dz is still zero, set to dx
    if (cube->dy == 0) {
        cube->dy = cube->dx;
    }
    if (cube->dz == 0) {
        cube->dz = cube->dx;
    }

    // Assign remaining members
    cube->cur_idx[0] = 0;
    cube->cur_idx[1] = 0;
    cube->cur_idx[2] = 0;
    cube->next_idx[0] = -1;
    cube->next_idx[1] = -1;
    cube->next_idx[2] = -1;
    cube->n_cubes = n_cubes;
    cube->n_stacks = n_stacks;
    cube->n_layers = n_layers;

    // Return success
    return(RUBIKS_CUBE_RETURN_SUCCESS);
}

// Destroy rubiks_cube struct
void rubiks_cube_destroy(struct rubiks_cube **cube_ptr){
    // Return if no valid rubiks_cube was provided
    if (cube_ptr == NULL || *cube_ptr == NULL) {
        return;
    }

    // Free memory of the 5 allocated arrays
    struct rubiks_cube *cube = *cube_ptr;
    free(cube->x);
    free(cube->y);
    free(cube->z);
    free(cube->density);
    free(cube->rock_id);

    // Free memory of rubiks_cube itself
    free(cube);

    // Assign NULL to cube pointer
    *cube_ptr = NULL;
}

// Find cube in which provided coordinates are located
enum rubiks_cube_return rubiks_cube_find_cube(double x, double y, double z,
                                              double dx, double dy, double dz,
                                              struct rubiks_cube *cube,
                                              double *x_ptr, double *y_ptr, double *z_ptr,
                                              double *density_ptr, int *rock_id_ptr,
                                              double *step_ptr){
    // Declare some variables
    unsigned long xi, yi, zi;
    _Bool x_flag, y_flag, z_flag;

    // Check if current cube holds the required X-coordinate
    if (dx >= 0) {
        x_flag = (cube->x[cube->cur_idx[0]]-cube->dx <= x) && (x < cube->x[cube->cur_idx[0]]+cube->dx);
    }
    else {
        x_flag = (cube->x[cube->cur_idx[0]]-cube->dx < x) && (x <= cube->x[cube->cur_idx[0]]+cube->dx);
    }

    // Check if current cube holds the required Y-coordinate
    if (dy >= 0) {
        y_flag = (cube->y[cube->cur_idx[1]].v-cube->dy <= y) && (y < cube->y[cube->cur_idx[1]].v+cube->dy);
    }
    else {
        y_flag = (cube->y[cube->cur_idx[1]].v-cube->dy < y) && (y <= cube->y[cube->cur_idx[1]].v+cube->dy);
    }

    // Check if current cube holds the required Z-coordinate
    if (dz >= 0) {
        z_flag = (cube->z[cube->cur_idx[2]].v-cube->dz <= z) && (z < cube->z[cube->cur_idx[2]].v+cube->dz);
    }
    else {
        z_flag = (cube->z[cube->cur_idx[2]].v-cube->dz < z) && (z <= cube->z[cube->cur_idx[2]].v+cube->dz);
    }

    // Check if all three flags are true
    if (x_flag && y_flag && z_flag) {
        // If so, the current cube is required
        xi = cube->cur_idx[0];
        yi = cube->cur_idx[1];
        zi = cube->cur_idx[2];
    }

    // Else, search for the cube normally
    else {
        // Check if the next predicted cube Z-coordinate is correct
        zi = cube->next_idx[2];
        if (dz >= 0) {
            // If particle is moving in positive direction, exclude upper limit of cube
            if (zi == -1 ||
                !((cube->z[zi].v-cube->dz <= z) && (z < cube->z[zi].v+cube->dz))) {
                // If not, find the cube that has the correct Z-coordinate
                zi = 0;
                while (zi < cube->n_layers &&
                       !((cube->z[zi].v-cube->dz <= z) && (z < cube->z[zi].v+cube->dz))) {
                    zi++;
                }
            }
        }
        else {
            // Else, exclude lower limit of cube
            if (zi == -1 ||
                !((cube->z[zi].v-cube->dz < z) && (z <= cube->z[zi].v+cube->dz))) {
                // If not, find the cube that has the correct Z-coordinate
                zi = 0;
                while (zi < cube->n_layers &&
                       !((cube->z[zi].v-cube->dz < z) && (z <= cube->z[zi].v+cube->dz))) {
                    zi++;
                }
            }
        }

        // If zi == n_layers, the requested cube does not exist
        if (zi == cube->n_layers) {
            return(RUBIKS_CUBE_RETURN_CUBE_NOT_FOUND);
        }

        // Determine how many Y-coordinates this Z-coordinate has
        unsigned long yj;
        if (zi < cube->n_layers-1) {
            yj = cube->z[zi+1].i;
        }
        else {
            yj = cube->n_stacks;
        }

        // Check if the next predicted cube Y-coordinate is correct
        yi = cube->next_idx[1];
        if (dy >= 0) {
            // If particle is moving in positive direction, exclude upper limit of cube
            if (yi == -1 ||
                !((cube->y[yi].v-cube->dy <= y) && (y < cube->y[yi].v+cube->dy))) {
                // If not, find the cube that has the correct Y-coordinate
                yi = cube->z[zi].i;
                while (yi < yj &&
                       !((cube->y[yi].v-cube->dy <= y) && (y < cube->y[yi].v+cube->dy))) {
                    yi++;
                }
            }
        }
        else {
            // Else, exclude lower limit of cube
            if (yi == -1 ||
                !((cube->y[yi].v-cube->dy < y) && (y <= cube->y[yi].v+cube->dy))) {
                // If not, find the cube that has the correct Y-coordinate
                yi = cube->z[zi].i;
                while (yi < yj &&
                       !((cube->y[yi].v-cube->dy < y) && (y <= cube->y[yi].v+cube->dy))) {
                    yi++;
                }
            }
        }

        // If yi == yj, the requested cube does not exist
        if (yi == yj) {
            return(RUBIKS_CUBE_RETURN_CUBE_NOT_FOUND);
        }

        // Determine how many X-coordinates this Y-coordinate has
        unsigned long xj;
        if (yi < cube->n_stacks-1) {
            xj = cube->y[yi+1].i;
        }
        else {
            xj = cube->n_cubes;
        }

        // Check if the next predicted cube X-coordinate is correct
        xi = cube->next_idx[0];
        if (dx >= 0) {
            // If particle is moving in positive direction, exclude upper limit of cube
            if (xi == -1 ||
                !((cube->x[xi]-cube->dx <= x) && (x < cube->x[xi]+cube->dx))) {
                // If not, find the cube that has the correct X-coordinate
                xi = cube->y[yi].i;
                while (xi < xj &&
                    !((cube->x[xi]-cube->dx <= x) && (x < cube->x[xi]+cube->dx))) {
                    xi++;
                }
            }
        }
        else {
            // Else, exclude lower limit of cube
            if (xi == -1 ||
                !((cube->x[xi]-cube->dx < x) && (x <= cube->x[xi]+cube->dx))) {
                // If not, find the cube that has the correct X-coordinate
                xi = cube->y[yi].i;
                while (xi < xj &&
                    !((cube->x[xi]-cube->dx < x) && (x <= cube->x[xi]+cube->dx))) {
                    xi++;
                }
            }
        }

        // If xi == xj, the requested cube does not exist
        if (xi == xj) {
            return(RUBIKS_CUBE_RETURN_CUBE_NOT_FOUND);
        }
    }

    // Save this index
    cube->cur_idx[0] = xi;
    cube->cur_idx[1] = yi;
    cube->cur_idx[2] = zi;

    // Obtain all requested values
    if (x_ptr != NULL) {
        *x_ptr = cube->x[xi];
    }
    if (y_ptr != NULL) {
        *y_ptr = cube->y[yi].v;
    }
    if (z_ptr != NULL) {
        *z_ptr = cube->z[zi].v;
    }
    if (density_ptr != NULL) {
        *density_ptr = cube->density[xi];
    }
    if (rock_id_ptr != NULL) {
        *rock_id_ptr = cube->rock_id[xi];
    }
    if (step_ptr != NULL) {
        double rx, ry, rz;

        // Determine how many times dx can fit in the remainder of this cube
        if (dx != 0) {
            rx = fabs(((cube->x[xi]+(dx/fabs(dx))*cube->dx)-x)/dx);
        }
        else {
            rx = DBL_MAX;
        }

        // Determine how many times dy can fit in the remainder of this cube
        if (dy != 0) {
            ry = fabs(((cube->y[yi].v+(dy/fabs(dy))*cube->dy)-y)/dy);
        }
        else {
            ry = DBL_MAX;
        }

        // Determine how many times dz can fit in the remainder of this cube
        if (dz != 0) {
            rz = fabs(((cube->z[zi].v+(dz/fabs(dz))*cube->dz)-z)/dz);
        }
        else {
            rz = DBL_MAX;
        }

        // Determine the smallest of these factors
        double min_r = min(rx, min(ry, rz));

        // Determine distance possible in each direction
        dx *= min_r;
        dy *= min_r;
        dz *= min_r;

        // Determine distance of vector made out of these distances and assign it
        *step_ptr = sqrt(pow(dx, 2)+pow(dy, 2)+pow(dz, 2));

        // Check if the next estimated Z-coordinate is different
        if (rz == min_r) {
            cube->next_idx[0] = -1;
            cube->next_idx[1] = -1;
            cube->next_idx[2] = max(0, min(zi+(int)(dz/fabs(dz)), cube->n_layers-1));
        }

        // Else, check if the next estimated Y-coordinate is different
        else if (ry == min_r) {
            cube->next_idx[0] = -1;
            cube->next_idx[1] = max(0, min(yi+(int)(dy/fabs(dy)), cube->n_stacks-1));
            cube->next_idx[2] = zi;
        }

        // Else, the next estimated X-coordinate is different
        else {
            cube->next_idx[0] = max(0, min(xi+(int)(dx/fabs(dx)), cube->n_cubes-1));
            cube->next_idx[1] = yi;
            cube->next_idx[2] = zi;
        }
    }

    // Return success
    return(RUBIKS_CUBE_RETURN_SUCCESS);
}
