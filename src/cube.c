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

// Declare struct for holding geometry data
struct rubiks_cube{
    // Dynamic array of cube center X-coordinates
    double *x;
    // Dynamic array of cube center Y-coordinate data
    double *y;
    // Dynamic array of cube center Z-coordinate data
    double *z;
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
    // Total number of x-values
    size_t nx;
    // Total number of y-values
    size_t ny;
    // Total number of z-values
    size_t nz;
};

// Private function for interpreting a line from a CSV file
static int interpret_CSV_line(char *line, double *x, double *y, double *z,
                              double *density, int *rock_id){
    // Initialize variables
    int ntok = 5;
    const char sep[2] = ",";
    char *field;

    // Obtain the x field
    field = strtok(line, sep);
    if (field == NULL) {
        return(ntok);
    }
    else {
        ntok--;
    }
    *x = strtod(field, NULL);

    // Obtain the y field
    field = strtok(NULL, sep);
    if (field == NULL) {
        return(ntok);
    }
    else {
        ntok--;
    }
    *y = strtod(field, NULL);

    // Obtain the z field
    field = strtok(NULL, sep);
    if (field == NULL) {
        return(ntok);
    }
    else {
        ntok--;
    }
    *z = strtod(field, NULL);

    // Obtain the density field
    field = strtok(NULL, sep);
    if (field == NULL) {
        return(ntok);
    }
    else {
        ntok--;
    }
    *density = strtod(field, NULL);

    // Obtain the rock_id field
    field = strtok(NULL, sep);
    if (field == NULL) {
        return(ntok);
    }
    else {
        ntok--;
    }
    *rock_id = (int)strtod(field, NULL);

    // Return number of not ok tokens
    return(ntok);
}

// Create and/or read the meta data file required for the rubiks_cube struct
static enum rubiks_cube_return read_meta_data(const char *file_path, size_t *nx_ptr,
                                              size_t *ny_ptr, size_t *nz_ptr){
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
        int rock_id_tmp, ntok;
        double x0, y0, z0, xi, yi, zi, xj, yj, zj, density_tmp;
        size_t nx, ny, nz;
        _Bool x_all, y_all, z_all;
        char * ret;

        // Initialize numbers to 1 (as the first line is not inspected)
        x_all = y_all = z_all = 0;
        nx = ny = nz = 1;

        // Read in first line of file
        ret = fgets(line, sizeof(line), file);
        if (ret == NULL) {
            // If this is the EOF, the file is empty
            fclose(file);
            fprintf(stderr, "ERROR: File is empty!\n");
            return(RUBIKS_CUBE_RETURN_EMPTY_FILE);
        }

        // Check if first character is part of a number
        if (!(line[0] >= '0' && line[0] <= '9')) {
            // If not, the file has a header
            fgets(line, sizeof(line), file);
        }

        // Interpret the first line
        ntok = interpret_CSV_line(line, &x0, &y0, &z0, &density_tmp, &rock_id_tmp);
        if (ntok) {
            goto format_error;
        }
        xi = x0;
        yi = y0;
        zi = z0;

        // Read new line
        ret = fgets(line, sizeof(line), file);

        // Read all remaining lines in the file and determine number of x, y and z values
        while (ret != NULL) {
            // Swap xj, yj, zj to xi, yi, zi
            xj = xi;
            yj = yi;
            zj = zi;

            // Interpret line
            ntok = interpret_CSV_line(line, &xi, &yi, &zi, &density_tmp, &rock_id_tmp);
            if (ntok) {
                goto format_error;
            }

            // Check if a new x-value was read in
            if (!x_all && (xi != xj)) {
                // Check if this x-value is not the first value read in
                if (xi == x0) {
                    // If so, set flag to True
                    x_all = 1;
                }
                else {
                    // Else, increase nx by 1
                    nx++;
                }
            }

            // Check if a new y-value was read in
            if (!y_all && (yi != yj)) {
                // Check if this y-value is not the first value read in
                if (yi == y0) {
                    // If so, set flag to True
                    y_all = 1;
                }
                else {
                    // Else, increase ny by 1
                    ny++;
                }
            }

            // Check if a new z-value was read in
            if (!z_all && (zi != zj)) {
                // Check if this z-value is not the first value read in
                if (zi == z0) {
                    // If so, set flag to True
                    z_all = 1;
                }
                else {
                    // Else, increase nz by 1
                    nz++;
                }
            }

            // Read new line
            ret = fgets(line, sizeof(line), file);
        }

        // Close the data file
        fclose(file);
        goto write_metadata;

        // If data reading went wrong in any way
        format_error:
            // Close file and raise error
            fclose(file);
            fprintf(stderr, "ERROR: Invalid file format (%s)!\n", file_path);
            return(RUBIKS_CUBE_RETURN_FORMAT_ERROR);

        // If data reading went right
        write_metadata:
            // Write the meta data to file
            file = fopen(file_path_meta, "w");
            fprintf(file, "%zu,%zu,%zu", nx, ny, nz);
            fclose(file);

        // Assign values to proper pointers
        *nx_ptr = nx;
        *ny_ptr = ny;
        *nz_ptr = nz;
    }
    else {
        // If so, read it
        FILE *file_meta = fopen(file_path_meta, "r");
        fgets(line, sizeof(line), file_meta);
        fclose(file_meta);

        // Initialize variables
        const char sep[2] = ",";
        char *field;

        // Obtain the x field
        field = strtok(line, sep);
        sscanf(field, "%zu", nx_ptr);

        // Obtain the y field
        field = strtok(NULL, sep);
        sscanf(field, "%zu", ny_ptr);

        // Obtain the z field
        field = strtok(NULL, sep);
        sscanf(field, "%zu", nz_ptr);
    }

    // Return that everything went okay
    return(RUBIKS_CUBE_RETURN_SUCCESS);
}

// Initialize rubiks_cube struct from file
enum rubiks_cube_return rubiks_cube_create(struct rubiks_cube **cube_ptr, 
                                           const char *file_path){
    // Initialize size variables
    size_t nx, ny, nz;

    // Obtain values for these variables
    enum rubiks_cube_return rc = read_meta_data(file_path, &nx, &ny, &nz);
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
    cube->x = (double *)malloc(sizeof(double)*nx);
    cube->y = (double *)malloc(sizeof(double)*ny);
    cube->z = (double *)malloc(sizeof(double)*nz);
    cube->density = (double *)malloc(sizeof(double)*nx*ny*nz);
    cube->rock_id = (int *)malloc(sizeof(int)*nx*ny*nz);

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
    int xi, yi, zi;
    double x_tmp, y_tmp, z_tmp;
    _Bool x_all, y_all, z_all;
    char line[80];

    // Initialize values of *_all
    xi = yi = zi = 0;
    x_all = y_all = z_all = 0;

    // Open the data file
    FILE *file = fopen(file_path, "r");

    // Read in first line of file
    fgets(line, sizeof(line), file);

    // Check if first character is part of a number
    if (!(line[0] >= '0' && line[0] <= '9')) {
        // If not, the file has a header, so read in the actual first line
        fgets(line, sizeof(line), file);
    }

    // Interpret the first line of the file outside of the following loop
    interpret_CSV_line(line, &cube->x[0], &cube->y[0], &cube->z[0],
                       &cube->density[0], &cube->rock_id[0]);

    // Set the values of dx, dy and dz for checking later
    cube->dx = cube->dy = cube->dz = 0;

    // Loop over all remaining lines in the file and read in their contents
    for (i=1; i<nx*ny*nz; i++) {
        // Read line and assign
        fgets(line, sizeof(line), file);
        interpret_CSV_line(line, &x_tmp, &y_tmp, &z_tmp,
                           &cube->density[i], &cube->rock_id[i]);

        // Check if a new x-value was read in
        if (!x_all && (x_tmp != cube->x[xi])) {
            // Check if this x-value is not the first value read in
            if (x_tmp == cube->x[0]) {
                // If so, set flag to True
                x_all = 1;
            }
            else {
                // Else, store value
                xi++;
                cube->x[xi] = x_tmp;

                // Attempt to determine dx if not done before
                if (cube->dx == 0) {
                    cube->dx = (cube->x[xi]-cube->x[xi-1])/2;
                }
            }
        }

        // Check if a new y-value was read in
        if (!y_all && (y_tmp != cube->y[yi])) {
            // Check if this y-value is not the first value read in
            if (y_tmp == cube->y[0]) {
                // If so, set flag to True
                y_all = 1;
            }
            else {
                // Else, store value
                yi++;
                cube->y[yi] = y_tmp;

                // Attempt to determine dy if not done before
                if (cube->dy == 0) {
                    cube->dy = (cube->y[yi]-cube->y[yi-1])/2;
                }
            }
        }

        // Check if a new z-value was read in
        if (!z_all && (z_tmp != cube->z[zi])) {
            // Check if this z-value is not the first value read in
            if (z_tmp == cube->z[0]) {
                // If so, set flag to True
                z_all = 1;
            }
            else {
                // Else, store value
                zi++;
                cube->z[zi] = z_tmp;

                // Attempt to determine dz if not done before
                if (cube->dz == 0) {
                    cube->dz = (cube->z[zi]-cube->z[zi-1])/2;
                }
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
    cube->nx = nx;
    cube->ny = ny;
    cube->nz = nz;

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

    // Calculate the index of the X-coordinate of the cube
    xi = (unsigned long)floor((x-(cube->x[0]-cube->dx))/(2*cube->dx));

    // If particle is moving in negative X-direction, decrease xi by 1
    if (dx < 0) {
        xi--;
    }

    // If xi < 0 or xi >= nx, the requested cube does not exist
    if (xi < 0 || xi >= cube->nx) {
        return(RUBIKS_CUBE_RETURN_CUBE_NOT_FOUND);
    }

    // Calculate the index of the Y-coordinate of the cube
    yi = (unsigned long)floor((y-(cube->y[0]-cube->dy))/(2*cube->dy));

    // If particle is moving in negative Y-direction, decrease yi by 1
    if (dy < 0) {
        yi--;
    }

    // If yi < 0 or yi >= ny, the requested cube does not exist
    if (yi < 0 || yi >= cube->ny) {
        return(RUBIKS_CUBE_RETURN_CUBE_NOT_FOUND);
    }

    // Calculate the index of the Z-coordinate of the cube
    zi = (unsigned long)floor((z-(cube->z[0]-cube->dz))/(2*cube->dz));

    // If particle is moving in negative Z-direction, decrease zi by 1
    if (dz < 0) {
        zi--;
    }

    // If zi < 0 or zi >= nz, the requested cube does not exist
    if (zi < 0 || zi >= cube->nz) {
        return(RUBIKS_CUBE_RETURN_CUBE_NOT_FOUND);
    }

    // Obtain all requested values
    if (x_ptr != NULL) {
        *x_ptr = cube->x[xi];
    }
    if (y_ptr != NULL) {
        *y_ptr = cube->y[yi];
    }
    if (z_ptr != NULL) {
        *z_ptr = cube->z[zi];
    }
    if (density_ptr != NULL) {
        *density_ptr = cube->density[zi*cube->nx*cube->ny+yi*cube->nx+xi];
    }
    if (rock_id_ptr != NULL) {
        *rock_id_ptr = cube->rock_id[zi*cube->nx*cube->ny+yi*cube->nx+xi];
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
            ry = fabs(((cube->y[yi]+(dy/fabs(dy))*cube->dy)-y)/dy);
        }
        else {
            ry = DBL_MAX;
        }

        // Determine how many times dz can fit in the remainder of this cube
        if (dz != 0) {
            rz = fabs(((cube->z[zi]+(dz/fabs(dz))*cube->dz)-z)/dz);
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
    }

    // Return success
    return(RUBIKS_CUBE_RETURN_SUCCESS);
}
