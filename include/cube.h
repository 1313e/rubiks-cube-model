/*
Copyright (C) 2021 Ellert van der Velden
All rights reserved.

This software is free software, distributed under the BSD-3 License.
You may redistribute and/or modify it without any restrictions, as long as the conditions specified in the terms of the BSD-3 license (included) are met.
*/

#ifndef cube_h
#define cube_h

// Versioning
#ifndef RUBIKS_CUBE_VERSION
#define RUBIKS_CUBE_VERSION 1.0
#endif

// Coordinate data struct declaration
struct coord_data;

// Rubiks cube struct declaration
struct rubiks_cube;

// Return codes
enum rubiks_cube_return{
    // Execution was successful
    RUBIKS_CUBE_RETURN_SUCCESS = 0,

    // Cube file could not be found
    RUBIKS_CUBE_RETURN_FILE_NOT_FOUND,

    // Cube memory could not be allocated
    RUBIKS_CUBE_RETURN_MALLOC_ERROR,

    // Cube file is empty
    RUBIKS_CUBE_RETURN_EMPTY_FILE,

    // Cube file has wrong format
    RUBIKS_CUBE_RETURN_FORMAT_ERROR,

    // Cube is not contained within rubiks cube
    RUBIKS_CUBE_RETURN_CUBE_NOT_FOUND
};

/**
 * Initialize rubiks_cube struct from file.
 *
 * @param cube      Pointer where to store the Rubik's cube data.
 * @param file_path Path to file to read in.
 * @return          On success, `RUBIKS_CUBE_RETURN_SUCCESS` is returned and an error code is returned otherwise.
 * __Notes__
 * 
 * The `rubiks_cube` struct assumes that all cubes have the same dimensions.
 * It is also assumed that the cubes in the data file are ordered such that the X-coordinate changes the fastest; then the Y-coordinate; and then the Z-coordinate.
 *
 * __Error codes__
 *
 *     RUBIKS_CUBE_RETURN_FILE_NOT_FOUND    Cube file could not be found
 *
 *     RUBIKS_CUBE_RETURN_MALLOC_ERROR      Cube memory could not be allocated
 *
 *     RUBIKS_CUBE_RETURN_EMPTY_FILE        Cube file is empty
 *
 *     RUBIKS_CUBE_RETURN_FORMAT_ERROR      Cube file has wrong format
 */
enum rubiks_cube_return rubiks_cube_create(
    struct rubiks_cube **cube, const char *file_path);

// Declaration of rubiks_cube_destroy
void rubiks_cube_destroy(
    struct rubiks_cube **cube);

/**
 * Get primary altitude of rubiks_cube struct.
 *
 * @param cube      The Rubik's cube.
 * @param alt_ptr   Pointer where to store the primary altitude.
 */
void rubiks_cube_primary_altitude(
    struct rubiks_cube *cube, double *alt_ptr);

/**
 * Find cube in which provided coordinates are located.
 *
 * @param x             X-coordinate of particle.
 * @param y             Y-coordinate of particle.
 * @param z             Z-coordinate of particle.
 * @param dx            X-component of direction vector of particle.
 * @param dy            Y-component of direction vector of particle.
 * @param dz            Z-component of direction vector of particle.
 * @param cube          The Rubik's cube.
 * @param x_ptr         Pointer where to store X-coordinate of center of current cube or `NULL`.
 * @param y_ptr         Pointer where to store Y-coordinate of center of current cube or `NULL`.
 * @param z_ptr         Pointer where to store Z-coordinate of center of current cube or `NULL`.
 * @param density_ptr   Pointer where to store density of current cube or `NULL`.
 * @param rock_id_ptr   Pointer where to store medium/rock ID of current cube or `NULL`.
 * @param step_ptr      Pointer where to store distance towards edge of current cube given its direction vector or `NULL`.
 * @return              On success, `RUBIKS_CUBE_RETURN_SUCCESS` is returned and an error code is returned otherwise.
 *
 * If an argument is not required to be retrieved, set it to `NULL`.
 *
 * Note that if one uses backward transport instead of forward, the direction vector should point in the direction the particle came from.
 *
 * __Error codes__
 * 
 *     RUBIKS_CUBE_RETURN_CUBE_NOT_FOUND    Cube containing provided coordinates does not exist.
 */
enum rubiks_cube_return rubiks_cube_find_cube(
    double x, double y, double z, double dx, double dy, double dz,
    struct rubiks_cube *cube, double *x_ptr, double *y_ptr, double *z_ptr,
    double *density_ptr, int *rock_id_ptr, double *step_ptr);
#endif
