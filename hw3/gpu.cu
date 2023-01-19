#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <cuda.h>
#include "common.h"

#define NUM_THREADS 256

extern double size;
//
//  benchmarking program
//

__device__ void apply_force_gpu(particle_t &particle, particle_t &neighbor)
{
    double dx = neighbor.x - particle.x;
    double dy = neighbor.y - particle.y;
    double r2 = dx * dx + dy * dy;
    if (r2 > cutoff * cutoff)
        return;
    // r2 = fmax( r2, min_r*min_r );
    r2 = (r2 > min_r * min_r) ? r2 : min_r * min_r;
    double r = sqrt(r2);

    //
    //  very simple short-range repulsive force
    //
    double coef = (1 - cutoff / r) / r2 / mass;
    particle.ax += coef * dx;
    particle.ay += coef * dy;
}

__global__ void compute_forces_gpu(particle_t *particles, int n, particle_t *bin_grid, int *bin_grid_count)
{
    // Get thread (particle) ID
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= n)
        return;

    particles[tid].ax = particles[tid].ay = 0;

    int sidelength = get_sidelength();

    // Compute what i and j values the particle is at
    double iIndex = 0;
    double jIndex = 0;

    double lowVal = 0;
    double highVal = cutoff;
    for (int j = 0; j < sidelength; j++)
    {
        if (particles[tid].x >= lowVal && particles[tid].x <= highVal)
        {
            iIndex = j;
        }
        if (particles[tid].y >= lowVal && particles[tid].y <= highVal)
        {
            jIndex = j;
        }

        lowVal += cutoff;
        highVal += cutoff;
    }

    // Top left bin
    if (
        (iIndex - 1) >= 0 &&
        (iIndex - 1) < sidelength &&
        (jIndex - 1) >= 0 &&
        (jIndex - 1) < sidelength)
    {
        int index2D = get2DIndex(iIndex - 1, jIndex - 1);
        int count = bin_grid_count[index2D];
        int index3D = get3DIndex(iIndex - 1, jIndex - 1, count);
        // Iterate through all other particles of that bin
        for (int j = 0; j < count; j++)
        {
            apply_force_gpu(
                particles[tid],
                bin_grid[index3D]);
        }
    }
    // Left bin
    if (
        (iIndex - 1) >= 0 &&
        (iIndex - 1) < sidelength &&
        (jIndex) >= 0 &&
        (jIndex) < sidelength)
    {
        int index2D = get2DIndex(iIndex - 1, jIndex);
        int count = bin_grid_count[index2D];
        int index3D = get3DIndex(iIndex - 1, jIndex, count);
        // Iterate through all other particles of that bin
        for (int j = 0; j < bin_grid.at(iIndex - 1).at(jIndex).size(); j++)
        {
            apply_force_gpu(
                particles[tid],
                bin_grid[index3D]);
        }
    }
    // Bot left bin
    if (
        (iIndex - 1) >= 0 &&
        (iIndex - 1) < sidelength &&
        (jIndex + 1) >= 0 &&
        (jIndex + 1) < sidelength)
    {
        int index2D = get2DIndex(iIndex - 1, jIndex + 1);
        int count = bin_grid_count[index2D];
        int index3D = get3DIndex(iIndex - 1, jIndex + 1, count);
        // Iterate through all other particles of that bin
        for (int j = 0; j < bin_grid.at(iIndex - 1).at(jIndex + 1).size(); j++)
        {
            apply_force_gpu(
                particles[tid],
                bin_grid[index3D]);
        }
    }
    // Top bin
    if (
        (iIndex) >= 0 &&
        (iIndex) < sidelength &&
        (jIndex - 1) >= 0 &&
        (jIndex - 1) < sidelength)
    {
        int index2D = get2DIndex(iIndex, jIndex - 1);
        int count = bin_grid_count[index2D];
        int index3D = get3DIndex(iIndex, jIndex - 1, count);
        // Iterate through all other particles of that bin
        for (int j = 0; j < bin_grid.at(iIndex).at(jIndex - 1).size(); j++)
        {
            apply_force_gpu(
                particles[tid],
                bin_grid[index3D]);
        }
    }
    // Top right bin
    if (
        (iIndex + 1) >= 0 &&
        (iIndex + 1) < sidelength &&
        (jIndex - 1) >= 0 &&
        (jIndex - 1) < sidelength)
    {
        int index2D = get2DIndex(iIndex + 1, jIndex - 1);
        int count = bin_grid_count[index2D];
        int index3D = get3DIndex(iIndex + 1, jIndex - 1, count);
        // Iterate through all other particles of that bin
        for (int j = 0; j < bin_grid.at(iIndex + 1).at(jIndex - 1).size(); j++)
        {
            apply_force_gpu(
                particles[tid],
                bin_grid[index3D]);
        }
    }
    // Right bin
    if (
        (iIndex + 1) >= 0 &&
        (iIndex + 1) < sidelength &&
        (jIndex) >= 0 &&
        (jIndex) < sidelength)
    {
        int index2D = get2DIndex(iIndex + 1, jIndex);
        int count = bin_grid_count[index2D];
        int index3D = get3DIndex(iIndex + 1, jIndex, count);
        // Iterate through all other particles of that bin
        for (int j = 0; j < bin_grid.at(iIndex + 1).at(jIndex).size(); j++)
        {
            apply_force_gpu(
                particles[tid],
                bin_grid[index3D]);
        }
    }
    // Bot right bin
    if (
        (iIndex + 1) >= 0 &&
        (iIndex + 1) < sidelength &&
        (jIndex + 1) >= 0 &&
        (jIndex + 1) < sidelength)
    {
        int index2D = get2DIndex(iIndex + 1, jIndex + 1);
        int count = bin_grid_count[index2D];
        int index3D = get3DIndex(iIndex + 1, jIndex + 1, count);
        // Iterate through all other particles of that bin
        for (int j = 0; j < bin_grid.at(iIndex + 1).at(jIndex + 1).size(); j++)
        {
            apply_force_gpu(
                particles[tid],
                bin_grid[index3D]);
        }
    }
    // Bot bin
    if (
        (iIndex) >= 0 &&
        (iIndex) < sidelength &&
        (jIndex + 1) >= 0 &&
        (jIndex + 1) < sidelength)
    {
        int index2D = get2DIndex(iIndex, jIndex + 1);
        int count = bin_grid_count[index2D];
        int index3D = get3DIndex(iIndex, jIndex + 1, count);
        // Iterate through all other particles of that bin
        for (int j = 0; j < bin_grid.at(iIndex).at(jIndex + 1).size(); j++)
        {
            apply_force_gpu(
                particles[tid],
                bin_grid[index3D]);
        }
    }
    // Middle bin
    if (
        (iIndex) >= 0 &&
        (iIndex) < sidelength &&
        (jIndex) >= 0 &&
        (jIndex) < sidelength)
    {
        int index2D = get2DIndex(iIndex, jIndex);
        int count = bin_grid_count[index2D];
        int index3D = get3DIndex(iIndex, jIndex, count);
        // Iterate through all other particles of that bin
        for (int j = 0; j < bin_grid.at(iIndex).at(jIndex).size(); j++)
        {
            apply_force_gpu(
                particles[tid],
                bin_grid[index3D]);
        }
    }

    // for(int j = 0 ; j < n ; j++) {
    //     apply_force_gpu(particles[tid], particles[j]);
    // }
}

__global__ void move_gpu(particle_t *particles, int n, double size)
{

    // Get thread (particle) ID
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid >= n)
        return;

    particle_t *p = &particles[tid];
    //
    //  slightly simplified Velocity Verlet integration
    //  conserves energy better than explicit Euler method
    //
    p->vx += p->ax * dt;
    p->vy += p->ay * dt;
    p->x += p->vx * dt;
    p->y += p->vy * dt;

    //
    //  bounce from walls
    //
    while (p->x < 0 || p->x > size)
    {
        p->x = p->x < 0 ? -(p->x) : 2 * size - p->x;
        p->vx = -(p->vx);
    }
    while (p->y < 0 || p->y > size)
    {
        p->y = p->y < 0 ? -(p->y) : 2 * size - p->y;
        p->vy = -(p->vy);
    }
}

int main(int argc, char **argv)
{
    // This takes a few seconds to initialize the runtime
    cudaDeviceSynchronize();

    if (find_option(argc, argv, "-h") >= 0)
    {
        printf("Options:\n");
        printf("-h to see this help\n");
        printf("-n <int> to set the number of particles\n");
        printf("-o <filename> to specify the output file name\n");
        return 0;
    }

    int n = read_int(argc, argv, "-n", 1000);

    char *savename = read_string(argc, argv, "-o", NULL);

    FILE *fsave = savename ? fopen(savename, "w") : NULL;
    particle_t *particles = (particle_t *)malloc(n * sizeof(particle_t));

    // GPU particle data structure
    particle_t *d_particles;
    cudaMalloc((void **)&d_particles, n * sizeof(particle_t));

    set_size(n);

    init_particles(n, particles);

    // Initialize our bins for O(n) shortcut checking
    int sidelength = get_sidelength();
    particle_t *d_bin_grid; // device bin grid
    cudaMalloc(
        (void **) % d_bin_grid,
        sidelength * sidelength * n); // A bin may have up to n particles, so we'll set the inner size to n
    // This may be horribly inefficient
    int *d_bin_grid_count; // how many particles are in each bin
    cudaMalloc(
        (void **) % d_bin_grid_count,
        sidelength * sidelength);

    cudaDeviceSynchronize();
    double copy_time = read_timer();

    // Copy the particles to the GPU
    cudaMemcpy(d_particles, particles, n * sizeof(particle_t), cudaMemcpyHostToDevice);

    cudaDeviceSynchronize();
    copy_time = read_timer() - copy_time;

    //
    //  simulate a number of time steps
    //
    cudaDeviceSynchronize();
    double simulation_time = read_timer();

    for (int step = 0; step < NSTEPS; step++)
    {
        // Zero out all the matrices
        for (int i = 0; i < sidelength; i++)
        {
            for (int j = 0; j < sidelength; j++)
            {
                for (int k = 0; k < n; k++)
                {
                    int index3D = get3DIndex(i, j, k);
                    d_bin_grid[index3D].x = -1;
                    d_bin_grid[index3D].y = -1;
                }

                int index2D = get2DIndex(i, j);
                d_bin_grid_count[index2D] = 0;
            }
        }

        // Add all particles to the matrix in their proper grid locations
        for (int i = 0; i < n; i++)
        {
            // Compute what i and j values the particle is at
            double iIndex = 0;
            double jIndex = 0;

            double lowVal = 0;
            double highVal = cutoff;
            for (int j = 0; j < sidelength; j++)
            {
                if (d_particles[i].x >= lowVal && d_particles[i].x <= highVal)
                {
                    iIndex = j;
                }
                if (d_particles[i].y >= lowVal && d_particles[i].y <= highVal)
                {
                    jIndex = j;
                }

                lowVal += cutoff;
                highVal += cutoff;
            }

            // Place the particle in the grid
            int index2D = get2DIndex(iIndex, jIndex);
            int countForBin = d_bin_grid_count[index2D];
            int index3D = get3DIndex(iIndex, jIndex, countForBin);
            d_bin_grid[index3D] = d_particles[i]; // Place this particle
            d_bin_grid_count[index2D] += 1;       // Count this particle
        }

        //
        //  compute forces
        //

        int blks = (n + NUM_THREADS - 1) / NUM_THREADS;
        compute_forces_gpu<<<blks, NUM_THREADS>>>(d_particles, n, d_bin_grid, d_bin_grid_count);

        //
        //  move particles
        //
        move_gpu<<<blks, NUM_THREADS>>>(d_particles, n, size);

        //
        //  save if necessary
        //
        if (fsave && (step % SAVEFREQ) == 0)
        {
            // Copy the particles back to the CPU
            cudaMemcpy(particles, d_particles, n * sizeof(particle_t), cudaMemcpyDeviceToHost);
            save(fsave, n, particles);
        }
    }
    cudaDeviceSynchronize();
    simulation_time = read_timer() - simulation_time;

    printf("CPU-GPU copy time = %g seconds\n", copy_time);
    printf("n = %d, simulation time = %g seconds\n", n, simulation_time);

    cudaFree(d_bin_grid);
    cudaFree(d_bin_grid_count);

    free(particles);
    cudaFree(d_particles);
    if (fsave)
        fclose(fsave);

    return 0;
}
