#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include <vector>
#include <iostream>

//
//  benchmarking program
//
int main(int argc, char **argv)
{
    int navg, nabsavg = 0;
    double davg, dmin, absmin = 1.0, absavg = 0.0;

    if (find_option(argc, argv, "-h") >= 0)
    {
        printf("Options:\n");
        printf("-h to see this help\n");
        printf("-n <int> to set the number of particles\n");
        printf("-o <filename> to specify the output file name\n");
        printf("-s <filename> to specify a summary file name\n");
        printf("-no turns off all correctness checks and particle output\n");
        return 0;
    }

    int n = read_int(argc, argv, "-n", 1000);

    char *savename = read_string(argc, argv, "-o", NULL);
    char *sumname = read_string(argc, argv, "-s", NULL);

    FILE *fsave = savename ? fopen(savename, "w") : NULL;
    FILE *fsum = sumname ? fopen(sumname, "a") : NULL;

    particle_t *particles = (particle_t *)malloc(n * sizeof(particle_t));
    set_size(n);
    init_particles(n, particles);

    double cutoff = get_cutoff();

    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer();

    for (int step = 0; step < NSTEPS; step++)
    {
        navg = 0;
        davg = 0.0;
        dmin = 1.0;

        // Bins of even spacing
        std::vector<std::vector<std::vector<particle_t> > > bin_grid;

        // Initialize the bin grid
        int sidelength = get_sidelength();
        for (int i = 0; i < sidelength; i++)
        {
            bin_grid.push_back(std::vector< std::vector<particle_t> >());
            for (int j = 0; j < sidelength; j++)
            {
                bin_grid.at(i).push_back(std::vector<particle_t>());
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
                if (particles[i].x >= lowVal && particles[i].x <= highVal)
                {
                    iIndex = j;
                }
                if (particles[i].y >= lowVal && particles[i].y <= highVal)
                {
                    jIndex = j;
                }

                lowVal += cutoff;
                highVal += cutoff;
            }

            // Place the particle in the grid
            bin_grid.at(iIndex).at(jIndex).push_back(particles[i]);
        }

        //
        //  compute forces
        //
        for (int i = 0; i < n; i++)
        {
            particles[i].ax = particles[i].ay = 0;

            // Compute what i and j values the particle is at
            double iIndex = 0;
            double jIndex = 0;

            double lowVal = 0;
            double highVal = cutoff;
            for (int j = 0; j < sidelength; j++)
            {
                if (particles[i].x >= lowVal && particles[i].x <= highVal)
                {
                    iIndex = j;
                }
                if (particles[i].y >= lowVal && particles[i].y <= highVal)
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
                // Iterate through all other particles of that bin
                for (int j = 0; j < bin_grid.at(iIndex - 1).at(jIndex - 1).size(); j++)
                {
                    apply_force(
                        particles[i],
                        bin_grid.at(iIndex - 1).at(jIndex - 1).at(j),
                        &dmin,
                        &davg,
                        &navg);
                }
            }
            // Left bin
            if (
                (iIndex - 1) >= 0 &&
                (iIndex - 1) < sidelength &&
                (jIndex) >= 0 &&
                (jIndex) < sidelength)
            {
                // Iterate through all other particles of that bin
                for (int j = 0; j < bin_grid.at(iIndex - 1).at(jIndex).size(); j++)
                {
                    apply_force(
                        particles[i],
                        bin_grid.at(iIndex - 1).at(jIndex).at(j),
                        &dmin,
                        &davg,
                        &navg);
                }
            }
            // Bot left bin
            if (
                (iIndex - 1) >= 0 &&
                (iIndex - 1) < sidelength &&
                (jIndex + 1) >= 0 &&
                (jIndex + 1) < sidelength)
            {
                // Iterate through all other particles of that bin
                for (int j = 0; j < bin_grid.at(iIndex - 1).at(jIndex + 1).size(); j++)
                {
                    apply_force(
                        particles[i],
                        bin_grid.at(iIndex - 1).at(jIndex + 1).at(j),
                        &dmin,
                        &davg,
                        &navg);
                }
            }
            // Top bin
            if (
                (iIndex) >= 0 &&
                (iIndex) < sidelength &&
                (jIndex - 1) >= 0 &&
                (jIndex - 1) < sidelength)
            {
                // Iterate through all other particles of that bin
                for (int j = 0; j < bin_grid.at(iIndex).at(jIndex - 1).size(); j++)
                {
                    apply_force(
                        particles[i],
                        bin_grid.at(iIndex).at(jIndex - 1).at(j),
                        &dmin,
                        &davg,
                        &navg);
                }
            }
            // Top right bin
            if (
                (iIndex + 1) >= 0 &&
                (iIndex + 1) < sidelength &&
                (jIndex - 1) >= 0 &&
                (jIndex - 1) < sidelength)
            {
                // Iterate through all other particles of that bin
                for (int j = 0; j < bin_grid.at(iIndex + 1).at(jIndex - 1).size(); j++)
                {
                    apply_force(
                        particles[i],
                        bin_grid.at(iIndex + 1).at(jIndex - 1).at(j),
                        &dmin,
                        &davg,
                        &navg);
                }
            }
            // Right bin
            if (
                (iIndex + 1) >= 0 &&
                (iIndex + 1) < sidelength &&
                (jIndex) >= 0 &&
                (jIndex) < sidelength)
            {
                // Iterate through all other particles of that bin
                for (int j = 0; j < bin_grid.at(iIndex + 1).at(jIndex).size(); j++)
                {
                    apply_force(
                        particles[i],
                        bin_grid.at(iIndex + 1).at(jIndex).at(j),
                        &dmin,
                        &davg,
                        &navg);
                }
            }
            // Bot right bin
            if (
                (iIndex + 1) >= 0 &&
                (iIndex + 1) < sidelength &&
                (jIndex + 1) >= 0 &&
                (jIndex + 1) < sidelength)
            {
                // Iterate through all other particles of that bin
                for (int j = 0; j < bin_grid.at(iIndex + 1).at(jIndex + 1).size(); j++)
                {
                    apply_force(
                        particles[i],
                        bin_grid.at(iIndex + 1).at(jIndex + 1).at(j),
                        &dmin,
                        &davg,
                        &navg);
                }
            }
            // Bot bin
            if (
                (iIndex) >= 0 &&
                (iIndex) < sidelength &&
                (jIndex + 1) >= 0 &&
                (jIndex + 1) < sidelength)
            {
                // Iterate through all other particles of that bin
                for (int j = 0; j < bin_grid.at(iIndex).at(jIndex + 1).size(); j++)
                {
                    apply_force(
                        particles[i],
                        bin_grid.at(iIndex).at(jIndex + 1).at(j),
                        &dmin,
                        &davg,
                        &navg);
                }
            }
            // Middle bin
                if (
                (iIndex) >= 0 &&
                (iIndex) < sidelength &&
                (jIndex) >= 0 &&
                (jIndex) < sidelength)
            {
                // Iterate through all other particles of that bin
                for (int j = 0; j < bin_grid.at(iIndex).at(jIndex).size(); j++)
                {
                    apply_force(
                        particles[i],
                        bin_grid.at(iIndex).at(jIndex).at(j),
                        &dmin,
                        &davg,
                        &navg);
                }
            }
            // for (int j = 0; j < n; j++)
            //     apply_force(particles[i], particles[j], &dmin, &davg, &navg);
        }

        //
        //  move particles
        //
        for (int i = 0; i < n; i++)
            move(particles[i]);

        if (find_option(argc, argv, "-no") == -1)
        {
            //
            // Computing statistical data
            //
            if (navg)
            {
                absavg += davg / navg;
                nabsavg++;
            }
            if (dmin < absmin)
                absmin = dmin;

            //
            //  save if necessary
            //
            if (fsave && (step % SAVEFREQ) == 0)
                save(fsave, n, particles);
        }
    }
    simulation_time = read_timer() - simulation_time;

    printf("n = %d, simulation time = %g seconds", n, simulation_time);

    if (find_option(argc, argv, "-no") == -1)
    {
        if (nabsavg)
            absavg /= nabsavg;
        //
        //  -The minimum distance absmin between 2 particles during the run of the simulation
        //  -A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
        //  -A simulation where particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
        //
        //  -The average distance absavg is ~.95 when most particles are interacting correctly and ~.66 when no particles are interacting
        //
        printf(", absmin = %lf, absavg = %lf", absmin, absavg);
        if (absmin < 0.4)
            printf("\nThe minimum distance is below 0.4 meaning that some particle is not interacting");
        if (absavg < 0.8)
            printf("\nThe average distance is below 0.8 meaning that most particles are not interacting");
    }
    printf("\n");

    //
    // Printing summary data
    //
    if (fsum)
        fprintf(fsum, "%d %g\n", n, simulation_time);

    //
    // Clearing space
    //
    if (fsum)
        fclose(fsum);
    free(particles);
    if (fsave)
        fclose(fsave);

    return 0;
}