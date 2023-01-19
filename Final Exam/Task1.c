// ----- This is a file for Task 1 -----

// --- Part A ---
/*
An example of a non-transcendental number that would not be represented accurately with bfloat16 is +5.55E-130.
This is because the number has an exponent that exceeds the bounds of the 8 bits used in bfloat16.
(The number is too small to be represented)
*/

// --- Part B ---
// Assuming we can't use math.h
#include "math.h"
#include <iostream>

typedef short int int16_t;

// Resolution of Taylor Series approximation
// Some relatively high value should suffice
#define ITERS 16

int16_t sin16(int16_t x)
{
    int16_t result = x;
    // Compute the Taylor series expansion of sine
    for (int16_t i = 1; i < ITERS; i++)
    {
        int16_t termVal = (2 * i) + 1; // Odd numbers for 3, 5, 7, 9, ...
        int16_t numerator = 1;
        int16_t denom = 1;

        // Compute numerator
        for (int16_t i = 0; i < termVal; i++)
            numerator *= x;
        // Compute denominator
        for (int16_t i = termVal; i > 0; i--)
            denom *= i;

        // Add to the result
        if (i % 2 == 0)
        {
            result += (numerator / denom);
        }
        else
        {
            result -= (numerator / denom);
        }
    }
    return result;
}

int16_t cos16(int16_t x)
{
    int16_t result = x;
    // Compute the Taylor series expansion of sine
    for (int16_t i = 1; i < ITERS; i++)
    {
        int16_t termVal = (2 * i); // Even numbers for 2, 4, 6, 8, ...
        int16_t numerator = 1;
        int16_t denom = 1;

        // Compute numerator
        for (int16_t i = 0; i < termVal; i++)
            numerator *= x;
        // Compute denominator
        for (int16_t i = termVal; i > 0; i--)
            denom *= i;
        // Add to the result
        if (i % 2 == 0)
        {
            result += (numerator / denom);
        }
        else
        {
            result -= (numerator / denom);
        }
    }
    return result;
}

int main()
{
    std::cout << cos16(45) << " -> " << cos16(45) << std::endl;
    return EXIT_SUCCESS;
}
