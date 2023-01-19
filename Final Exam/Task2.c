// ----- This is a file for Task 2 -----

#define N 1024
typedef struct
{
    float real;
    float imag;
} complex_t;
// Shared mem
static complex_t FFT_input[N];

float cabs(complex_t x)
{
    return (sqrt(x.real * x.real + x.imag * x.imag));
}

void processorA()
{
    int sum = 0;
    for (int i = 0; i < N; i++)
    {
        sum += cabs(FFT_input[i]);
    }
}

void processorB()
{
    int max = 0;
    for (int i = 0; i < N; i++)
    {
        if (cabs(FFT_input[i]) > max)
            max = FFT_input[i];
    }
}

void processorC()
{
    for (int i = 0; i < N; i++)
    {
        if (FFT_input[i].real < 0.0f)
            FFT_input[i].real = 0.0f;
    }
}

// --- Part A ---
/*
Assuming each processor has its own local cache (of layers L1 or more),
it will need to update values of FFT_input[i] anytime there is a cache miss due to the updating in Processor C.
In the case of Processor B, this could prove much worse than in the other two due to reading from FFT_input[i] two separate times during each iteration of the loop.
*/

// --- Part B ---
/*
First possible solution: have processor C complete its work before the work of processor A and B complete.
In this scenario, A and B should not run into any cache misses during runtime.
Second possible solution: have the total work divided equally amongst the three processors (iterate through only N/3 elements of FFT_input, copied to local cache).
Each processor will also do the sum, max, and .real updating, but only to its own third of the array and only to its own local instance of the array.
Then, when complete, each processor updates shared memory.
*/