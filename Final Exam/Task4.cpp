#include <vector>
#include <iostream>
#include <cstdlib>
#include <cmath>

typedef struct
{
    float real;
    float imag;
} complex_t;

// Size of the recording file
#define N sizeOfInput
#define WINDOW_SIZE 1024

std::vector<float> anechoicRecording = new std::vector<float>(N);
std::vector<float> impulseResponse = new std::vector<float>(N);

std::vector<float> computeArtificialReverb()
{
    // Signal stored in std::vector<float> anechoicRecording;
    // Impulse stored in std::vector<float> impulseRecording;

    // This is where we store windows
    std::vector<std::vector<complex_t>> fftAnechoicWindows = new std::vector<std::vector<complex_t>>();
    std::vector<std::vector<complex_t>> fftImpulseWindows = new std::vector<std::vector<complex_t>>();
    std::vector<std::vector<complex_t>> fftSums = new std::vector<std::vector<complex_t>>();

    // Split into windows
    for (int i = 0; i < N; i += WINDOW_SIZE)
    {
        // Save this window of the recording
        std::vector<float> anechoicRecordingWindow = new std::vector<float>(2 * WINDOW_SIZE);
        std::vector<float> impulseResponseWindow = new std::vector<float>(2 * WINDOW_SIZE);
        for (int j = i; j < i + WINDOW_SIZE; j++)
        {
            if (j < anechoicRecording.size())
            {
                anechoicRecordingWindow(j - i) = anechoicRecording(j);
            }
            else
            {
                anechoicRecordingWindow(j - i) = 0;
            }

            if (j < impulseResponseWindow)
            {
                impulseResponseWindow(j - i) = impulseResponse(j);
            }
            else
            {
                impulseResponseWindow(j - i) = 0;
            }
        }

        // Now pad with 1024 zeroes
        for (int j = 1024; j < (2 * WINDOW_SIZE); j++)
        {
            anechoicRecordingWindow(j) = 0;
            impulseResponseWindow(j) = 0;
        }

        // Now that we have this, FFT both
        std::vector<complex_t> fftAnechoic = fft2048(anechoicRecordingWindow);
        fftAnechoicWindows.insert(fftAnechoic);
        std::vector<complex_t> fftImpulse = fft2048(impulseResponseWindow);
        fftImpulseWindows.insert(fftImpulse);

        fftSums.insert(new std::vector<complex_t>(WINDOW_SIZE * 2));
        for (int j = 0; j < fftAnechoicWindows.size(); j++)
        {
            for (int k = fftImpulseWindows.size(); k >= 0; k--)
            {
                std::vector<complex_t> summedWindow = new std::vector<complex_t>();
                for (int l = 0; l < (2 * WINDOW_SIZE); l++)
                {
                    // Compute the sum of two complex numbers here
                    complex_t value = {
                        fftImpulseWindows.at(i).at(l).real + fftAnechoicWindows.at(i).at(l),
                        fftImpulseWindows.at(i).at(l).imag + fftAnechoicWindows.at(i).at(l)};
                    summedWindow.insert(value);
                }

                // Now insert into the sum
                fftSums.insert(summedWindow);
            }
        }
    }

    // Now that we have the windows, compute the IFFT of each to get the time domain
    std::vector<std::vector<float>> ifftSignals = new std::vector<std::vector<float>>();
    for (int i = 0; i < fftSums.size(); i++)
    {
        ifftSignals.insert(ifft2048(fftSums.at(i)));
    }

    // Then, simply concatenate these time domain signals together
    std::vector<float> timeDomainSignal = new std::vector<float>();
    for (int i = 0; i < ifftSignals.size(); i++)
    {
        for (int j = 0; j < ifftSignals.at(0).size(); j++)
        {
            timeDomainSignal.insert(ifftSignals.at(i).at(j));
        }
    }

    return timeDomainSignal;
}