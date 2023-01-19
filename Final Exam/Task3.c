// ----- This is a file for Task 3 -----
#include "lib/image.h" // Library from Data Structures and Algorithms class
#include <string>
#include <iostream>

// Enable an extreme amount of printing. Useful for debugging
#define VERBOSE false

// ===== Image functions =======================================================

/**
 * Function to load an image from a file and return as an Image object
 * @param filename - the name of the file to load
 * @return - the loaded image object as an Image<Pixel>
 */
Image<Pixel> loadImage(const std::string &filename)
{
    try
    {
        // Try to load image given filename
        return readFromFile(filename);
    }
    catch (std::exception &e)
    {
        // The file could not be loaded
        std::cerr << "Error loading image: " << e.what() << std::endl;
        Image<Pixel> emptyImage;
        return emptyImage;
    }
}

/**
 * Writes an image to a file given the image and filename
 * @param image - the image to write
 * @param filename - the filename to write to
 */
void writeImageToFile(const Image<Pixel> &image, const std::string &filename)
{
    try
    {
        writeToFile(image, filename);
    }
    catch (std::exception &e)
    {
        // The file could not be written
        std::cerr << "Error writing image: " << e.what() << std::endl;
    }
}

// ===== Image processing functions ============================================

// Take in a pixel, and convert it to black and white
Pixel toBlackAndWhite(Pixel pixel)
{
    // Get the average of the RGB values
    uint8_t average = (pixel.red + pixel.green + pixel.blue) / 3;
    // Set the RGB values to the average
    Pixel newPixel = {average, average, average, 255};
    return newPixel;
}

// Take in an image, and convert all pixels to black and white
Image<Pixel> imageToBlackAndWhite(Image<Pixel> image)
{
    // Create a new image with the same dimensions as the original
    Image<Pixel> newImage(image.width(), image.height());
    // Loop through all pixels in the image
    for (int i = 0; i < image.width(); i++)
    {
        for (int j = 0; j < image.height(); j++)
        {
            // Get the pixel at the current position
            // For some reason, the image library uses (y, x) instead???
            Pixel pixel = image(j, i);

#if VERBOSE
            std::cout << "Pixel (" << i << ", " << j << ")" << std::endl;
#endif

            // Convert the pixel to black and white
            Pixel newPixel = toBlackAndWhite(pixel);
            // Set the pixel in the new image
            newImage(j, i) = newPixel;
        }
    }
    return newImage;
}

// Take in a pixel, and find the closest palette color to it
Pixel findClosestPaletteColor(Pixel pixel)
{
    // Start by extracting the RGB values as an average
    uint8_t average = (pixel.red + pixel.green + pixel.blue) / 3;

    // Compute a rounding of this average
    uint8_t newColor = 0;
    if (average < 128)
    {
        newColor = 255;
    }

    // Return the pixel
    Pixel newPixel = {newColor, newColor, newColor, 255};
    return newPixel;
}

// Find the quantitative error of a new pixel compared to an old pixel
uint8_t findQuantitativeError(Pixel oldPixel, Pixel newPixel)
{
    // Get the average of the RGB values
    uint8_t oldAverage = (oldPixel.red + oldPixel.green + oldPixel.blue) / 3;
    uint8_t newAverage = (newPixel.red + newPixel.green + newPixel.blue) / 3;

    // Compute the error
    uint8_t error = oldAverage - newAverage;
    return error;
}

// Update a pixel based on a passed in error
Pixel updatePixel(Pixel pixel, uint8_t error)
{
    // Get the average value of the pixel
    uint8_t average = (pixel.red + pixel.green + pixel.blue) / 3;

    // Compute the new average
    uint8_t newAverage = average + error;

    // Return the new pixel
    Pixel newPixel = {newAverage, newAverage, newAverage, 255};
    return newPixel;
}

// Take in an image, and apply Floyd-Steinberg dithering to it
Image<Pixel> applyDitheringToImage(Image<Pixel> image)
{
    Image<Pixel> newImage(image.width(), image.height());
// Copy the old image to the new image
#pragma omp for
    for (int y = 0; y < image.height(); y++)
    {
        for (int x = 0; x < image.width(); x++)
        {
            newImage(y, x) = image(y, x);
        }
    }
#pragma omp barrier

// Now loop through and apply dithering
#pragma omp for
    for (int y = 0; y < newImage.height(); y++)
    {
        for (int x = 0; x < newImage.width(); x++)
        {
            // Save the old pixel
            Pixel oldPixel = newImage(y, x);
            // Create a new pixel with the nearest palette color
            Pixel newPixel = findClosestPaletteColor(oldPixel);
            // Compute the error
            uint8_t quantError = findQuantitativeError(oldPixel, newPixel);

            // Set the new pixel into the image
            newImage(y, x) = newPixel;

            // Propagate the error to the surrounding pixels (if they exist)
            if (x + 1 < newImage.width())
                newImage(y, x + 1) =
                    updatePixel(newImage(y, x + 1), quantError * 7 / 16);
            if (x - 1 >= 0 && y + 1 < newImage.height())
                newImage(y + 1, x - 1) =
                    updatePixel(newImage(y + 1, x - 1), quantError * 3 / 16);
            if (y + 1 < newImage.height())
                newImage(y + 1, x) =
                    updatePixel(newImage(y + 1, x), quantError * 5 / 16);
            if (x + 1 < newImage.width() && y + 1 < newImage.height())
                newImage(y + 1, x + 1) =
                    updatePixel(newImage(y + 1, x + 1), quantError * 1 / 16);
        }
    }
#pragma omp barrier

    // Now we've (hopefully) dithered an image
    return newImage;
}

// ===== Driver code ===========================================================

int main()
{
    std::cout << "Convert an image to black and white!" << std::endl;

    // Image name - Must be a PNG file
    const std::string *filename = new std::string("lena.png");

    // Create the input filename
    const std::string *inFilepath = new std::string("in/" + *filename);
    // Load the image
    Image<Pixel> image = loadImage(*inFilepath);
    // Check if the image was loaded
    if (image.width() == 0 || image.height() == 0)
    {
        std::cerr << "Could not load image!" << std::endl;
        return EXIT_FAILURE;
    }
    else
    {
        std::cout << "Image successfully loaded. (" << image.width() << "x"
                  << image.height() << ")" << std::endl;
    }

    // Convert the image to black and white
    Image<Pixel> imageBW = imageToBlackAndWhite(image);
    std::cout << "Image converted to black and white." << std::endl;

    // Dither the image
    Image<Pixel> imageDither = applyDitheringToImage(imageBW);
    std::cout << "Image dithered." << std::endl;

    // Create the output filename
    const std::string *outFilepath = new std::string("out/" + *filename);
    // Write the image to a file
    writeImageToFile(imageDither, *outFilepath);
    std::cout << "Image successfully written to file." << std::endl;

    return EXIT_SUCCESS;
}