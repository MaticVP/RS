#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>


#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image.h"
#include "stb_image_write.h"

#define IMAGE_WIDTH 512
#define LEVELS 4


// Define image structure
// struct Image {
//     uint8_t pixels[IMAGE_WIDTH][IMAGE_HEIGHT];
// };

void copy_submatrix1D(const double *source,double *dest, int start_row, int start_col, int end_row, int end_col,int width) {
    int dest_rows = end_row - start_row;
    int dest_cols = end_col - start_col;
    int i, j;
    for (i = 0; i < dest_rows; i++) {
        for (j = 0; j < dest_cols; j++) {
            dest[j+i*width] = source[(start_col + j)+width*(start_row + i)];
        }
    }
}

// improve for showing
void improve(double *image, double* output) {
    int i, j;
    int height = IMAGE_WIDTH;
    int width = IMAGE_WIDTH;
    // Approximation coefficient
    for(int level=LEVELS;level>0;level--){
        for (i = 0; i < height/(1<<level); i++) {
            for (j = 0; j < width/(1<<level); j++) {
                // calculate average of 4 pixels in the neighborhood
                output[i * width + j] = image[i * width + j]/4;
                output[(i + height/2) * width + j] = image[(i + height/2) * width + j]+0.5;
                output[i * width + j + height/2] = image[i * width + j + height/2]+0.5;
                output[(i + height/2) * width + j + height/2] = image[(i + height/2) * width + j + height/2]+0.5;
            }
        }
    }
}



// calculate DWT
void forward_wavelet_transform(const double *image, double* output, int levels) {
    int i, j, l;
    double temp;
    int height = IMAGE_WIDTH;
    int width = IMAGE_WIDTH;

    double* g = (double *)malloc(sizeof(double ) * width * height);

    for (l = levels; l > 0; l--) {
        int n2 = IMAGE_WIDTH / (1 << (levels-l));
        for (i = 0; i < n2; i+=2) {
            for (j = 0; j < n2; j+=2) {
                double a, b, c, d;

                if(levels==l) {
                    a = image[j + i * width];
                    b = image[(j + 1) + i * width];
                    c = image[j + (i + 1) * width];
                    d = image[(j + 1) + (i + 1) * width];
                } else{
                    a = output[j + i * width];
                    b = output[(j + 1) + i * width];
                    c = output[j + (i + 1) * width];
                    d = output[(j + 1) + (i + 1) * width];
                }

                double A = (a + b + c + d); // Approximation coefficient
                double H = (-a - b + c + d); // Horizontal coefficients
                double V = (-a + b - c + d); // Vertical coefficients
                double D = (a - b - c + d); // Diagonal coefficients

                int halfN2 = n2 / 2;
                int indexA = (j / 2) + (i / 2) * width;
                int indexH = (j / 2) + (i / 2 + halfN2) * width;
                int indexV = (j / 2 + halfN2) + (i / 2) * width;
                int indexD = (j / 2 + halfN2) + (i / 2 + halfN2) * width;

                g[indexA] = A;
                g[indexH] = H;
                g[indexV] = V;
                g[indexD] = D;
            }
        }

        for (i = 0; i < n2; i++) {
            for (j = 0; j < n2; j++) {
                output[i*width + j] = g[i*width + j];
            }
        }
    }
}


// calculate Inverse DWT (IDWT) 
void inverse_wavelet_transform(double *image, double* output, int levels) {
    int i, j, l;
    double temp;
    int height = IMAGE_WIDTH;
    int width = IMAGE_WIDTH;

    int max_n2 = height / 2;
    double* A = (double*)malloc(IMAGE_WIDTH * IMAGE_WIDTH * sizeof(double*));
    double* H = (double*)malloc(IMAGE_WIDTH * IMAGE_WIDTH * sizeof(double*));
    double* V = (double*)malloc(IMAGE_WIDTH * IMAGE_WIDTH * sizeof(double*));
    double* D = (double*)malloc(IMAGE_WIDTH * IMAGE_WIDTH * sizeof(double*));

    double *g = (double *)malloc(IMAGE_WIDTH * IMAGE_WIDTH * sizeof(double *));

    for (i = 0; i < IMAGE_WIDTH * IMAGE_WIDTH; i++) {
        g[i] = 0;
    }

    for (l=levels;l>0;l--){
        int n2 = IMAGE_WIDTH / (1 << (l));

        if(l==levels) {
            copy_submatrix1D(output, A, 0, 0, n2, n2, width);
            copy_submatrix1D(output, H, n2, 0, 2 * n2, n2, width);
            copy_submatrix1D(output, V, 0, n2, n2, 2 * n2, width);
            copy_submatrix1D(output, D, n2, n2, 2 * n2, 2 * n2, width);
        } else{
            copy_submatrix1D(g, A, 0, 0, n2, n2,width);
            copy_submatrix1D(g, H, n2, 0, 2 * n2, n2,width);
            copy_submatrix1D(g, V, 0, n2, n2, 2 * n2,width);
            copy_submatrix1D(g, D, n2, n2, 2 * n2, 2 * n2,width);
        }

        for (i = 0; i < n2 ; i++) {
            for (j = 0; j < n2; j++) {
                g[(2*i)*width+(2*j)] = (A[i*width+j] - H[i*width+j] - V[i*width+j] + D[i*width+j]) / 4;
                g[(2*i)*width+(2*j+1)] = (A[i*width+j] - H[i*width+j] + V[i*width+j] + D[i*width+j]) / 4;
                g[(2*i+1)*width+(2*j)] = (A[i*width+j] + H[i*width+j] - V[i*width+j] - D[i*width+j]) / 4;
                g[(2*i+1)*width+(2*j+1)] = (A[i*width+j] + H[i*width+j] + V[i*width+j] + D[i*width+j]) /4;
            }
        }
    }

    for (i = 0; i < IMAGE_WIDTH; i++) {
        for (j = 0; j < IMAGE_WIDTH; j++) {
            output[i * IMAGE_WIDTH + j] = g[i*width+j];
        }
    }
}

int main() {
    // Define sample grayscale image

    // Load image from file and allocate space for the output image
    char image_name[] = "./elaine.jpg";
    int width, height, cpp;
    // load only gray scale image
    unsigned char *h_imageIn;
    h_imageIn = stbi_load(image_name, &width, &height, &cpp, STBI_grey);
    if (h_imageIn == NULL)
    {
        printf("Error reading loading image %s!\n", image_name);
        exit(EXIT_FAILURE);
    }
    printf("Loaded image %s of size %dx%d.\n", image_name, width, height);
    printf("Image is %d bytes per pixel.\n", cpp);
    // Save grayscale image to file
    printf("Size of image is %ld, %ld\n", sizeof(unsigned char), sizeof(h_imageIn));
    


    double *image_pixels = (double*)malloc(sizeof(double) * width * height);
    double *output = (double*)malloc(sizeof(double) * width * height);
    double *show_output = (double*)malloc(sizeof(double) * width * height);

    // convert to grayscale 
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            image_pixels[i*width + j] = h_imageIn[i * width + j]/255.0;
        }
    }


    forward_wavelet_transform(image_pixels, output, LEVELS);

    // used only for showing purposes, not needed for IDWT
    improve(output, show_output);
    
    inverse_wavelet_transform(image_pixels, output, LEVELS);

    // Save image to file

    int image_size = width * height;


    // Save image to file
    for(int i = 0; i < height; i++) {
        for(int j = 0; j < width; j++) {
            h_imageIn[i*width + j] = (char)(output[i*width + j]*255);
        }
    }
    // Free memory

    // free(image.pixels[0]);
    // free(image.pixels);
    //free(clusters->num_points);
    stbi_write_jpg("dwt_output.jpg", width, height,STBI_grey, h_imageIn, 100);
    free(image_pixels);

    return 0;
}
