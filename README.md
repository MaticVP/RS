# Discrete Wavelet Transform


## Introduction 

The Discrete Wavelet Transform (DWT) is a well-known transformation used in signal processing and image compression that transforms a discrete signal into a set of wavelet coefficients. Wavelets, as basis functions, offer a powerful alternative to traditional methods like Fourier series in signal and image analysis. A wavelet is a waveform of effectively limited duration that has an average value of zero and nonzero norm. Unlike Fourier basis functions, which are sinusoids of varying frequencies, wavelets are localized in both time and frequency domains, enabling them to capture localized features and transient behavior more effectively. This localization property makes wavelets particularly well-suited for analyzing signals with discontinuities or sharp changes, as they can provide both time and frequency information simultaneously. Moreover, wavelets offer a multi-resolution analysis, allowing for the representation of signals at different levels of detail. This property makes DWT a powerful tool in various applications such as image compression (e.g., JPEG2000), denoising, and feature extraction in pattern recognition. Its ability to represent data at multiple resolutions simultaneously makes it a versatile and widely used method in modern signal processing.

### Haar Wavelet Transform

The Haar Wavelet Transform represents one of the simplest transformations in wavelet theory due to its simplicity and computational efficiency. It operates by recursively decomposing a signal into pairs of averages and differences, which correspond to the approximation and detail coefficients, respectively. This binary splitting approach makes the Haar transform particularly intuitive, as it uses a step function that switches between 1 and -1. The transform's effectiveness lies in its ability to localize both time and frequency information, making it adept at analyzing signals with abrupt changes or discontinuities. Due to its straightforward implementation and minimal computational requirements, the Haar Wavelet Transform is widely used in various fields, including signal processing, image compression, and data analysis.

Let's take a 1D signal represented by the array [4, 6, 10, 12, 14, 16, 18, 20] and apply the Haar Wavelet Transform to it.

1. **First Level Decomposition:**
   - Pair the numbers and compute the averages and differences:
     - (4, 6): Average = (4 + 6) / 2 = 5, Difference = (4 - 6) / 2 = -1
     - (10, 12): Average = (10 + 12) / 2 = 11, Difference = (10 - 12) / 2 = -1
     - (14, 16): Average = (14 + 16) / 2 = 15, Difference = (14 - 16) / 2 = -1
     - (18, 20): Average = (18 + 20) / 2 = 19, Difference = (18 - 20) / 2 = -1
   - This results in:
     - Averages: [5, 11, 15, 19]
     - Differences: [-1, -1, -1, -1]

2. **Second Level Decomposition:**
   - Now, decompose the averages array [5, 11, 15, 19]:
     - (5, 11): Average = (5 + 11) / 2 = 8, Difference = (5 - 11) / 2 = -3
     - (15, 19): Average = (15 + 19) / 2 = 17, Difference = (15 - 19) / 2 = -2
   - This results in:
     - Averages: [8, 17]
     - Differences: [-3, -2]

3. **Third Level Decomposition:**
   - Finally, decompose the averages array [8, 17]:
     - (8, 17): Average = (8 + 17) / 2 = 12.5, Difference = (8 - 17) / 2 = -4.5
   - This results in:
     - Averages: [12.5]
     - Differences: [-4.5]

4. **Combining Results:**
   - The final output consists of the last average coefficient and all the different coefficients from each level:
     - Final Approximation Coefficient: [12.5]
     - Detail Coefficients: [-4.5, -3, -2, -1, -1, -1, -1]

So, after applying the Haar Wavelet Transform to the array [4, 6, 10, 12, 14, 16, 18, 20], we get the transformed signal:
[ 12.5, -4.5, -3, -2, -1, -1, -1, -1 ]

This transformation effectively condenses the original signal into a set of coefficients that highlight the overall trend (approximation) and the detailed fluctuations (differences) at various levels of resolution.

### 2D Haar Wavelet 

The recursive algorithm for the 2D Haar Wavelet Transform using 2x2 blocks involves iteratively applying the Haar transform to progressively smaller submatrices of the image. This approach ensures that the image is decomposed into multiple levels of resolution, capturing both coarse and fine details. We can describe it with the following algorithm: 


```
Function DWTHaar2D(matrix, levels):
    N = size(matrix)  # Get the size of the matrix
        # Apply the 2D Haar Transform to the NxN submatrix
        if(levels > 0):
            for i = 0 to N-1 by 2:
                for j = 0 to N-1 by 2:
                    # Get the 2x2 block
                    a = matrix[i][j]
                    b = matrix[i][j+1]
                    c = matrix[i+1][j]
                    d = matrix[i+1][j+1]
                    
                    # Compute the averages and differences
                    A = a + b + c + d  # Approximation
                    H = (-a - b + c + d)  # Horizontal detail
                    V = (-a + b - c + d)  # Vertical detail
                    D = (a - b - c + d)  # Diagonal detail
                    
                    # Place the computed values back into the matrix
                    matrix[i/2][j/2] = A  # Top-left quadrant
                    matrix[i/2 + N/2][j/2 ] = H  # Top-right quadrant
                    matrix[i/2][j/2 + N/2] = V  # Bottom-left quadrant
                    matrix[(i/2) + N/2][(j/2) + N/2] = D  # Bottom-right quadrant

            DWTHaar2D(A, level-1)


        # Reduce the size for the next level of decomposition
    return matrix
```

Where `matrix` represents the square matrix dimension `N x N,` while the `levels` represent the numbers of decomposition/resolution levels that we want to calculate. In the folder images, you can image `elaine.jpg` and its Wavelet transformation `dwt_output.jpg` on the first level.  

The inverse Haar wavelet transform can be described with the following algorithm: 


```
Function InverseHaarWavelet(matrix, levels):
    N = size(matrix)  # Get the size of the matrix

    n2 = n // 2  # Half of the size 
    g = copy(matrix)  # Copy the input matrix to preserve original data
    

    if level >= 1:
        # Recursive call to process the top-left quadrant
        A = InverseHaarWavelet(matrix[:n2, :n2], levels-1)
        
        # Invert the filter using the Haar wavelet coefficients
        H = matrix[n2:2*n2, :n2]
        V = matrix[:n2, n2:2*n2]
        D = matrix[n2:2*n2, n2:2*n2]

        # Reconstruct the original values from the coefficients
        for i = 0 to n2-1:
            for j = 0 to n2-1:
                g[2*i, 2*j] = (A[i, j] - H[i, j] - V[i, j] + D[i, j]) / 4
                g[2*i, 2*j + 1] = (A[i, j] - H[i, j] + V[i, j] - D[i, j]) / 4
                g[2*i + 1, 2*j] = (A[i, j] + H[i, j] - V[i, j] - D[i, j]) / 4
                g[2*i + 1, 2*j + 1] = (A[i, j] + H[i, j] + V[i, j] + D[i, j]) / 4

    return g

```

## Task - HLS synthesis of wavelet transformation 

- (3 pts) In the first step, you will write functions for the forward and inverse wavelet transform and test their functionality using a sample image from the folder images. A code snippet can be found in the folder `src/sw.`

- (12 pts) In the second step, you will use the Vitis HLS tool to develop and optimize kernels for the forward and inverse DWT:

  - (3/12 pts) Analyze the delay and resource utilization of the base version for DWT and IDWT with one level of decomposition for matrix sizes 256 x 256, 512 x 512, and 1024 x 1024. 
  - (3/12 pts) Analyze the delay and resource utilization of the base version for DWT and IDWT for a matrix size of 512 x 512 with the following number of decomposition levels: 4, 5, 6, and 7.
  - (4/12 pts) Using pragmas, optimize the kernel. Report your findings for a matrix size of 512 x 512 with two levels of decomposition.
  - (2/12 pts) Starting from the previously optimized kernel, try to employ the dataflow pragma. Report the resource utilization and delay.

The code for the kernels can be found in the folders `hls/src/dwt.c` and `hls/src/idwt.c`. To run HLS synthesis, you will use the sbatch script `sbatch_script.sh`, which runs the script `run_hls.sh`.

- (5 pts) The final five points will be awarded based on your project presentation.

## Connecting to SLING NSC cluster

- Open the terminal and connect to the cluster using SSH

```console
$ ssh username-ID@nsc-login2.ijs.si
```

## Literature 

Jeff Calder, [The Discrete Wavelet Transform](https://www-users.cse.umn.edu/~jwcalder/5467/lec_dwt.pdf), Mathematics of Image and Data Analysis Math 546 