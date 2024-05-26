//
//  Created by Ratko Pilipovic on 22/11/2023.
//  Copyright © 023 Ratko Pilipovic  UL FRI. All rights reserved.
//

#include "matv_mult_opt.h"
#define MAT_SIZE 64
#include <iostream>

using namespace std;

int main() {

    int a[MAT_SIZE][MAT_SIZE];
    int b[MAT_SIZE];
    

    // initialize matrix and vector
    for (int i = 0; i < MAT_SIZE; i++) {
        for (int j = 0; j < MAT_SIZE; j++) {
            a[i][j] = i+j;
        }
        b[i] = i;
    }

    // perform matrix vector multiplication on CPU
    int c_sw[MAT_SIZE];
    for (int i = 0; i < MAT_SIZE; i++) {
        c_sw[i] = 0;
        for (int j = 0; j < MAT_SIZE; j++) {
            c_sw[i] += a[i][j] * b[j];
        }
    }

    // perform matrix vector multiplication on FPGA
    int c_hw[MAT_SIZE];
    matv_mult_opt(c_hw, (int*)a, (int*)b);

    // compare results
    for (int i = 0; i < MAT_SIZE; i++) {
        if (c_sw[i] != c_hw[i]) {
            cout << "Error: c_sw[" << i << "] = " << c_sw[i] << ", c_hw[" << i << "] = " << c_hw[i] << endl;
            return 1;
        }
    }
}