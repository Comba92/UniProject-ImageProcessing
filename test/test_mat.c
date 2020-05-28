#include <stdio.h>
#include <string.h>
#include "ip_lib.h"
#include "bmp.h"

int main(int argc, char * argv[]) {
    char * in1 = "flower2.bmp";
    char * in2 = "caf.bmp";
    char * out = "out.bmp";

    int fun, var;

    if(argc == 1) {
        printf("arg 1 - function to test\n \
            arg 2 - variable\n");
        return 0;
    }

    if(argc > 1) {
        fun = atoi(argv[1]);
        if(argc > 2)
            var = atoi(argv[2]);
    }

    Bitmap * bmpimg = bm_load(in1);
    Bitmap * bmpimg2 = bm_load(in2);
    ip_mat * matimg = bitmap_to_ip_mat(bmpimg);
    ip_mat * matimg2 = bitmap_to_ip_mat(bmpimg2);
    bm_free(bmpimg);
    bm_free(bmpimg2);

    ip_mat * result;

    switch(fun) {
        case 1:
            printf("Testing corruption...\n");
            result = ip_mat_corrupt(matimg, var);
            break;
        case 2:
            printf("Testing brighten...\n");
            result = ip_mat_brighten(matimg, var);
            clamp(result, 0, 255);
            break;
        case 3:
            printf("Testing gray...\n");
            result = ip_mat_to_gray_scale(matimg);
            break;
        case 4:
            printf("Testing blend...\n");
            result = ip_mat_blend(matimg, matimg2, var);
            break;
        case 5:
            printf("Testing concat...\n");
            result = ip_mat_concat(matimg, matimg2, var);
            break;
        case 6: 
            printf("Testing subset...\n");
            result = ip_mat_subset(matimg, 2, 4, 2, 4);
            break;
        case 7:
            printf("Testing padding...\n");
            result = ip_mat_padding(matimg, 50, 50);
            break;
        case 8:
            printf("Testing sharpen...\n");
            result = ip_mat_convolve(matimg, create_sharpen_filter());
            clamp(result, 0, 255);
            break;
        case 9:
            printf("Testing edge...\n");
            result = ip_mat_convolve(matimg, create_edge_filter());
            clamp(result, 0, 255);
            break;
        case 10:
            printf("Testing emboss...\n");
            result = ip_mat_convolve(matimg, create_emboss_filter());
            clamp(result, 0, 255);
            break;
        case 11:
            printf("Testing average...\n");
            result = ip_mat_convolve(matimg, create_average_filter(3, 3, 2));
            break;
        case 12:
            printf("Testing gaussian...\n");
            result = ip_mat_convolve(matimg, create_gaussian_filter(6, 6, 3, 1));
            break;
        case 13: 
            printf("Testing free...\n");
            ip_mat_free(matimg);
            return 0;
    }

    bmpimg = ip_mat_to_bitmap(result);
    ip_mat_free(result);
    bm_save(bmpimg, out);
    bm_free(bmpimg);

    return(0);
}
