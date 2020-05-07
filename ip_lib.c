/*
 Created by Sebastiano Vascon on 23/03/20.
*/

#include <stdio.h>
#include "ip_lib.h"
#include "bmp.h"

void ip_mat_show(ip_mat * t){
    unsigned int i,l,j;
    printf("Matrix of size %d x %d x %d (hxwxk)\n",t->h,t->w,t->k);
    for (l = 0; l < t->k; l++) {
        printf("Slice %d\n", l);
        for(i=0;i<t->h;i++) {
            for (j = 0; j < t->w; j++) {
                printf("%f ", get_val(t,i,j,l));
            }
            printf("\n");
        }
        printf("\n");
    }
}

void ip_mat_show_stats(ip_mat * t){
    unsigned int k;

    compute_stats(t);

    for(k=0;k<t->k;k++){
        printf("Channel %d:\n", k);
        printf("\t Min: %f\n", t->stat[k].min);
        printf("\t Max: %f\n", t->stat[k].max);
        printf("\t Mean: %f\n", t->stat[k].mean);
    }
}

ip_mat * bitmap_to_ip_mat(Bitmap * img){
    unsigned int i=0,j=0;

    unsigned char R,G,B;

    unsigned int h = img->h;
    unsigned int w = img->w;

    ip_mat * out = ip_mat_create(h, w,3,0);

    for (i = 0; i < h; i++)              /* rows */
    {
        for (j = 0; j < w; j++)          /* columns */
        {
            bm_get_pixel(img, j,i,&R, &G, &B);
            set_val(out,i,j,0,(float) R);
            set_val(out,i,j,1,(float) G);
            set_val(out,i,j,2,(float) B);
        }
    }

    return out;
}

Bitmap * ip_mat_to_bitmap(ip_mat * t){

    Bitmap *b = bm_create(t->w,t->h);

    unsigned int i, j;
    for (i = 0; i < t->h; i++)              /* rows */
    {
        for (j = 0; j < t->w; j++)          /* columns */
        {
            bm_set_pixel(b, j,i, (unsigned char) get_val(t,i,j,0),
                    (unsigned char) get_val(t,i,j,1),
                    (unsigned char) get_val(t,i,j,2));
        }
    }
    return b;
}

float get_val(ip_mat * a, unsigned int i,unsigned int j,unsigned int k){
    if(i<a->h && j<a->w &&k<a->k){  /* j>=0 and k>=0 and i>=0 is non sense*/
        return a->data[i][j][k];
    }else{
        printf("Errore get_val!!!");
        exit(1);
    }
}

void set_val(ip_mat * a, unsigned int i,unsigned int j,unsigned int k, float v){
    if(i<a->h && j<a->w &&k<a->k){
        a->data[i][j][k]=v;
    }else{
        printf("Errore set_val!!!");
        exit(1);
    }
}

float get_normal_random(){
    float y1 = ( (float)(rand()) + 1. )/( (float)(RAND_MAX) + 1. );
    float y2 = ( (float)(rand()) + 1. )/( (float)(RAND_MAX) + 1. );
    return cos(2*PI*y2)*sqrt(-2.*log(y1));

}

ip_mat * ip_mat_create(unsigned int h, unsigned int w, unsigned int k, float v) {
    ip_mat * out = (ip_mat *) malloc(sizeof(ip_mat));
    out->h = h;
    out->w = w;
    out->k = k;

    out->stat = (stats *) malloc(sizeof(stats));
    out->stat->min = v;
    out->stat->max = v;
    out->stat->mean = v;

    int i,j,l;
    out->data = (float ***) malloc(sizeof(float **) * h);
    for(i=0; i<h; i++) {
        out->data[i] = (float **) malloc(sizeof(float *) * w);
        for(j=0; j<w; j++){
            out->data[i][j] = (float *) malloc(sizeof(float) * k);
            for(l=0; l<k; l++)
                out->data[i][j][l] = v;
        }
    }

    return out;
}

void ip_mat_free(ip_mat *a) {
    free(a->stat);
    
    /* AD OGNI MALLOC UN CORRISPETTIVO FREE!! */
    int i,j;
    for(i=0; i<a->h; i++) {
        for(j=0; j<a->w; j++)
            free(a->data[i][j]);
        free(a->data[i]);
    }
    free(a->data);

    free(a);
}

void compute_stats(ip_mat * t) {
    /* Max e min definite come macro, sono più efficienti */
    #define max(a,b) (a > b ? a : b);
    #define min(a,b) (a < b ? a : b);
    float minVal, maxVal, sum=0;
    int i,j,k;

    /*Setto al primo valore di data */
    minVal = get_val(t, 0,0,0);
    maxVal = get_val(t,0,0,0);

    for(i=0; i < t->h; i++)
        for(j=0; j < t->w; j++)
            for(k=0; k < t->k; k++) {
                minVal = min(minVal, get_val(t, i,j,k));
                maxVal = max(maxVal, get_val(t, i,j,k));
                sum += get_val(t, i,j,k);
            }
    t->stat->min = minVal;
    t->stat->max = maxVal;
    t->stat->mean = sum / (float)(t->w * t->h * t->k);
}

void ip_mat_init_random(ip_mat * t, float mean, float var) {
    ip_mat * out = (ip_mat *) malloc(sizeof(ip_mat));
    int i, j, k;
    for(i=0; i<out->h; i++)
        for(j=0; j<out->w; j++)
            for(k=0; k < out->k; k++)
                set_val(out, i,j,k, 
                    get_normal_random() * var + mean);
}

ip_mat * ip_mat_copy(ip_mat * in) {
    ip_mat * out = (ip_mat *) malloc(sizeof(ip_mat));
    out->h = in->h;
    out->w = in->w;
    out->k = in->k;
    out->stat = in->stat;
    out->data = in->data;

    return out;
}

ip_mat * ip_mat_subset(ip_mat * t, unsigned int row_start, unsigned int row_end, unsigned int col_start, unsigned int col_end) {
    ip_mat * out = ip_mat_create(row_end - row_start,
        col_end - col_start, t->k, 0);
    
    int i, j, k;
    for(i=0; row_start != row_end; i++, row_start++)
        for(j=0; col_start != col_end; j++, col_start++)
            for(k=0; k < out->k; k++)
                set_val(out, i,j,k, 
                    get_val(t, row_start, col_start, k));

    return out;
}

ip_mat * ip_mat_concat(ip_mat * a, ip_mat *b, int dimensione) {
    ip_mat * out;
    switch(dimensione) {
        case 0:
            out = ip_mat_create(a->h + b->h, a->w, a->k, 0);
            break;
        case 1:
            out = ip_mat_create(a->h, a->w + b->w, a->k, 0);
            break;
        case 2:
            out = ip_mat_create(a->h, a->w, a->k + b->k, 0);
            break;
    }
    return out;
}

ip_mat * ip_mat_sum(ip_mat *a, ip_mat *b) {
    ip_mat * out = ip_mat_create(a->h, a->w, a->k, 0);

    int i, j, k;
    for(i=0; i<out->h; i++)
        for(j=0; j<out->w; j++)
            for(k=0; k < out->k; k++)
                set_val(out, i,j,k, 
                    get_val(a,i,j,k) + get_val(b,i,j,k));
    return out;
}

ip_mat * ip_mat_sub(ip_mat *a, ip_mat *b) {
    ip_mat * out = ip_mat_create(a->h, a->w, a->k, 0);

    int i, j, k;
    for(i=0; i<out->h; i++)
        for(j=0; j<out->w; j++)
            for(k=0; k < out->k; k++)
                set_val(out, i,j,k, 
                    get_val(a,i,j,k) - get_val(b,i,j,k));
    return out;
}

ip_mat * ip_mat_mul_scalar(ip_mat *a, float c) {
    ip_mat * out = ip_mat_create(a->h, a->w, a->k, 0);

    int i, j, k;
    for(i=0; i<out->h; i++)
        for(j=0; j<out->w; j++)
            for(k=0; k < out->k; k++)
                set_val(out, i,j,k, 
                    get_val(a,i,j,k) * c);
    return out;
}

ip_mat *  ip_mat_add_scalar(ip_mat *a, float c) {
    ip_mat * out = ip_mat_create(a->h, a->w, a->k, 0);

    int i, j, k;
    for(i=0; i<out->h; i++)
        for(j=0; j<out->w; j++)
            for(k=0; k < out->k; k++)
                set_val(out, i,j,k, 
                    get_val(a,i,j,k) + c);
    return out;
}

ip_mat * ip_mat_mean(ip_mat * a, ip_mat * b) {
    ip_mat * out = ip_mat_create(a->h, a->w, a->k, 0);

    int i, j, k;
    for(i=0; i<out->h; i++)
        for(j=0; j<out->w; j++)
            for(k=0; k < out->k; k++)
                set_val(out, i,j,k, 
                    (get_val(a,i,j,k) + get_val(b,i,j,k)) / 2);
    return out;
}

ip_mat * ip_mat_to_gray_scale(ip_mat * in) {
    ip_mat * out = ip_mat_create(in->h, in->w, in->k, 0);

    float mean;
    int i, j, k;
    for(i=0; i<out->h; i++)
        for(j=0; j<out->w; j++) {
            mean = (get_val(in,i,j,0) + get_val(in,i,j,1) +
                get_val(in,i,j,2)) / 3;
            for(k=0; k < out->k; k++)
                set_val(out, i,j,k,mean);
        }
    return out;
}

ip_mat * ip_mat_blend(ip_mat * a, ip_mat * b, float alpha) {
    ip_mat * out = ip_mat_create(a->h, a->w, a->k, 0);

    float blend;
    int i, j, k;
    for(i=0; i<out->h; i++)
        for(j=0; j<out->w; j++)
            for(k=0; k < out->k; k++) {
                blend = alpha*get_val(a,i,j,k) + 
                    (1-alpha)*get_val(b,i,j,k);
                set_val(out,i,j,k, blend);
            }
                
    return out;
}

ip_mat * ip_mat_brighten(ip_mat * a, float bright) {
    return ip_mat_add_scalar(a, bright);
}

ip_mat * ip_mat_corrupt(ip_mat * a, float amount) {
    return ip_mat_add_scalar(a, get_normal_random()*amount);
}


