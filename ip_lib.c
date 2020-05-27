/*
 Created by Sebastiano Vascon on 23/03/20.
*/

#include <stdio.h>
#include "ip_lib.h"
#include "bmp.h"

#define max(a,b) (a > b ? a : b);
#define min(a,b) (a < b ? a : b);
#define mat_iterate(m,i,j,ch) for(ch=0;ch<m->k;ch++) for(i=0;i<m->h;i++) for(j=0;j<m->w;j++)

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

float get_normal_random(float media, float std){

    float y1 = ( (float)(rand()) + 1. )/( (float)(RAND_MAX) + 1. );
    float y2 = ( (float)(rand()) + 1. )/( (float)(RAND_MAX) + 1. );
    float num = cos(2*PI*y2)*sqrt(-2.*log(y1));

    return media + num*std;
}

ip_mat * ip_mat_create(unsigned int h, unsigned int w, unsigned int k, float v) {
    unsigned int i,j,l;
    ip_mat * out = (ip_mat *) malloc(sizeof(ip_mat));
    out->h = h;
    out->w = w;
    out->k = k;

    out->stat = (stats *) malloc(sizeof(stats) * k);

    out->data = (float ***) malloc(sizeof(float **) * h);
    for(i=0; i<h; i++) {
        out->data[i] = (float **) malloc(sizeof(float *) * w);
        for(j=0; j<w; j++){
            out->data[i][j] = (float *) malloc(sizeof(float) * k);
            for(l=0; l<k; l++)
                set_val(out, i,j,l, v);
        }
    }

    return out;
}

void ip_mat_free(ip_mat *a) {
    if(a) {
        unsigned int i,j;
        for(i=0; i<a->h; i++) {
            for(j=0; j<a->w; j++)
                free(a->data[i][j]);
            free(a->data[i]);
        }
        free(a->data);

        free(a->stat);
        free(a);
    }
}

void compute_stats(ip_mat * t) {
    float minVal, maxVal, sum=0;
    unsigned int i,j,k;

    /*Setto al primo valore di data */
    minVal = get_val(t,0,0,0);
    maxVal = get_val(t,0,0,0);

    for(k=0; k < t->k; k++) {
        for(i=0; i< t->h; i++)
            for(j=0; j < t->w; j++) {
                minVal = min(minVal, get_val(t, i,j,k));
                maxVal = max(maxVal, get_val(t, i,j,k));
                sum += get_val(t, i,j,k);
            }
        t->stat[k].min = minVal;
        t->stat[k].max = maxVal;
        t->stat[k].mean = sum / (t->h * t->w);
    }
}

void ip_mat_init_random(ip_mat * t, float mean, float var) {
    unsigned int i, j, k;
    for(i=0; i<t->h; i++)
        for(j=0; j<t->w; j++)
            for(k=0; k < t->k; k++)
                set_val(t, i,j,k, 
                    get_normal_random(mean, var));
}

ip_mat * ip_mat_copy(ip_mat * in) {
    ip_mat * out = ip_mat_create(in->h, in->w, in->k, 0);
    unsigned int i,j,k;
    for(i=0; i<out->h; i++)
        for(j=0; j<out->w; j++)
            for(k=0; k < out->k; k++)
            	set_val(out, i,j,k, 
            		get_val(in, i,j,k));

    return out;
}

ip_mat * ip_mat_subset(ip_mat * t, unsigned int row_start, unsigned int row_end, unsigned int col_start, unsigned int col_end) {
    ip_mat * out = ip_mat_create(row_end - row_start + 1,
        col_end - col_start + 1, t->k, 0);
    
    unsigned int i, j, k;
    for(i=0; i + row_start <= row_end; i++)
        for(j=0; j + col_start <= col_end; j++) 
            for(k=0; k < out->k; k++)
                set_val(out, i,j,k, 
                    get_val(t, i + row_start, j + col_start, k));

    return out;
}


ip_mat * ip_mat_concat(ip_mat * a, ip_mat *b, int dimensione) {
    ip_mat * out = NULL;
    unsigned int i, j, k;
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

    for(k=0; k<out->k; k++) {
        for(i=0; i<a->h; i++)
            for(j=0; j<a->w; j++)
                set_val(out, i,j,k,
                    get_val(a,i,j,k));
        for(i=0; i<b->h; i++)
            for(j=0; j<b->w; j++)
                switch(dimensione) {
                    case 0:
                        set_val(out, a->h+i,j,k,
                            get_val(b,i,j,k));
                        break;
                    case 1:
                        set_val(out, i,a->w+j,k,
                            get_val(b,i,j,k));
                        break;
                    case 2:
                        set_val(out, i,j,a->k+k,
                            get_val(b,i,j,k));
                        break;
                }
    }

    return out;
}

ip_mat * ip_mat_sum(ip_mat *a, ip_mat *b) {
    ip_mat * out = ip_mat_copy(a);

    unsigned int i, j, k;
    for(i=0; i<out->h; i++)
        for(j=0; j<out->w; j++)
            for(k=0; k < out->k; k++)
                set_val(out, i,j,k, 
                    (get_val(a,i,j,k) + get_val(b,i,j,k)));
    
    return out;
}

ip_mat * ip_mat_sub(ip_mat *a, ip_mat *b) {
    ip_mat * out = ip_mat_copy(a);

    unsigned int i, j, k;
    for(i=0; i<out->h; i++)
        for(j=0; j<out->w; j++)
            for(k=0; k < out->k; k++)
                set_val(out, i,j,k, 
                    (get_val(a,i,j,k) - get_val(b,i,j,k)));

    return out;
}

ip_mat * ip_mat_mul_scalar(ip_mat *a, float c) {
    ip_mat * out = ip_mat_copy(a);

    unsigned int i, j, k;
    for(i=0; i<out->h; i++)
        for(j=0; j<out->w; j++)
            for(k=0; k < out->k; k++)
                set_val(out, i,j,k, 
                    (get_val(a,i,j,k) * c));

    return out;
}

ip_mat *  ip_mat_add_scalar(ip_mat *a, float c) {
    ip_mat * out = ip_mat_copy(a);

    unsigned int i, j, k;
    for(i=0; i<out->h; i++)
        for(j=0; j<out->w; j++)
            for(k=0; k < out->k; k++)
                set_val(out, i,j,k, 
                    (get_val(a,i,j,k) + c));

    return out;
}

ip_mat * ip_mat_mean(ip_mat * a, ip_mat * b) {
    ip_mat * out = ip_mat_copy(a);

    unsigned int i, j, k;
    for(i=0; i<out->h; i++)
        for(j=0; j<out->w; j++)
            for(k=0; k < out->k; k++)
                set_val(out, i,j,k, 
                    (get_val(a,i,j,k) + get_val(b,i,j,k)) / 2);

    return out;
}

ip_mat * ip_mat_to_gray_scale(ip_mat * in) {
    ip_mat * out = ip_mat_copy(in);

    float mean, sum;
    unsigned int i, j, k;
    for(i=0; i<out->h; i++)
        for(j=0; j<out->w; j++) {
            sum = 0;
            for(k=0; k < out->k; k++)
                sum += get_val(in, i,j,k);
            mean = sum / (out->k);
            for(k=0; k < out->k; k++)
                set_val(out, i,j,k, mean);
        }

    return out;
}

ip_mat * ip_mat_blend(ip_mat * a, ip_mat * b, float alpha) {
    ip_mat * out = ip_mat_copy(a);

    float blend;
    unsigned int i, j, k;
    for(i=0; i<a->h; i++)
        for(j=0; j<a->w; j++)
            for(k=0; k < a->k; k++) {
                blend = (alpha * get_val(a,i,j,k))
                     + ((1.0-alpha) * get_val(b,i,j,k));
                set_val(out,i,j,k, blend);
            }

    return out;
}

ip_mat * ip_mat_brighten(ip_mat * a, float bright) {
    ip_mat * out = ip_mat_add_scalar(a, bright);
    clamp(out, 0, 255);
    return out;
}

ip_mat * ip_mat_corrupt(ip_mat * a, float amount) {
    ip_mat * out = NULL;
    ip_mat * temp1 = ip_mat_copy(a);
    ip_mat * temp2;
    ip_mat_init_random(temp1, 0, 0.5);
    temp2 = ip_mat_mul_scalar(temp1, amount);
    out = ip_mat_sum(a, temp2);
    ip_mat_free(temp1);
    ip_mat_free(temp2);
    return out;
}

ip_mat * ip_mat_padding(ip_mat * a, unsigned int pad_h, unsigned int pad_w) {
    ip_mat * out = ip_mat_create(a->h + (pad_h * 2), a->w + (pad_w * 2), 3, 0);

    unsigned int i, j, k;
    for(i=0; i<a->h; i++)
        for(j=0; j<a->w; j++)
            for(k=0; k < a->k; k++)
                set_val(out, i+pad_h,j+pad_w,k, 
                    get_val(a, i,j,k));
    return out;
}

ip_mat * ip_mat_convolve(ip_mat * a, ip_mat * f) {
    unsigned int i,j,k,fi,fj;
    float sum;
    unsigned int pad_h = (f->h - 1)/2;
    unsigned int pad_w = (f->w - 1)/2;
    ip_mat * out = ip_mat_create(a->h, a->w, a->k, 0); 
    ip_mat * padded = ip_mat_padding(a, pad_h, pad_w);

    for(k=0; k < out->k; k++)
        for(i=0; i < out->h; i++)
            for(j=0; j < out->w; j++) {
                sum = 0;
                for(fi=0; fi < f->h; fi++) 
                    for(fj=0; fj < f->w; fj++)
                        sum += get_val(padded, i+fi,j+fj,k) * 
                            get_val(f, fi, fj, 0);
                set_val(out, i,j,k, sum);
            }
    
    ip_mat_free(padded);
    clamp(out, 0, 255);
    return out;
}

/* 0 -1 0 
 * -1 5 -1
 * 0 -1 0 */
ip_mat * create_sharpen_filter() {
    ip_mat * out = ip_mat_create(3,3,3,0);
    
    unsigned int i,j,k;
    for(i=0; i<out->h; i++)
        for(j=0; j<out->w; j++)
            for(k=0; k<out->k; k++)
                if(i == 1 && j == 1)
                    set_val(out, i,j,k, 5);
                else if ((i == 1 || j == 1) && 
                    (i == 0 || i == 2 || j == 0 || j == 2))
                    set_val(out, i,j,k, -1);
                
    return out;
}

/* -1 -1 -1 
 * -1 8 -1
 * -1 -1 -1 */
ip_mat * create_edge_filter() {
    ip_mat * out = ip_mat_create(3,3,3,0);
    unsigned int i,j,k;
    for(i=0; i<out->h; i++)
        for(j=0; j<out->w; j++)
            for(k=0; k<out->k; k++)
                if(i == 1 && j == 1) 
                    set_val(out, i,j,k, 8);
                else set_val(out,i,j,k, (-1));
    return out;
}

/* -2 -1 0 
 * -1 1 1
 * 0 1 2 */
ip_mat * create_emboss_filter() {
    ip_mat * out = ip_mat_create(3,3,3,0);

    unsigned int i,j,k;
    for(i=0; i<out->h; i++)
        for(j=0; j<out->w; j++)
            for(k=0; k<out->k; k++) {
                set_val(out, 0,0,k, -2);
                set_val(out, 0,1,k, -1);
                set_val(out, 1,0,k, -1);
                set_val(out, 1,1,k, 1);
                set_val(out, 1,2,k, 1);
                set_val(out, 2,1,k, 1);
                set_val(out, 2,2,k, 2);
            }

    return out;
}

ip_mat * create_average_filter(unsigned int w, unsigned int h, unsigned int k) {
    ip_mat * out = ip_mat_create(h,w,k,0);

    unsigned int i,j,l;
    for(i=0; i<out->h; i++)
        for(j=0; j<out->w; j++)
            for(l=0; l<out->k; l++)
                set_val(out, i,j,l, 
                    (1.0 / ((w * h)))
                );

    return out;
}

float gauss_noise(int x, int y, float sigma) {
    return 1 / (2 * PI * sigma * exp((pow(x, 2) + pow(y, 2)) / (2 * sigma)));
}

ip_mat * create_gaussian_filter(unsigned int w, unsigned int h, unsigned int l, float sigma) {
    ip_mat * out = NULL;
    ip_mat * temp = ip_mat_create(h, w, l, 0);
    
    int cx = h/2, cy = w/2;
    int dx, dy;
    float sum = 0;

    unsigned int i,j;
    for(i=0; i<temp->h; i++)
        for(j=0; j<temp->w; j++) {
            dx = i - cx, dy = j - cy;
            set_val(temp, i,j,0, gauss_noise(dx, dy, sigma));
            sum += get_val(temp, i,j,0);
        }
            

    out = ip_mat_mul_scalar(temp, (1/sum));
    ip_mat_free(temp);
    return out;
}

void rescale(ip_mat * t, float new_max) {
	unsigned int i,j,k;
	float result;
	compute_stats(t);

    for(k=0; k<t->k; k++)
        for(i=0; i<t->h; i++)
            for(j=0; j<t->w; j++) {
                result = ((get_val(t,i,j,k) - t->stat[k].min) / 
                	(t->stat[k].max - t->stat[k].min)) * new_max;
                set_val(t, i,j,k, result);
            }
}

void clamp(ip_mat * t, float low, float high) {
	unsigned int i, j, k;
	float pixel;
    for(i=0; i<t->h; i++)
        for(j=0; j<t->w; j++)
        	for(k=0; k<t->k; k++) {
        		pixel = get_val(t, i,j,k);
        		pixel = max(low, pixel);
        		pixel = min(high, pixel);
        		set_val(t, i,j,k, pixel);
        	}
}
