#include <stdio.h> 
#include <inttypes.h>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <string.h>


# define TAP126 ((uint64_t)1<<62)

#define TAP99 ((uint64_t)1<<35)

#define BUFF_SIZE 150000

#define ALLZERO64 0
#define ALLONE64 0xFFFFFFFFFFFFFFFF

/*
This is an LFSR based PRNG with 
Irreducible polynomial
p(x) = x^127 + x^126 + X^99 + 1
and 128 bit state register. 

*/


void init_rand_register(uint64_t * x){
    int i;
    int r;
    (*x) = 0;
    for (i=0;i<8;i++){
        
        r = rand() % 128;
        //printf("%d\n",r);
        (*x) <<= 8;
        (*x) = (*x) | (uint64_t) r;
    }

}



int main(){


    int sweep_count;
    int i;
    
    uint64_t out_bit;
    uint64_t in_bit;
    uint64_t carry;
    FILE * fp_in;
    FILE * fp_out;
    uint8_t *s;

    s = malloc(BUFF_SIZE*sizeof(uint8_t));

    
    /*  a and b combined creates the state register of length 128 bits. 
        the bits are indexed b_127,b_126, ..., a_63,a_62, ..., a_1, a_0.
        b_127 is the output bit out_bit. This if the msb of b. 
        a_0 is the input bit a_0. This is the lsb of a. 
    */
    uint64_t b;
    uint64_t a;

    uint64_t ra;
    uint64_t rb;
   
    
    /* initialize random seed: */
    srand (time(NULL));
    
    // initialize the lfsr state a random seed. In real application this will be done using a pregenerated random seed. 
    init_rand_register(&a);
    init_rand_register(&b);
    ra = 0;
    rb = 0;

    
    printf("%" PRIu64 " \n", a);
    printf("%" PRIu64 " \n", b);
    fp_in = fopen("error_corrected_key.txt","r");
    fscanf(fp_in,"%s",s);

    printf("*** %d\n",(int)strlen(s));

    fp_out = fopen("finalKey.bin","wb");

  

    for(sweep_count=0;sweep_count<51;sweep_count++){
        i = 0;
        ra = 0;
        rb = 0;
        while(s[i]){
            //tap all the bits. This part implements the polynomial p(x)
            out_bit = b>>63;
            in_bit = out_bit ^ ((TAP126 & b) >> 62) ^ ((TAP99 & b) >> 35); 
            carry = a >> 63;

            // shift the registers. 
            b = (b<<1) | carry;
            a = (a<<1) | in_bit; 

            if (s[i]=='1'){
                ra = ra^a;
                rb = rb^b;
            }
            i++;
        }

        fwrite(&ra, sizeof(uint64_t), 1, fp_out);
        fwrite(&rb, sizeof(uint64_t), 1, fp_out);
    }        

}


