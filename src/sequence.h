//
// Created by sachetto on 31/03/2022.
//

#ifndef BIOC_SEQUENCE_H
#define BIOC_SEQUENCE_H

#include <stdlib.h>
#include <stdint.h>

typedef struct bioc_seq_t {
    char *nucleotides;
    size_t len;
    //TODO: more
} bioc_seq;

bioc_seq * seq(char *nucleotides);
bioc_seq * bioc_complement(bioc_seq *s);
bioc_seq * bioc_reverse(bioc_seq *s);
bioc_seq * bioc_reverse_complement(bioc_seq *s);
void       bioc_add_nucleotide(bioc_seq *seq, char nucleotide);
uint64_t   bioc_count_nucleotide(bioc_seq *seq, char nucleotide);
double     bioc_GC_content(bioc_seq *seq);
#endif//BIOC_SEQUENCE_H
