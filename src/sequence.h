//
// Created by sachetto on 31/03/2022.
//

#ifndef BIOC_SEQUENCE_H
#define BIOC_SEQUENCE_H

#include <stdlib.h>
#include <stdint.h>

typedef struct bioc_seq_t {
    char *nucleotides;
    int64_t len;
    //TODO: more
} bioc_seq;

bioc_seq * seq(char *nucleotides);
bioc_seq * bioc_complement(bioc_seq *s);
bioc_seq * bioc_reverse(bioc_seq *s);
bioc_seq * bioc_reverse_complement(bioc_seq *s);
bioc_seq * bioc_transcribe(bioc_seq *s);
bioc_seq * bioc_back_transcribe(bioc_seq *s);
int64_t    bioc_find(bioc_seq *seq, char *sub);
int64_t    bioc_find_with_bounds(bioc_seq *seq, char *sub, int64_t start, int64_t end);
int64_t    bioc_count(bioc_seq *seq, char *sub);
int64_t    bioc_count_overlap(bioc_seq *seq, char *sub);
int64_t    bioc_count_with_bounds(bioc_seq *seq, char *sub, int64_t start, int64_t end);
int64_t    bioc_count_with_bounds_overlap(bioc_seq *seq, char *sub, int64_t start, int64_t end);
void       bioc_add_nucleotide(bioc_seq *seq, char nucleotide);
uint64_t   bioc_count_nucleotide(bioc_seq *seq, char nucleotide);
double     bioc_GC_content(bioc_seq *seq);
void bioc_seq_free(bioc_seq *seq);
#endif//BIOC_SEQUENCE_H
