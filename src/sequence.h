//
// Created by sachetto on 31/03/2022.
//

#ifndef BIOC_SEQUENCE_H
#define BIOC_SEQUENCE_H

#include <stdlib.h>
#include <stdint.h>

typedef struct seq {
    char *nucleotides;
    size_t len;
    //TODO: more
} sequence;

sequence * seq(char *nucleotides);
sequence * complement(sequence *s);
sequence * reverse(sequence *s);
sequence * reverse_complement(sequence *s);
void add_nucleotide(sequence *seq, char nucleotide);
uint64_t count_nucleotide(sequence *seq, char nucleotide);
double GC_content(sequence *seq);
#endif//BIOC_SEQUENCE_H
