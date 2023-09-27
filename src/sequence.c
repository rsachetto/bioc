//
// Created by sachetto on 31/03/2022.
//

#define _GNU_SOURCE

#include "sequence.h"
#include <string.h>
#include <printf.h>

#include <stdio.h>
#include <stdbool.h>

#define STB_DS_IMPLEMENTATION
#include "3d_party/stb/stb_ds.h"

inline static char nucleotide_complement(char n) {
    switch (n) {
        case 'A': return 'T';
        case 'a': return 't';
        case 'T': return 'A';
        case 't': return 'a';
        case 'U': return 'A';
        case 'u': return 'a';
        case 'G': return 'C';
        case 'g': return 'c';
        case 'C': return 'G';
        case 'c': return 'g';
        default : return 'U';
    }
}

static int bioc_print_sequence(FILE *stream, const struct printf_info *info,  const void *const *args) {
    const bioc_seq *s;
    char *buffer;
    int len;

    /* Format the output into a string. */
    s = *((const bioc_seq **) (args[0]));
    len = asprintf (&buffer, "%s", s->nucleotides);

    if (len == -1)
        return -1;

    /* Pad to the minimum field width and print to the stream. */
    len = fprintf (stream, "%*s", (info->left ? -info->width : info->width), buffer);

    /* Clean up and return. */
    free (buffer);
    return len;
}


int print_widget_arginfo (const struct printf_info *info, size_t n, int *argtypes, int *size) {
    if (n > 0)
        argtypes[0] = PA_POINTER;
    return 1;
}

bioc_seq * seq(char *nucleotides) {
    static bool registred = false;
    if(!registred) {
        //TODO: maybe this is not the place for this initialization.
        register_printf_specifier('Q', bioc_print_sequence, print_widget_arginfo);
        registred = true;
    }

    bioc_seq *new_seq = (bioc_seq *) calloc(1, sizeof(bioc_seq));

    if(nucleotides == NULL || nucleotides[0] == '\0') {
        return new_seq;
    }

    size_t n_len = strlen(nucleotides) + 1;

    arrsetlen(new_seq->nucleotides, n_len);

    memcpy(new_seq->nucleotides, nucleotides, n_len);
    new_seq->len = n_len - 1;

    return new_seq;
}

bioc_seq *bioc_complement(bioc_seq *s) {
    if(s == NULL) {
        return NULL;
    }

    bioc_seq *comp = seq(NULL);
    size_t s_len = s->len;

    arrsetlen(comp->nucleotides, s_len + 1);
    comp->nucleotides[s_len] = 0;

    for(size_t i = 0; i < s_len; i++) {
        comp->nucleotides[i] = nucleotide_complement(s->nucleotides[i]);
    }

    return comp;
}

bioc_seq *bioc_reverse(bioc_seq *s) {

    if(s == NULL) {
        return NULL;
    }

    bioc_seq *rev = seq(NULL);
    int s_len = (int)s->len;

    arrsetlen(rev->nucleotides, s_len + 1);
    rev->nucleotides[s_len] = 0;

    for(int i = s_len-1; i >= 0; i--) {
        rev->nucleotides[s_len - 1 - i] = s->nucleotides[i];
    }

    return rev;
}

bioc_seq *bioc_reverse_complement(bioc_seq *s) {

    if(s == NULL) {
        return NULL;
    }

    bioc_seq *rev = seq(NULL);
    int s_len = (int)s->len;

    arrsetlen(rev->nucleotides, s_len + 1);
    rev->nucleotides[s_len] = 0;

    for(int i = s_len-1; i >= 0; i--) {
        rev->nucleotides[s_len - 1 - i] = nucleotide_complement(s->nucleotides[i]);
    }

    return rev;
}

uint64_t bioc_count_nucleotide(bioc_seq *seq, char nucleotide) {
    int n = seq->len;
    uint64_t count = 0;
    //TODO: omp?
    for(int i = 0; i < n; i++) {
        if(seq->nucleotides[i] == nucleotide) {
            count++;
        }
    }

    return count;
}

//TODO: macro maybe?
double bioc_GC_content(bioc_seq *seq) {
    return 100.0 * (double)(bioc_count_nucleotide(seq, 'G') + bioc_count_nucleotide(seq, 'C') + bioc_count_nucleotide(seq, 'g') + bioc_count_nucleotide(seq, 'c'))/seq->len;
}

void bioc_add_nucleotide(bioc_seq *seq, char nucleotide) {
    arrput(seq->nucleotides, nucleotide);
    seq->len++;
}

