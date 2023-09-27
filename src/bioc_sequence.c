//
// Created by sachetto on 31/03/2022.
//

#include <bits/stdint-intn.h>
#define _GNU_SOURCE

#include "bioc_sequence.h"
#include <printf.h>
#include <string.h>

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

    int64_t n_len = (int64_t) strlen(nucleotides) + 1;

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
    comp->len = s_len;

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
    rev->len = s_len;

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
    int64_t s_len = (int)s->len;
    rev->len = s_len;

    arrsetlen(rev->nucleotides, s_len + 1);
    rev->nucleotides[s_len] = 0;

    for(int64_t i = s_len-1; i >= 0; i--) {
        rev->nucleotides[s_len - 1 - i] = nucleotide_complement(s->nucleotides[i]);
    }

    return rev;
}

bioc_seq *bioc_transcribe(bioc_seq *s) {

    if(s == NULL) {
        return NULL;
    }

    //TODO: this is a common pattern, make a function or a macro
    bioc_seq *trans = seq(NULL);
    int64_t s_len = (int)s->len;
    trans->len = s_len;

    arrsetlen(trans->nucleotides, s_len + 1);
    trans->nucleotides[s_len] = 0;

    for(int64_t i = 0; i < s_len; i++) {
        switch (s->nucleotides[i]) {
            case 'T':
                trans->nucleotides[i] = 'U';
                break;
            case 't':
                trans->nucleotides[i] = 'u';
                break;
           default:
                trans->nucleotides[i] = s->nucleotides[i];
        }
    }

    return trans;
}

bioc_seq *bioc_back_transcribe(bioc_seq *s) {

    if(s == NULL) {
        return NULL;
    }
    //TODO: this is a common pattern, make a function or a macro
    bioc_seq *trans = seq(NULL);
    int64_t s_len = (int)s->len;
    trans->len = s_len - 1;

    arrsetlen(trans->nucleotides, s_len + 1);
    trans->nucleotides[s_len] = 0;

    for(int64_t i = 0; i < s_len; i++) {
        switch (s->nucleotides[i]) {
            case 'U':
                trans->nucleotides[i] = 'T';
                break;
            case 'u':
                trans->nucleotides[i] = 't';
                break;
           default:
                trans->nucleotides[i] = s->nucleotides[i];
        }
    }

    return trans;
}

uint64_t bioc_count_nucleotide(bioc_seq *seq, char nucleotide) {
    int64_t n = seq->len;
    uint64_t count = 0;
    //TODO: omp?
    for(int64_t i = 0; i < n; i++) {
        if(seq->nucleotides[i] == nucleotide) {
            count++;
        }
    }

    return count;
}

//TODO: macro maybe?
double bioc_GC_content(bioc_seq *seq) {
    return 100.0 * (double)(bioc_count_nucleotide(seq, 'G') + bioc_count_nucleotide(seq, 'C') + bioc_count_nucleotide(seq, 'g') + bioc_count_nucleotide(seq, 'c'))/(double)seq->len;
}

void bioc_add_nucleotide(bioc_seq *seq, char nucleotide) {
    arrput(seq->nucleotides, nucleotide);
    seq->len++;
}


int64_t bioc_find_with_bounds(bioc_seq *seq, char *sub, int64_t start, int64_t end) {
    //TODO: check bounds
    if(end > seq->len) end = seq->len; //TODO: warning??

    char *tmp = seq->nucleotides + start;

    char *find = strstr(tmp, sub);

    if(find == NULL) {
        return -1;
    }

    int64_t index = (long) (find - seq->nucleotides);

    if(index > end) {
        return -1;
    }

    return index;

}


int64_t bioc_find(bioc_seq *seq, char *sub) {
    return bioc_find_with_bounds(seq, sub, 0, seq->len);
}

static int64_t bioc_count_common(bioc_seq *seq, char *sub, int64_t start, int64_t end, bool overlap) {
    //TODO: check bounds
    if(end > seq->len) end = seq->len; //TODO: warning??

    char *tmp = seq->nucleotides + start;

    char *find = strstr(tmp, sub);

    int64_t count = 0;
    size_t n_sub = 1;

    if(!overlap) n_sub = strlen(sub);

    while(find != NULL) {
        count++;
        find = strstr(find + n_sub, sub);
        int64_t index = (long) (find - seq->nucleotides);

        if(index > end) {
            return -1;
        }
    }

    return count;
}

int64_t bioc_count(bioc_seq *seq, char *sub) {
    return  bioc_count_common(seq, sub, 0, seq->len, false);
}

int64_t bioc_count_overlap(bioc_seq *seq, char *sub) {
    return  bioc_count_common(seq, sub, 0, seq->len, true);
}

int64_t bioc_count_bounded(bioc_seq *seq, char *sub, int64_t start, int64_t end) {
    return  bioc_count_common(seq, sub, start, end, false);
}

int64_t bioc_count_with_overlap_bounded(bioc_seq *seq, char *sub, int64_t start, int64_t end) {
    return  bioc_count_common(seq, sub, start, end, true);
}

bool bioc_seq_starts_with_bounded(bioc_seq *seq, const char *prefix, int64_t start, int64_t end) {
    //TODO: check bound properly
    int64_t s = (int64_t) strlen(seq->nucleotides + start);

    if(seq->nucleotides + start + s > seq->nucleotides + end) {
        //TODO: check
        return false;
    }
    return strncmp(prefix, seq->nucleotides + start, strlen(prefix)) == 0;
}

bool bioc_seq_starts_with(bioc_seq *seq, const char *prefix) {
    return bioc_seq_starts_with_bounded(seq, prefix, 0, seq->len);
}

bool bioc_seq_ends_with_bounded(bioc_seq *seq, char* suffix, int64_t start, int64_t end) {
    //TODO: check bound properly
    int64_t seq_l = seq->len;
    int64_t suffix_l = (int64_t)strlen(suffix);

    if (seq_l < suffix_l) {
        return false;
    }

    return strncmp(seq->nucleotides + start + (end - suffix_l), suffix, suffix_l) == 0;
}

bool bioc_seq_ends_with(bioc_seq *seq, char* suffix) {
    return bioc_seq_ends_with_bounded(seq, suffix, 0, seq->len);
}

void bioc_seq_free(bioc_seq *seq) {
    arrfree(seq->nucleotides);
    free(seq);
    seq = NULL;
}
