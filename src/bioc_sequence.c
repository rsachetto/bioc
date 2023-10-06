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
#include "3d_party/hash/ht.h"
#include "3d_party/stb/stb_ds.h"

typedef struct bioc_translation_table_entry_t {
    char *key;
    char value;
} bioc_translation_table_entry;

static bool registred = false;

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

bioc_seq * empty_seq_with_len(int64_t len) {

    if(!registred) {
        //TODO: maybe this is not the place for this initialization.
        register_printf_specifier('Q', bioc_print_sequence, print_widget_arginfo);
        registred = true;
    }

    bioc_seq *new_seq = (bioc_seq *) calloc(1, sizeof(bioc_seq));

    arrsetlen(new_seq->nucleotides, len + 1);

    new_seq->len = len;

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
}


static void make_table(bioc_translation_table_entry **table, bioc_genetic_code code) {

    assert(*table);

    static const char *base1 = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG";
    static const char *base2 = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG";
    static const char *base3 = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG";

    const char * translation_data;

    switch (code) {
        case BIOC_STANDART:
            translation_data = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
            break;
        case BIOC_VERTEBRATE_MITOCHONDRIAL:
            translation_data = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG";
            break;
        case BIOC_INVERTEBRATE_MITOCHONDRIAL:
            translation_data = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG";
            break;
    }

    size_t l = 64;

    for(size_t i = 0; i < l; i++) {
        char *combination = calloc(4, 1);
        combination[0] = base1[i];
        combination[1] = base2[i];
        combination[2] = base3[i];
        combination[3] = '\0';
        shput(*table, combination, translation_data[i]);

        if(combination[0] == 'T') {
            combination[0] = 'U';
        }

        if(combination[1] == 'T') {
            combination[1] = 'U';
        }

        if(combination[2] == 'T') {
            combination[2] = 'U';
        }

        shput(*table, combination, translation_data[i]);
    }

}

bioc_seq * bioc_translate(bioc_seq *s, bioc_genetic_code code, bool to_stop) {

    int64_t s_len = s->len;

    if (s_len < 3) {
        return seq(NULL);
    }

    static bioc_translation_table_entry *standard = NULL;
    static bioc_translation_table_entry *vertebrate = NULL;
    static bioc_translation_table_entry *invertebrate = NULL;

    bioc_translation_table_entry *tab = NULL;

    switch (code) {
        case BIOC_STANDART:
        {
            if(standard == NULL) {
                sh_new_strdup(standard);
                shdefault(standard, '?');
                make_table(&standard, code);
            }
            tab = standard;
            break;
        }
        case BIOC_VERTEBRATE_MITOCHONDRIAL:
        {
            if(vertebrate == NULL) {
                sh_new_strdup(vertebrate);
                shdefault(standard, '?');
                make_table(&vertebrate, code);
            }
            tab = vertebrate;
            break;

        }
        case BIOC_INVERTEBRATE_MITOCHONDRIAL:
        {
            if(invertebrate == NULL) {
                sh_new_strdup(invertebrate);
                shdefault(standard, '?');
                make_table(&invertebrate, code);
            }
            tab = invertebrate;
            break;

        }
    }

    int64_t rest = s_len % 3;
    if(rest != 0) {
        s_len -= rest;
    }

    bioc_seq *translated = seq(NULL);

    int64_t i = 0;
    while (i < s_len) {

        char b[4] = {0,0,0,0}; //Do we need this?

        b[0] = s->nucleotides[i];
        b[1] = s->nucleotides[i + 1];
        b[2] = s->nucleotides[i + 2];

        i +=3;

        char tr = shget(tab, b);

        if(to_stop) {
            if (tr == '*') {
                break;
            }
        }

        bioc_add_nucleotide(translated, tr);
    }
    arrput(translated->nucleotides, '\0');

    return translated;
}

bioc_seq *bioc_concatenate(bioc_seq *s1, bioc_seq *s2) {

    int64_t s1_l = s1->len;
    int64_t s2_l = s2->len;

    bioc_seq *concat = empty_seq_with_len(s1_l + s2_l);

    //TODO: memcpy
    for(int i = 0; i < s1_l; i++) {
        concat->nucleotides[i] = s1->nucleotides[i];
    }

    for(int i = 0; i < s2_l; i++) {
        concat->nucleotides[i + s1_l] = s2->nucleotides[i];
    }

    concat->nucleotides[s1_l + s2_l] = '\0';

    return concat;

}