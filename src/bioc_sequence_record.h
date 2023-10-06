//
// Created by sachetto on 31/03/2022.
//

#ifndef BIOC_BIOC_SEQUENCE_RECORD_H
#define BIOC_BIOC_SEQUENCE_RECORD_H

#include "3d_party/hash/ht.h"
#include "3d_party/stb/stb_ds.h"
#include "bioc_sequence.h"

#define FOR_EACH_RECORD(record, records) \
    for(int it_count = 0; it_count < arrlen(records); it_count++)\
        if((record = &(records[it_count])))

enum bioc_sequence_file_format {
    FASTA,
    GENBANK,
    FASTQ_ILLUMINA
};

typedef struct bioc_sequence_record_t {
    bioc_seq *seq;
    char id[256];
    ht *indexes;
} bioc_sequence_record;

typedef bioc_sequence_record *bioc_sequence_record_list;

bioc_sequence_record_list bioc_parse_sequence_file(char *sequence_file_path, enum bioc_sequence_file_format sequence_file_f);
void bioc_sequence_record_list_free(bioc_sequence_record_list l);
void bioc_sequence_record_free(bioc_sequence_record *r);

#endif // BIOC_BIOC_SEQUENCE_RECORD_H
