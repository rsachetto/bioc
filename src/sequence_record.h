//
// Created by sachetto on 31/03/2022.
//

#ifndef BIOC_SEQUENCE_RECORD_H
#define BIOC_SEQUENCE_RECORD_H

#include "sequence.h"
#include "3d_party/hash/ht.h"
#include "3d_party/stb/stb_ds.h"

#define FOR_EACH_RECORD(record, records) \
    for(int it_count = 0; it_count < arrlen(records); it_count++)\
        if((record = &(records[it_count])))

enum sequence_file_format {
    FASTA,
    GENBANK,
    FASTQ_ILLUMINA
};

typedef struct sequence_record_t {
    sequence *seq;
    char id[256];
    ht *indexes;
} sequence_record;

typedef sequence_record *sequence_record_list;

sequence_record_list parse_sequence_file(char *sequence_file_path, enum sequence_file_format sequence_file_f);

#endif//BIOC_SEQUENCE_RECORD_H
