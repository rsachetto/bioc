//
// Created by sachetto on 31/03/2022.
//

#include "sequence_record.h"
#include "3d_party/stb/stb_ds.h"
#include "file_utils/file_utils.h"
#include "parsing_utils.h"
#include "sequence.h"
#include <stdio.h>
#include <ctype.h>

bioc_sequence_record new_sequence_record() {
    bioc_sequence_record s;
    s.seq = seq("");
    s.seq = calloc(1, sizeof(bioc_seq));
    s.indexes = ht_create();
    return s;
}

#define ADD_OR_BREAK(content) \
    if (!*(content)) break;   \
    (content)++;

static bioc_sequence_record_list parse_fasta_file(char *content, size_t size) {

    bioc_sequence_record_list seqs = NULL;

    while (*content) {

        while (*content == ';') {
            skip_line(&content);
        }

        while (*content == '>') {

            ADD_OR_BREAK(content);
            
            bioc_sequence_record s = new_sequence_record();

            int c = 0;

            while (*content != '\n') {
                s.id[c] = *content;
                c++;
                ADD_OR_BREAK(content);
            }

            s.id[c] = 0;

            ADD_OR_BREAK(content);

            while (*content != '>') {
                if(*content != '\n' && !isspace(*content)) {
                    //TODO: check for error
                    bioc_add_nucleotide(s.seq, *content);
                }
                ADD_OR_BREAK(content);
            }

            arrput(seqs, s);
        }
    }

    return seqs;
}

static bioc_sequence_record_list parse_fastq_file(char *content, size_t size) {

    bioc_sequence_record_list seqs = NULL;

    while (*content) {
 
        while(isspace(*content)) {
            content++;
        }
        
        if (*content == '@') {

            ADD_OR_BREAK(content);
            
            bioc_sequence_record s = new_sequence_record();

            int c = 0;

            while (*content != '\n') {
                s.id[c] = *content;
                c++;
                ADD_OR_BREAK(content);
            }

            s.id[c] = 0;

            ADD_OR_BREAK(content);

            while (*content != '+') {
                if(*content != '\n' && !isspace(*content)) {
                    //TODO: check for error
                    bioc_add_nucleotide(s.seq, *content);
                }
                ADD_OR_BREAK(content);
            }

            arrput(seqs, s);

            while(*content != '@' ) {
                content++;
            }
        }
    }

    return seqs;
}
bioc_sequence_record_list bioc_parse_sequence_file(char *sequence_file_path, enum bioc_sequence_file_format sequence_file_f) {

    size_t sequence_file_size = 0;

    char *file_content = read_entire_file_with_mmap(sequence_file_path, &sequence_file_size);

    if (file_content == NULL) {
        fprintf(stderr, "Error opening file %s. Exiting!\n", sequence_file_path);
        return NULL;
    }

    bioc_sequence_record_list l = NULL;

    switch (sequence_file_f) {
        case FASTA:
            l = parse_fasta_file(file_content, sequence_file_size);
            break;
        case GENBANK:
            fprintf(stderr, "GenBank format not implemented yet!\n");
            break;        
        case FASTQ_ILLUMINA:
            l = parse_fastq_file(file_content, sequence_file_size);
            break;         
        default:
            fprintf(stderr, "Invalid format!\n");
            break;
    }

    munmap(file_content, sequence_file_size);
    return l;
}
