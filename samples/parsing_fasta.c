//
// Created by sachetto on 31/03/2022.
//

#include "../src/bioc.h"
#include <stdio.h>

int main() {
    sequence_record_list seqs;
    seqs = parse_sequence_file("sample.fasta", FASTA);

    sequence_record *s;

    FOR_EACH_RECORD(s, seqs) {
        printf("%s\n", s->id);
        printf("%s\n", s->seq->nucleotides);
    }
    
    sequence_record_list seqs_fq;
    seqs_fq = parse_sequence_file("sample.fq", FASTQ_ILLUMINA);

    FOR_EACH_RECORD(s, seqs_fq) {
        printf("%s\n", s->id);
        printf("%s\n", s->seq->nucleotides);
    }    
}
