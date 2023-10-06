//

#include "../src/bioc.h"
#include <stdio.h>

int main() {
    bioc_sequence_record_list seqs;
    seqs = bioc_parse_sequence_file("sample.fasta", FASTA);

    bioc_sequence_record *s;

    FOR_EACH_RECORD(s, seqs) {
        printf("%s\n", s->id);
        printf("%s\n", s->seq->nucleotides);
    }
    
    bioc_sequence_record_list seqs_fq;
    seqs_fq = bioc_parse_sequence_file("sample.fq", FASTQ_ILLUMINA);

    FOR_EACH_RECORD(s, seqs_fq) {
        printf("%s\n", s->id);
        printf("%s\n", s->seq->nucleotides);
    }

    bioc_sequence_record_list_free(seqs);
    bioc_sequence_record_list_free(seqs_fq);
}
