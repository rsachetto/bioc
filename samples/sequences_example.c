//
// Created by sachetto on 31/03/2022.
//

#include "../src/bioc.h"
#include <stdio.h>

int main() {
    bioc_seq *s = seq("GATCGATGGGCCTATATAGGATCGAAAATCGc");
    printf("SEQUENCE:           %Q\n", s);

    bioc_seq *tmp = bioc_complement(s);
    printf("COMPLEMENT:         %Q\n", tmp);
    bioc_seq_free(tmp);

    tmp = bioc_reverse(s);
    printf("REVERSE:            %Q\n", tmp);
    bioc_seq_free(tmp);

    tmp = bioc_reverse_complement(s);
    printf("REVERSE COMPLEMENT: %Q\n", tmp);
    bioc_seq_free(tmp);

    printf("Gs: %ld\n", bioc_count_nucleotide(s, 'G'));
    printf("GC content: %lf\n", bioc_GC_content(s));

    printf("find AUG: %ld\n", bioc_find(s, "AUG"));
    printf("find GGC: %ld\n", bioc_find(s, "GGC"));

    bioc_seq *s2 = seq("GUCAUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAGUUG");
    printf("find AUG: %ld\n", bioc_find(s2, "AUG"));

    bioc_seq *s3 =  seq("AGTACACTGGT");
    printf("find ACT: %ld\n", bioc_find(s3, "ACT"));
    printf("find TAG: %ld\n", bioc_find(s3, "TAG"));
    printf("count A: %ld\n", bioc_count(s3, "A"));
    printf("count ACT: %ld\n", bioc_count(s3, "ACT"));

    bioc_seq *s4 =  seq("AAAA");
    printf("count AA: %ld\n", bioc_count(s4, "AA"));
    printf("count overlap AA: %ld\n", bioc_count_overlap(s4, "AA"));

    tmp = bioc_transcribe(s3);
    printf("TRANSCRIPTION:      %Q\n", tmp);

    bioc_seq *tmp2 = bioc_back_transcribe(tmp);

    printf("BACK TRANSCRIPTION: %Q\n", tmp2);
    bioc_seq_free(tmp2);
    bioc_seq_free(tmp);


    bioc_seq_free(s);
    bioc_seq_free(s2);
    bioc_seq_free(s3);
    bioc_seq_free(s4);
}
