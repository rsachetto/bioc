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
    printf("%Q starts_with GUC: %s\n", s2, bioc_seq_starts_with(s2, "GUC") ? "Yes": "No");
    printf("%Q starts_with AUG: %s\n", s2, bioc_seq_starts_with(s2, "AUG") ? "Yes": "No");
    printf("%Q starts_with AUG: %s\n", s2, bioc_seq_starts_with_bounded(s2, "AUG", 3, s2->len) ? "Yes": "No");


    printf("%Q ends_with UUG: %s\n", s2, bioc_seq_ends_with(s2, "UUG") ? "Yes": "No");
    printf("%Q ends_with AUG: %s\n", s2, bioc_seq_ends_with(s2, "AUG") ? "Yes": "No");
    printf("%.*s ends_with AUG: %s\n", 18, s2->nucleotides, bioc_seq_ends_with_bounded(s2, "AUG", 0, 18) ? "Yes": "No");


    //TODO:
    //add a function to work on a list:
    //my_rna.startswith(("UCC", "UCA", "UCG"), 1)
    // True

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

    bioc_seq *s5 = seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG");
    tmp = bioc_translate(s5, BIOC_STANDART, false);
    printf("Translation: %Q\n", tmp);
    bioc_seq_free(tmp);

    tmp = bioc_translate(s5, BIOC_STANDART, true);
    printf("Translation: %Q\n", tmp);
    bioc_seq_free(tmp);

    tmp = bioc_translate(s5, BIOC_VERTEBRATE_MITOCHONDRIAL, false);
    printf("Translation vertebrate: %Q\n", tmp);
    bioc_seq_free(tmp);

    tmp = bioc_translate(s5, BIOC_VERTEBRATE_MITOCHONDRIAL, true);
    printf("Translation vertebrate stop codon: %Q\n", tmp);
    bioc_seq_free(tmp);

    bioc_seq *s6 = seq("ACGT");
    bioc_seq  *s7 = seq("AACCGG");
    tmp = bioc_concatenate(s6, s7);
    printf("Concat: %Q\n", tmp);
    bioc_seq_free(tmp);

    bioc_seq_free(s);
    bioc_seq_free(s2);
    bioc_seq_free(s3);
    bioc_seq_free(s4);
    bioc_seq_free(s5);
    bioc_seq_free(s6);
    bioc_seq_free(s7);
}
