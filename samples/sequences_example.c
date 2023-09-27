//
// Created by sachetto on 31/03/2022.
//

#include "../src/bioc.h"
#include <stdio.h>

int main() {
    bioc_seq *s = seq("GATCGATGGGCCTATATAGGATCGAAAATCGc");
    printf("SEQUENCE:           %Q\n", s);
    printf("COMPLEMENT:         %Q\n", bioc_complement(s));
    printf("REVERSE:            %Q\n", bioc_reverse(s));
    printf("REVERSE COMPLEMENT: %Q\n", bioc_reverse_complement(s));
    printf("Gs: %ld\n", bioc_count_nucleotide(s, 'G'));
    printf("GC content: %lf\n", bioc_GC_content(s));

    printf("find AUG: %ld\n", bioc_find(s, "AUG", 0, s->len));
    printf("find GGC: %ld\n", bioc_find(s, "GGC", 0, s->len));

    bioc_seq *s2 = seq("GUCAUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAGUUG");
    printf("find AUG: %ld\n", bioc_find(s2, "AUG", 0, s2->len));
}
