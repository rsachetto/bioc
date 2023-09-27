//
// Created by sachetto on 31/03/2022.
//

#include "../src/bioc.h"
#include <stdio.h>

int main() {
    sequence *s = seq("GATCGATGGGCCTATATAGGATCGAAAATCGc");
    printf("SEQUENCE:           %Q\n", s);
    printf("COMPLEMENT:         %Q\n", complement(s));
    printf("REVERSE:            %Q\n", reverse(s));
    printf("REVERSE COMPLEMENT: %Q\n", reverse_complement(s));
    printf("Gs: %ld\n", count_nucleotide(s, 'G'));
    printf("GC content: %lf\n", GC_content(s));
}
