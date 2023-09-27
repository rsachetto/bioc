//
// Created by sachetto on 31/03/2022.
//

#include "../src/bioc.h"
#include <stdio.h>

int main() {

    sequence *s = seq("AGTACACTGGT");
    for(int  i = 0; i < 10000; i++) {
        printf("%Q\n", s);
        printf("%Q\n", bioc_complement(s));
        printf("%Q\n", bioc_reverse_complement(s));

    }
}
