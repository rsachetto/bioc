#include "parsing_utils.h"

void skip_line(char **src) {

    while(**src != '\n') {
        (*src)++;
    }

    (*src)++;
}

