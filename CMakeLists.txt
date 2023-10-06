cmake_minimum_required(VERSION 3.21)
project(bioc C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -fsanitize=address,undefined")

add_library(bioc src/bioc.h src/bioc_sequence.c src/bioc_sequence.h
            src/bioc_sequence_record.c src/bioc_sequence_record.h
            src/file_utils/file_utils.c src/file_utils/file_utils.h
            src/3d_party/hash/ht.c src/3d_party/hash/ht.h
            src/bioc_parsing_utils.c src/bioc_parsing_utils.h)

add_subdirectory(samples)
